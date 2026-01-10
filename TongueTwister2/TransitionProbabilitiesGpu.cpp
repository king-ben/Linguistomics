#include <iostream>
#include <algorithm>
#include <vector>
#include "Msg.hpp"
#include "Ctmc.hpp"
#include "Node.hpp"
#include "ParameterTree.hpp"
#include "RateMatrix.hpp"
#include "Threads.hpp"
#include "TransitionMatrixMap.hpp"
#include "TransitionProbabilitiesGpu.hpp"
#include "TransitionProbabilityCalculator.hpp"
#include "Tree.hpp"
#include "GPUMatrixExponentialBatch.hpp"
#include "CtmcAccelerated.hpp"

static const size_t MAX_BATCH_SIZE = 512;



TransitionProbabilitiesGpu::TransitionProbabilitiesGpu(ThreadPool* p, ParameterTree* t, Ctmc* sm, size_t cleanupFreq) : 
    pool(p),
    myTree(t),
    subModel(sm),
    map(nullptr),
    calculatorPool(nullptr),
    calculatorPoolCapacity(0),
    calculatorPoolSize(0),
    maxCalculatorsUsedRecently(0),
    calculatorShrinkCounter(0),
    cleanupFrequency(cleanupFreq),
    cleanupCounter(0),
    matrixPool(sm->getNumStates(), sm->getNumStates(), 128) {

    numStates = subModel->getNumStates();
    map = new TransitionMatrixMap(numStates, 128, &matrixPool);
    
    // initialize compute backend 
    computeBackend = ComputeBackend::Auto;
    batchThreshold = 8;
    gpuBatcher = &GPUMatrixExponentialBatch::getInstance();
    
    // for large matrices, prefer batched approach with lower threshold 
    if (numStates >= 64)
        batchThreshold = 4;
    
    // pre-allocate buffers 
    usedKeysBuffer.reserve(256);
    keysToEraseBuffer.reserve(64);
    batchBranchLengths.reserve(256);
    batchOutputs.reserve(256);
    
    if (cleanupFrequency > 100)
        cleanupFrequency = 100;
    
    rebuildFromTree();
}

TransitionProbabilitiesGpu::~TransitionProbabilitiesGpu(void) {

    delete map;
    delete subModel;
    for (size_t i=0; i<calculatorPoolSize; i++)
        delete calculatorPool[i];
    delete [] calculatorPool;
}

void TransitionProbabilitiesGpu::computeDirtyMatrices(void) {

    const size_t n = map->getDirtyCount();
    if (n == 0)
        return;
    
    if (n > maxCalculatorsUsedRecently)
        maxCalculatorsUsedRecently = n;
    
    // select compute path based on backend setting
    bool useBatched = false;
    
    switch (computeBackend)
        {
        case ComputeBackend::Auto:
            useBatched = (n >= batchThreshold);
            break;
        case ComputeBackend::CPUThreaded:
            useBatched = false;
            break;
        case ComputeBackend::CPUBatched:
        case ComputeBackend::GPUBatched:
            useBatched = true;
            break;
        }
    
    if (useBatched)
        computeDirtyMatricesBatched();
    else
        computeDirtyMatricesThreaded();
}

void TransitionProbabilitiesGpu::computeDirtyMatricesThreaded(void) {

    const size_t n = map->getDirtyCount();
    
    ensureCalculatorPoolCapacity(n);
    
    size_t processed = 0;
    while (processed < n)
        {
        size_t batchSize = n - processed;
        if (batchSize > MAX_BATCH_SIZE)
            batchSize = MAX_BATCH_SIZE;
        
        for (size_t i=0; i<batchSize; i++)
            {
            size_t idx = processed + i;
            double branchLength = map->getDirtyBranchLength(idx);
            DoubleMatrix* matrix = map->getDirtyMatrix(idx);
            calculatorPool[idx]->set(branchLength, matrix);
            pool->pushTask(calculatorPool[idx]);
            }
        
        pool->wait();
        processed += batchSize;
        }
}

void TransitionProbabilitiesGpu::computeDirtyMatricesBatched(void) {

    const size_t n = map->getDirtyCount();
    
    // collect all dirty entries 
    batchBranchLengths.clear();
    batchOutputs.clear();
    batchBranchLengths.reserve(n);
    batchOutputs.reserve(n);
    
    for (size_t i=0; i<n; i++)
        {
        batchBranchLengths.push_back(map->getDirtyBranchLength(i));
        batchOutputs.push_back(map->getDirtyMatrix(i));
        }
    
    // get rate matrix Q from the substitution model
    const double* qData = getRateMatrixData();
    
    if (qData != nullptr)
        {
        // create temporary DoubleMatrix wrapper for Q 
        DoubleMatrix Q(numStates, numStates);
        memcpy(Q.begin(), qData, numStates * numStates * sizeof(double));
        
        // batch compute all matrix exponentials using ThreadPool 
        gpuBatcher->computeBatch(Q, batchBranchLengths, batchOutputs, pool);
        }
    else
        {
        // fallback to threaded approach if we can't get Q directly 
        computeDirtyMatricesThreaded();
        }
}

const double* TransitionProbabilitiesGpu::getRateMatrixData(void) const {

    // try to get Q data from CtmcAccelerated 
    CtmcAccelerated* accModel = dynamic_cast<CtmcAccelerated*>(subModel);
    
    if (accModel != nullptr)
        return accModel->getQData();
    
    return nullptr;
}

void TransitionProbabilitiesGpu::ensureCalculatorPoolCapacity(size_t n) {

    if (n <= calculatorPoolSize)
        return;
    
    if (n > calculatorPoolCapacity)
        {
        size_t newCapacity = calculatorPoolCapacity * 2;
        if (newCapacity < n)
            newCapacity = n;
        if (newCapacity < 16)
            newCapacity = 16;
        
        TransitionProbabilityCalculator** newPool = new TransitionProbabilityCalculator*[newCapacity];
        
        for (size_t i=0; i<calculatorPoolSize; i++)
            newPool[i] = calculatorPool[i];
        
        delete [] calculatorPool;
        calculatorPool = newPool;
        calculatorPoolCapacity = newCapacity;
        }
    
    for (size_t i=calculatorPoolSize; i<n; i++)
        calculatorPool[i] = new TransitionProbabilityCalculator(subModel);
    
    calculatorPoolSize = n;
}

void TransitionProbabilitiesGpu::shrinkCalculatorPoolIfNeeded(void) {

    ++calculatorShrinkCounter;
    if (calculatorShrinkCounter < CALCULATOR_SHRINK_FREQUENCY)
        return;
    
    calculatorShrinkCounter = 0;
    
    size_t targetSize = maxCalculatorsUsedRecently * 2;
    if (targetSize < CALCULATOR_MIN_POOL_SIZE)
        targetSize = CALCULATOR_MIN_POOL_SIZE;
    
    maxCalculatorsUsedRecently = 0;
    
    if (calculatorPoolSize > targetSize && calculatorPoolSize > targetSize + targetSize / 4)
        {
        for (size_t i=targetSize; i<calculatorPoolSize; i++)
            {
            delete calculatorPool[i];
            calculatorPool[i] = nullptr;
            }
        calculatorPoolSize = targetSize;
        }
}

DoubleMatrix* TransitionProbabilitiesGpu::find(double branchLength) {

    return map->find(branchLength);
}

DoubleMatrix& TransitionProbabilitiesGpu::getTransitionProbability(double v) {

    DoubleMatrix* m = map->find(v);
    if (m == nullptr)
        Msg::error("Could not find transition probabilities for branch length");
    return *m;
}

void TransitionProbabilitiesGpu::keep(void) {

    map->keep();
    newEntriesThisProposal.clear();
    
    if (cleanupFrequency > 0)
        {
        ++cleanupCounter;
        if (cleanupCounter >= cleanupFrequency)
            {
            cleanupOrphanedEntries();
            cleanupCounter = 0;
            }
        }
    
    shrinkCalculatorPoolIfNeeded();
}

void TransitionProbabilitiesGpu::restore(void) {

    map->restore();
    
    // remove entries that were created during this proposal
    for (uint64_t key : newEntriesThisProposal)
        map->erase(key);
    newEntriesThisProposal.clear();
}

void TransitionProbabilitiesGpu::updateAllBranches(void) {

    // clear tracking at start of new proposal
    newEntriesThisProposal.clear();
    
    // Scan all trees for branch lengths.
    // For each branch length:
    //   - If entry doesn't exist, create it and track for restore
    //   - Mark the entry for update (backs up and adds to dirty list)
    
    std::vector<Tree*> allTrees;
    myTree->getTrees(allTrees);
    
    for (Tree* t : allTrees)
        {
        const std::vector<Node*>& dp = t->getPostOrder();
        for (Node* p : dp)
            {
            if (p != t->getRoot())
                {
                double bl = p->getBranchLength();
                
                // check if this branch length already has an entry
                DoubleMatrix* matrix = map->find(bl);
                if (matrix == nullptr)
                    {
                    // new branch length - create entry and track for restore
                    uint64_t key = branchLengthToKey(bl);
                    newEntriesThisProposal.push_back(key);
                    map->getOrCreate(bl);
                    }
                
                // mark for update (entry now guaranteed to exist)
                map->markForUpdate(bl);
                }
            }
        }
    
    computeDirtyMatrices();
}

void TransitionProbabilitiesGpu::updateBranch(void) {

    newEntriesThisProposal.clear();
    
    std::vector<Tree*> allTrees;
    myTree->getTrees(allTrees);
    
    for (Tree* t : allTrees)
        {
        const std::vector<Node*>& dp = t->getPostOrder();
        for (Node* p : dp)
            {
            if (p != t->getRoot())
                {
                double bl = p->getBranchLength();
                DoubleMatrix* matrix = map->find(bl);
                if (matrix == nullptr)
                    {
                    uint64_t key = branchLengthToKey(bl);
                    newEntriesThisProposal.push_back(key);
                    
                    matrix = map->getOrCreate(bl);
                    subModel->computeTransitionProbs(bl, matrix, &singleCache);
                    }
                }
            }
        }
}

void TransitionProbabilitiesGpu::rebuildFromTree(void) {

    map->clear();
    newEntriesThisProposal.clear();
    
    std::vector<Tree*> allTrees;
    myTree->getTrees(allTrees);
    
    for (Tree* t : allTrees)
        {
        const std::vector<Node*>& dp = t->getPostOrder();
        for (Node* p : dp)
            {
            if (p != t->getRoot())
                map->getOrCreate(p->getBranchLength());
            }
        }
    
    map->markAllForUpdate();
    computeDirtyMatrices();
    map->keep();
}

void TransitionProbabilitiesGpu::cleanupOrphanedEntries(void) {

    usedKeysBuffer.clear();
    keysToEraseBuffer.clear();
    
    std::vector<Tree*> allTrees;
    myTree->getTrees(allTrees);
    
    for (Tree* t : allTrees)
        {
        const std::vector<Node*>& dp = t->getPostOrder();
        for (Node* p : dp)
            {
            if (p != t->getRoot())
                usedKeysBuffer.push_back(branchLengthToKey(p->getBranchLength()));
            }
        }
    
    std::sort(usedKeysBuffer.begin(), usedKeysBuffer.end());
    auto last = std::unique(usedKeysBuffer.begin(), usedKeysBuffer.end());
    usedKeysBuffer.erase(last, usedKeysBuffer.end());
    
    size_t numEntries = map->size();
    for (size_t i=0; i<numEntries; ++i)
        {
        uint64_t key = map->getEntryKey(i);
        if (!std::binary_search(usedKeysBuffer.begin(), usedKeysBuffer.end(), key))
            keysToEraseBuffer.push_back(key);
        }
    
    for (uint64_t key : keysToEraseBuffer)
        map->erase(key);
}

void TransitionProbabilitiesGpu::setComputeBackend(ComputeBackend backend) {

    computeBackend = backend;
}

bool TransitionProbabilitiesGpu::isGPUAvailable(void) const {

    return gpuBatcher != nullptr && gpuBatcher->isAvailable();
}

const char* TransitionProbabilitiesGpu::getGPUDeviceName(void) const {

    if (gpuBatcher != nullptr)
        return gpuBatcher->getDeviceName();
    return "N/A";
}

void TransitionProbabilitiesGpu::print(void) const {

    std::cout << "TransitionProbabilitiesGPU:" << std::endl;
    std::cout << "  States: " << numStates << std::endl;
    std::cout << "  Compute backend: ";
    switch (computeBackend)
        {
        case ComputeBackend::Auto: std::cout << "Auto"; break;
        case ComputeBackend::CPUThreaded: std::cout << "CPU Threaded"; break;
        case ComputeBackend::CPUBatched: std::cout << "CPU Batched"; break;
        case ComputeBackend::GPUBatched: std::cout << "GPU Batched"; break;
        }
    std::cout << std::endl;
    std::cout << "  Batch threshold: " << batchThreshold << std::endl;
    std::cout << "  GPU available: " << (isGPUAvailable() ? "Yes" : "No") << std::endl;
    if (isGPUAvailable())
        std::cout << "  GPU device: " << getGPUDeviceName() << std::endl;
    std::cout << "  Calculator pool: " << calculatorPoolSize << " (capacity " << calculatorPoolCapacity << ")" << std::endl;
    matrixPool.printStatistics();
    map->print();
}
