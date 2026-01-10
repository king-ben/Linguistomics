#include <iostream>
#include <unordered_set>
#include <vector>
#include "Msg.hpp"
#include "Ctmc.hpp"
#include "Node.hpp"
#include "ParameterTree.hpp"
#include "Threads.hpp"
#include "TransitionMatrixMap.hpp"
#include "TransitionProbabilitiesCpu.hpp"
#include "TransitionProbabilityCalculator.hpp"
#include "Tree.hpp"

// Maximum batch size for parallel task processing.
// Must be less than ThreadPool::queueCapacity (1024) to avoid overflow.
static const size_t MAX_BATCH_SIZE = 512;



TransitionProbabilitiesCpu::TransitionProbabilitiesCpu(ThreadPool* p, ParameterTree* t, Ctmc* sm, size_t cleanupFreq) : 
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
    
    // pre-allocate buffers for cleanup to avoid repeated allocations
    usedKeysBuffer.reserve(256);
    keysToEraseBuffer.reserve(64);
    
    // reduce default cleanup frequency for faster memory reclamation
    if (cleanupFrequency > 100)
        cleanupFrequency = 100;
    
    rebuildFromTree();
}

TransitionProbabilitiesCpu::~TransitionProbabilitiesCpu(void) {

    delete map;
    delete subModel;
    for (size_t i=0; i<calculatorPoolCapacity; i++)
        delete calculatorPool[i];
    delete [] calculatorPool;
}

void TransitionProbabilitiesCpu::computeDirtyMatrices(void) {

    const size_t n = map->getDirtyCount();
    if (n == 0)
        return;
    
    ensureCalculatorPoolCapacity(n);
    
    // track max usage for shrinking decisions
    if (n > maxCalculatorsUsedRecently)
        maxCalculatorsUsedRecently = n;
    
    // process tasks in batches to avoid ThreadPool queue overflow.
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

void TransitionProbabilitiesCpu::ensureCalculatorPoolCapacity(size_t n) {

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
        
        for (size_t i = 0; i < calculatorPoolSize; i++)
            newPool[i] = calculatorPool[i];
        
        delete [] calculatorPool;
        calculatorPool = newPool;
        calculatorPoolCapacity = newCapacity;
        }
    
    for (size_t i=calculatorPoolSize; i<n; i++)
        calculatorPool[i] = new TransitionProbabilityCalculator(subModel);
    
    calculatorPoolSize = n;
}

// Each TransitionProbabilityCalculator contains a MathCache with 16 matrices,
// so this can reclaim significant memory when dirty counts are low.
void TransitionProbabilitiesCpu::shrinkCalculatorPoolIfNeeded(void) {

    ++calculatorShrinkCounter;
    if (calculatorShrinkCounter < CALCULATOR_SHRINK_FREQUENCY)
        return;
    
    calculatorShrinkCounter = 0;
    
    // determine target size based on recent usage with 2x headroom
    size_t targetSize = maxCalculatorsUsedRecently * 2;
    if (targetSize < CALCULATOR_MIN_POOL_SIZE)
        targetSize = CALCULATOR_MIN_POOL_SIZE;
    
    // Reset usage tracking for next period
    maxCalculatorsUsedRecently = 0;
    
    // only shrink if we can reclaim significant memory (at least 25%)
    if (calculatorPoolSize > targetSize && calculatorPoolSize > targetSize + targetSize / 4)
        {
        // delete excess calculators (each has a MathCache with 16 matrices)
        for (size_t i=targetSize; i<calculatorPoolSize; i++)
            {
            delete calculatorPool[i];
            calculatorPool[i] = nullptr;
            }
        calculatorPoolSize = targetSize;
        }
}

DoubleMatrix* TransitionProbabilitiesCpu::find(double branchLength) {

    return map->find(branchLength);
}

DoubleMatrix& TransitionProbabilitiesCpu::getTransitionProbability(double v) {

    DoubleMatrix* m = map->find(v);
    if (m == nullptr)
        {
        // print what's in the map
        int count = 0;
        std::vector<Tree*> allTrees;
        myTree->getTrees(allTrees);
        for (Tree* t : allTrees)
            {
            const std::vector<Node*>& dp = t->getPostOrder();
            for (Node* p : dp)
                {
                if (p != t->getRoot() && count < 20)
                    {
                    count++;
                    }
                }
            }
        Msg::error("Could not find transition probabilities for branch length");
        }
    return *m;
}

void TransitionProbabilitiesCpu::keep(void) {

    map->keep();
    
    newEntriesThisProposal.clear();
    
    // This prevents unbounded growth of the map from accumulated
    // branch lengths that are no longer used by any tree.
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

void TransitionProbabilitiesCpu::print(void) const {

    std::cout << "TransitionProbabilities:" << std::endl;
    std::cout << "  States: " << numStates << std::endl;
    std::cout << "  Calculator pool: " << calculatorPoolSize << " (capacity " << calculatorPoolCapacity << ")" << std::endl;
    std::cout << "  Max calculators used recently: " << maxCalculatorsUsedRecently << std::endl;
    std::cout << "  Cleanup frequency: " << cleanupFrequency << std::endl;
    std::cout << "  Cleanup counter: " << cleanupCounter << std::endl;
    std::cout << "  New entries this proposal: " << newEntriesThisProposal.size() << std::endl;
    matrixPool.printStatistics();
    map->print();
}

void TransitionProbabilitiesCpu::rebuildFromTree(void) {

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
    
    // compute all transition probabilities
    map->markAllForUpdate();
    computeDirtyMatrices();
    
    // accept initial state
    map->keep();
}

void TransitionProbabilitiesCpu::restore(void) {

    map->restore();
    
    // remove entries that were created during this proposal.
    for (uint64_t key : newEntriesThisProposal)
        map->erase(key);
    newEntriesThisProposal.clear();
}

void TransitionProbabilitiesCpu::updateAllBranches(void) {

    // clear tracking at start of new proposal
    newEntriesThisProposal.clear();
    
    // Scan all trees for branch lengths.
    // For each branch length:
    //   - If entry doesn't exist, create it and track for restore
    //   - Mark the entry for update (backs up and adds to dirty list)
    
    std::vector<Tree*> allTrees;
    myTree->getTrees(allTrees);
        
    for (size_t treeIdx=0; treeIdx<allTrees.size(); treeIdx++)
        {
        Tree* t = allTrees[treeIdx];
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
        
    // compute only the entries we marked
    computeDirtyMatrices();
}

void TransitionProbabilitiesCpu::updateBranch(void) {

    // clear tracking at start of new proposal
    newEntriesThisProposal.clear();
    
    // Scan all subtrees for any missing branch lengths.
    // This handles the case where full tree branch changes propagate
    // to subtrees via BranchMapping, resulting in new subtree branch lengths.
    
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
                    // track this new entry so we can remove it on restore.
                    uint64_t key = branchLengthToKey(bl);
                    newEntriesThisProposal.push_back(key);
                    
                    // new branch length in subtree - create and compute
                    matrix = map->getOrCreate(bl);
                    subModel->computeTransitionProbs(bl, matrix, &singleCache);
                    }
                }
            }
        }
}

void TransitionProbabilitiesCpu::cleanupOrphanedEntries(void) {

    // clear and reuse pre-allocated buffers
    usedKeysBuffer.clear();
    keysToEraseBuffer.clear();
    
    // collect all branch lengths currently in use
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
                usedKeysBuffer.push_back(branchLengthToKey(bl));
                }
            }
        }
    
    // sort for binary search (faster than unordered_set for small sizes)
    std::sort(usedKeysBuffer.begin(), usedKeysBuffer.end());
    
    // remove duplicates
    auto last = std::unique(usedKeysBuffer.begin(), usedKeysBuffer.end());
    usedKeysBuffer.erase(last, usedKeysBuffer.end());
    
    // collect keys to erase (can't erase while iterating)
    size_t numEntries = map->size();
    for (size_t i=0; i<numEntries; ++i)
        {
        uint64_t key = map->getEntryKey(i);
        if (!std::binary_search(usedKeysBuffer.begin(), usedKeysBuffer.end(), key))
            keysToEraseBuffer.push_back(key);
        }
    
    // now erase the orphaned entries
    for (uint64_t key : keysToEraseBuffer)
        map->erase(key);
}
