#include <algorithm>
#include <iomanip>
#include <map>
#include "Ctmc.hpp"
#include "CtmcAccelerated.hpp"
#include "JsonData.hpp"
#include "LikelihoodCalculator.hpp"
#include "Model.hpp"
#include "Msg.hpp"
#include "ParameterAlignment.hpp"
#include "ParameterExchangeabilities.hpp"
#include "ParameterFrequencies.hpp"
#include "ParameterIndelRates.hpp"
#include "ParameterTree.hpp"
#include "Partition.hpp"
#include "RandomVariable.hpp"
#include "RateMatrix.hpp"
#include "States.hpp"
#include "Threads.hpp"
#include "TransitionProbabilitiesCpu.hpp"
#include "TransitionProbabilitiesGpu.hpp"



Model::Model(RandomVariable* r, SubstitutionModel mt, ThreadPool* p) : 
    rng(r), modelType(mt), pool(p),
    states(nullptr), rateMatrix(nullptr), tiProbs(nullptr), numStates(0),
    cachedTotalLnL(0.0), cacheInitialized(false) {
    
    // initialize the model states
    if (initializeStates() == false)
        Msg::error("Problem initializing the model states");
        
    // initialize the insertion/deletion rates
    if (initializeIndelRates() == false)
        Msg::error("Problem initializing the indel rates");
        
    // initialize the alignment parameters (indel rates must be instantiated first)
    if (initializeAlignments() == false)
        Msg::error("Problem initializing the alignments");
        
    // initialize the tree parameter
    if (initializeTree() == false)
        Msg::error("Problem initalizing the tree");
        
    // initialize the substitution parameters and machinery for calculating transition probabilities
    if (initializeSubstitutionModel() == false)
        Msg::error("Problem initializing transition probabilities");
        
    // initialize the likelihood calculators
    if (initializeCalculators() == false)
        Msg::error("Problem initalizing the likelihood calculators");
        
    /* Initialize likelihood caching infrastructure
       Note: actual cache values are set after first lnLikelihood() call */
    size_t numCognates = calculators.size();
    cachedLnL.resize(numCognates, 0.0);
    dirtyFlags.resize(numCognates, true);  // all dirty initially
    proposedLnL.resize(numCognates, 0.0);
    proposedTotalLnL = 0.0;
    
    // reserve space for temporary storage
    dirtyIndices.reserve(numCognates);
    dirtyCalcs.reserve(numCognates);
}

Model::~Model(void) {

    delete states;
    for (size_t i=0; i<parameters.size(); i++)
        delete parameters[i];
    for (size_t i=0; i<calculators.size(); i++)
        delete calculators[i];
    delete tiProbs;
    if (rateMatrix != nullptr)
        delete rateMatrix;
}

bool Model::alignmentsAreConsistent(std::vector<ParameterAlignment*>& alignments) {

    JsonData& jd = JsonData::jsonInstance();
    if (jd.hasKey("Taxa") == false)
        Msg::error("Could not find the \"Taxa\" in JSON object");
    std::vector<std::string> canonicalTaxonList = jd.getValue<std::vector<std::string>>("Taxa");

    // check that the taxa in each alignment are in the same order as in the canonical taxon list
    for (size_t i=0; i<alignments.size(); i++)
        {
        const std::vector<std::string>& tNames = alignments[i]->getTaxonNames();
        size_t k = 0;
        for (size_t j=0; j<tNames.size(); j++)
            {
            size_t idx = 0;
            for (idx=0; idx<canonicalTaxonList.size(); idx++)
                {
                if (tNames[j] == canonicalTaxonList[idx])
                    break;
                }
            
            if (idx < k)
                Msg::error("Scrambled alignment");
            if (idx == canonicalTaxonList.size())
                Msg::error("Could not find taxon in canonical taxon list");
            k = idx;
            }
        }

    return true;
}

size_t Model::getCognateIndex(LikelihoodCalculator* calc) const {

    for (size_t i = 0; i < calculators.size(); i++)
        {
        if (calculators[i] == calc)
            return i;
        }
    return SIZE_MAX;
}

size_t Model::getNumDirtyCognates(void) const {

    size_t count = 0;
    for (bool dirty : dirtyFlags)
        {
        if (dirty)
            count++;
        }
    return count;
}

bool Model::hasAnyCognateDirty(void) const {

    for (bool dirty : dirtyFlags)
        {
        if (dirty)
            return true;
        }
    return false;
}

bool Model::initializeAlignments(void) {
    
    // decide whether we only sample complete alignments
    bool onlyCompleteAlignments = false;
    UserSettings& settings = UserSettings::userSettings();
    if (settings.getUseOnlyCompleteWords() == true)
        onlyCompleteAlignments = true;
    
    // check that the necessary keys are in the JSON object
    JsonData& jd = JsonData::jsonInstance();
    if (jd.hasKey("Taxa") == false)
        Msg::error("Could not find taxa list in JSON object");
    if (jd.hasKey("Words") == false)
        Msg::error("Could not find word information in the JSON file");
    if (jd.hasKey("NumberOfStates") == false)
        Msg::error("Could not find the number of states in JSON object");

    // read the number of states
    numStates = jd.getValue<size_t>("NumberOfStates");

    // read the taxon names
    std::vector<std::string> canonicalTaxonList = jd.getValue<std::vector<std::string>>("Taxa");
    if (canonicalTaxonList.size() > sizeof(unsigned)*8) // cannot represent taxon composition using unsigned
        Msg::error("This program only allows for the analysis of up to " + std::to_string(sizeof(unsigned)*8) + " taxa");

    // read the individual cognates
    nlohmann::json jWords = *jd.getJsonPointer("Words");
    size_t numWords = jWords.size();
    std::vector<ParameterAlignment*> alignments;
    for (size_t i=0; i<numWords; i++)
        {
        nlohmann::json jw = jWords[i];
        std::string alnName = jw["Name"];
        ParameterAlignment* aln = new ParameterAlignment(this, rng, alnName);
        aln->initializeFromJson(jw, numStates, canonicalTaxonList);
        if (aln->getNumTaxa() != canonicalTaxonList.size() && onlyCompleteAlignments == true)
            {
            delete aln;
            }
        else
            {
            parameters.push_back(aln);
            alignments.push_back(aln);
            //aln->print();
            }
        }
        
    // find the unique taxon combinations
    for (ParameterAlignment* p : alignments)
        {
        const unsigned& mask = p->getTaxonMask();
        std::set<unsigned>::iterator it = uniqueTaxonCombinations.find(mask);
        if (it == uniqueTaxonCombinations.end())
            uniqueTaxonCombinations.insert(mask);
        }
        
    // check the consistency of the alignments
    if (alignmentsAreConsistent(alignments) == false)
        Msg::error("Alignments are inconsistent");
        
    // print a summary of the alignments
    std::cout << "   Cognate Alignment Information" << std::endl;
    std::cout << "   * Number of cognate alignments: " << alignments.size() << std::endl;
    std::cout << "   * Taxon completeness:" << std::endl;
    std::map<std::string, int> taxonCompleteness;
    for (size_t i=0; i<canonicalTaxonList.size(); i++)
        taxonCompleteness.insert( std::make_pair(canonicalTaxonList[i],0) );
    for (ParameterAlignment* p : alignments)
        {
        const std::vector<std::string>& tNames = p->getTaxonNames();
        for (size_t i = 0; i<tNames.size(); i++)
            {
            std::map<std::string, int>::iterator it = taxonCompleteness.find(tNames[i]);
            if (it == taxonCompleteness.end())
                Msg::error("Could not find taxon \"" + tNames[i] + "\" in canonical taxon list");
            it->second++;
            }
        }
    size_t maxNameLen = 0;
    for (size_t i=0; i<canonicalTaxonList.size(); i++)
        {
        if (canonicalTaxonList[i].length() > maxNameLen)
            maxNameLen = canonicalTaxonList[i].length();
        }
    for (auto [key,val] : taxonCompleteness)   
        {
        std::cout << "     " << key << " ";
        for (size_t i=0; i<maxNameLen-key.length(); i++)
            std::cout << " ";
        std::cout << std::fixed << std::setprecision(2);
        std::cout << std::setw(6) << ((double)val/alignments.size()) * 100.0 << "%" << std::endl;
        }
    std::vector<int> alignmentSize(canonicalTaxonList.size()+1, 0);
    for (ParameterAlignment* p : alignments)
        alignmentSize[p->getNumTaxa()]++;
    std::cout << "   * Number of unique taxon combinations: " << uniqueTaxonCombinations.size() << std::endl;
    std::cout << "   * Cognate size distribution:" << std::endl;
    std::cout << "     x-axis: Percentage of alignments" << std::endl;
    std::cout << "     y-axis: Number of taxa" << std::endl;
    std::map<int,int> plotData;
    for (size_t i=0; i<alignmentSize.size(); i++)
        plotData.insert( std::make_pair(i,alignmentSize[i]) );
    plotAscii(plotData);
    std::cout << std::endl;
    
    return true;
}

bool Model::initializeCalculators(void) {

    // find the tree parameter
    ParameterTree* t = findParameter<ParameterTree>();
    if (t == nullptr)
        Msg::error("Could not find tree parameter");
        
    // find the indel rates
    ParameterIndelRates* r = findParameter<ParameterIndelRates>();
    if (r == nullptr)
        Msg::error("Could not find indel rates parameter");
 
     // find the equilibirum frequencies
    ParameterFrequencies* f = findParameter<ParameterFrequencies>();
       
    // assign a likelihood calculator to each alignment
    size_t cognateIdx = 0;
    for (size_t i=0; i<parameters.size(); i++)
        {
        if (isType<ParameterAlignment>(parameters[i]))
            {
            ParameterAlignment* aln = dynamic_cast<ParameterAlignment*>(parameters[i]);
            LikelihoodCalculator* calc = new LikelihoodCalculator(tiProbs, aln, t, r, f);
            calculators.push_back(calc);
            
            // store the cognate index in the alignment parameter for quick lookup
            aln->setCognateIndex(cognateIdx);
            cognateIdx++;
            }
        }

    return true;
}

bool Model::initializeIndelRates(void) {

    ParameterIndelRates* rates = new ParameterIndelRates(this, rng, "Insertion/Deletion rates", 2.0, 2.0);
    parameters.push_back(rates);
    return true;
}

void Model::initializeLikelihoodCache(void) {

    /* This should be called after the first lnLikelihood() computation
       to populate the cache with valid values */
    if (cacheInitialized)
        return;
        
    /* The cache is already populated by the first lnLikelihood() call
       so mark it as initialized */
    cacheInitialized = true;
    
    // clear all dirty flags since we just computed everything
    for (size_t i = 0; i < dirtyFlags.size(); i++)
        dirtyFlags[i] = false;
}

bool Model::initializeStates(void) {

    states = new States;
    numStates = states->getNumStates();
    if (!states)
        return false;
    states->print();
    return true;
}

bool Model::initializeSubstitutionModel(void) {

    // find the tree parameter
    ParameterTree* t = findParameter<ParameterTree>();
    if (t == nullptr)
        Msg::error("Could not find tree parameter");
        
    std::cout << "   Substitution Model" << std::endl;
    
    // get compute backend preference from settings
    UserSettings& settings = UserSettings::userSettings();
    MatrixExpBackend matExpBackend = settings.getMatrixExpBackend();
        
    // initialize the substitution model
    Ctmc* subModel = nullptr; // note that subModel is not owned by Model (see below, near end of function)
    rateMatrix = nullptr;
    if (modelType == jc69)
        {
        std::cout << "   * Model Type: Jukes-Cantor" << std::endl;
        subModel = new CtmcJC69(numStates);
        }
    else if (modelType == f81)
        {
        std::cout << "   * Model Type: Felsenstein 1981" << std::endl;
        ParameterFrequencies* freqs = new ParameterFrequencies(this, rng, "Segment frequencies", numStates);
        parameters.push_back(freqs);
        subModel = new CtmcF81(numStates, freqs);
        }
    else if (modelType == gtr)
        {
        std::cout << "   * Model Type: GTR" << std::endl;
        ParameterFrequencies* freqs = new ParameterFrequencies(this, rng, "Segment frequencies", numStates);
        ParameterExchangeabilities* rates = new ParameterExchangeabilities(this, rng, "Exchangeability rates", numStates);
        parameters.push_back(freqs);
        parameters.push_back(rates);
        rateMatrix = new RateMatrixGTR(numStates, rates, freqs);
        rateMatrix->updateRateMatrix();
        subModel = new CtmcAccelerated(numStates, rateMatrix);
        if (rateMatrix == nullptr)
            Msg::error("Failed to instantiate rate matrix");
        }
    else if (modelType == custom)
        {
        std::cout << "   * Model Type: Natural Class" << std::endl;
        Partition* part = states->getPartition();
        ParameterFrequencies* freqs = new ParameterFrequencies(this, rng, "Segment frequencies", numStates);
        ParameterExchangeabilities* rates = new ParameterExchangeabilities(this, rng, "Exchangeability rates", numStates, part->numSubsets());
        parameters.push_back(freqs);
        parameters.push_back(rates);
        rateMatrix = new RateMatrixCustom(numStates, rates, freqs, states->getPartition());
        rateMatrix->updateRateMatrix();
        subModel = new CtmcAccelerated(numStates, rateMatrix);
        if (rateMatrix == nullptr)
            Msg::error("Failed to instantiate rate matrix");
        }
    else
        Msg::error("Unknown model type");
    if (subModel == nullptr)
        Msg::error("Failed to create substitution model");

    // instantiate the object that manages the transition probabilities
    // based on the selected compute backend
    switch (matExpBackend)
        {
        case cpuThreaded:
            std::cout << "   * Compute Backend: CPU Threaded (per-branch threading)" << std::endl;
            tiProbs = new TransitionProbabilitiesCpu(pool, t, subModel);
            break;
            
        case cpuBatched:
            std::cout << "   * Compute Backend: CPU Batched (Accelerate BLAS)" << std::endl;
            {
            TransitionProbabilitiesGpu* gpuTiProbs = new TransitionProbabilitiesGpu(pool, t, subModel);
            gpuTiProbs->setComputeBackend(ComputeBackend::CPUBatched);
            tiProbs = gpuTiProbs;
            }
            break;
            
        case gpuBatched:
            std::cout << "   * Compute Backend: GPU Batched" << std::endl;
            {
            TransitionProbabilitiesGpu* gpuTiProbs = new TransitionProbabilitiesGpu(pool, t, subModel);
            gpuTiProbs->setComputeBackend(ComputeBackend::GPUBatched);
            if (gpuTiProbs->isGPUAvailable())
                std::cout << "   * GPU Device: " << gpuTiProbs->getGPUDeviceName() << std::endl;
            else
                std::cout << "   * GPU not available, using CPU fallback" << std::endl;
            tiProbs = gpuTiProbs;
            }
            break;
            
        case autoBackend:
        default:
            std::cout << "   * Compute Backend: Auto (adaptive)" << std::endl;
            {
            TransitionProbabilitiesGpu* gpuTiProbs = new TransitionProbabilitiesGpu(pool, t, subModel);
            gpuTiProbs->setComputeBackend(ComputeBackend::Auto);
            tiProbs = gpuTiProbs;
            }
            break;
        }
    
    tiProbs->updateAllBranches();
    
    std::cout << std::endl;
    
    return true;
}

bool Model::initializeTree(void) {

    ParameterTree* t = new ParameterTree(this, rng, "Tree");
    t->initialize(uniqueTaxonCombinations);
    parameters.push_back(t);
    std::cout << "   Tree information" << std::endl;
    std::cout << "   * Number of taxa: " << t->getNumTaxa() << std::endl;
    std::cout << std::endl;
    //t->print();
    return true;
}

void Model::keepLikelihoodCache(void) {

    /* Called when MCMC proposal is accepted.
       Copy the proposed new values into the permanent cache.
       The proposedLnL values become the new cached values. */
    for (size_t i = 0; i < dirtyFlags.size(); i++)
        {
        if (dirtyFlags[i])
            {
            cachedLnL[i] = proposedLnL[i];
            }
        dirtyFlags[i] = false;
        }
    cachedTotalLnL = proposedTotalLnL;
}

#if 1
double Model::lnLikelihood(void) {

    /* Instead of recomputing ALL cognate likelihoods every MCMC cycle,
       we only recompute those that are marked dirty.

       This method does NOT modify cachedLnL or cachedTotalLnL.
       It stores proposed values in proposedLnL and proposedTotalLnL.
       The cache is only updated when keepLikelihoodCache() is called (accept). */
    
    // if cache not initialized, compute everything (first call)
    if (!cacheInitialized)
        {
        // first call: compute all likelihoods and populate cache
        for (LikelihoodCalculator* like : calculators)
            pool->pushTask(like);
        pool->wait();
        
        cachedTotalLnL = 0.0;
        for (size_t i = 0; i < calculators.size(); i++)
            {
            cachedLnL[i] = calculators[i]->getResult();
            cachedTotalLnL += cachedLnL[i];
            dirtyFlags[i] = false;
            }
        
        // initialize proposed values to match cached values
        proposedLnL = cachedLnL;
        proposedTotalLnL = cachedTotalLnL;
        
        cacheInitialized = true;
        return cachedTotalLnL;
        }
    
    // find which calculators are dirty and need recomputation
    dirtyCalcs.clear();
    dirtyIndices.clear();
    
    for (size_t i = 0; i < calculators.size(); i++)
        {
        if (dirtyFlags[i])
            {
            dirtyCalcs.push_back(calculators[i]);
            dirtyIndices.push_back(i);
            }
        }
    
    // if nothing is dirty, return cached total (no proposal active)
    if (dirtyCalcs.empty())
        return cachedTotalLnL;
    
    // push only dirty calculators to thread pool
    for (LikelihoodCalculator* calc : dirtyCalcs)
        pool->pushTask(calc);
    pool->wait();
    
    /* Compute proposed total by starting from cached total
       and adjusting for the cognates that were recomputed */
    proposedTotalLnL = cachedTotalLnL;
    
    for (size_t j = 0; j < dirtyCalcs.size(); j++)
        {
        size_t i = dirtyIndices[j];
        double oldVal = cachedLnL[i];
        double newVal = dirtyCalcs[j]->getResult();
        
        // store proposed value (don't modify cached value yet!)
        proposedLnL[i] = newVal;
        
        // update proposed total
        proposedTotalLnL += (newVal - oldVal);
        }
    
    return proposedTotalLnL;
}
#else
double Model::lnLikelihood(void) {

    // Bypass caching - compute everything every time
    for (LikelihoodCalculator* like : calculators)
        pool->pushTask(like);
    pool->wait();
    
    double lnL = 0.0;
    for (LikelihoodCalculator* like : calculators)
        lnL += like->getResult();
    return lnL;
}
#endif

double Model::lnPrior(void) {

    double lnP = 0.0;
    for (Parameter* parm : parameters)
        lnP += parm->lnPriorProbability();
    return lnP;
}

void Model::markAllCognatesDirty(void) {

    for (size_t i = 0; i < dirtyFlags.size(); i++)
        dirtyFlags[i] = true;
}

void Model::markCognateDirty(size_t cognateIdx) {

    if (cognateIdx < dirtyFlags.size())
        dirtyFlags[cognateIdx] = true;
}

void Model::plotAscii(const std::map<int,int>& data) {

    if (data.empty()) 
        Msg::error("No data to plot");

    int sum = 0;
    for (const auto& [x, y] : data) 
        sum += y;
        
    // determine bounds
    int minY = data.begin()->second;
    int maxY = data.begin()->second;
    for (const auto& [x, y] : data) 
        {
        minY = std::min(minY, y);
        maxY = std::max(maxY, y);
        }
        
    for (const auto& [key,val] : data)
        {
        std::cout << "     " << std::setw(2) << key << " |";
        int nTicks = static_cast<int>(std::round((((double)val/sum) * 100.0) / 5.0));
        
        for (int j=0; j<nTicks; j++)
            std::cout << "*";
        for (int j=0; j<20-nTicks; j++)
            std::cout << " ";
        std::cout << " (" << std::fixed << std::setprecision(1) << ((double)val/sum) * 100.0  << "%)" << std::endl;
        }
    std::cout << "        +";
    for (int j=0; j<20; j++)
        std::cout << "-";
    std::cout << std::endl;

}

void Model::restoreLikelihoodCache(void) {

    /* Called when MCMC proposal is rejected

       Since lnLikelihood() does not modify cachedLnL or cachedTotalLnL,
       there's nothing to restore. The cached values are still valid.

       We just need to:
       1. Clear the dirty flags
       2. Reset proposedLnL to match cachedLnL for dirty entries
          (so they're ready for the next proposal) */
    
    for (size_t i = 0; i < dirtyFlags.size(); i++)
        {
        if (dirtyFlags[i])
            {
            proposedLnL[i] = cachedLnL[i];  // reset proposed to cached
            dirtyFlags[i] = false;
            }
        }
    
    proposedTotalLnL = cachedTotalLnL;
}
