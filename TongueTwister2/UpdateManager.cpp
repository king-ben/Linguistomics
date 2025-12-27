#include <iostream>
#include <unordered_map>
#include "Model.hpp"
#include "Msg.hpp"
#include "ParameterAlignment.hpp"
#include "ParameterExchangeabilities.hpp"
#include "ParameterFrequencies.hpp"
#include "ParameterIndelRates.hpp"
#include "ParameterTree.hpp"
#include "RandomVariable.hpp"
#include "RateMatrix.hpp"
#include "TransitionProbabilityManager.hpp"
#include "UpdateAlignment.hpp"
#include "UpdateBranchLength.hpp"
#include "UpdateExchangeabilities.hpp"
#include "UpdateFrequencies.hpp"
#include "UpdateIndelRates.hpp"
#include "UpdateManager.hpp"
#include "UpdateTopology.hpp"



UpdateManager::UpdateManager(Model* m, RandomVariable* r) : model(m), rng(r) {
    
    tiProbs = model->getTiProbs();
    rateMatrix = model->getRateMatrix();
    const std::vector<Parameter*> parameters = model->getParameters();
        
    for (Parameter* parm : parameters)
        {
        if (dynamic_cast<ParameterAlignment*>(parm) != nullptr)
            {
            UpdateAlignment* updater = new UpdateAlignment(model, rng, dynamic_cast<ParameterAlignment*>(parm), 1.5, -5.0, 0.5);
            updates.push_back(updater);
            alignmentUpdates.push_back(updater);
            }
        else if (dynamic_cast<ParameterExchangeabilities*>(parm) != nullptr)
            {
            UpdateExchangeabilities* updater = new UpdateExchangeabilities(model, rng, dynamic_cast<ParameterExchangeabilities*>(parm));
            updates.push_back(updater);
            otherUpdates.push_back(updater);
            }
        else if (dynamic_cast<ParameterFrequencies*>(parm) != nullptr)
            {
            UpdateFrequencies* updater = new UpdateFrequencies(model, rng, dynamic_cast<ParameterFrequencies*>(parm));
            updates.push_back(updater);
            otherUpdates.push_back(updater);
            }
        else if (dynamic_cast<ParameterIndelRates*>(parm) != nullptr)
            {
            UpdateIndelRates* updater = new UpdateIndelRates(model, rng, dynamic_cast<ParameterIndelRates*>(parm));
            updates.push_back(updater);
            otherUpdates.push_back(updater);
            }
        else if (dynamic_cast<ParameterTree*>(parm) != nullptr)
            {
            // branch length updater
            UpdateBranchLength* updater1 = new UpdateBranchLength(model, rng, dynamic_cast<ParameterTree*>(parm));
            updates.push_back(updater1);
            otherUpdates.push_back(updater1);
            
            // tree topology updater
            const std::vector<ParameterAlignment*> alnVec = model->getParametersOfType<ParameterAlignment>();
            UpdateTopology* updater2 = new UpdateTopology(model, rng, dynamic_cast<ParameterTree*>(parm), alnVec);
            updates.push_back(updater2);
            otherUpdates.push_back(updater2);
            }
        else 
            Msg::error("Could not find update in update manager object");
        }
        
    setProposalProbabilities();
    buildAliasTable();
}

UpdateManager::~UpdateManager(void) {

    for (size_t i=0; i<updates.size(); i++)
        delete updates[i];
}

void UpdateManager::accept(Update* u) {

    u->accept();

    Parameter* parm = u->getUpdatedParameter();
    parm->keep();
    
    if (u->getRateMatrixNeedsUpdate() == true && rateMatrix != nullptr)
        rateMatrix->keep();
    
    // transition probabilities: just clear dirty flags
    if (u->getAllTiprobsNeedUpdate() == true || u->getSingleBranchChanged() == true)
        tiProbs->keep();
}

void UpdateManager::buildAliasTable(void) {

    // Walker, A. J. 1977. An efficient method for generating 
    //    discrete random variables with general distributions.
    //    ACM Transactions on Mathematical Soware, 3(3):253-256
    
    const size_t n = proposalProbabilities.size();
    if (n == 0)
        return;
    
    // initialize tables
    aliasProbability.resize(n);
    aliasIndex.resize(n);
    
    // scale probabilities by n (so average is 1.0)
    std::vector<double> scaledProb(n);
    for (size_t i=0; i<n; i++)
        scaledProb[i] = proposalProbabilities[i] * static_cast<double>(n);
    
    // partition into small (<1) and large (>=1) groups
    // using indices rather than separate vectors to avoid allocations
    std::vector<size_t> small, large;
    small.reserve(n);
    large.reserve(n);
    
    for (size_t i=0; i<n; i++)
        {
        if (scaledProb[i] < 1.0)
            small.push_back(i);
        else
            large.push_back(i);
        }
    
    // build the alias table by pairing small and large probabilities
    while (!small.empty() && !large.empty())
        {
        size_t s = small.back();
        small.pop_back();
        size_t l = large.back();
        large.pop_back();
        
        aliasProbability[s] = scaledProb[s];
        aliasIndex[s] = l;
        
        // transfer probability mass from large to small
        scaledProb[l] = (scaledProb[l] + scaledProb[s]) - 1.0;
        
        if (scaledProb[l] < 1.0)
            small.push_back(l);
        else
            large.push_back(l);
        }
    
    // handle remaining entries (should be ~1.0 due to floating point)
    while (!large.empty())
        {
        size_t l = large.back();
        large.pop_back();
        aliasProbability[l] = 1.0;
        aliasIndex[l] = l;
        }
    
    while (!small.empty())
        {
        size_t s = small.back();
        small.pop_back();
        aliasProbability[s] = 1.0;
        aliasIndex[s] = s;
        }
}

void UpdateManager::print(void) {

    size_t i = 0;
    for (Update* u : updates)
        std::cout << ++i << " -- " << u << " " << u->parameterType() << std::endl;
}

Update* UpdateManager::randomlyChooseUpdate(void) {

    // Walker's alias method: O(1) selection
    // Generate two random numbers: one for bin selection, one for alias decision
    
    const size_t n = updates.size();
    if (n == 0)
        {
        Msg::error("No updates available");
        return nullptr;
        }
    
    // pick a random bin uniformly
    double u = rng->uniformRv();
    size_t bin = static_cast<size_t>(u * n);
    if (bin >= n)  // guard against u == 1.0
        bin = n - 1;
    
    // flip biased coin to decide: original or alias
    double v = rng->uniformRv();
    if (v < aliasProbability[bin])
        return updates[bin];
    else
        return updates[aliasIndex[bin]];
}

void UpdateManager::reject(Update* u) {

    Parameter* parm = u->getUpdatedParameter();
    parm->restore();
    
    if (u->getRateMatrixNeedsUpdate() == true && rateMatrix != nullptr)
        rateMatrix->restore();
    
    // transition probabilities: restore from backup
    if (u->getAllTiprobsNeedUpdate() == true || u->getSingleBranchChanged() == true)
        tiProbs->restore();
}

void UpdateManager::setProposalProbabilities(void) {

    // initialize probability vector (all managed here, not in Update objects)
    proposalProbabilities.resize(updates.size(), 0.0);
    
    // map updates to their indices for quick lookup
    std::unordered_map<Update*, size_t> updateIndex;
    for (size_t i=0; i<updates.size(); i++)
        updateIndex[updates[i]] = i;
    
    // alignments get probability 1.0 each (will be normalized)
    for (Update* u : alignmentUpdates)
        proposalProbabilities[updateIndex[u]] = 1.0;
    
    // calculate sum of alignment probabilities
    double alignmentSum = 0.0;
    for (Update* u : alignmentUpdates)
        alignmentSum += proposalProbabilities[updateIndex[u]];
    
    // other updates share the alignment sum equally
    double otherProb = (otherUpdates.size() > 0) ? alignmentSum / static_cast<double>(otherUpdates.size()) : 0.0;
    for (Update* u : otherUpdates)
        proposalProbabilities[updateIndex[u]] = otherProb;
    
    // normalize to sum to 1.0
    double sum = 0.0;
    for (double p : proposalProbabilities)
        sum += p;
    
    if (sum > 0.0)
        {
        double invSum = 1.0 / sum;
        for (double& p : proposalProbabilities)
            p *= invSum;
        }
}

void UpdateManager::summary(void) {

    size_t longestNameLength = 0;
    for (Update* u : updates)
        {
        std::string name = u->getUpdateName();
        if (name.length() > longestNameLength)
            longestNameLength = name.length();
        }
    
    // map for probability lookup
    std::unordered_map<Update*, size_t> updateIndex;
    for (size_t i = 0; i < updates.size(); i++)
        updateIndex[updates[i]] = i;
        
    std::cout << "   MCMC Update Summary" << std::endl;
    for (Update* u : updates)
        {
        std::cout << "   * ";
        std::string name = u->getUpdateName();
        std::cout << name << " ";
        for (size_t i=0; i<longestNameLength-name.length(); i++)
            std::cout << " ";
        std::cout << std::setw(4) << u->getNumTries() << " ";
        std::cout << std::setw(4) << u->getNumAcceptances() << " ";
        if (u->getNumTries() > 0)
            std::cout << std::fixed << std::setprecision(1) << std::setw(5) << (double)(u->getNumAcceptances() * 100.0) / u->getNumTries() << "% ";
        else 
            std::cout << std::setw(5) << "N/A" << " ";
        std::cout << std::fixed << std::setprecision(4) << std::setw(6) << proposalProbabilities[updateIndex[u]];
        std::cout << std::endl;
        }
    std::cout << std::endl;
}

void UpdateManager::updateDependants(Update* u) {

    // update rate matrix if needed (must happen before transition probabilities)
    if (u->getRateMatrixNeedsUpdate() == true && rateMatrix != nullptr)
        rateMatrix->updateRateMatrix();
    
    // update transition probabilities based on what changed
    if (u->getAllTiprobsNeedUpdate() == true)
        {
        // model parameters changed -> recompute ALL matrices
        tiProbs->updateAllBranches();
        }
    else if (u->getSingleBranchChanged() == true)
        {
        // single branch length changed -> scan subtrees for new branch lengths
        tiProbs->updateBranch();
        }
    // else: no transition probability update needed (e.g., alignment or indel rates)
}
