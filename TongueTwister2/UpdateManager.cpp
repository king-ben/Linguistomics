#include <iostream>
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
#include "UpdateExchangeabilities.hpp"
#include "UpdateFrequencies.hpp"
#include "UpdateIndelRates.hpp"
#include "UpdateManager.hpp"
#include "UpdateTree.hpp"



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
            UpdateTree* updater = new UpdateTree(model, rng, dynamic_cast<ParameterTree*>(parm));
            updates.push_back(updater);
            otherUpdates.push_back(updater);
            }
        else 
            Msg::error("Could not find update in update manager object");
        }
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

void UpdateManager::print(void) {

    size_t i = 0;
    for (Update* u : updates)
        std::cout << ++i << " -- " << u << " " << u->parameterType() << std::endl;
}

Update* UpdateManager::randomlyChooseUpdate(void) {

    double alignmentUpdateProb = 0.5;

    Update* u = nullptr;
    if (rng->uniformRv() < alignmentUpdateProb)
        u = alignmentUpdates[static_cast<int>(rng->uniformRv()*alignmentUpdates.size())];
    else 
        u = otherUpdates[static_cast<int>(rng->uniformRv()*otherUpdates.size())];
    if (u == nullptr)
        Msg::error("Failed to randomly choose a parameter to update");
        
    return u;
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

void UpdateManager::summary(void) {

    size_t longestNameLength = 0;
    for (Update* u : updates)
        {
        std::string name = u->getUpdateName();
        if (name.length() > longestNameLength)
            longestNameLength = name.length();
        }
        
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
