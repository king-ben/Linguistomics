#include <iomanip>
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
#include "TransitionProbabilities.hpp"
#include "UpdateAlignment.hpp"
#include "UpdateBranchLength.hpp"
#include "UpdateExchangeabilities.hpp"
#include "UpdateFrequencies.hpp"
#include "UpdateIndelRates.hpp"
#include "UpdateManager.hpp"
#include "UpdateTopology.hpp"
#include "UserSettings.hpp"
#include "SubModel.hpp"



UpdateManager::UpdateManager(Model* m, RandomVariable* r) : model(m), rng(r) {

    UserSettings& settings = UserSettings::userSettings();
    bool isMultiFamily = settings.getIsMultiFamily();
    bool linked = (settings.getMultiTreeModel() == linkedFamilies);
    
    if (isMultiFamily)
        {
        // collect alignment updates from all families
        const std::vector<Parameter*>& parameters = model->getParameters();
        for (Parameter* parm : parameters)
            {
            if (ParameterAlignment* p = dynamic_cast<ParameterAlignment*>(parm))
                {
                UpdateAlignment* updater = new UpdateAlignment(model, rng, p);
                updates.push_back(updater);
                alignmentUpdates.push_back(updater);
                }
            }
        
        // one set of tree updates per family
        int numFamilies = model->getNumFamilies();
        for (int i=0; i<numFamilies; i++)
            {
            SubModel* sm = model->getSubModel(i);
            ParameterTree* tree = sm->getTree();
            
            UpdateBranchLength* upd = new UpdateBranchLength(model, rng, tree);
            updates.push_back(upd);
            otherUpdates.push_back(upd);
            }
        
        // substitution parameter updates: one set (linked) or one per family (unlinked)
        int numParamSets = linked ? 1 : model->getNumFamilies();
        for (Parameter* parm : model->getParameters())
            {
            if (ParameterExchangeabilities* p =
                    dynamic_cast<ParameterExchangeabilities*>(parm))
                {
                UpdateExchangeabilities* upd =
                    new UpdateExchangeabilities(model, rng, p);
                updates.push_back(upd);
                otherUpdates.push_back(upd);
                }
            else if (ParameterFrequencies* p =
                    dynamic_cast<ParameterFrequencies*>(parm))
                {
                UpdateFrequencies* upd = new UpdateFrequencies(model, rng, p);
                updates.push_back(upd);
                otherUpdates.push_back(upd);
                }
            else if (ParameterIndelRates* p =
                    dynamic_cast<ParameterIndelRates*>(parm))
                {
                UpdateIndelRates* upd = new UpdateIndelRates(model, rng, p);
                updates.push_back(upd);
                otherUpdates.push_back(upd);
                }
            }
        }
    else
        {
        // existing single-family logic unchanged
        const std::vector<Parameter*> parameters = model->getParameters();
        for (Parameter* parm : parameters)
            {
            if (ParameterAlignment* p = dynamic_cast<ParameterAlignment*>(parm))
                {
                UpdateAlignment* updater = new UpdateAlignment(model, rng, p);
                updates.push_back(updater);
                alignmentUpdates.push_back(updater);
                }
            }
        for (Parameter* parm : parameters)
            {
            if (ParameterExchangeabilities* p =
                    dynamic_cast<ParameterExchangeabilities*>(parm))
                {
                UpdateExchangeabilities* updater =
                    new UpdateExchangeabilities(model, rng, p);
                updates.push_back(updater);
                otherUpdates.push_back(updater);
                }
            else if (ParameterFrequencies* p =
                    dynamic_cast<ParameterFrequencies*>(parm))
                {
                UpdateFrequencies* updater = new UpdateFrequencies(model, rng, p);
                updates.push_back(updater);
                otherUpdates.push_back(updater);
                }
            else if (ParameterIndelRates* p =
                    dynamic_cast<ParameterIndelRates*>(parm))
                {
                UpdateIndelRates* updater = new UpdateIndelRates(model, rng, p);
                updates.push_back(updater);
                otherUpdates.push_back(updater);
                }
            else if (ParameterTree* p = dynamic_cast<ParameterTree*>(parm))
                {
                UpdateBranchLength* updater1 = new UpdateBranchLength(model, rng, p);
                updates.push_back(updater1);
                otherUpdates.push_back(updater1);
                }
            }
        }

    for (size_t i=0; i<updates.size(); i++)
        updateIndex[updates[i]] = i;
    
    setProposalProbabilities();
    buildAliasTable();
}

UpdateManager::~UpdateManager(void) {

    for (size_t i=0; i<updates.size(); i++)
        delete updates[i];
}

void UpdateManager::accept(Update* u) {

    // increment acceptance count
    updateInfo.accept(u);

    // keep the primary updated parameter
    Parameter* parm = u->getUpdatedParameter();
    parm->keep();
    
    // keep any additional parameters that were modified (e.g., alignments during topology update)
    for (Parameter* p : u->getAdditionalModifiedParameters())
        p->keep();
    
    if (u->getRateMatrixNeedsUpdate() == true)
        {
        RateMatrix* rm = rateMatrixForUpdate(u);
        if (rm != nullptr)
            rm->keep();
        }
    
    // transition probabilities: just clear dirty flags
    if (u->getAllTiprobsNeedUpdate() == true || u->getSingleBranchChanged() == true)
        {
        if (model->getNumFamilies() > 0)
            {
            if (u->getAllTiprobsNeedUpdate() == true)
                {
                for (int i=0; i<model->getNumFamilies(); i++)
                    model->getSubModel(i)->getTiProbs()->keep();
                }
            else if (u->getSingleBranchChanged() == true)
                {
                Parameter* parm = u->getUpdatedParameter();

                if (ParameterTree* treeParm = dynamic_cast<ParameterTree*>(parm))
                    {
                    for (int i=0; i<model->getNumFamilies(); i++)
                        {
                        if (model->getSubModel(i)->getTree() == treeParm)
                            {
                            model->getSubModel(i)->getTiProbs()->keep();
                            break;
                            }
                        }
                    }
                }
            }
        else
            {
            tiProbs->keep();
            }
        }
        
    // keep the likelihood cache
    model->keepLikelihoodCache();
}

void UpdateManager::buildAliasTable(void) {

    /* Walker, A. J. 1977. An efficient method for generating 
          discrete random variables with general distributions.
          ACM Transactions on Mathematical Software, 3(3):253-256 */
    
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

void UpdateManager::markCognatesDirty(Update* u) {

    /* Mark cognates dirty based on the update type

       This determines which cognate likelihoods need to be recalculated.

       Rules:
       - UpdateAlignment: Only the specific cognate being updated
       - UpdateIndelRates: ALL cognates (TKF91 parameters are shared)
       - UpdateFrequencies: ALL cognates (stationary frequencies are shared)
       - UpdateExchangeabilities: ALL cognates (rate matrix is shared)
       - UpdateBranchLength: ALL cognates (tree is shared)
       - UpdateTopology: ALL cognates (tree is shared)

       This is called AFTER the update is applied but BEFORE lnLikelihood(). */
    
    Parameter* parm = u->getUpdatedParameter();
    
    // check if this is an alignment update (single cognate affected)
    if (ParameterAlignment* alnParm = dynamic_cast<ParameterAlignment*>(parm))
        {
        // only mark the specific cognate that was updated
        size_t cognateIdx = alnParm->getCognateIndex();
       if (cognateIdx != SIZE_MAX)
            {
            model->markCognateDirty(cognateIdx);
            }
        else
            {
            // fallback: if cognateIndex not set, mark all dirty
            model->markAllCognatesDirty();
            }
        }
    else
        {
        /* All other updates affect all cognates
           This includes:
           - ParameterIndelRates
           - ParameterFrequencies
           - ParameterExchangeabilities
           - ParameterTree (branch lengths and topology) */
        model->markAllCognatesDirty();
        }
}

void UpdateManager::print(void) {

    size_t i = 0;
    for (Update* u : updates)
        std::cout << ++i << " -- " << u << " " << u->parameterType() << std::endl;
}

Update* UpdateManager::randomlyChooseUpdate(void) {

    // Walker's alias method: O(1) selection
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
    size_t chosenIdx = (v < aliasProbability[bin]) ? bin : aliasIndex[bin];
        
    return updates[chosenIdx];
}

void UpdateManager::reject(Update* u) {

    updateInfo.reject(u);

    // restore the primary updated parameter
    Parameter* parm = u->getUpdatedParameter();
    parm->restore();
    
    // restore any additional parameters that were modified (e.g., alignments during topology update)
    for (Parameter* p : u->getAdditionalModifiedParameters())
        p->restore();
    
    if (u->getRateMatrixNeedsUpdate() == true)
        {
        RateMatrix* rm = rateMatrixForUpdate(u);
        if (rm != nullptr)
            rm->restore();
        }
    
    // transition probabilities: restore from backup
    if (u->getAllTiprobsNeedUpdate() == true || u->getSingleBranchChanged() == true)
        {
        if (model->getNumFamilies() > 0)
            {
            if (u->getAllTiprobsNeedUpdate() == true)
                {
                for (int i=0; i<model->getNumFamilies(); i++)
                    model->getSubModel(i)->getTiProbs()->restore();
                }
            else if (u->getSingleBranchChanged() == true)
                {
                Parameter* parm = u->getUpdatedParameter();

                if (ParameterTree* treeParm = dynamic_cast<ParameterTree*>(parm))
                    {
                    for (int i=0; i<model->getNumFamilies(); i++)
                        {
                        if (model->getSubModel(i)->getTree() == treeParm)
                            {
                            model->getSubModel(i)->getTiProbs()->restore();
                            break;
                            }
                        }
                    }
                }
            }
        else
            {
            tiProbs->restore();
            }
        }
        
    // restore the likelihood cache
    model->restoreLikelihoodCache();
}

void UpdateManager::setProposalProbabilities(void) {

    proposalProbabilities.resize(updates.size(), 0.0);

    // check user setting - same pattern as Model::initializeAlignments()
    UserSettings& settings = UserSettings::userSettings();
    bool sampleAlignments = settings.getSampleAlignments();

    if (sampleAlignments)
        {
        for (Update* u : alignmentUpdates)
            proposalProbabilities[updateIndex[u]] = 1.0;
        }
    // else: alignment probabilities remain 0.0

    // sum of alignment weights (zero if sampling is off)
    double alignmentSum = 0.0;
    for (Update* u : alignmentUpdates)
        alignmentSum += proposalProbabilities[updateIndex[u]];

    // other updates share equal weight relative to alignments
    // if alignments are off, give other updates equal weight using their count
    double referenceSum = (sampleAlignments && alignmentSum > 0.0)
                          ? alignmentSum
                          : static_cast<double>(otherUpdates.size());

    double otherProb = (otherUpdates.size() > 0)
                       ? referenceSum / static_cast<double>(otherUpdates.size())
                       : 0.0;
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

    updateInfo.print();
}

void UpdateManager::tune(void) {

    updateInfo.tune();
}

void UpdateManager::updateDependants(Update* u) {

    UserSettings& settings = UserSettings::userSettings();
    bool multiFamily = (model->getNumFamilies() > 0);
    bool linked = (settings.getMultiTreeModel() == linkedFamilies);

    Parameter* parm = u->getUpdatedParameter();

    // Update rate matrix if needed
    if (u->getRateMatrixNeedsUpdate() == true)
        {
        RateMatrix* rm = rateMatrixForUpdate(u);
        if (rm != nullptr)
            rm->updateRateMatrix();
        }

    // Update transition probabilities
    if (u->getAllTiprobsNeedUpdate() == true)
        {
        if (multiFamily == false)
            {
            tiProbs->updateAllBranches();
            }
        else if (linked)
            {
            // Shared substitution parameters affect all families
            for (int i=0; i<model->getNumFamilies(); i++)
                model->getSubModel(i)->getTiProbs()->updateAllBranches();
            }
        else
            {
            // Unlinked model: update only the family whose substitution
            // parameter changed
            bool found = false;

            for (int i=0; i<model->getNumFamilies(); i++)
                {
                SubModel* sm = model->getSubModel(i);

                if (sm->getFrequencies() == parm ||
                    sm->getExchangeabilities() == parm)
                    {
                    sm->getTiProbs()->updateAllBranches();
                    found = true;
                    break;
                    }
                }

            if (found == false)
                {
                // Fallback: if we cannot determine ownership, be safe
                // and update all families.
                for (int i=0; i<model->getNumFamilies(); i++)
                    model->getSubModel(i)->getTiProbs()->updateAllBranches();
                }
            }
        }
    else if (u->getSingleBranchChanged() == true)
        {
        if (multiFamily == false)
            {
            tiProbs->updateBranch();
            }
        else
            {
            if (ParameterTree* treeParm = dynamic_cast<ParameterTree*>(parm))
                {
                bool found = false;

                for (int i=0; i<model->getNumFamilies(); i++)
                    {
                    if (model->getSubModel(i)->getTree() == treeParm)
                        {
                        model->getSubModel(i)->getTiProbs()->updateBranch();
                        found = true;
                        break;
                        }
                    }

                if (found == false)
                    Msg::error("Could not find family for updated tree parameter");
                }
            else
                {
                Msg::error("Single branch changed but updated parameter is not a ParameterTree");
                }
            }
        }

    markCognatesDirty(u);
}

RateMatrix* UpdateManager::rateMatrixForUpdate(Update* u) {

    Parameter* parm = u->getUpdatedParameter();

    if (model->getNumFamilies() == 0)
        return rateMatrix;

    UserSettings& settings = UserSettings::userSettings();

    // Linked model: one shared rate matrix
    if (settings.getMultiTreeModel() == linkedFamilies)
        return model->getRateMatrix();

    // Unlinked model: find which family owns the updated parameter
    for (int i=0; i<model->getNumFamilies(); i++)
        {
        SubModel* sm = model->getSubModel(i);

        if (sm->getFrequencies() == parm ||
            sm->getExchangeabilities() == parm)
            {
            return sm->getRateMatrix();
            }
        }

    return nullptr;
}

void UpdateManager::zeroOut(void) {

    updateInfo.zeroOut();
}
