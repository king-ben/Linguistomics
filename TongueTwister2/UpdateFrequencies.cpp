#include <map>
#include "Model.hpp"
#include "Msg.hpp"
#include "ParameterFrequencies.hpp"
#include "Probability.hpp"
#include "RandomVariable.hpp"
#include "UpdateFrequencies.hpp"

double UpdateFrequencies::minVal = 0.001;



UpdateFrequencies::UpdateFrequencies(Model* m, RandomVariable* r, ParameterFrequencies* p) : Update(m, r), myParm(p) {

    numStates = myParm->getNumStates();
}

void UpdateFrequencies::setDependants(void) {

    clearDependencyFlags();
    
    updatedParameter      = myParm;
    allTiprobsNeedUpdate  = true;   // Model changed → all matrices need recomputation
    singleBranchChanged   = false;
    
    // Rate matrix needs update for GTR, but not for F81
    if (model->getRateMatrix() == nullptr)
        rateMatrixNeedsUpdate = false;
    else
        rateMatrixNeedsUpdate = true;
}

double UpdateFrequencies::update(void) {

    return update(1);
}

double UpdateFrequencies::update(int k) {
    
    std::vector<double>& freqs = myParm->getFrequencies();
    double bounce = 1.0;
    double lnP = 0.0;
    if (k == 1)
        {
        // resize vectors for move
        oldValues.resize(2, 0.0);
        newValues.resize(2, 0.0);
        alphaForward.resize(2, 0.0);
        alphaReverse.resize(2, 0.0);
        
        // parameterize dirichlet for forward move
        double alpha0 = 100.0;
        int indexToUpdate = static_cast<int>(rng->uniformRv()*numStates);
        oldValues[0] = freqs[indexToUpdate];
        oldValues[1] = 1.0 - oldValues[0];
        alphaForward[0] = oldValues[0] * alpha0 + bounce;
        alphaForward[1] = oldValues[1] * alpha0 + bounce;
        if (alphaForward[0] < 0.0)
            Msg::error("Negative alpha[0] in k=1 proposal for frequencies");
        if (alphaForward[1] < 0.0)
            Msg::error("Negative alpha[1] in k=1 proposal for frequencies");
        
        // propose new values
        bool err = Probability::Dirichlet::rv(rng, alphaForward, newValues);
        if (err == true)
            Msg::error("Problem proposing dirichlet values in k=1 for frequencies");
        double factor = newValues[1] / oldValues[1];
        for (size_t i=0; i<numStates; i++)
            freqs[i] *= factor;
        freqs[indexToUpdate] = newValues[0];
        Probability::Helper::normalize(newValues, minVal);

        // parameterize dirichlet for reverse move
        alphaReverse[0] = newValues[0] * alpha0 + bounce;
        alphaReverse[1] = newValues[1] * alpha0 + bounce;

        // calculate proposal probability
        lnP = Probability::Dirichlet::lnPdf(alphaReverse, oldValues) - Probability::Dirichlet::lnPdf(alphaForward, newValues);
        lnP += (numStates - 2) * log(factor); // Jacobian
        }
    else if (k > 1 && k < static_cast<int>(numStates))
        {
        double alpha0 = 1000.0;
        
        std::vector<int> indicesToUpdate(numStates);
        for (size_t i=0; i<numStates; i++)
            indicesToUpdate[i] = static_cast<int>(i);
        for (int i=0; i<k; i++)
            {
            int idx = static_cast<int>(rng->uniformRv() * (numStates-i));
            indicesToUpdate[idx] = indicesToUpdate[numStates-i-1];
            indicesToUpdate.pop_back();
            }
        if (indicesToUpdate.size() != (size_t)k)
            Msg::error("Problem selecting indices to update");
        
        std::map<size_t,size_t> mapper;
        for (size_t i=0; i<indicesToUpdate.size(); i++)
            mapper.insert( std::make_pair(indicesToUpdate[i], i) );
            
        oldValues.resize(k+1, 0.0);
        newValues.resize(k+1, 0.0);
        alphaForward.resize(k+1, 0.0);
        alphaReverse.resize(k+1, 0.0);
        
        for (size_t i=0; i<numStates; i++)
            {
            std::map<size_t,size_t>::iterator it = mapper.find(i);
            if (it != mapper.end())
                oldValues[it->second] += freqs[it->first];
            else
                oldValues[k] += freqs[i];
            }
        
        for (int i=0; i<k+1; i++)
            alphaForward[i] = oldValues[i] * alpha0 + bounce;
        
        // draw a new value for the reduced vector
        bool err = false;
        do
            {
            err = Probability::Dirichlet::rv(rng, alphaForward, newValues);
            if (err == true)
                Msg::warning("Trying to recover from Dirichlet random variable error in state frequencies update (k=intermediate)");
            } while (err == true);
        Probability::Helper::normalize(newValues, minVal);
        
        // fill in the Dirichlet parameters for the reverse probability calculations
        for (int i=0; i<k+1; i++)
            alphaReverse[i] = newValues[i] * alpha0 + bounce;
        
        // fill in the full vector
        double factor = newValues[k] / oldValues[k];
        for (size_t i=0; i<numStates; i++)
            {
            std::map<size_t,size_t>::iterator it = mapper.find(i);
            if (it != mapper.end())
                freqs[i] = newValues[it->second];
            else
                freqs[i] *= factor;
            }
        
        // Hastings ratio
        lnP = Probability::Dirichlet::lnPdf(alphaReverse, oldValues) - Probability::Dirichlet::lnPdf(alphaForward, newValues);
        lnP += (numStates - k - 1) * log(factor); // Jacobian
        }
    else
        {
        // update all of the rates
        double alpha0 = 100000.0;
        oldValues.resize(numStates, 0.0);
        newValues.resize(numStates, 0.0);
        alphaForward.resize(numStates, 0.0);
        alphaReverse.resize(numStates, 0.0);

        for (size_t i=0; i<numStates; i++)
            {
            oldValues[i] = freqs[i];
            alphaForward[i] = oldValues[i] * alpha0 + bounce;
            }
        
        bool err = false;
        do
            {
            err = Probability::Dirichlet::rv(rng, alphaForward, newValues);
            if (err == true)
                Msg::warning("Trying to recover from Dirichlet random variable error in state frequencies update (k=N)");
            } while (err == true);
        Probability::Helper::normalize(newValues, minVal);
        
        for (size_t i=0; i<numStates; i++)
            alphaReverse[i] = newValues[i] * alpha0 + bounce;

        lnP = Probability::Dirichlet::lnPdf(alphaReverse, oldValues) - Probability::Dirichlet::lnPdf(alphaForward, newValues);
        
        freqs = newValues;
        }

    setDependants();
    
    return lnP;
}

double UpdateFrequencies::update(double power) {

    if (rng->uniformRv() < priorSampleProb(power))
        return updateFromPrior();
    return update();
}

double UpdateFrequencies::updateFromPrior(void) {

    // draw from the prior distribution, which is a flat Dirichlet distribution
    std::vector<double>& freqs = myParm->getFrequencies();
    std::vector<double>& alpha = myParm->getAlpha();
    Probability::Dirichlet::rv(rng, alpha, freqs);

    setDependants();

    return 0.0;
}
