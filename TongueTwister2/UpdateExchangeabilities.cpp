#include <cmath>
#include <map>
#include "Msg.hpp"
#include "ParameterExchangeabilities.hpp"
#include "Probability.hpp"
#include "RandomVariable.hpp"
#include "UpdateExchangeabilities.hpp"

double UpdateExchangeabilities::minVal = 0.000001;



UpdateExchangeabilities::UpdateExchangeabilities(Model* m, RandomVariable* r, ParameterExchangeabilities* p) : Update(m, r), myParm(p) {

    numRates = myParm->getNumRates();
}

void UpdateExchangeabilities::notifyDependants(void) {

}

void UpdateExchangeabilities::setDependants(void) {

    clearDependencyFlags();
    
    updatedParameter      = myParm;
    rateMatrixNeedsUpdate = true;
    allTiprobsNeedUpdate  = true;   // Model changed → all matrices need recomputation
    singleBranchChanged   = false;
}

double UpdateExchangeabilities::update(void) {

    numTries++;
    
    double bounce = 1.0;
    std::vector<double>& rates = myParm->getRates();
    
    // choose the number of rates to update
    int k = 1;
    double u = rng->uniformRv();
    if (u <= 0.5)
        k = 1;
    else if (u > 0.5 && u <= 0.9)
        k = 10;
    else
        k = static_cast<int>(numRates);
        
    k = 1;
            
    // update the rates
    double lnP = 0.0;
    if (k == 1)
        {
        // resize vectors for move
        oldValues.resize(2, 0.0);
        newValues.resize(2, 0.0);
        alphaForward.resize(2, 0.0);
        alphaReverse.resize(2, 0.0);

        // parameterize dirichlet for forward move
        double alpha0 = 200.0;
        int indexToUpdate = static_cast<int>(rng->uniformRv()*numRates);
        oldValues[0] = rates[indexToUpdate];
        oldValues[1] = 1.0 - oldValues[0];
        alphaForward[0] = oldValues[0] * alpha0 + bounce;
        alphaForward[1] = oldValues[1] * alpha0 + bounce;
        if (alphaForward[0] < 0.0)
            Msg::error("Negative alpha[0] in k=1 proposal for exchangeability rates");
        if (alphaForward[1] < 0.0)
            Msg::error("Negative alpha[1] in k=1 proposal for exchangeability rates");
                
        // propose new values
        bool err = Probability::Dirichlet::rv(rng, alphaForward, newValues);
        if (err == true)
            Msg::error("Problem proposing dirichlet values in k=1 for exchangeability rates");
        double factor = newValues[1] / oldValues[1];
        for (size_t i=0; i<numRates; i++)
            rates[i] *= factor;
        rates[indexToUpdate] = newValues[0];
        Probability::Helper::normalize(rates, minVal);

        // parameterize dirichlet for reverse move
        alphaReverse[0] = newValues[0] * alpha0 + bounce;
        alphaReverse[1] = newValues[1] * alpha0 + bounce;

        // calculate proposal probability
        lnP = Probability::Dirichlet::lnPdf(alphaReverse, oldValues) - Probability::Dirichlet::lnPdf(alphaForward, newValues);
        lnP += (numRates - 2) * log(factor); // Jacobian
        }
    else if (k > 1 && k < static_cast<int>(numRates))
        {
        // resize vectors for move
        oldValues.resize(k+1, 0.0);
        newValues.resize(k+1, 0.0);
        alphaForward.resize(k+1, 0.0);
        alphaReverse.resize(k+1, 0.0);

        // parameterize dirichlet for forward move
        double alpha0 = 1000.0;
        std::vector<int> indicesToUpdate(numRates);
        for (size_t i=0; i<numRates; i++)
            indicesToUpdate[i] = static_cast<int>(i);
        for (int i=0; i<k; i++)
            {
            int idx = static_cast<int>(rng->uniformRv() * (numRates-i));
            indicesToUpdate[idx] = indicesToUpdate[numRates-i-1];
            indicesToUpdate.pop_back();
            }
        if (indicesToUpdate.size() != (size_t)k)
            Msg::error("Problem selecting indices to update");

        std::map<size_t,size_t> mapper;
        for (size_t i=0; i<indicesToUpdate.size(); i++)
            mapper.insert( std::make_pair(indicesToUpdate[i], i) );
        for (size_t i=0; i<numRates; i++)
            {
            std::map<size_t,size_t>::iterator it = mapper.find(i);
            if (it != mapper.end())
                oldValues[it->second] += rates[it->first];
            else
                oldValues[k] += rates[i];
            }
        for (int i=0; i<k+1; i++)
            {
            alphaForward[i] = oldValues[i] * alpha0 + bounce;
            if (alphaForward[i] < 0.0)
                Msg::error("Negative alpha[" + std::to_string(i)  + ")] in k=" + std::to_string(k) + " proposal for exchangeability rates");
            }
        
        // propose new values
        bool err = Probability::Dirichlet::rv(rng, alphaForward, newValues);
        if (err == true)
            Msg::error("Problem proposing dirichlet values in k=" + std::to_string(k) + " for exchangeability rates");
        double factor = newValues[k] / oldValues[k];
        for (size_t i=0; i<numRates; i++)
            {
            std::map<size_t,size_t>::iterator it = mapper.find(i);
            if (it != mapper.end())
                rates[i] = newValues[it->second];
            else
                rates[i] *= factor;
            }
        Probability::Helper::normalize(rates, minVal);

        // parameterize dirichlet for reverse move
        for (int i=0; i<k+1; i++)
            alphaReverse[i] = newValues[i] * alpha0 + bounce;

        // calculate proposal probability
        lnP = Probability::Dirichlet::lnPdf(alphaReverse, oldValues) - Probability::Dirichlet::lnPdf(alphaForward, newValues);
        lnP += (numRates - k - 1) * log(factor); // Jacobian
        }
    else
        {
        // resize vectors for move
        oldValues.resize(numRates);
        newValues.resize(numRates);
        alphaForward.resize(numRates);
        alphaReverse.resize(numRates);

        // parameterize dirichlet for forward move
        double alpha0 = 4000.0;
        for (size_t i=0; i<numRates; i++)
            {
            oldValues[i] = rates[i];
            alphaForward[i] = oldValues[i] * alpha0 + bounce;
            if (alphaForward[i] < 0.0)
                Msg::error("Negative alpha[" + std::to_string(i)  + ")] in k=" + std::to_string(numRates) + " proposal for exchangeability rates");
            }
        
        // propose new values
        bool err = Probability::Dirichlet::rv(rng, alphaForward, newValues);
        if (err == true)
            Msg::error("Problem proposing dirichlet values in k=" + std::to_string(numRates) + " for exchangeability rates");
        Probability::Helper::normalize(newValues, minVal);
        
        // parameterize dirichlet for reverse move
        for (size_t i=0; i<numRates; i++)
            {
            rates[i] = newValues[i];
            alphaReverse[i] = newValues[i] * alpha0 + bounce;
            }

        // calculate proposal probability
        lnP = Probability::Dirichlet::lnPdf(alphaReverse, oldValues) - Probability::Dirichlet::lnPdf(alphaForward, newValues);
        }
        
    // quick check that the exchangeability rates are still valid
    double sum = 0.0;
    bool isBad = false;
    for (size_t i=0; i<rates.size(); i++)
        {
        sum += rates[i];
        if (rates[i] < 0.0)
            isBad = true;
        }
    if (isBad == true || fabs(1.0-sum) > 0.000001)
        lnP = -1e100;
        
    setDependants();

    return lnP;
}

double UpdateExchangeabilities::update(double power) {

    if (rng->uniformRv() < priorSampleProb(power))
        return updateFromPrior();
    return update();
}

double UpdateExchangeabilities::updateFromPrior(void) {

    // draw from the prior distribution, which is a flat Dirichlet distribution
    std::vector<double>& rates = myParm->getRates();
    std::vector<double>& alpha = myParm->getAlpha();
    Probability::Dirichlet::rv(rng, alpha, rates);

    setDependants();

    return 0.0;
}
