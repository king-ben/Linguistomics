#include <cmath>
#include <map>
#include "Msg.hpp"
#include "ParameterExchangeabilities.hpp"
#include "Probability.hpp"
#include "RandomVariable.hpp"
#include "Simplex.hpp"
#include "UpdateExchangeabilities.hpp"

double UpdateExchangeabilities::minVal = 0.000001;



UpdateExchangeabilities::UpdateExchangeabilities(Model* m, RandomVariable* r, ParameterExchangeabilities* p) : Update(m, r), myParm(p) {

    numRates = myParm->getNumRates();

    updateNames[0] = "Exchangeabilities: Dirichlet k=1";
    updateNames[1] = "Exchangeabilities: Mass Transfer";
    updateNames[2] = "Exchangeabilities: ALR MVN";
    updateNames[3] = "Exchangeabilities: Prior";
    tuningValues[0] = 200.0;
    tuningValues[1] = 200.0;
    tuningValues[2] = 0.005;
    tuningValues[3] = 0.0;
    for (size_t i=0; i<4; i++)
        updateHashes[i] = hashUpdateName(updateNames[i]);
}

uint64_t UpdateExchangeabilities::getUpdateId(void) {

    return updateHashes[lastUpdate];
}

void UpdateExchangeabilities::setDependants(void) {

    clearDependencyFlags();
    
    updatedParameter      = myParm;
    rateMatrixNeedsUpdate = true;
    allTiprobsNeedUpdate  = true;   // model changed -> all matrices need recomputation
    singleBranchChanged   = false;
}

double UpdateExchangeabilities::update(void) {
    
    std::vector<double>& oldRates = myParm->getRates(1);
    std::vector<double>& newRates = myParm->getRates(); // 0
    
    double lnHastings = 0.0;
    double u = rng->uniformRv();
    if (u < 0.33)
        {
        lastUpdate = 1;
        lnHastings = Simplex::updateMassTransfer(rng, oldRates, newRates, tuningValues[1], 0.00001);
        }
    else if (u < 0.66)
        {
        lastUpdate = 0;
        lnHastings = Simplex::updateCenteredDirichlet(rng, oldRates, newRates, tuningValues[0], 1, 0.00001);
        }
    else 
        {
        lastUpdate = 2;
        lnHastings = Simplex::updateALRMVN(rng, oldRates, newRates, tuningValues[2], 0.00001);
        }
        
    setDependants();
        
    return lnHastings;
}

double UpdateExchangeabilities::update(double power) {

    if (rng->uniformRv() < priorSampleProb(power))
        {
        std::vector<double>& oldRates = myParm->getRates(1);
        std::vector<double>& newRates = myParm->getRates(); // 0
        std::vector<double>& alpha = myParm->getAlpha();
        lastUpdate = 3;
        double lnHastings = Simplex::updateFromPrior(rng, oldRates, newRates, alpha);
        setDependants();
        return lnHastings;
        }
    return update();
}
