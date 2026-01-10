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

    updateInfo.resize(4);
    updateInfo[0].updateIdx = 0;
    updateInfo[0].updateName = "Exchangeabilities: Dirichlet k=1";
    updateInfo[0].updateHash = hashUpdateName("Exchangeabilities: Dirichlet k=1");
    updateInfo[0].updateType = simplex;
    updateInfo[0].tuningParameter = 200.0;

    updateInfo[1].updateIdx = 1;
    updateInfo[1].updateName = "Exchangeabilities: Mass Transfer";
    updateInfo[1].updateHash = hashUpdateName("Exchangeabilities: Mass Transfer");
    updateInfo[1].updateType = undef;
    updateInfo[1].tuningParameter = 200.0;

    updateInfo[2].updateIdx = 2;
    updateInfo[2].updateName = "Exchangeabilities: ALR MVN";
    updateInfo[2].updateHash = hashUpdateName("Exchangeabilities: ALR MVN");
    updateInfo[2].updateType = factor;
    updateInfo[2].tuningParameter = 0.005;

    updateInfo[3].updateIdx = 3;
    updateInfo[3].updateName = "Exchangeabilities: Prior";
    updateInfo[3].updateHash = hashUpdateName("Exchangeabilities: Prior");
    updateInfo[3].updateType = undef;
    updateInfo[3].tuningParameter = 0.0;
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
        lnHastings = Simplex::updateMassTransfer(rng, oldRates, newRates, updateInfo[1].tuningParameter, 0.00001);
        }
    else if (u < 0.66)
        {
        lastUpdate = 0;
        lnHastings = Simplex::updateCenteredDirichlet(rng, oldRates, newRates, updateInfo[0].tuningParameter, 1, 0.00001);
        }
    else 
        {
        lastUpdate = 2;
        lnHastings = Simplex::updateALRMVN(rng, oldRates, newRates, updateInfo[2].tuningParameter, 0.00001);
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
