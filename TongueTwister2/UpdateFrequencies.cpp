#include <map>
#include "Model.hpp"
#include "Msg.hpp"
#include "ParameterFrequencies.hpp"
#include "Probability.hpp"
#include "RandomVariable.hpp"
#include "Simplex.hpp"
#include "UpdateFrequencies.hpp"

double UpdateFrequencies::minVal = 0.00001;



UpdateFrequencies::UpdateFrequencies(Model* m, RandomVariable* r, ParameterFrequencies* p) : Update(m, r), myParm(p) {

    numStates = myParm->getNumStates();

    updateNames[0] = "Stationary Frequencies: Dirichlet k=1";
    updateNames[1] = "Stationary Frequencies: Mass Transfer";
    updateNames[2] = "Stationary Frequencies: ALR MVN";
    updateNames[3] = "Stationary Frequencies: Prior";
    tuningValues[0] = 100.0;
    tuningValues[1] = 300.0;
    tuningValues[2] = 0.005;
    tuningValues[3] = 0.0;
    for (size_t i=0; i<4; i++)
        updateHashes[i] = hashUpdateName(updateNames[i]);
}

uint64_t UpdateFrequencies::getUpdateId(void) {

    return updateHashes[lastUpdate];
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

    std::vector<double>& oldFreqs = myParm->getFrequencies(1);
    std::vector<double>& newFreqs = myParm->getFrequencies(); // 0
    
    double lnHastings = 0.0;
    double u = rng->uniformRv();
    if (u < 0.33)
        {
        lastUpdate = 1;
        lnHastings = Simplex::updateMassTransfer(rng, oldFreqs, newFreqs, tuningValues[1], 0.00001);
        }
    else if (u < 0.66)
        {
        lastUpdate = 0;
        lnHastings = Simplex::updateCenteredDirichlet(rng, oldFreqs, newFreqs, tuningValues[0], 1, 0.00001);
        }
    else 
        {
        lastUpdate = 2;
        lnHastings = Simplex::updateALRMVN(rng, oldFreqs, newFreqs, tuningValues[2], 0.00001);
        }
        
    setDependants();
        
    return lnHastings;
}

double UpdateFrequencies::update(double power) {

    if (rng->uniformRv() < priorSampleProb(power))
        {
        std::vector<double>& oldFreqs = myParm->getFrequencies(1);
        std::vector<double>& newFreqs = myParm->getFrequencies(); // 0
        std::vector<double>& alpha = myParm->getAlpha();
        lastUpdate = 3;
        double lnHastings = Simplex::updateFromPrior(rng, oldFreqs, newFreqs, alpha);
        setDependants();
        return lnHastings;
        }
    return update();
}

