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

    updateInfo.resize(4);
    updateInfo[0].updateIdx = 0;
    updateInfo[0].updateName = "Stationary Frequencies: Dirichlet k=1";
    updateInfo[0].updateHash = hashUpdateName("Stationary Frequencies: Dirichlet k=1");
    updateInfo[0].updateType = simplex;
    updateInfo[0].tuningParameter = 100.0;

    updateInfo[1].updateIdx = 1;
    updateInfo[1].updateName = "Stationary Frequencies: Mass Transfer";
    updateInfo[1].updateHash = hashUpdateName("Stationary Frequencies: Mass Transfer");
    updateInfo[1].updateType = undef;
    updateInfo[1].tuningParameter = 300.0;

    updateInfo[2].updateIdx = 2;
    updateInfo[2].updateName = "Stationary Frequencies: ALR MVN";
    updateInfo[2].updateHash = hashUpdateName("Stationary Frequencies: ALR MVN");
    updateInfo[2].updateType = factor;
    updateInfo[2].tuningParameter = 0.005;

    updateInfo[3].updateIdx = 3;
    updateInfo[3].updateName = "Stationary Frequencies: Prior";
    updateInfo[3].updateHash = hashUpdateName("Stationary Frequencies: Prior");
    updateInfo[3].updateType = undef;
    updateInfo[3].tuningParameter = 0.0;
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
        lnHastings = Simplex::updateMassTransfer(rng, oldFreqs, newFreqs, updateInfo[1].tuningParameter, 0.00001);
        }
    else if (u < 0.66)
        {
        lastUpdate = 0;
        lnHastings = Simplex::updateCenteredDirichlet(rng, oldFreqs, newFreqs, updateInfo[0].tuningParameter, 1, 0.00001);
        }
    else 
        {
        lastUpdate = 2;
        lnHastings = Simplex::updateALRMVN(rng, oldFreqs, newFreqs, updateInfo[2].tuningParameter, 0.00001);
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

