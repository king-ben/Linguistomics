#include <algorithm>
#include "RandomVariable.hpp"
#include "Update.hpp"



Update::Update(Model* m, RandomVariable* r) : 
    model(m), rng(r), updatedParameter(nullptr), 
    rateMatrixNeedsUpdate(false), allTiprobsNeedUpdate(false), 
    singleBranchChanged(false), changedBranchLength(0.0) {
    
}

void Update::clearDependencyFlags(void) {

    rateMatrixNeedsUpdate = false;
    allTiprobsNeedUpdate = false;
    singleBranchChanged = false;
    changedBranchLength = 0.0;
}

double Update::priorSampleProb(double power) {

    // Hill, A. V. 1910. The possible effects of the aggregation of the molecules of 
    //   hæmoglobin on its dissociation curves. Journal of Physiology, 40, iv–vii.

    double h = 4.0;      // cliffiness (large h -> steeper cliff)
    double beta50 = 0.02; // cliff position (at what power is the probability 50/50)
    power = std::clamp(power, 0.0, 1.0); // keeps power in bounds

    // beta50 must be > 0, h must be > 0 for sane behavior
    if (beta50 <= 0.0) 
        return 0.0;
    if (h <= 0.0)      
        return 0.0;

    if (power == 0.0) 
        return 1.0;

    const double x = power / beta50;
    const double xh = std::pow(x, h);
    return 1.0 / (1.0 + xh);
}
