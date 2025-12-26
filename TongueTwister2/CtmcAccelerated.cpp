#include "CtmcAccelerated.hpp"
#include "MathCacheAccelerated.hpp"
#include "RateMatrix.hpp"
#include "Msg.hpp"
#include <cmath>



CtmcAccelerated::CtmcAccelerated(size_t ns, RateMatrix* r) : Ctmc(ns), rateMatrix(r) {

}

CtmcAccelerated::~CtmcAccelerated(void) {

}

MathCache* CtmcAccelerated::allocateMathCache(void) {

    return new MathCacheAccelerated;
}

void CtmcAccelerated::computeTransitionProbs(double branchLength, DoubleMatrix* output, MathCache* cache) const {

    if (branchLength < 0.0)
        Msg::error("Cannot compute transition probabilities for negative branch length: " + std::to_string(branchLength));
    
    MathCacheAccelerated* accCache = static_cast<MathCacheAccelerated*>(cache);
    DoubleMatrix& Q = rateMatrix->getQ();
    
    accCache->computeMatrixExponential(Q, branchLength, *output);
    
#   if defined(CHECK_TIPROBS)
    if (checkTransitionProbabilities(output, 0.00001) == false)
        Msg::error("Problem calculating transition probabilities");
#   endif
}

const double* CtmcAccelerated::getQData(void) const {

    return rateMatrix->getQ().begin();
}
