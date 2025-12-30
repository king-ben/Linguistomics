#include "Ctmc.hpp"
#include "TransitionProbabilityCalculator.hpp"



TransitionProbabilityCalculator::TransitionProbabilityCalculator(Ctmc* m) : ctmc(m) {

    cache = ctmc->allocateMathCache();
}

TransitionProbabilityCalculator::~TransitionProbabilityCalculator(void) {

    if (cache != nullptr)
        delete cache;
}

void TransitionProbabilityCalculator::computeTransitionProbabilities(void) {

    ctmc->computeTransitionProbs(branchLength, output, cache);
}

void TransitionProbabilityCalculator::set(double bl, DoubleMatrix* out) {

    branchLength = bl;
    output = out;
}
