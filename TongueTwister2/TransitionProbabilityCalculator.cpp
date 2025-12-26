#include "Ctmc.hpp"
#include "TransitionProbabilityCalculator.hpp"



TransitionProbabilityCalculator::TransitionProbabilityCalculator(Ctmc* m) : model(m) {

    cache = model->allocateMathCache();
}

TransitionProbabilityCalculator::~TransitionProbabilityCalculator(void) {

    if (cache != nullptr)
        delete cache;
}

void TransitionProbabilityCalculator::computeTransitionProbabilities(void) {

    model->computeTransitionProbs(branchLength, output, cache);
}

void TransitionProbabilityCalculator::set(double bl, DoubleMatrix* out) {

    branchLength = bl;
    output = out;
}
