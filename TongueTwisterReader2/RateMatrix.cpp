#include "Container.hpp"
#include "Msg.hpp"
#include "RateMatrix.hpp"



RateMatrix::RateMatrix(size_t ns) : numStates(ns), mean(nullptr), variance(nullptr) {

}

RateMatrix::~RateMatrix(void) {

    for (size_t i=0; i<matrices.size(); i++)
        delete matrices[i];
    if (mean != nullptr)
        delete mean;
    if (variance != nullptr)
        delete variance;
}

void RateMatrix::calculateMeanAndVariance(void) {

    if (mean != nullptr)
        delete mean;
    if (variance != nullptr)
        delete variance;
    
    mean = new DoubleMatrix(numStates,numStates);
    variance = new DoubleMatrix(numStates,numStates);
    
    // Using Welford's online algorithm for numerical stability
    DoubleMatrix runningMean(numStates, numStates);
    DoubleMatrix M2(numStates, numStates);
    
    for (size_t n = 0; n < matrices.size(); n++)
        {
        DoubleMatrix& m = *matrices[n];
        for (size_t i = 0; i < numStates; i++)
            {
            for (size_t j = 0; j < numStates; j++)
                {
                double x = m(i, j);
                double delta = x - runningMean(i, j);
                runningMean(i, j) += delta / (n + 1);
                double delta2 = x - runningMean(i, j);
                M2(i, j) += delta * delta2;
                }
            }
        }

    size_t nSamples = matrices.size();
    for (size_t i = 0; i < numStates; i++)
        {
        for (size_t j = 0; j < numStates; j++)
            {
            (*mean)(i, j) = runningMean(i, j);
            if (nSamples > 1)
                (*variance)(i, j) = M2(i, j) / (nSamples - 1);
            else 
                (*variance)(i, j) = 0.0;
            }
        }
}
