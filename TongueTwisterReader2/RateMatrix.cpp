#include "Container.hpp"
#include "Msg.hpp"
#include "RateMatrix.hpp"
#include "Statistics.hpp"



RateMatrix::RateMatrix(size_t ns) : numStates(ns), mean(nullptr), variance(nullptr) {

}

RateMatrix::~RateMatrix(void) {

    for (size_t i=0; i<matrices.size(); i++)
        delete matrices[i];
    if (mean != nullptr)
        delete mean;
    if (variance != nullptr)
        delete variance;
    if (lowerCI != nullptr)
        delete lowerCI;
    if (upperCI != nullptr)
        delete upperCI;
}

void RateMatrix::calculateMeanAndVariance(void) {

    if (mean != nullptr)
        delete mean;
    if (variance != nullptr)
        delete variance;
    if (lowerCI != nullptr)
        delete lowerCI;
    if (upperCI != nullptr)
        delete upperCI;
    
    mean = new DoubleMatrix(numStates,numStates);
    variance = new DoubleMatrix(numStates,numStates);
    lowerCI = new DoubleMatrix(numStates,numStates);
    upperCI = new DoubleMatrix(numStates,numStates);
    
    size_t numSamples = matrices.size();
    std::vector<float> x(numSamples);
    for (size_t i = 0; i < numStates; i++)
        {
        for (size_t j = 0; j < numStates; j++)
            {
            for (size_t n = 0; n < matrices.size(); n++)
                {
                DoubleMatrix& m = *matrices[n];
                x[n] = m(i,j);
                }
            std::pair<float,float> mv = Statistics::getMeanAndVariance(x);
            CredibleInterval ci = Statistics::getCredibleInterval(x);
            (*mean)(i,j) = mv.first;
            (*variance)(i,j) = mv.second;
            (*lowerCI)(i,j) = ci.lower;
            (*upperCI)(i,j) = ci.upper;
            }        
        }
}
