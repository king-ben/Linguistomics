#include "Msg.hpp"
#include "RateMatrixAverage.hpp"
#include "StateFrequencies.hpp"



RateMatrixAverage::RateMatrixAverage(size_t ns) : RateMatrix(ns) {

    DoubleMatrix* m = new DoubleMatrix(numStates,numStates);
    
    // fill in off-diagonal entries
    double x = (double)1.0 / numStates;
    for (size_t i=0; i<numStates; i++)
        {
        for (size_t j=i+1; j<numStates; j++)
            {
            (*m)(i,j) = x;
            (*m)(j,i) = x;
            }
        }
        
    // fill in the diagonal elements
    for (size_t i=0; i<numStates; i++)
        (*m)(i,i) = 0.0;
        
    matrices.push_back(m);
    
    mean = new DoubleMatrix(numStates,numStates);
    variance = new DoubleMatrix(numStates,numStates);
    lowerCI = new DoubleMatrix(numStates,numStates);
    upperCI = new DoubleMatrix(numStates,numStates);
    for (size_t i=0; i<numStates; i++)
        {
        for (size_t j=0; j<numStates; j++)
            {
            (*mean)(i,j) = (*m)(i,j);
            (*variance)(i,j) = 0.0;
            (*lowerCI)(i,j) = (*m)(i,j);
            (*upperCI)(i,j) = (*m)(i,j);
            }
        }
}

RateMatrixAverage::RateMatrixAverage(const RateMatrix& rf, const StateFrequencies& freqs) : RateMatrix(freqs.size()) {

    const std::vector<DoubleMatrix*>& q = rf.getRateMatrices();
    const SampleVector& f = freqs.getValues();
    if (rf.getNumSamples() != freqs.numSamples())
        Msg::error("Unequal sample sizes when constructing average rate matrix");
    size_t numSamples = q.size();
    
    for (size_t n=0; n<numSamples; n++)
        {
        DoubleMatrix* m = new DoubleMatrix(numStates,numStates);
        
        // fill in off-diagonal entries
        for (size_t i=0; i<numStates; i++)
            {
            for (size_t j=i+1; j<numStates; j++)
                {
                (*m)(i,j) = f[i][n] * (*q[n])(i,j);
                (*m)(j,i) = f[j][n] * (*q[n])(j,i);
                }
            }
            
        // fill in the diagonal elements
        double averageRate = 0.0;
        for (size_t i=0; i<numStates; i++)
            {
            double sum = 0.0;
            for (size_t j=0; j<numStates; j++)
                {
                if (i != j)
                    sum += (*m)(i,j);
                }
            (*m)(i,i) = 0.0;
            averageRate += sum;
            }
        if (fabs(1.0-averageRate) > 0.001)
            Msg::error("Average rate should be 1.0");
                
        matrices.push_back(m);
        }
        
    calculateMeanAndVariance();
}
