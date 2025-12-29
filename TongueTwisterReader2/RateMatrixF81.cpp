#include "RateMatrixF81.hpp"
#include "StateFrequencies.hpp"



RateMatrixF81::RateMatrixF81(const StateFrequencies& f) : RateMatrix(f.size()) {
            
    const SampleVector& freqs = f.getValues();
    size_t numSamples = freqs[0].size();
    
    for (size_t n=0; n<numSamples; n++)
        {
        DoubleMatrix* m = new DoubleMatrix(numStates,numStates);
        
        // fill in off-diagonal entries
        for (size_t i=0; i<numStates; i++)
            {
            for (size_t j=i+1; j<numStates; j++)
                {
                (*m)(i,j) = freqs[j][n];
                (*m)(j,i) = freqs[i][n];
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
            (*m)(i,i) = -sum;
            averageRate += freqs[i][n] * sum;
            }
            
        // rescale the rates so the average rate is one
        double scaleFactor = 1.0 / averageRate;
        for (size_t i=0; i<numStates; i++)
            for (size_t j=0; j<numStates; j++)
                (*m)(i,j) *= scaleFactor;
                
        matrices.push_back(m);
        }
        
    calculateMeanAndVariance();
}
