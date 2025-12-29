#include "Exchangeabilities.hpp"
#include "Msg.hpp"
#include "RateMatrixGTR.hpp"
#include "StateFrequencies.hpp"



RateMatrixGTR::RateMatrixGTR(const StateFrequencies& f, const Exchangeabilities& r) : RateMatrix(f.size()) {
            
    const SampleVector& freqs = f.getValues();
    const SampleVector& rates = r.getValues();
    size_t numSamples = freqs[0].size();
    
    for (size_t n=0; n<numSamples; n++)
        {
        DoubleMatrix* m = new DoubleMatrix(numStates,numStates);
        
        // fill in off-diagonal entries
        for (size_t i=0, k=0; i<numStates; i++)
            {
            for (size_t j=i+1; j<numStates; j++)
                {
                (*m)(i,j) = rates[k][n] * freqs[j][n];
                (*m)(j,i) = rates[k][n] * freqs[i][n];
                k++;
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
