#include "Exchangeabilities.hpp"
#include "Msg.hpp"
#include "RateMatrixHelper.hpp"
#include "RateMatrixNaturalClass.hpp"
#include "StateFrequencies.hpp"



RateMatrixNaturalClass::RateMatrixNaturalClass(const StateFrequencies& f, const Exchangeabilities& r, Partition* part) : RateMatrix(f.size()) {

    RateMatrixHelper helper(numStates, part);
    if (helper.getNumRates() != static_cast<int>(r.getNumRates()))
        Msg::error("Mismatch in number of exchangeability rates for Custom model");
    
    const int* rateMap = helper.getRateIndexMapPtr();

    const SampleVector& freqs = f.getValues();
    const SampleVector& rates = r.getValues();
    size_t numSamples = freqs[0].size();
    const size_t ns = numStates;  // local copy for clarity
    
    for (size_t n = 0; n < numSamples; n++)
        {
        DoubleMatrix* m = new DoubleMatrix(ns, ns);
        double* qData = m->begin();
        
        double averageRate = 0.0;
        
        for (size_t i = 0; i < ns; i++)
            {
            double* qRow_i = qData + i * ns;
            const int* mapRow_i = rateMap + i * ns;
            double rowSum = 0.0;
            const double fi = freqs[i][n];
            
            // upper triangle: j > i
            for (size_t j = i + 1; j < ns; j++)
                {
                int rateIdx = mapRow_i[j];
                double rk = rates[rateIdx][n];
                double qij = rk * freqs[j][n];
                double qji = rk * fi;
                
                qRow_i[j] = qij;
                qData[j * ns + i] = qji;
                
                rowSum += qij;
                }
            
            // lower triangle
            for (size_t j = 0; j < i; j++)
                rowSum += qRow_i[j];
            
            // diagonal
            qRow_i[i] = -rowSum;
            
            // accumulate average rate
            averageRate += fi * rowSum;
            }
        
        // rescale
        if (averageRate > 0.0)
            {
            double scaleFactor = 1.0 / averageRate;
            size_t totalElements = ns * ns;
            for (size_t i = 0; i < totalElements; i++)
                qData[i] *= scaleFactor;
            }
            
        matrices.push_back(m);
        }
    
    calculateMeanAndVariance();
}

