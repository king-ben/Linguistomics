#include <cmath>
#include <iomanip>
#include <iostream>
#include "Msg.hpp"
#include "ParameterExchangeabilities.hpp"
#include "ParameterFrequencies.hpp"
#include "RateMatrix.hpp"
#include "RateMatrixHelper.hpp"



RateMatrix::RateMatrix(size_t ns) : numStates(ns), activeMatrix(0) {

    Q[0] = new DoubleMatrix(numStates, numStates);
    Q[1] = new DoubleMatrix(numStates, numStates);
}

RateMatrix::~RateMatrix(void) {

    delete Q[0];
    delete Q[1];
}

void RateMatrix::flipActiveValues(void) {

    activeMatrix = 1 - activeMatrix;
}

void RateMatrix::keep(void) {

    // no need to flip the active rate matrix
}

void RateMatrix::print(void) {

    std::cout << std::fixed << std::setprecision(5);
    DoubleMatrix& QM = *Q[activeMatrix];
    for (size_t i=0; i<numStates; i++)
        {
        for (size_t j=0; j<numStates; j++)
            std::cout << QM(i, j) << " ";
        std::cout << std::endl;
        }
}

void RateMatrix::restore(void) {

    flipActiveValues();
}

RateMatrixGTR::RateMatrixGTR(size_t ns, ParameterExchangeabilities* e, ParameterFrequencies* f) : 
    RateMatrix(ns), exchParm(e), freqParm(f) {

}

void RateMatrixGTR::updateRateMatrix(void) {

    // flip to work on inactive matrix
    flipActiveValues();

    // get raw pointers for inner loop performance
    DoubleMatrix& QM = *Q[activeMatrix];
    double* qData = QM.begin();
    const double* r = exchParm->getRates().data();
    const double* f = freqParm->getFrequencies().data();
    const size_t n = numStates;
    
    // fill off-diagonal and accumulate row sums simultaneously
    double averageRate = 0.0;
    for (size_t i=0; i<n; i++)
        {
        double* qRow_i = qData + i * n;
        double rowSum = 0.0;
        const double fi = f[i];
        
        // upper triangle: j > i
        size_t k = (i * (2 * n - i - 1)) / 2;  // index into symmetric rate array
        for (size_t j=i+1; j<n; j++, k++)
            {
            double rk = r[k];
            double qij = rk * f[j];
            double qji = rk * fi;
            
            qRow_i[j] = qij;
            qData[j * n + i] = qji;  // symmetric element
            
            rowSum += qij;
            }
        
        // lower triangle was filled when we processed earlier rows
        // Just sum those values
        for (size_t j=0; j<i; j++)
            rowSum += qRow_i[j];
        
        // set diagonal
        qRow_i[i] = -rowSum;
        
        // accumulate average rate
        averageRate += fi * rowSum;
        }
    
    // rescale so average rate = 1.0
    if (averageRate > 0.0)
        {
        double scaleFactor = 1.0 / averageRate;
        size_t totalElements = n * n;
        for (size_t i = 0; i < totalElements; i++)
            qData[i] *= scaleFactor;
        }
}

RateMatrixCustom::RateMatrixCustom(size_t ns, ParameterExchangeabilities* e, ParameterFrequencies* f, Partition* p) : 
    RateMatrix(ns), exchParm(e), freqParm(f), rateIndexMap(nullptr) {

    rateMatrixHelper = new RateMatrixHelper(ns, p);
    
    if (rateMatrixHelper->getNumRates() != static_cast<int>(e->getNumRates()))
        Msg::error("Mismatch in number of exchangeability rates for Custom model");
    
    // cache pointer to rate index map for inner loop
    rateIndexMap = rateMatrixHelper->getRateIndexMapPtr();
}

RateMatrixCustom::~RateMatrixCustom(void) {

    delete rateMatrixHelper;
}

void RateMatrixCustom::updateRateMatrix(void) {

    // flip to work on inactive matrix
    flipActiveValues();

    // get raw pointers for maximum inner loop performance
    DoubleMatrix& QM = *Q[activeMatrix];
    double* qData = QM.begin();
    const double* r = exchParm->getRates().data();
    const double* f = freqParm->getFrequencies().data();
    const int* rateMap = rateIndexMap;
    const size_t n = numStates;
    
    // fill off-diagonal and accumulate row sums simultaneously
    double averageRate = 0.0;
    
    for (size_t i=0; i<n; i++)
        {
        double* qRow_i = qData + i * n;
        const int* mapRow_i = rateMap + i * n;
        double rowSum = 0.0;
        const double fi = f[i];
        
        // upper diagonal
        for (size_t j=i+1; j<n; j++)
            {
            int rateIdx = mapRow_i[j];
            double rk = r[rateIdx];
            double qij = rk * f[j];
            double qji = rk * fi;
            
            qRow_i[j] = qij;
            qData[j * n + i] = qji;  // symmetric element
            
            rowSum += qij;
            }
        
        // lower triangle: sum values filled from earlier rows
        for (size_t j=0; j<i; j++)
            rowSum += qRow_i[j];
        
        // set diagonal
        qRow_i[i] = -rowSum;
        
        // accumulate average rate
        averageRate += fi * rowSum;
        }
    
    // rescale so average rate = 1.0
    if (averageRate > 0.0)
        {
        double scaleFactor = 1.0 / averageRate;
        size_t totalElements = n * n;
        for (size_t i = 0; i < totalElements; i++)
            qData[i] *= scaleFactor;
        }
}
