#include <cmath>
#include <iomanip>
#include <iostream>
#include "MathCache.hpp"
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

    if (activeMatrix == 0)
        activeMatrix = 1;
    else
        activeMatrix = 0;
}

void RateMatrix::keep(void) {

    // no need to flip the active rate matrix
}

void RateMatrix::print(void) {

    std::cout << std::fixed << std::setprecision(5);
    for (size_t i=0; i<numStates; i++)
        {
        for (size_t j=0; j<numStates; j++)
            std::cout << (*Q[activeMatrix])(i,j) << " ";
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

    // flip the active rate matrix 
    flipActiveValues();

    // fill in off-diagonal elements of the rate matrix
    DoubleMatrix& QM = *Q[activeMatrix];
    std::vector<double>& r = exchParm->getRates();
    std::vector<double>& f = freqParm->getFrequencies();
    for (size_t i=0, k=0; i<numStates; i++)
        {
        for (size_t j=i+1; j<numStates; j++)
            {
            QM(i,j) = r[k] * f[j];
            QM(j,i) = r[k] * f[i];
            k++;
            }
        }

    // fill in the diagonal elements of the rate matrix
    for (size_t i=0; i<numStates; i++)
        {
        double sum = 0.0;
        for (size_t j=0; j<numStates; j++)
            {
            if (i != j)
                sum += QM(i,j);
            }
        QM(i,i) = -sum;
        }
    
    // rescale the rate matrix
    double averageRate = 0.0;
    for (size_t i=0; i<numStates; i++)
        averageRate += -f[i] * QM(i,i);
    double scaleFactor = 1.0 / averageRate;
    QM.multiply(scaleFactor); // QM *= scaleFactor
}

RateMatrixCustom::RateMatrixCustom(size_t ns, ParameterExchangeabilities* e, ParameterFrequencies* f, Partition* p) : 
    RateMatrix(ns), exchParm(e), freqParm(f) {

    rateMatrixHelper = new RateMatrixHelper(ns, p);
    if (rateMatrixHelper->getLabels().size() != e->getNumRates())
        Msg::error("Mismatch in number of exchangeability rates for Custom model");
}

RateMatrixCustom::~RateMatrixCustom(void) {

    delete rateMatrixHelper;
}

void RateMatrixCustom::updateRateMatrix(void) {

    // flip the active rate matrix 
    flipActiveValues();

    // fill in off-diagonal elements of the rate matrix
    DoubleMatrix& QM = *Q[activeMatrix];
    std::vector<double>& r = exchParm->getRates();
    std::vector<double>& f = freqParm->getFrequencies();
    IntMatrix& map = rateMatrixHelper->getMap();
    for (size_t i=0; i<numStates; i++)
        {
        for (size_t j=i+1; j<numStates; j++)
            {
            int changeType = map(i,j);
            QM(i,j) = r[changeType] * f[j];
            QM(j,i) = r[changeType] * f[i];
            }
        }

    // fill in the diagonal elements of the rate matrix
    for (size_t i=0; i<numStates; i++)
        {
        double sum = 0.0;
        for (size_t j=0; j<numStates; j++)
            {
            if (i != j)
                sum += QM(i,j);
            }
        QM(i,i) = -sum;
        }
    
    // rescale the rate matrix
    double averageRate = 0.0;
    for (size_t i=0; i<numStates; i++)
        averageRate += -f[i] * QM(i,i);
    double scaleFactor = 1.0 / averageRate;
    QM.multiply(scaleFactor); // QM *= scaleFactor
}
