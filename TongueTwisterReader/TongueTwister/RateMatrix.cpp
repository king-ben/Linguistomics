#include <cmath>
#include <iomanip>
#include <iostream>
#include "MatrixMath.hpp"
#include "RateMatrix.hpp"
#include "RateMatrixHelper.hpp"



RateMatrix::RateMatrix(void) {
    Q[0] = NULL;
    Q[1] = NULL;
    numStates = 0;
    activeMatrix = 0;
    isInitialized = false;
    rateMatrixHelper = NULL;
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

void RateMatrix::initialize(int d, RateMatrixHelper* h) {

    if (isInitialized == true)
        return;
        
    rateMatrixHelper = h;
    numStates = d;
    activeMatrix = 0;

    Q[0] = new DoubleMatrix(numStates, numStates);
    Q[1] = new DoubleMatrix(numStates, numStates);
        
    for (int s=0; s<2; s++)
        equilibriumFrequencies[s].resize(numStates);
        
    isInitialized = true;
}

void  RateMatrix::updateRateMatrix(std::vector<double>& rates, std::vector<double>& f) {

    // set the stationary frequencies
    for (int i=0; i<numStates; i++)
        equilibriumFrequencies[activeMatrix][i] = f[i];
    
    // initialize the rate matrix
    DoubleMatrix& QM = *Q[activeMatrix];
    
    // fill in off diagonal components of the rate matrix in
    // a model-dependent manner
    if ((int)rates.size() == numStates * (numStates-1) / 2)
        {
        // gtr model
        for (int i=0, k=0; i<numStates; i++)
            {
            for (int j=i+1; j<numStates; j++)
                {
                QM(i,j) = rates[k] * f[j];
                QM(j,i) = rates[k] * f[i];
                k++;
                }
            }
        }
    else
        {
        // custom model
        int** map = rateMatrixHelper->getMap();
        for (int i=0; i<numStates; i++)
            {
            for (int j=i+1; j<numStates; j++)
                {
                int changeType = map[i][j];
                QM(i,j) = rates[changeType] * f[j];
                QM(j,i) = rates[changeType] * f[i];
                }
            }
        }
                
    // fill in the diagonal elements of the rate matrix
    for (int i=0; i<numStates; i++)
        {
        double sum = 0.0;
        for (int j=0; j<numStates; j++)
            {
            if (i != j)
                sum += QM(i,j);
            }
        QM(i,i) = -sum;
        }
    
    // rescale the rate matrix
    double averageRate = 0.0;
    for (int i=0; i<numStates; i++)
        averageRate += -f[i] * QM(i,i);
    double scaleFactor = 1.0 / averageRate;
    QM.multiply(scaleFactor); // QM *= scaleFactor
            
#   if 0
    std::cout << std::fixed << std::setprecision(5);
    for (int i=0; i<numStates; i++)
        {
        for (int j=0; j<numStates; j++)
            std::cout << QM[i][j] << " ";
        std::cout << std::endl;
        }
    for (int i=0; i<numStates; i++)
        std::cout << equilibriumFrequencies[activeMatrix][i] << " ";
    std::cout << std::endl;
#   endif
}
