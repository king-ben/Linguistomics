#include "RateMatrixJC69.hpp"



RateMatrixJC69::RateMatrixJC69(size_t ns) : RateMatrix(ns) {

    DoubleMatrix* m = new DoubleMatrix(numStates,numStates);
    
    // fill in off-diagonal entries
    for (size_t i=0; i<numStates; i++)
        {
        for (size_t j=i+1; j<numStates; j++)
            {
            (*m)(i,j) = 1.0;
            (*m)(j,i) = 1.0;
            }
        }
        
    // fill in the diagonal elements
    double freq = (double)1.0 / numStates;
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
        averageRate += freq * sum;
        }
        
    // rescale the rates so the average rate is one
    double scaleFactor = 1.0 / averageRate;
    for (size_t i=0; i<numStates; i++)
        for (size_t j=0; j<numStates; j++)
            (*m)(i,j) *= scaleFactor;
            
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
