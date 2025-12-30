#include "Analysis.hpp"
#include "AnalysisComparison.hpp"
#include "Container.hpp"



void AnalysisComparison::compare(Analysis* a1, Analysis* a2, size_t numStates, size_t nSamples) {

    std::vector<float> f1;
    std::vector<float> f2;
    double sumW = 0.0, sumF = 0.0;
    
    for (size_t n = 0; n < nSamples; n++)
        {
        DoubleMatrix* m1 = a1->randomlyChooseRateMatrixAndFreqs(f1);
        DoubleMatrix* m2 = a2->randomlyChooseRateMatrixAndFreqs(f2);
        numStates = m1->getNumRows();
        
        double ss = 0.0;
        for (size_t i=0; i<numStates; i++)
            {
            for (size_t j=i+1; j<numStates; j++)
                {
                double x = f1[i] * (*m1)(i,j) - f2[i] * (*m2)(i,j);
                double y = f1[j] * (*m1)(j,i) - f2[j] * (*m2)(j,i);
                ss += (x*x + y*y);
                }
            }
        sumW += ss;
        
        ss = 0.0;
        for (size_t i=0; i<numStates; i++)
            ss += (f1[i] - f2[i]) * (f1[i] - f2[i]);
        sumF += ss;
        }
    double dW = sumW / nSamples;
    double dF = sumF / nSamples;
        
    
    // compare rates
    
    std::cout << "Model Comparison:   " << a1->modelName() << "-" << a2->modelName() << std::endl;
    std::cout << "Average W distance: " << dW << std::endl;
    std::cout << "Average F distance: " << dF << std::endl;
    std::cout << "Average W distance: " << sqrt(dW) << std::endl;
    std::cout << "Average F distance: " << sqrt(dF) << std::endl;
}
