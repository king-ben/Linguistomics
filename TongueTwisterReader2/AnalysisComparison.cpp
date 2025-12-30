#include "Analysis.hpp"
#include "AnalysisComparison.hpp"
#include "Container.hpp"



void AnalysisComparison::compare(Analysis* a1, Analysis* a2, size_t nSamples) {

    size_t numStates = 0;
    
    // compare the rate matrices
    double sum = 0.0;
    for (size_t n=0; n<nSamples; n++)
        {
        DoubleMatrix* m1 = a1->randomlyChooseRateMatrix();
        DoubleMatrix* m2 = a2->randomlyChooseRateMatrix();
        numStates = m1->getNumRows();
        
        double ss = 0.0;
        for (size_t i=0; i<numStates; i++)
            {
            for (size_t j=i+1; j<numStates; j++)
                {
                ss += std::pow((*m1)(i,j) - (*m2)(i,j), 2.0);
                ss += std::pow((*m1)(j,i) - (*m2)(j,i), 2.0);
                }
            }
        sum += ss;
        }
    double dQ = sum / nSamples;
        
    // compare the state frequencies
    std::vector<float> f1(numStates);
    std::vector<float> f2(numStates);
    sum = 0.0;
    for (size_t n=0; n<nSamples; n++)
        {
        a1->randomlyChooseFreqs(f1);
        a2->randomlyChooseFreqs(f2);
        double ss = 0.0;
        for (size_t i=0; i<numStates; i++)
            ss += std::pow(f1[i] - f2[i], 2.0);
        sum += ss;
        }
    double dF = sum / nSamples;
    
    // compare rates, if possible
    
    std::cout << "Model Comparison:   " << a1->modelName() << "-" << a2->modelName() << std::endl;
    std::cout << "Average Q distance: " << dQ << std::endl;
    std::cout << "Average F distance: " << dF << std::endl;
}
