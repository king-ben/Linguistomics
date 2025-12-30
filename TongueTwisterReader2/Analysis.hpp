#ifndef Analysis_hpp
#define Analysis_hpp

#include <fstream>
#include <string>
#include <vector>
#include "Container.hpp"
class AlignmentDistribution;
class Exchangeabilities;
class Partition;
class RandomVariable;
class RateMatrix;
class StateFrequencies;
class ThreadPool;
class TreeSamples;


class Analysis;
namespace AnalysisComparison {

    void compare(Analysis* a1, Analysis* a2, size_t nSamples=100);
};



class Analysis {

    friend void AnalysisComparison::compare(Analysis* a1, Analysis* a2, size_t nSamples);

    public:
                                            Analysis(void) = delete;
                                            Analysis(RandomVariable* r, ThreadPool* tp, std::string directoryName, double burnFraction);
                                           ~Analysis(void);
        std::string                         modelName(void);
        void                                print(void);
        void                                printSorted(void);
        void                                randomlyChooseFreqs(std::vector<float>& f);
        DoubleMatrix*                       randomlyChooseRateMatrix(void);
        DoubleMatrix*                       randomlyChooseRateMatrixAndFreqs(std::vector<float>& f);
        void                                writeNytril(std::string pathName);
    
    private:
        void                                nytrilOutput(std::ofstream& file, int maxAlignment);
        RandomVariable*                     rng;
        ThreadPool*                         pool;
        Exchangeabilities*                  rates;
        Partition*                          part;
        RateMatrix*                         Q;
        StateFrequencies*                   freqs;
        TreeSamples*                        trees;
        std::vector<AlignmentDistribution*> alignments;
};

#endif
