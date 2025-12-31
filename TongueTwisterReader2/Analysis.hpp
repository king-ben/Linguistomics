#ifndef Analysis_hpp
#define Analysis_hpp

#include <fstream>
#include <string>
#include <vector>
#include "Container.hpp"
class AlignmentDistribution;
class Exchangeabilities;
class IndelRates;
class Partition;
class RandomVariable;
class RateMatrix;
class StateFrequencies;
class ThreadPool;
class TreeSamples;


class Analysis;
namespace AnalysisComparison {

    void compare(Analysis* a1, Analysis* a2, size_t numStates, size_t nSamples=100);
};



class Analysis {

    friend void AnalysisComparison::compare(Analysis* a1, Analysis* a2, size_t numStates, size_t nSamples);

    public:
                                            Analysis(void) = delete;
                                            Analysis(RandomVariable* r, ThreadPool* tp, std::string directoryName, double burnFraction);
                                           ~Analysis(void);
        size_t                              getNumStates(void) { return numStates; }
        std::string                         modelName(void);
        void                                print(void);
        void                                printSorted(void);
        void                                randomlyChooseFreqs(std::vector<float>& f);
        DoubleMatrix*                       randomlyChooseRateMatrix(void);
        DoubleMatrix*                       randomlyChooseRateMatrixAndFreqs(std::vector<float>& f);
        void                                writeNytril(std::string pathName);
    
    private:
        void                                nytrilOutput(std::ofstream& file, int maxAlignment);
        void                                writeMatrix(std::ofstream& file, DoubleMatrix &m, std::string name);
        size_t                              numStates;
        RandomVariable*                     rng;
        ThreadPool*                         pool;
        IndelRates*                         indelRates;
        Exchangeabilities*                  rates;
        Partition*                          part;
        RateMatrix*                         Q;
        RateMatrix*                         aveQ;
        StateFrequencies*                   freqs;
        StateFrequencies*                   ncFreqs;
        TreeSamples*                        trees;
        std::vector<AlignmentDistribution*> alignments;
};

#endif
