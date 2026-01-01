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


struct Distances {

    double wD;
    double qD;
    double aD;
};

class Analysis;
namespace AnalysisComparison {

    void                        compareAnalyses(std::string fileName, std::vector<Analysis*>& analyses);
    Distances                   compare(Analysis* a1, Analysis* a2, size_t numStates, size_t nSamples=100);
    double                      compareAlignments(Analysis* a1, Analysis* a2);
    std::pair<double,double>    compareQ(Analysis* a1, Analysis* a2, size_t numStates, size_t nSamples);
};



class Analysis {

    friend Distances AnalysisComparison::compare(Analysis* a1, Analysis* a2, size_t numStates, size_t nSamples);

    public:
                                                    Analysis(void) = delete;
                                                    Analysis(RandomVariable* r, ThreadPool* tp, std::string directoryName, double burnFraction);
                                                   ~Analysis(void);
        std::string                                 getName(void) { return analysisName; }
        size_t                                      getNumStates(void) { return numStates; }
        const std::vector<AlignmentDistribution*>   getAlignments(void) const { return alignments; }
        RateMatrix*                                 getRateMatrix(void) { return Q; }
        RateMatrix*                                 getWeightMatrix(void) { return aveQ; }
        std::string                                 modelName(void);
        void                                        print(void);
        void                                        printSorted(void);
        void                                        randomlyChooseFreqs(std::vector<float>& f);
        DoubleMatrix*                               randomlyChooseRateMatrix(void);
        DoubleMatrix*                               randomlyChooseRateMatrixAndFreqs(std::vector<float>& f);
        void                                        setName(std::string s) { analysisName = s; }
        void                                        writeNytril(std::string pathName);
        void                                        writeR(std::string pathName, std::string fileName, std::string modelName, int idx);
    
    private:
        void                                        nytrilOutput(std::ofstream& file, int maxAlignment);
        void                                        writeMatrix(std::ofstream& file, DoubleMatrix &m, std::string name);
        size_t                                      numStates;
        RandomVariable*                             rng;
        ThreadPool*                                 pool;
        IndelRates*                                 indelRates;
        Exchangeabilities*                          rates;
        Partition*                                  part;
        RateMatrix*                                 Q;
        RateMatrix*                                 aveQ;
        StateFrequencies*                           freqs;
        StateFrequencies*                           ncFreqs;
        TreeSamples*                                trees;
        std::string                                 analysisName;
        std::vector<AlignmentDistribution*>         alignments;
};

#endif
