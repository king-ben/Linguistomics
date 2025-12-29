#ifndef Analysis_hpp
#define Analysis_hpp

#include <fstream>
#include <string>
#include <vector>
class AlignmentDistribution;
class Exchangeabilities;
class Partition;
class RandomVariable;
class RateMatrix;
class StateFrequencies;
class ThreadPool;



class Analysis {

    public:
                                            Analysis(void) = delete;
                                            Analysis(RandomVariable* r, ThreadPool* tp, std::string directoryName, double burnFraction);
                                           ~Analysis(void);
        void                                print(void);
        void                                printSorted(void);
        void                                writeNytril(std::string pathName);
    
    private:
        void                                nytrilOutput(std::ofstream& file, int maxAlignment);
        RandomVariable*                     rng;
        ThreadPool*                         pool;
        Exchangeabilities*                  rates;
        Partition*                          part;
        RateMatrix*                         Q;
        StateFrequencies*                   freqs;
        std::vector<AlignmentDistribution*> alignments;
};

#endif
