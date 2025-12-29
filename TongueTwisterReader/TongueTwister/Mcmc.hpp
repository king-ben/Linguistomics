#ifndef Mcmc_hpp
#define Mcmc_hpp

#include <chrono>
#include <fstream>
#include <string>
#include <vector>
class Model;
class RandomVariable;

typedef std::chrono::high_resolution_clock::time_point Timer;

class Mcmc {

    public:
                            Mcmc(void) = delete;
                            Mcmc(Model** m, RandomVariable* r);
        void                run(void);
    
    private:
        std::vector<double> calculatePowers(int numStones, double alpha, double beta);
        void                chooseModelsToSwap(int& idx1, int& idx2);
        void                closeOutputFiles(void);
        std::string         formattedTime(Timer& t1, Timer& t2);
        void                initialize(void);
        Model*              getColdModel(void);
        int                 numDigits(double lnX);
        void                openOutputFiles(void);
        double              power(int idx);
        void                print(int gen, double* curLnL, double* curLnP, bool accept, Timer& t1, Timer& t2);
        void                print(int powIdx, int numPowers, std::string phase, int gen, double* curLnL, double* curLnP);
        double              safeExponentiation(double lnX);
        void                sample(int gen, double lnL, double lnP);
        void                runPathSampling(void);
        void                runPosterior(void);
        Model**             modelPtr;
        RandomVariable*     rv;
        int                 numChains;
        int                 numMcmcCycles;
        int                 printFrequency;
        int                 sampleFrequency;
        double              temperature;
        std::ofstream       parmStrm;
        std::ofstream       treeStrm;
        double*             parmValues;
        int                 numParmValues;
        int                 maxPriorPrint;
        int                 maxLikePrint;
        int                 maxGenPrint;
        int                 preburninLength;
        int                 numTunes;
        int                 tuneLength;
        int                 burninLength;
        int                 sampleLength;
        int                 stoneSampleFrequency;
};

#endif
