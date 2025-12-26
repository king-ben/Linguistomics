#ifndef McmcMarginal_hpp
#define McmcMarginal_hpp

#include "Mcmc.hpp"
#include "McmcPhase.hpp"
class McmcPhase;
class McmcTimer;
class Model;
class RandomVariable;
class SteppingStones;



class McmcMarginal : public Mcmc {

    public:
                            McmcMarginal(void) = delete;
                            McmcMarginal(RandomVariable* r, Model* m);
                           ~McmcMarginal(void);
        void                run(void);
            
    private:
        int                 calculateNumCycles(void);
        std::pair<int,int>  calculatePrintDigits(void);
        void                printStatus(size_t stoneIdx, size_t cycle, size_t cnt, std::string& phase, double lnL1, double lnL2, McmcTimer* timer);
        McmcPhase*          mcmcPhases;
        SteppingStones*     samples;
        int                 numCycles;
        int                 numDigitsPhase;
};

#endif
