#ifndef McmcPosterior_hpp
#define McmcPosterior_hpp

#include <cstddef>
#include "Mcmc.hpp"
class McmcTimer;
class Model;
class RandomVariable;



class McmcPosterior : public Mcmc {

    public:
                        McmcPosterior(void) = delete;
                        McmcPosterior(RandomVariable* r, Model* m);
                       ~McmcPosterior(void);
        void            run(void);
    
    protected:
        void            printStatus(int cycle, int nCycles, double lnL1, double lnL2, McmcTimer* timer);
};

#endif
