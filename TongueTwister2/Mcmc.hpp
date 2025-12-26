#ifndef Mcmc_hpp
#define Mcmc_hpp

#include <cstddef>
#include <string>
class Model;
class RandomVariable;
class UpdateManager;



class Mcmc {

    public:
                        Mcmc(void) = delete;
                        Mcmc(RandomVariable* r, Model* m);
        virtual        ~Mcmc(void);
        virtual void    run(void) = 0;
    
    protected:
        RandomVariable* rng;
        Model*          myModel;
        UpdateManager*  updateMngr;
        std::string     baseOutputFileName;
        int             numCycles;
        int             printFrequency;
        int             sampleFrequency;
        int             numDigits;
};

#endif
