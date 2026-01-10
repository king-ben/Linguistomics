#ifndef ParameterExchangeabilities_hpp
#define ParameterExchangeabilities_hpp

#include <vector>
#include "Parameter.hpp"



class ParameterExchangeabilities : public Parameter {

    public:
                                    ParameterExchangeabilities(void) = delete;
                                    ParameterExchangeabilities(Model* m, RandomVariable* r, std::string n, size_t ns);
                                    ParameterExchangeabilities(Model* m, RandomVariable* r, std::string n, size_t ns, size_t ng);
                                   ~ParameterExchangeabilities(void);
        std::vector<double>&        getAlpha(void) { return alpha; }
        size_t                      getNumRates(void) { return numRates; }
        size_t                      getNumStates(void) { return numStates; }
        std::vector<double>&        getRates(void) { return rates[0]; }
        const std::vector<double>&  getRates(void) const { return rates[0]; }
        std::vector<double>&        getRates(size_t idx) { return rates[idx]; }
        const std::vector<double>&  getRates(size_t idx) const { return rates[idx]; }
        void                        keep(void);
        double                      lnPriorProbability(void);
        void                        print(void);
        void                        restore(void);
    
    private:
                                    // accessed during rate matrix construction
        size_t                      numStates;         // 8 bytes
        size_t                      numRates;          // 8 bytes
        
                                    // vectors (24 bytes each)
        std::vector<double>         rates[2];          // 48 bytes
        std::vector<double>         alpha;             // 24 bytes
};

#endif
