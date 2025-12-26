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
        std::vector<double>&        getRates(void) { return rates[0]; }
        const std::vector<double>&  getRates(void) const { return rates[0]; }
        void                        keep(void);
        double                      lnPriorProbability(void);
        void                        print(void);
        void                        restore(void);
    
    private:
        size_t                      numStates;
        size_t                      numRates;
        std::vector<double>         rates[2];
        std::vector<double>         alpha;
        bool                        isPriorFlat;
};

#endif
