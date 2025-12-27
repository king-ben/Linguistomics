#ifndef UpdateExchangeabilities_hpp
#define UpdateExchangeabilities_hpp

#include "Update.hpp"
class Model;
class ParameterExchangeabilities;



class UpdateExchangeabilities : public Update {

    public:
                                    UpdateExchangeabilities(void) = delete;
                                    UpdateExchangeabilities(Model* m, RandomVariable* r, ParameterExchangeabilities* p);
        std::string                 getUpdateName(void) { return "Exchangeability Rates Update"; }
        std::string                 parameterType(void) { return "ParameterExchangeabilities"; }
        void                        setDependants(void);
        double                      update(void);
        double                      update(double power);
        double                      updateFromPrior(void);
    
    private:
        ParameterExchangeabilities* myParm;
        static double               minVal;
        size_t                      numRates;
        std::vector<double>         oldValues;
        std::vector<double>         newValues;
        std::vector<double>         alphaForward;
        std::vector<double>         alphaReverse;
};

#endif
