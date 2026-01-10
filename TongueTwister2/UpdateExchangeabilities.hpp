#ifndef UpdateExchangeabilities_hpp
#define UpdateExchangeabilities_hpp

#include "Update.hpp"
class Model;
class ParameterExchangeabilities;



class UpdateExchangeabilities : public Update {

    public:
                                    UpdateExchangeabilities(void) = delete;
                                    UpdateExchangeabilities(Model* m, RandomVariable* r, ParameterExchangeabilities* p);
        double                      getTuningParameter(void) { return tuningValues[lastUpdate]; }
        uint64_t                    getUpdateId(void);
        std::string                 getUpdateName(void) { return updateNames[lastUpdate]; }
        std::string                 parameterType(void) { return "ParameterExchangeabilities"; }
        void                        setDependants(void);
        double                      update(void);
        double                      update(double power);
    
    private:
        ParameterExchangeabilities* myParm;
        static double               minVal;
        size_t                      numRates;
        std::vector<double>         oldValues;
        std::vector<double>         newValues;
        std::vector<double>         alphaForward;
        std::vector<double>         alphaReverse;
        size_t                      lastUpdate;
        std::string                 updateNames[4];
        uint64_t                    updateHashes[4];
        double                      tuningValues[4];
};

#endif
