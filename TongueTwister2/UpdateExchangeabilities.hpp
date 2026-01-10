#ifndef UpdateExchangeabilities_hpp
#define UpdateExchangeabilities_hpp

#include "Update.hpp"
class Model;
class ParameterExchangeabilities;



class UpdateExchangeabilities : public Update {

    public:
                                    UpdateExchangeabilities(void) = delete;
                                    UpdateExchangeabilities(Model* m, RandomVariable* r, ParameterExchangeabilities* p);
        double                      getTuningParameter(void) { return updateInfo[lastUpdate].tuningParameter; }
        size_t                      getUpdateIdx(void) { return updateInfo[lastUpdate].updateIdx; }
        uint64_t                    getUpdateId(void) { return updateInfo[lastUpdate].updateHash; }
        UpdateType                  getUpdateType(void) { return updateInfo[lastUpdate].updateType;}
        std::string                 getUpdateName(void) { return updateInfo[lastUpdate].updateName; }
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
//        std::vector<std::string>    updateNames;
//        std::vector<uint64_t>       updateHashes;
};

#endif
