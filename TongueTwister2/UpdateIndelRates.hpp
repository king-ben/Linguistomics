#ifndef UpdateIndelRates_hpp
#define UpdateIndelRates_hpp

#include "Update.hpp"
class Model;
class ParameterIndelRates;



class UpdateIndelRates : public Update {

    public:
                                    UpdateIndelRates(void) = delete;
                                    UpdateIndelRates(Model* m, RandomVariable* r, ParameterIndelRates* p);
        double                      getTuningParameter(void) { return updateInfo[0].tuningParameter; }
        size_t                      getUpdateIdx(void) { return 0; }
        uint64_t                    getUpdateId(void) { return updateInfo[0].updateHash; }
        UpdateType                  getUpdateType(void) { return updateInfo[0].updateType;}
        std::string                 getUpdateName(void) { return "Indel Rates Update"; }
        std::string                 parameterType(void) { return "ParameterIndelRates"; }
        void                        setDependants(void);
        double                      update(void);
        double                      update(double power);
        double                      updateFromPrior(void);
    
    private:
        ParameterIndelRates*        myParm;
};

#endif
