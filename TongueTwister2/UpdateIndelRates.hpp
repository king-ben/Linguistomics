#ifndef UpdateIndelRates_hpp
#define UpdateIndelRates_hpp

#include "Update.hpp"
class Model;
class ParameterIndelRates;



class UpdateIndelRates : public Update {

    public:
                                    UpdateIndelRates(void) = delete;
                                    UpdateIndelRates(Model* m, RandomVariable* r, ParameterIndelRates* p);
        std::string                 getUpdateName(void) { return "Indel Rates Update"; }
        void                        notifyDependants(void);
        std::string                 parameterType(void) { return "ParameterIndelRates"; }
        void                        setDependants(void);
        double                      update(void);
        double                      update(double power);
        double                      updateFromPrior(void);
    
    private:
        ParameterIndelRates*        myParm;
};

#endif
