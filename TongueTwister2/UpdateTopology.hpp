#ifndef UpdateTopology_hpp
#define UpdateTopology_hpp


#include "Update.hpp"
class Model;
class ParameterTree;
class TransitionProbabilityManager;



class UpdateTopology : public Update {

    public:
                                        UpdateTopology(void) = delete;
                                        UpdateTopology(Model* m, RandomVariable* r, ParameterTree* p);
        std::string                     getUpdateName(void) { return "Branch Length Update"; }
        void                            notifyDependants(void);
        std::string                     parameterType(void) { return "ParameterTree"; }
        void                            setDependants(void);
        double                          update(void);
        double                          update(double power);
        double                          updateFromPrior(void);
    
    private:
        ParameterTree*                  myParm;
        TransitionProbabilityManager*   tiProbs;
        double                          maxBrlen;
};

#endif
