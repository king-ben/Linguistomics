#ifndef UpdateBranchLength_hpp
#define UpdateBranchLength_hpp

#include "Update.hpp"
class Model;
class ParameterTree;
class TransitionProbabilityManager;



class UpdateBranchLength : public Update {

    public:
                                        UpdateBranchLength(void) = delete;
                                        UpdateBranchLength(Model* m, RandomVariable* r, ParameterTree* p);
        std::string                     getUpdateName(void) { return "Branch Length Update"; }
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
