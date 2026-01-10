#ifndef UpdateBranchLength_hpp
#define UpdateBranchLength_hpp

#include "Update.hpp"
class Model;
class ParameterTree;
class TransitionProbabilities;



class UpdateBranchLength : public Update {

    public:
                                        UpdateBranchLength(void) = delete;
                                        UpdateBranchLength(Model* m, RandomVariable* r, ParameterTree* p);
        double                          getTuningParameter(void) { return updateInfo[0].tuningParameter; }
        size_t                          getUpdateIdx(void) { return 0; }
        uint64_t                        getUpdateId(void) { return updateInfo[0].updateHash; }
        UpdateType                      getUpdateType(void) { return updateInfo[0].updateType;}
        std::string                     getUpdateName(void) { return "Branch Length Update"; }
        std::string                     parameterType(void) { return "ParameterTree"; }
        void                            setDependants(void);
        double                          update(void);
        double                          update(double power);
        double                          updateFromPrior(void);
    
    private:
        ParameterTree*                  myParm;
        TransitionProbabilities*        tiProbs;
        double                          maxBrlen;
};

#endif
