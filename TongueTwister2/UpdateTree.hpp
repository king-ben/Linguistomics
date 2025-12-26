#ifndef UpdateTree_hpp
#define UpdateTree_hpp

#include "Update.hpp"
class Model;
class ParameterTree;
class TransitionProbabilityManager;



class UpdateTree : public Update {

    public:
                                        UpdateTree(void) = delete;
                                        UpdateTree(Model* m, RandomVariable* r, ParameterTree* p);
        std::string                     getUpdateName(void) { return "Tree Update"; }
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
