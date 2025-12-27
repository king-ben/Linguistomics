#ifndef UpdateTopology_hpp
#define UpdateTopology_hpp

#include <vector>
#include "Update.hpp"
class Model;
class ParameterAlignment;
class ParameterTree;
class TransitionProbabilityManager;



class UpdateTopology : public Update {

    public:
                                            UpdateTopology(void) = delete;
                                            UpdateTopology(Model* m, RandomVariable* r, ParameterTree* p, const std::vector<ParameterAlignment*>& alnVec);
        std::string                         getUpdateName(void) { return "Branch Length Update"; }
        std::string                         parameterType(void) { return "ParameterTree"; }
        void                                setDependants(void);
        double                              update(void);
        double                              update(double power);
        double                              updateFromPrior(void);
    
    private:
        ParameterTree*                      myParm;
        std::vector<ParameterAlignment*>    myAlignments;
        TransitionProbabilityManager*       tiProbs;
        double                              maxBrlen;
};

#endif
