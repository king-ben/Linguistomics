#ifndef UpdateTopology_hpp
#define UpdateTopology_hpp

#include <vector>
#include "Update.hpp"
class Model;
class Node;
class ParameterAlignment;
class ParameterTree;
class TransitionProbabilityManager;
class Tree;
class UpdateAlignment;



class UpdateTopology : public Update {

    public:
                                            UpdateTopology(void) = delete;
                                            UpdateTopology(Model* m, RandomVariable* r, ParameterTree* p, const std::vector<UpdateAlignment*>& alnVec);
        std::string                         getUpdateName(void) { return "LOCAL Topology Update"; }
        std::string                         parameterType(void) { return "ParameterTree"; }
        void                                setDependants(void);
        double                              update(void);
        double                              update(double power);
        double                              updateFromPrior(void);
        
                                            // Override: alignments are also modified during topology updates
        std::vector<Parameter*>             getAdditionalModifiedParameters(void) override;
    
    private:
        void                                applyNni(Node* u, Node* v, Node* a, Node* c);
        Node*                               randomlyChooseInternalBranch(Tree* t);

        ParameterTree*                      myParm;
        std::vector<UpdateAlignment*>       myAlignmentUpdates;
        TransitionProbabilityManager*       tiProbs;
        double                              maxBrlen;
        double                              tuning;
};

#endif
