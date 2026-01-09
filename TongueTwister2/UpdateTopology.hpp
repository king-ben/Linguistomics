#ifndef UpdateTopology_hpp
#define UpdateTopology_hpp

#include <vector>
#include "Update.hpp"
class Model;
class Node;
class ParameterAlignment;
class ParameterTree;
class TransitionProbabilities;
class Tree;
class UpdateAlignment;



class UpdateTopology : public Update {

    public:
                                            UpdateTopology(void) = delete;
                                            UpdateTopology(Model* m, RandomVariable* r, ParameterTree* p, const std::vector<UpdateAlignment*>& alnVec);
        double                              getTuningParameter(void) { return tuningParameter; }
        uint64_t                            getUpdateId(void) { return updateId; }
        std::string                         getUpdateName(void) override { return "LOCAL Topology Update"; }
        std::string                         parameterType(void) override { return "ParameterTree"; }
        void                                setDependants(void) override;
        double                              update(void) override;
        double                              update(double power) override;
        double                              updateFromPrior(void) override;
        
                                            // alignments are also modified during topology updates
        std::vector<Parameter*>             getAdditionalModifiedParameters(void) override;
    
    private:
        void                                applyNni(Node* u, Node* v, Node* a, Node* c);
        Node*                               randomlyChooseInternalBranch(Tree* t);

        ParameterTree*                      myParm;
        std::vector<UpdateAlignment*>       myAlignmentUpdates;
        TransitionProbabilities*            tiProbs;
        double                              maxBrlen;
};

#endif
