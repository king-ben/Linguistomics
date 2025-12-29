#ifndef ParameterTree_hpp
#define ParameterTree_hpp

#include <map>
#include <string>
#include "Parameter.hpp"
#include "RbBitSet.h"
class Alignment;
class RandomVariable;
class Tree;

struct TreePair {
            TreePair(void) { trees[0] = NULL; trees[1] = NULL; }
            TreePair(Tree* t1, Tree* t2) { trees[0] = t1; trees[1] = t2; }
    Tree*   trees[2];
};



class ParameterTree : public Parameter {

    public:
                                        ParameterTree(void) = delete;
                                        ParameterTree(const ParameterTree& pt) = delete;
                                        ParameterTree(RandomVariable* r, Model* m, std::string treeStr, std::vector<std::string> tNames, std::vector<Alignment*>& wordAlignments, double itl);
                                       ~ParameterTree(void);
        void                            accept(void);
        void                            clearSubtrees(void);
        void                            fillParameterValues(double* x, int& start, int maxNumValues);
        Tree*                           getActiveTree(void) { return fullTree.trees[0]; }
        Tree*                           getActiveTree(RbBitSet& mask);
        Tree*                           getActiveTree(const RbBitSet& mask);
        std::map<RbBitSet,TreePair>&    getSubtrees(void) { return subTrees; }
        std::string                     getJsonString(void);
        std::string                     getHeader(void) { return ""; }
        int                             getNumValues(void) { return 1; }
        std::string                     getUpdateName(int idx);
        std::string                     getString(void);
        char*                           getCString(void);
        double                          lnPriorProbability(void);
        void                            print(void);
        void                            printNewick(void);
        void                            reject(void);
        double                          update(int iter);
        void                            updateSubtrees(void);
                
    private:
        bool                            checkSubtreeCompatibility(double tolerance);
        int                             countMaskBits(std::vector<bool>& m);
        void                            initializeSubtrees(std::vector<Alignment*>& alns);
        void                            nniArea(std::vector<Node*>& backbone, Node*& incidentNode);
        double                          updateBrlen(void);
        double                          updateBranchlengthsFromPrior(void);
        double                          updateTopologyFromPrior(void);
        double                          updateNni(void);
        double                          updateSpr(void);
        double                          updateTreeLength(void);
        double                          lambda;
        TreePair                        fullTree;
        std::map<RbBitSet,TreePair>     subTrees;
        bool                            isClockConstrained;
};

#endif
