#ifndef ParameterTree_hpp
#define ParameterTree_hpp

#include <set>
#include <string>
#include <unordered_map>
#include <vector>
#include "BranchMapping.hpp"
#include "json.hpp"
#include "Parameter.hpp"
class Node;
class Tree;



struct TreePair {

            TreePair(void) { trees[0] = nullptr; trees[1] = nullptr; }
            TreePair(Tree* t1, Tree* t2) { trees[0] = t1; trees[1] = t2; }
    Tree*   trees[2];
};



class ParameterTree : public Parameter {

    public:
                                                    ParameterTree(void) = delete;
                                                    ParameterTree(Model* m, RandomVariable* r, std::string n);
                                                   ~ParameterTree(void);
        void                                        applyNniToSubtrees(void);
        double                                      getBrlenLambda(void) { return brlenLambda; }
        size_t                                      getNumTaxa(void);
        const std::vector<std::string>&             getCanonicalTaxonList(void) const { return canonicalTaxonList; }
        double                                      getMaximumBrlen(void) { return MAX_BRLEN; }
        Tree*                                       getTree(void) { return fullTree.trees[0]; }
        Tree*                                       getTree(unsigned& mask);
        Tree*                                       getTree(const unsigned& mask);
        void                                        getTrees(std::vector<Tree*>& allTrees);
        void                                        initialize(const std::set<unsigned>& uniqueTaxonCombinations);
        void                                        initializeBranchMappings(void);
        void                                        keep(void);
        double                                      lnPriorProbability(void);
        void                                        print(void);
        void                                        restore(void);
        void                                        setBranchLength(Node* p, double v);
        void                                        setTopologyChanged(bool tf) { topologyChanged = tf; }
    
    private:
        bool                                        checkTipToTipDistances(double threshhold);
        void                                        makeSubtree(Tree& t, const unsigned& taxonMask);
        void                                        rebuildSubtrees(void);
        void                                        restoreTopology(void);
        void                                        saveTopology(void);
        size_t                                      taxonIndex(const std::string& tName);
        constexpr static const double               MAX_BRLEN = 2.0;

                                                    // accessed frequently during MCMC updates
        TreePair                                    fullTree;         // 16 bytes
        double                                      brlenLambda;      // 8 bytes
        double                                      lnProbLessMax;    // 8 bytes
        
                                                    // containers accessed during tree operations
        std::unordered_map<unsigned,TreePair>       subTrees;
        std::unordered_map<unsigned,BranchMapping>  branchMappings;
        std::vector<std::string>                    canonicalTaxonList;
        
                                                    // flag checked only for LOCAL move, which is currently off
        bool                                        topologyChanged;
};

#endif
