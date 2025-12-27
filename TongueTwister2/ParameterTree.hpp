#ifndef ParameterTree_hpp
#define ParameterTree_hpp

#include <set>
#include <string>
#include <unordered_map>
#include <vector>
#include "BranchMapping.hpp"
#include "json.hpp"
#include "Parameter.hpp"
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
    
    private:
        bool                                        checkTipToTipDistances(double threshhold);
        void                                        makeSubtree(Tree& t, const unsigned& taxonMask);
        size_t                                      taxonIndex(const std::string& tName);
        constexpr static const double               MAX_BRLEN = 2.0; // 2 is a very long branch and we don't want to mess up the branch length hashing
        double                                      lnProbLessMax;
        double                                      brlenLambda;
        TreePair                                    fullTree;
        std::unordered_map<unsigned,TreePair>       subTrees;
        std::vector<std::string>                    canonicalTaxonList;
        std::unordered_map<unsigned,BranchMapping>  branchMappings;
};

#endif
