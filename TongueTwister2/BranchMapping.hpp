#ifndef BranchMapping_hpp
#define BranchMapping_hpp

#include <string>
#include <unordered_map>
#include <vector>
class Node;
class Tree;



class BranchMapping {

    friend class ParameterTree;

    public:
                                                    BranchMapping(void) = default;
        bool                                        affectsSubtree(int fullTreeNodeOffset);
        void                                        buildMapping(Tree* fullTree, Tree* subTree);
        void                                        copyBranchLengths(Tree* src, Tree* dst);
        void                                        print(void);
        void                                        propagateChange(Node* fullTreeNode, double delta, Tree* subTree);
        double                                      recomputeBranchLength(int subTreeNodeOffset, Tree* fullTree);
        void                                        recomputeAllBranchLengths(Tree* fullTree, Tree* subTree);

    private:
        Node*                                       findMRCA(Tree* t, std::vector<std::string>& taxa);
        Node*                                       findNodeByOffset(Tree* t, int offset);
        void                                        getDescendantLeaves(Node* p, std::vector<std::string>& leaves);
        std::unordered_map<int, std::vector<int>>   subToFullBranches;
        std::unordered_map<int, int>                fullToSubBranch;
};

#endif
