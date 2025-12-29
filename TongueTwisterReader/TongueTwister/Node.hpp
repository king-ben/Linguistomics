#ifndef Node_H
#define Node_H

#include <set>
#include <string>
#include <vector>
#include "Container.hpp"
class NodeSet;
class RbBitSet;
class Tree;



class Node {

    public:
                            Node(void);
                            Node(int idx);
                           ~Node(void);
        void                addDescendant(Node* p);
        void                clean(void);
        Node*               getAncestor(void) { return ancestor; }
        double              getBranchLength(void) { return length; }
        NodeSet*&           getDescendants(void) { return descendants; }
        std::vector<Node*>  getDescendantsVector(void);
        void                getDescendantsVector(std::vector<Node*>& dv);
        bool                getFlag(void) { return flag; }
        int                 getIndex(void) { return index; }
        bool                getIsLeaf(void) { return isLeaf; }
        std::string         getName(void) { return name; }
        int                 getOffset(void) { return offset; }
        int                 getOffset(void) const { return offset; }
        Node*               getSisterNode(void);
        DoubleMatrix*       getTransitionProbability(void) { return transitionProbability; }
        double              getNodeTime(void) { return nodeTime; }
        bool                isDescendant(Node* p);
        size_t              numDescendants(void);
        double              oldestDescendantTime(void);
        void                print(void);
        void                removeDescendant(Node* p);
        void                removeDescendants(void);
        void                setAncestor(Node* p) { ancestor = p; }
        void                setBranchLength(double x) { length = x; }
        void                setFlag(bool tf) { flag = tf; }
        void                setIndex(int x) { index = x; }
        void                setIsLeaf(bool tf) { isLeaf = tf; }
        void                setMyTree(Tree* t) { myTree = t; }
        void                setName(std::string s) { name = s; }
        void                setOffset(int x) { offset = x; }
        void                setTransitionProbability(DoubleMatrix* p) { transitionProbability = p; }
        void                setNodeTime(double x) { nodeTime = x; }

    protected:
        NodeSet*            descendants;
        Node*               ancestor;
        Tree*               myTree;
        DoubleMatrix*       transitionProbability;
        std::string         name;
        double              length;
        double              nodeTime;
        int                 offset;
        int                 index;
        bool                isLeaf;
        bool                flag;
};

struct CompNode {

    bool operator()(const Node* n1, const Node* n2) const {
        
        if (n1->getOffset() > n2->getOffset())
            return true;
        return false;
        }
};

#endif
