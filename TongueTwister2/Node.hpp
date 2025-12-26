#ifndef Node_H
#define Node_H

#include <set>
#include <string>
#include <vector>
#include "NodeSet.hpp"

#define MAX_NAME_LENGTH 64



class Node {

    friend class Tree;
    
    public:
                            Node(void);
                            Node(int idx);
                           ~Node(void);
        void                addDescendant(Node* p);
        void                clean(void);
        Node*               getAncestor(void) { return ancestor; }
        double              getBranchLength(void) { return length; }
        NodeSet&            getDescendants(void) { return descendants; }
        const NodeSet&      getDescendants(void) const { return descendants; }
        bool                getFlag(void) { return flag; }
        int                 getIndex(void) { return index; }
        bool                getIsLeaf(void) { return isLeaf; }
        bool                getIsOutgroup(void) { return isOutgroup; }
        char*               getName(void) { return name; }
        int                 getOffset(void) { return offset; }
        int                 getOffset(void) const { return offset; }
        Node*               getSisterNode(void);
        bool                isDescendant(Node* p);
        size_t              numDescendants(void);
        void                print(void);
        void                removeDescendant(Node* p);
        void                removeDescendants(void);
        void                setAncestor(Node* p) { ancestor = p; }
        void                setBranchLength(double x) { length = x; }
        void                setFlag(bool tf) { flag = tf; }
        void                setIndex(int x) { index = x; }
        void                setIsLeaf(bool tf) { isLeaf = tf; }
        void                setIsOutgroup(bool tf) { isOutgroup = tf; }
        void                setName(char* c) { strcpy(name, c); }
        void                setName(const char* c) { strcpy(name, c); }
        void                setOffset(int x) { offset = x; }

    private:
        NodeSet             descendants;
        Node*               ancestor;
        char                name[MAX_NAME_LENGTH];
        double              length;
        int                 offset;
        int                 index;
        bool                isLeaf;
        bool                flag;
        bool                isOutgroup;
};

#endif

