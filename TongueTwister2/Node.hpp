#ifndef Node_H
#define Node_H

#include <cstring>

#define MAX_NAME_LENGTH 64



class Node {

    friend class Tree;
    
    public:
                            Node(void);
                            Node(int idx);
                           ~Node(void);
        
                            // descendant management 
        void                addDescendant(Node* p);
        bool                isDescendant(Node* p) const;
        size_t              numDescendants(void) const;
        void                orderDescendantsByOffset(void);
        void                removeDescendant(Node* p);
        void                removeDescendants(void);
        
                            // fast inline access to children
        Node*               getLeft(void) const { return left; }
        Node*               getRight(void) const { return right; }
        Node*               getDescendant(size_t idx) const { return (idx == 0) ? left : right; }
        
                            // iteration support without allocation
                            // Usage: for (int i = 0; i < p->numDescendants(); i++) { Node* child = p->getDescendant(i); }
        void                clean(void);
        Node*               getAncestor(void) const { return ancestor; }
        double              getBranchLength(void) const { return length; }
        bool                getFlag(void) const { return flag; }
        int                 getIndex(void) const { return index; }
        bool                getIsLeaf(void) const { return isLeaf; }
        bool                getIsOutgroup(void) const { return isOutgroup; }
        const char*         getName(void) const { return name; }
        char*               getName(void) { return name; }
        int                 getOffset(void) const { return offset; }
        Node*               getSisterNode(void) const;
        void                print(void) const;
        void                setAncestor(Node* p) { ancestor = p; }
        void                setBranchLength(double x) { length = x; }
        void                setFlag(bool tf) { flag = tf; }
        void                setIndex(int x) { index = x; }
        void                setIsLeaf(bool tf) { isLeaf = tf; }
        void                setIsOutgroup(bool tf) { isOutgroup = tf; }
        void                setName(const char* c) { std::strncpy(name, c, MAX_NAME_LENGTH - 1); name[MAX_NAME_LENGTH - 1] = '\0'; }
        void                setOffset(int x) { offset = x; }

    private:
                            // Note: The total size should be 112 bytes, with the most frequently
                            // used data in the first 48 bytes, which shold fit into a single
                            // cache line with some room to spare.
        Node*               ancestor;       // 8 bytes - parent pointer
        Node*               left;           // 8 bytes - first child (replaces NodeSet)
        Node*               right;          // 8 bytes - second child
        double              length;         // 8 bytes - branch length (used in likelihood calc)
        int                 offset;         // 4 bytes - position in node array
        int                 index;          // 4 bytes - logical index
        bool                isLeaf;         // 1 byte  - frequently checked
        bool                flag;           // 1 byte  - traversal flag
        bool                isOutgroup;     // 1 byte
        char                padding_[5];    // 5 bytes - pad to 48 bytes
        
                            // name is rarely accessed during computation
        char                name[MAX_NAME_LENGTH];  // 64 bytes
};

#endif
