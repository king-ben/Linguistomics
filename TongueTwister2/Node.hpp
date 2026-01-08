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
        void                    addDescendant(Node* p);
        void                    clean(void);
        Node*                   getAncestor(void) const { return ancestor; }
        Node*                   getLeft(void) const { return left; }
        Node*                   getRight(void) const { return right; }
        double                  getBranchLength(void) const { return length; }
        Node*                   getDescendant(size_t idx) const { return (idx == 0) ? left : right; }
        bool                    getFlag(void) const { return flag; }
        int                     getIndex(void) const { return index; }
        bool                    getIsLeaf(void) const { return isLeaf; }
        bool                    getIsOutgroup(void) const { return isOutgroup; }
        const char*             getName(void) const { return name; }
        char*                   getName(void) { return name; }
        int                     getOffset(void) const { return offset; }
        Node*                   getSisterNode(void) const;
        bool                    isDescendant(Node* p) const;
        size_t                  numDescendants(void) const;
        void                    orderDescendantsByOffset(void);
        void                    print(void) const;
        void                    removeDescendant(Node* p);
        void                    removeDescendants(void);
        void                    setAncestor(Node* p) { ancestor = p; }
        void                    setBranchLength(double x) { length = x; }
        void                    setFlag(bool tf) { flag = tf; }
        void                    setIndex(int x) { index = x; }
        void                    setIsLeaf(bool tf) { isLeaf = tf; }
        void                    setIsOutgroup(bool tf) { isOutgroup = tf; }
        void                    setName(const char* c) { std::strncpy(name, c, MAX_NAME_LENGTH - 1); name[MAX_NAME_LENGTH - 1] = '\0'; }
        void                    setOffset(int x) { offset = x; }

    private:
                                // The total size should be 112 bytes, with the most frequently
                                // used data in the first 48 bytes, which shold fit into a single
                                // cache line with some room to spare.
        Node*                   ancestor;               //  8 bytes -   8 bytes - ancestor pointer
        Node*                   left;                   //  8 bytes -  16 bytes - first descendant
        Node*                   right;                  //  8 bytes -  24 bytes - second descendant
        double                  length;                 //  8 bytes -  32 bytes - branch length (used in likelihood calculation)
        int                     offset;                 //  4 bytes -  36 bytes - position in node array
        int                     index;                  //  4 bytes -  40 bytes - logical index
        bool                    isLeaf;                 //  1 byte  -  41 bytes - frequently checked
        bool                    flag;                   //  1 byte  -  42 bytes - traversal flag
        bool                    isOutgroup;             //  1 byte  -  43 bytes 
        char                    name[MAX_NAME_LENGTH];  // 64 bytes - 112 bytes
};

#endif
