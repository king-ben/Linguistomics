#ifndef Node_hpp
#define Node_hpp

#include <set>
#include <string>
#include <vector>

struct BranchStatistics {

    float prob;
    float lowerCi;
    float upperCi;
};


class Node {

    public:
                            Node(void);
        virtual            ~Node(void);
        void                addNeighbor(Node* p) { neighbors.insert(p); }
        Node*               getAncestor(void) { return ancestor; }
        float               getBrlen(void) { return brlen; }
        std::vector<Node*>  getDescendants(void);
        int                 getIndex(void) { return index; }
        bool                getIsLeaf(void) { return isLeaf; }
        std::string         getName(void) { return name; }
        std::set<Node*>&    getNeighbors(void) { return neighbors; }
        size_t              getOffset(void) { return offset; }
        void                removeAllNeighbors(void) { neighbors.clear(); }
        void                removeNeighbor(Node* p) { neighbors.erase(p); }
        void                setAncestor(Node* p) { ancestor = p; }
        void                setBrlen(float x) { brlen = x; }
        void                setIndex(int x) { index = x; }
        void                setIsLeaf(bool tf) { isLeaf = tf; }
        void                setName(std::string s) { name = s; }
        void                setOffset(size_t x) { offset = x; }
    
    protected:
        std::set<Node*>     neighbors;
        std::string         name;
        Node*               ancestor;
        float               brlen;
        int                 index;
        size_t              offset;
        bool                isLeaf;
};

class NodeConTree : public Node {

    public:
                            NodeConTree(void);
        float               getProb(void) { return branchStats.prob; }
        float               getLowerCi(void) { return branchStats.lowerCi; }
        float               getUpperCi(void) { return branchStats.upperCi; }
        void                setLowerCi(float x) { branchStats.lowerCi = x; }
        void                setUpperCi(float x) { branchStats.upperCi = x; }
        void                setProb(float x) { branchStats.prob = x; }
    
    private:
        BranchStatistics    branchStats;
};

#endif
