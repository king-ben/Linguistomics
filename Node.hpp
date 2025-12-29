#ifndef Node_hpp
#define Node_hpp

#include <set>
#include <string>
#include <vector>

struct BranchStatistics {

    double brlen;
    double prob;
    double lowerCi;
    double upperCi;
};


class Node {

    public:
                            Node(void);
        void                addNeighbor(Node* p) { neighbors.insert(p); }
        Node*               getAncestor(void) { return ancestor; }
        double              getBrlen(void) { return branchStats.brlen; }
        double              getCladeProbability(void) { return cladeProbability; }
        std::vector<Node*>  getDescendants(void);
        int                 getIndex(void) { return index; }
        bool                getIsLeaf(void) { return isLeaf; }
        double              getLowerCi(void) { return branchStats.lowerCi; }
        std::string         getName(void) { return name; }
        std::set<Node*>&    getNeighbors(void) { return neighbors; }
        double              getProb(void) { return branchStats.prob; }
        double              getUpperCi(void) { return branchStats.upperCi; }
        void                removeAllNeighbors(void) { neighbors.clear(); }
        void                removeNeighbor(Node* p) { neighbors.erase(p); }
        void                setAncestor(Node* p) { ancestor = p; }
        void                setBrlen(double x) { branchStats.brlen = x; }
        void                setCladeProbability(double x) { cladeProbability = x; }
        void                setIndex(int x) { index = x; }
        void                setIsLeaf(bool tf) { isLeaf = tf; }
        void                setLowerCi(double x) { branchStats.lowerCi = x; }
        void                setName(std::string s) { name = s; }
        void                setProb(double x) { branchStats.prob = x; }
        void                setUpperCi(double x) { branchStats.upperCi = x; }
    
    private:
        std::set<Node*>     neighbors;
        std::string         name;
        Node*               ancestor;
        double              cladeProbability;
        BranchStatistics    branchStats;
        int                 index;
        bool                isLeaf;
};

#endif
