#ifndef Node_hpp
#define Node_hpp

#include <set>
#include <string>
#include <vector>



class Node {

    public:
                            Node(void);
        virtual            ~Node(void);
        virtual Node*       clone(void) const;
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
        virtual            ~NodeConTree(void) {}
        Node*               clone(void) const override;
        float               getProb(void) { return prob; }
        float               getLowerCi(void) { return lowerCi; }
        float               getUpperCi(void) { return upperCi; }
        void                setLowerCi(float x) { lowerCi = x; }
        void                setUpperCi(float x) { upperCi = x; }
        void                setProb(float x) { prob = x; }
    
    private:
        float               prob;
        float               lowerCi;
        float               upperCi;
};

#endif
