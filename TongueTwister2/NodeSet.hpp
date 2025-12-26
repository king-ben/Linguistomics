#ifndef NodeSet_hpp
#define NodeSet_hpp

#include <set>
#include "NodeComparator.hpp"
class Node;


class NodeSet {

    public:
        Node*                           operator[](size_t idx);
        Node*                           operator[](size_t idx) const;
        void                            addNode(Node* p) { nodes.insert(p); }
        void                            deleteNode(Node* p) { nodes.erase(p); }
        void                            deleteAllNodes(void) { nodes.clear(); }
        std::set<Node*,NodeComparator>& getNodes(void) { return nodes; }
        size_t                          size(void) { return nodes.size(); }
        size_t                          size(void) const { return nodes.size(); }
    
    private:
        std::set<Node*,NodeComparator>  nodes;
};

#endif
