#ifndef NodeSet_hpp
#define NodeSet_hpp

#include <set>
#include "Node.hpp"



class NodeSet {

    public:
        void                        addNode(Node* p) { nodes.insert(p); }
        void                        deleteNode(Node* p) { nodes.erase(p); }
        void                        deleteAllNodes(void) { nodes.clear(); }
        std::set<Node*,CompNode>&   getNodes(void) { return nodes; }
        size_t                      size(void) { return nodes.size(); }
    
    private:
        std::set<Node*,CompNode>    nodes;
};

#endif
