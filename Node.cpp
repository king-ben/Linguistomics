#include "Node.hpp"



Node::Node(void) {

    index = 0;
    ancestor = NULL;
    name = "";
    isLeaf = false;
    cladeProbability    = 0.0;
    branchStats.brlen   = 0.0;
    branchStats.prob    = 0.0;
    branchStats.lowerCi = 0.0;
    branchStats.upperCi = 0.0;
}

std::vector<Node*> Node::getDescendants(void) {

    std::vector<Node*> des;
    for (Node* n : neighbors)
        {
        if (n != ancestor)
            des.push_back(n);
        }
    return des;
}
