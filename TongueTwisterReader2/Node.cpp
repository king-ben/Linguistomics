#include "Node.hpp"



Node::Node(void) {

    index = 0;
    ancestor = nullptr;
    name = "";
    isLeaf = false;
    offset = 0;
    brlen = 0.0f;
}

Node::~Node(void) {

}

Node* Node::clone(void) const {

    Node* copy = new Node;
    copy->name = name;
    copy->brlen = brlen;
    copy->index = index;
    copy->offset = offset;
    copy->isLeaf = isLeaf;
    // Note: ancestor and neighbors are set by Tree::clone()
    return copy;
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

NodeConTree::NodeConTree(void) : Node() {

    prob = 0.0f;
    lowerCi = 0.0f;
    upperCi = 0.0f;
}

Node* NodeConTree::clone(void) const {

    NodeConTree* copy = new NodeConTree;
    copy->name = name;
    copy->brlen = brlen;
    copy->index = index;
    copy->offset = offset;
    copy->isLeaf = isLeaf;
    copy->prob = prob;
    copy->lowerCi = lowerCi;
    copy->upperCi = upperCi;
    return copy;
}
