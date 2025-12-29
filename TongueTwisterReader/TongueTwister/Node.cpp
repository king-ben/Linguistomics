#include <iostream>
#include "Node.hpp"
#include "NodeSet.hpp"
#include "RbBitSet.h"
#include "Tree.hpp"



Node::Node(void) {

    index                 = 0;
    isLeaf                = false;
    name                  = "";
    ancestor              = NULL;
    length                = 0.0;
    nodeTime              = 0.0;
    offset                = 0;
    myTree                = NULL;
    descendants           = new NodeSet;
    transitionProbability = NULL;
}

Node::Node(int idx) {

    index                 = 0;
    isLeaf                = false;
    name                  = "";
    ancestor              = NULL;
    length                = 0.0;
    nodeTime              = 0.0;
    offset                = idx;
    myTree                = NULL;
    descendants           = new NodeSet;
    transitionProbability = NULL;
}

Node::~Node(void) {

    delete descendants;
}

void Node::addDescendant(Node* p) {

    descendants->addNode(p);
}

void Node::clean(void) {

    index      = 0;
    isLeaf     = false;
    name       = "";
    ancestor   = NULL;
    length     = 0.0;
    nodeTime   = 0.0;
    myTree     = NULL;
    descendants->deleteAllNodes();
    transitionProbability = NULL;
}

std::vector<Node*> Node::getDescendantsVector(void) {

    std::set<Node*,CompNode>& des = descendants->getNodes();
    std::vector<Node*> dv;
    for (Node* n : des)
        dv.push_back( n );
    return dv;
}

void Node::getDescendantsVector(std::vector<Node*>& dv) {

    std::set<Node*,CompNode>& des = descendants->getNodes();
    dv.clear();
    for (Node* n : des)
        dv.push_back( n );
}

Node* Node::getSisterNode(void) {

    if (ancestor == NULL)
        return NULL;
    std::set<Node*,CompNode>& des = ancestor->getDescendants()->getNodes();
    if (des.size() != 2)
        return NULL;
    for (Node* n : des)
        {
        if (n != this)
            return n;
        }
    return NULL;
}

bool Node::isDescendant(Node* p) {

    std::set<Node*,CompNode>& des = descendants->getNodes();
    std::set<Node*,CompNode>::iterator it = des.find(p);
    if (it != des.end())
        return true;
    return false;
}

size_t Node::numDescendants(void) {

    return descendants->size();
}

double Node::oldestDescendantTime(void) {

    double oldestDescendantTime = 0.0;
    for (Node* n : descendants->getNodes())
        {
        if (n->getNodeTime() > oldestDescendantTime)
            oldestDescendantTime = n->getNodeTime();
        }
    return oldestDescendantTime;
}

void Node::print(void) {

    std::cout << "Node (" << this << ")" << std::endl;
    std::cout << "    Ancestor: " << ancestor->getIndex() << std::endl;
    std::cout << "    Descendants: ";
    std::set<Node*,CompNode>& des = descendants->getNodes();
    for(Node* n : des)
        {
        std::cout << n->getIndex() << " ";
        }
    std::cout << std::endl;
    std::cout << "        index: " << index  << std::endl;
    std::cout << "       isLeaf: " << isLeaf << std::endl;
    std::cout << "         name: \"" << name << "\""  << std::endl;
}

void Node::removeDescendant(Node* p) {

    descendants->deleteNode(p);
}

void Node::removeDescendants(void) {

    descendants->deleteAllNodes();
}
