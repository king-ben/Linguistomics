#include <iostream>
#include "Node.hpp"



Node::Node(void) {

    clean();
}

Node::Node(int idx) {

    clean();
    offset = idx;
}

Node::~Node(void) {

}

void Node::addDescendant(Node* p) {

    descendants.addNode(p);
}

void Node::clean(void) {

    index      = 0;
    offset     = 0;
    name[0]    = '\0';
    ancestor   = nullptr;
    length     = 0.0;
    isLeaf     = false;
    flag       = false;
    isOutgroup = false;
    descendants.deleteAllNodes();
}

Node* Node::getSisterNode(void) {

    if (ancestor == nullptr)
        return nullptr;
    std::set<Node*,NodeComparator>& des = ancestor->getDescendants().getNodes();
    if (des.size() != 2)
        return nullptr;
    for (Node* n : des)
        {
        if (n != this)
            return n;
        }
    return nullptr;
}

bool Node::isDescendant(Node* p) {

    std::set<Node*,NodeComparator>& des = descendants.getNodes();
    std::set<Node*,NodeComparator>::iterator it = des.find(p);
    if (it != des.end())
        return true;
    return false;
}

size_t Node::numDescendants(void) {

    return descendants.size();
}

void Node::print(void) {

    std::cout << "Node (" << this << ")" << std::endl;
    std::cout << "    Ancestor: " << ancestor->getIndex() << std::endl;
    std::cout << "    Descendants: ";
    std::set<Node*,NodeComparator>& des = descendants.getNodes();
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

    descendants.deleteNode(p);
}

void Node::removeDescendants(void) {

    descendants.deleteAllNodes();
}
