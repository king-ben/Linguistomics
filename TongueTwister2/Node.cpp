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

    if (left == nullptr)
        left = p;
    else if (right == nullptr)
        right = p;
    else
        {
        // Error: trying to add more than 2 children to a binary tree node
        std::cerr << "Error: Node::addDescendant - node already has 2 children" << std::endl;
        }
}

void Node::clean(void) {

    ancestor   = nullptr;
    left       = nullptr;
    right      = nullptr;
    length     = 0.0;
    offset     = 0;
    index      = 0;
    isLeaf     = false;
    flag       = false;
    isOutgroup = false;
    name[0]    = '\0';
}

bool Node::isDescendant(Node* p) const {

    return (left == p || right == p);
}

size_t Node::numDescendants(void) const {

    if (left == nullptr)
        return 0;
    if (right == nullptr)
        return 1;
    return 2;
}

Node* Node::getSisterNode(void) const {

    if (ancestor == nullptr)
        return nullptr;
        
    Node* ancLeft = ancestor->getLeft();
    Node* ancRight = ancestor->getRight();
    
    if (ancLeft == this)
        return ancRight;
    if (ancRight == this)
        return ancLeft;
        
    return nullptr;
}

void Node::orderDescendantsByOffset(void) {

    if (left != nullptr && right != nullptr)
        {
        if (left->offset > right->offset)
            {
            Node* temp = left;
            left = right;
            right = temp;
            }
        }
}

void Node::print(void) const {

    std::cout << "Node (" << this << ")" << std::endl;
    std::cout << "    Ancestor: ";
    if (ancestor)
        std::cout << ancestor->getIndex();
    else
        std::cout << "nullptr";
    std::cout << std::endl;
    
    std::cout << "    Descendants: ";
    if (left)
        std::cout << left->getIndex() << " ";
    if (right)
        std::cout << right->getIndex() << " ";
    if (!left && !right)
        std::cout << "(none)";
    std::cout << std::endl;
    
    std::cout << "        index: " << index  << std::endl;
    std::cout << "       offset: " << offset << std::endl;
    std::cout << "       isLeaf: " << isLeaf << std::endl;
    std::cout << "       length: " << length << std::endl;
    std::cout << "         name: \"" << name << "\""  << std::endl;
}

void Node::removeDescendant(Node* p) {

    if (left == p)
        {
        left = right;  // shift right to left position
        right = nullptr;
        }
    else if (right == p)
        {
        right = nullptr;
        }
}

void Node::removeDescendants(void) {

    left = nullptr;
    right = nullptr;
}
