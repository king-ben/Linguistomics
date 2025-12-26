#include "Node.hpp"
#include "NodeComparator.hpp"



bool NodeComparator::operator()(const Node* n1, const Node* n2) const {
    
    if (n1->getOffset() > n2->getOffset())
        return true;
    return false;
}
