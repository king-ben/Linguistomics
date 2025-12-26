#include "NodeSet.hpp"



Node* NodeSet::operator[](size_t idx) {

    if (idx == 0)
        return *nodes.begin();
    else if (idx == 1 && nodes.size() == 2)
        return *nodes.rbegin();

    size_t i = 0;
    for (auto p : nodes)
        {
        if (i == idx)
            return p;
        i++;
        }
        
    return nullptr;
}

Node* NodeSet::operator[](size_t idx) const {

    if (idx == 0)
        return *nodes.begin();
    else if (idx == 1 && nodes.size() == 2)
        return *nodes.rbegin();

    size_t i = 0;
    for (auto p : nodes)
        {
        if (i == idx)
            return p;
        i++;
        }
        
    return nullptr;
}
