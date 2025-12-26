#ifndef NodeComparator_hpp
#define NodeComparator_hpp

class Node;



struct NodeComparator {

    bool operator()(const Node* n1, const Node* n2) const;
};

#endif
