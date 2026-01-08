#include "Node.hpp"
#include "NodeFactory.hpp"



NodeFactory::NodeFactory(void) : initialSize(4096) {

    Node* initialNodeVector = new Node[initialSize]; // makes certain that most of the nodes are contiguous
    for (size_t i=0; i<initialSize; i++)
        {
        Node* n = &initialNodeVector[i];
        allocated.push_back(n);
        pool.push_back(n);
        }
}

NodeFactory::~NodeFactory(void) {

    delete [] allocated[0];
    for (size_t i=initialSize; i<allocated.size(); i++)
        delete allocated[i];
}

Node* NodeFactory::getNode(void) {

    if (pool.empty() == true)
        {
        /* If the node pool is empty, we allocate a new node and return it. We
           do not need to add it to the node pool. */
        Node* n = new Node;
        allocated.push_back(n);
        return n;
        }
    
    // Return a node from the node pool, remembering to remove it from the pool.
    Node* n = pool.back();
    pool.pop_back();
    return n;
}

void NodeFactory::returnToPool(Node* n) {

    n->clean();
    pool.push_back( n );
}
