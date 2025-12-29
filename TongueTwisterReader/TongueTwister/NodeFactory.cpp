#include "Node.hpp"
#include "NodeFactory.hpp"



NodeFactory::NodeFactory(void) {

}

NodeFactory::~NodeFactory(void) {

    for (std::set<Node*>::iterator n=allocated.begin(); n != allocated.end(); n++)
        delete (*n);
}

void NodeFactory::drainPool(void) {

    for (std::vector<Node*>::iterator n=pool.begin(); n != pool.end(); n++)
        {
        allocated.erase(*n);
        delete (*n);
        }
}

Node* NodeFactory::getNode(void) {

    if ( pool.empty() == true )
        {
        /* If the node pool is empty, we allocate a new node and return it. We
           do not need to add it to the node pool. */
        Node* n = new Node;
        allocated.insert( n );
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
