#include "Msg.hpp"
#include "StateMatrixFactory.hpp"



StateMatrixFactory::StateMatrixFactory(void) {
    numStates = 0;

    // get the number of states from the configuration file
//    auto it = j.find("NumberOfStates");
//    if (it == j.end())
//        Msg::error("Could not the number of states in the JSON file");
//    int numStates = j["NumberOfStates"];
//    if (numStates <= 1)
//        Msg::error("There must be at least two states in the model");
}

StateMatrixFactory::~StateMatrixFactory(void) {

    for (std::set<StateMatrix_t*>::iterator m=allocated.begin(); m != allocated.end(); m++)
        delete (*m);
}

void StateMatrixFactory::drainPool(void) {

    for (std::vector<StateMatrix_t*>::iterator m=pool.begin(); m != pool.end(); m++)
        {
        allocated.erase(*m);
        delete (*m);
        }
}

StateMatrix_t* StateMatrixFactory::getNode(void) {

    if ( pool.empty() == true )
        {
        /* If the node pool is empty, we allocate a new node and return it. We
           do not need to add it to the node pool. */
        StateMatrix_t* m = new StateMatrix_t;
        allocated.insert( m );
        return m;
        }
    
    // Return a node from the node pool, remembering to remove it from the pool.
    StateMatrix_t* m = pool.back();
    pool.pop_back();
    return m;
}

void StateMatrixFactory::returnToPool(StateMatrix_t* m) {

    //n->clean();
    pool.push_back( m );
}
