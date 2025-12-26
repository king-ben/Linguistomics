#include "Msg.hpp"
#include "TransitionMatrixManager.hpp"



TransitionMatrixManager::TransitionMatrixManager(void) : numStates(0), isInitialized(false) {

}

TransitionMatrixManager::~TransitionMatrixManager(void) {

    for (DoubleMatrix* m : allocated)
        delete m;
}

DoubleMatrix* TransitionMatrixManager::getMatrix(void) {

    if (pool.empty() == true)
        {
        /* If the matrix pool is empty, we allocate a new matrix and return it. We
           do not need to add it to the matrix pool. */
        DoubleMatrix* m = new DoubleMatrix(numStates, numStates);
        allocated.insert(m);
        return m;
        }
    
    // return a matrix from the matrix pool, remembering to remove it from the pool.
    DoubleMatrix* m = pool.back();
    pool.pop_back();
    return m;
}

void TransitionMatrixManager::initialize(size_t ns, size_t initialCapacity) {

    // if already initialized with the same numStates, optionally grow the pool
    if (isInitialized)
        {
        if (numStates != ns)
            Msg::error("TransitionMatrixManager: Cannot reinitialize with different numStates");
        // Already initialized with correct size - optionally add more matrices
        // For now, just return without adding more
        return;
        }

    numStates = ns;
    
    for (size_t i=0; i<initialCapacity; i++)
        {
        DoubleMatrix* m = new DoubleMatrix(numStates, numStates);
        allocated.insert(m);
        pool.push_back(m);
        }

    isInitialized = true;
}

void TransitionMatrixManager::returnMatrix(DoubleMatrix* m) {

    pool.push_back(m);
}
