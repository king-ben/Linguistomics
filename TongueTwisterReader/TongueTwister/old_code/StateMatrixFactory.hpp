#ifndef StateMatrixFactory_hpp
#define StateMatrixFactory_hpp

#include <set>
#include <vector>
#include "EigenSystem.hpp"



class StateMatrixFactory {

    public:
        static StateMatrixFactory&  nodeFactory(void) {
                                        static StateMatrixFactory smf;
                                        return smf;
                                    }
        void                        drainPool(void);
        StateMatrix_t*              getNode(void);
        int                         getNumAllocated(void) { return (int) allocated.size(); }
        int                         getNumInPool(void) { return (int)pool.size(); }
        void                        returnToPool(StateMatrix_t* m);
    
    private:
                                    StateMatrixFactory(void);
                                   ~StateMatrixFactory(void);
                                    StateMatrixFactory(const StateMatrixFactory& tp) = delete;
        std::vector<StateMatrix_t*> pool;
        std::set<StateMatrix_t*>    allocated;
        int                         numStates;
};

#endif
