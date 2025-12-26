#ifndef TransitionMatrixManager_hpp
#define TransitionMatrixManager_hpp

#include <set>
#include <vector>
#include "Container.hpp"



class TransitionMatrixManager {

    public:
        static TransitionMatrixManager& getInstance(void) {
                                            static TransitionMatrixManager mngr;
                                            return mngr;
                                        }
        DoubleMatrix*                   getMatrix(void);
        void                            initialize(size_t ns, size_t initialCapacity);
        void                            returnMatrix(DoubleMatrix* m);
    
    private:
                                        TransitionMatrixManager(void);
                                       ~TransitionMatrixManager(void);
                                        TransitionMatrixManager(const TransitionMatrixManager& s) = delete;
        std::vector<DoubleMatrix*>      pool;
        std::set<DoubleMatrix*>         allocated;
        size_t                          numStates;
        bool                            isInitialized;
                        
};

#endif
