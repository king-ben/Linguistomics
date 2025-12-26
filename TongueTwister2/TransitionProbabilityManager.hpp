#ifndef TransitionProbabilityManager_hpp
#define TransitionProbabilityManager_hpp

#include "Container.hpp"



/* Abstract base class for transition probability management.
   Allows switching between different compute backends (CPU threaded,
   CPU batched with Accelerate, GPU batched) at runtime. */
class TransitionProbabilityManager {

    public:
        virtual                            ~TransitionProbabilityManager(void) = default;
        
                                            // lookup transition matrix by branch length
        virtual DoubleMatrix*               find(double branchLength) = 0;
        virtual DoubleMatrix&               getTransitionProbability(double v) = 0;
        
                                            // MCMC cycle management
        virtual void                        keep(void) = 0;
        virtual void                        restore(void) = 0;
        
                                            // update methods 
        virtual void                        updateBranch(void) = 0;
        virtual void                        updateAllBranches(void) = 0;
        virtual void                        rebuildFromTree(void) = 0;
        
                                            // diagnostics
        virtual void                        print(void) const = 0;
};

#endif
