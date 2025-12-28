#ifndef TransitionProbabilities_hpp
#define TransitionProbabilities_hpp

#include <vector>
#include "Container.hpp"
#include "MathCache.hpp"
#include "MatrixPool.hpp"
#include "TransitionProbabilityManager.hpp"

class Ctmc;
class ParameterTree;
class ThreadPool;
class TransitionMatrixMap;
class TransitionProbabilityCalculator;



class TransitionProbabilities : public TransitionProbabilityManager {

    public:
                                            TransitionProbabilities(void) = delete;
                                            TransitionProbabilities(ThreadPool* p, ParameterTree* t, Ctmc* sm, size_t cleanupFrequency = 100);
                                           ~TransitionProbabilities(void) override;
        
                                            // lookup transition matrix by branch length
        DoubleMatrix*                       find(double branchLength) override;
        DoubleMatrix&                       getTransitionProbability(double v) override;
        
                                            // MCMC cycle management
        void                                keep(void) override;
        void                                restore(void) override;
        
                                            // update after a single branch length change (scans subtrees for new values) 
        void                                updateBranch(void) override;
        
                                            // update all transition probabilities (e.g., when rate matrix parameters change) 
        void                                updateAllBranches(void) override;
        
                                            // rebuild map from tree (call after topology changes)
        void                                rebuildFromTree(void) override;
        
        void                                print(void) const override;
        
    private:
        void                                ensureCalculatorPoolCapacity(size_t n);
        void                                shrinkCalculatorPoolIfNeeded(void);
        void                                computeDirtyMatrices(void);
        void                                cleanupOrphanedEntries(void);
        
        size_t                              numStates;
        ThreadPool*                         pool;
        ParameterTree*                      myTree;
        Ctmc*                               subModel;
        
                                            // matrix pool (owned) - must be declared before map
        MatrixPool                          matrixPool;
        
                                            // transition matrix map (uses matrixPool)
        TransitionMatrixMap*                map;
        
                                            // calculator pool with shrinking support
        TransitionProbabilityCalculator**   calculatorPool;
        size_t                              calculatorPoolCapacity;
        size_t                              calculatorPoolSize;
        size_t                              maxCalculatorsUsedRecently;
        size_t                              calculatorShrinkCounter;
        static const size_t                 CALCULATOR_SHRINK_FREQUENCY = 100;
        static const size_t                 CALCULATOR_MIN_POOL_SIZE = 16;
        
                                            // single-threaded cache for single branch updates
        MathCache                           singleCache;
        
                                            // automatic cleanup of orphaned entries
        size_t                              cleanupFrequency;
        size_t                              cleanupCounter;
        std::vector<uint64_t>               newEntriesThisProposal;
        
                                            // avoids repeated allocations during cleanup
        std::vector<uint64_t>               usedKeysBuffer;
        std::vector<uint64_t>               keysToEraseBuffer;
};

#endif
