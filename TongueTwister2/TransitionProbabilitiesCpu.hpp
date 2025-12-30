#ifndef TransitionProbabilitiesCpu_hpp
#define TransitionProbabilitiesCpu_hpp

#include <vector>
#include "Container.hpp"
#include "MathCache.hpp"
#include "MatrixPool.hpp"
#include "TransitionProbabilities.hpp"

class Ctmc;
class ParameterTree;
class ThreadPool;
class TransitionMatrixMap;
class TransitionProbabilityCalculator;



class TransitionProbabilitiesCpu : public TransitionProbabilities {

    public:
                                            TransitionProbabilitiesCpu(void) = delete;
                                            TransitionProbabilitiesCpu(ThreadPool* p, ParameterTree* t, Ctmc* sm, size_t cleanupFrequency = 100);
                                           ~TransitionProbabilitiesCpu(void) override;
        
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
        
                                            // pointers accessed every update cycle
        ThreadPool*                         pool;
        ParameterTree*                      myTree;
        Ctmc*                               subModel;
        TransitionMatrixMap*                map;
        TransitionProbabilityCalculator**   calculatorPool;
        
                                            // size values used in computation
        size_t                              numStates;
        size_t                              calculatorPoolCapacity;
        size_t                              calculatorPoolSize;
        size_t                              maxCalculatorsUsedRecently;
        size_t                              calculatorShrinkCounter;
        size_t                              cleanupFrequency;
        size_t                              cleanupCounter;
        
                                            // matrix pool (large object)
        MatrixPool                          matrixPool;
        
                                            // single-threaded cache for single branch updates
        MathCache                           singleCache;
        
                                            // vectors for cleanup operations
        std::vector<uint64_t>               newEntriesThisProposal;
        std::vector<uint64_t>               usedKeysBuffer;
        std::vector<uint64_t>               keysToEraseBuffer;
        
                                            // static constants
        static const size_t                 CALCULATOR_SHRINK_FREQUENCY = 100;
        static const size_t                 CALCULATOR_MIN_POOL_SIZE = 16;
};

#endif
