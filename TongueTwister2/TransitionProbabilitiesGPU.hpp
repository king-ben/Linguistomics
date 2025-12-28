#ifndef TransitionProbabilitiesGPU_hpp
#define TransitionProbabilitiesGPU_hpp

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
class GPUMatrixExponentialBatch;



enum class ComputeBackend {
    Auto,           // automatically choose based on batch size 
    CPUThreaded,    // original thread pool approach 
    CPUBatched,     // batched Accelerate without thread pool 
    GPUBatched      // GPU batch (future) 
};



class TransitionProbabilitiesGPU : public TransitionProbabilityManager {

    public:
                                            TransitionProbabilitiesGPU(void) = delete;
                                            TransitionProbabilitiesGPU(ThreadPool* p, ParameterTree* t, Ctmc* sm, size_t cleanupFrequency = 100);
                                           ~TransitionProbabilitiesGPU(void) override;
        
                                            // lookup transition matrix by branch length
        DoubleMatrix*                       find(double branchLength) override;
        DoubleMatrix&                       getTransitionProbability(double v) override;
        
                                            // MCMC cycle management
        void                                keep(void) override;
        void                                restore(void) override;
        
                                            // update methods
        void                                updateBranch(void) override;
        void                                updateAllBranches(void) override;
        void                                rebuildFromTree(void) override;
        
                                            // backend selection
        void                                setComputeBackend(ComputeBackend backend);
        ComputeBackend                      getComputeBackend(void) const { return computeBackend; }
        void                                setBatchThreshold(size_t threshold) { batchThreshold = threshold; }
        size_t                              getBatchThreshold(void) const { return batchThreshold; }
        
                                            // diagnostics
        void                                print(void) const override;
        bool                                isGPUAvailable(void) const;
        const char*                         getGPUDeviceName(void) const;
        
    private:
        void                                ensureCalculatorPoolCapacity(size_t n);
        void                                shrinkCalculatorPoolIfNeeded(void);
        void                                computeDirtyMatrices(void);
        void                                computeDirtyMatricesThreaded(void);
        void                                computeDirtyMatricesBatched(void);
        void                                cleanupOrphanedEntries(void);
        const double*                       getRateMatrixData(void) const;
        
        size_t                              numStates;
        ThreadPool*                         pool;
        ParameterTree*                      myTree;
        Ctmc*                               subModel;
        
                                            // matrix pool (owned) - must be declared before map
        MatrixPool                          matrixPool;
        
                                            // transition matrix map (uses matrixPool)
        TransitionMatrixMap*                map;
        
                                            // compute backend selection
        ComputeBackend                      computeBackend;
        size_t                              batchThreshold;
        GPUMatrixExponentialBatch*          gpuBatcher;
        
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
        
                                            // avoids repeated allocations
        std::vector<uint64_t>               usedKeysBuffer;
        std::vector<uint64_t>               keysToEraseBuffer;
        std::vector<double>                 batchBranchLengths;
        std::vector<DoubleMatrix*>          batchOutputs;
};

#endif
