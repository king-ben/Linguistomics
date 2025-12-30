#ifndef TransitionProbabilitiesGpu_hpp
#define TransitionProbabilitiesGpu_hpp

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
class GPUMatrixExponentialBatch;



enum class ComputeBackend {
    Auto,           // automatically choose based on batch size 
    CPUThreaded,    // original thread pool approach 
    CPUBatched,     // batched Accelerate without thread pool 
    GPUBatched      // GPU batch (future) 
};



class TransitionProbabilitiesGpu : public TransitionProbabilities {

    public:
                                            TransitionProbabilitiesGpu(void) = delete;
                                            TransitionProbabilitiesGpu(ThreadPool* p, ParameterTree* t, Ctmc* sm, size_t cleanupFrequency = 100);
                                           ~TransitionProbabilitiesGpu(void) override;
        
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
        
                                            // accessed frequently during MCMC
        ThreadPool*                         pool;
        ParameterTree*                      myTree;
        Ctmc*                               subModel;
        TransitionMatrixMap*                map;
        GPUMatrixExponentialBatch*          gpuBatcher;
        TransitionProbabilityCalculator**   calculatorPool;
        
                                            // size values used in loops
        size_t                              numStates;
        size_t                              batchThreshold;
        size_t                              calculatorPoolCapacity;
        size_t                              calculatorPoolSize;
        size_t                              maxCalculatorsUsedRecently;
        size_t                              calculatorShrinkCounter;
        size_t                              cleanupFrequency;
        size_t                              cleanupCounter;
        
                                            // matrix pool (large object, keep together)
        MatrixPool                          matrixPool;
        
                                            // single-threaded cache
        MathCache                           singleCache;
        
                                            // vectors for batch operations
        std::vector<uint64_t>               newEntriesThisProposal;
        std::vector<uint64_t>               usedKeysBuffer;
        std::vector<uint64_t>               keysToEraseBuffer;
        std::vector<double>                 batchBranchLengths;
        std::vector<DoubleMatrix*>          batchOutputs;
        
                                            // enum at end (4 bytes, padding acceptable at end)
        ComputeBackend                      computeBackend;
        
                                            // static constants
        static const size_t                 CALCULATOR_SHRINK_FREQUENCY = 100;
        static const size_t                 CALCULATOR_MIN_POOL_SIZE = 16;
};

#endif
