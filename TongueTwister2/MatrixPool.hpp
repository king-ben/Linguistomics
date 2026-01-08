#ifndef MatrixPool_hpp
#define MatrixPool_hpp

#include <vector>
#include "Container.hpp"



class MatrixPool {

    public:
                                        MatrixPool(void) = delete;
                                        MatrixPool(size_t nRows, size_t nCols, size_t initialCapacity = 64);
                                       ~MatrixPool(void);
        
                                        // non-copyable
                                        MatrixPool(const MatrixPool&) = delete;
                                        MatrixPool& operator=(const MatrixPool&) = delete;
                                        MatrixPool(MatrixPool&&) = delete;
                                        MatrixPool& operator=(MatrixPool&&) = delete;
        
                                        // core operations
        DoubleMatrix*                   acquire(void);
        void                            release(DoubleMatrix* m);
        
                                        // batch operations - more efficient for bulk acquire/release
        void                            acquireBatch(std::vector<DoubleMatrix*>& out, size_t count);
        void                            releaseBatch(DoubleMatrix* const* matrices, size_t count);
        
                                        // pool management
        void                            reserve(size_t capacity);
        void                            shrinkToFit(size_t targetFreeCount);
        void                            clear(void);
        
                                        // accessors
        size_t                          getNumRows(void) const { return numRows; }
        size_t                          getNumCols(void) const { return numCols; }
        size_t                          totalAllocated(void) const { return allMatrices.size(); }
        size_t                          numFree(void) const { return freeList.size(); }
        size_t                          numInUse(void) const { return allMatrices.size() - freeList.size(); }
        
                                        // diagnostics
        void                            printStatistics(void) const;
    
    private:
        void                            grow(size_t count);
        size_t                          numRows;
        size_t                          numCols;
        std::vector<DoubleMatrix*>      freeList;    // freeList is accessed most frequently (push_back/pop_back)
        std::vector<DoubleMatrix*>      allMatrices; // only accessed for diagnostics or cleanup
        static constexpr size_t         minGrowth = 16;
        static constexpr double         growthFactor = 1.5;
};

#endif
