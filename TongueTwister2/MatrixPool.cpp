#include <algorithm>
#include <iostream>
#include "Container.hpp"
#include "MatrixPool.hpp"



MatrixPool::MatrixPool(size_t nRows, size_t nCols, size_t initialCapacity) :
    numRows(nRows), numCols(nCols) {

    // pre-allocate storage for pointers
    allMatrices.reserve(initialCapacity);
    freeList.reserve(initialCapacity);
        
    // allocate initial matrices
    if (initialCapacity > 0)
        grow(initialCapacity);
}

MatrixPool::~MatrixPool(void) {

    // delete all matrices owned by the pool
    for (DoubleMatrix* m : allMatrices)
        delete m;
}

DoubleMatrix* MatrixPool::acquire(void) {

    if (freeList.empty())
        {
        // grow by growth factor, but at least minGrowth
        size_t currentSize = allMatrices.size();
        size_t growBy = static_cast<size_t>(currentSize * growthFactor);
        if (growBy < minGrowth)
            growBy = minGrowth;
        grow(growBy);
        }
    
    DoubleMatrix* m = freeList.back();
    freeList.pop_back();
    return m;
}

void MatrixPool::acquireBatch(std::vector<DoubleMatrix*>& out, size_t count) {

    // ensure we have enough matrices
    if (freeList.size() < count)
        {
        size_t needed = count - freeList.size();
        size_t growBy = static_cast<size_t>(allMatrices.size() * growthFactor);
        if (growBy < needed)
            growBy = needed;
        if (growBy < minGrowth)
            growBy = minGrowth;
        grow(growBy);
        }
    
    // acquire from end of free list (cache friendly)
    out.reserve(out.size() + count);
    for (size_t i = 0; i < count; ++i)
        {
        out.push_back(freeList.back());
        freeList.pop_back();
        }
}

void MatrixPool::clear(void) {

    // delete all matrices and reset
    for (DoubleMatrix* m : allMatrices)
        delete m;
    allMatrices.clear();
    freeList.clear();
}

void MatrixPool::grow(size_t count) {

    allMatrices.reserve(allMatrices.size() + count);
    freeList.reserve(freeList.size() + count);
    for (size_t i = 0; i < count; ++i)
        {
        DoubleMatrix* m = new DoubleMatrix(numRows, numCols);
        allMatrices.push_back(m);
        freeList.push_back(m);
        }
}

void MatrixPool::printStatistics(void) const {

    size_t matrixBytes = numRows * numCols * sizeof(double);
    size_t totalBytes = allMatrices.size() * matrixBytes;
    size_t freeBytes = freeList.size() * matrixBytes;
    size_t inUseBytes = totalBytes - freeBytes;
    
    std::cout << "MatrixPool Statistics:" << std::endl;
    std::cout << "  Matrix dimensions: " << numRows << " x " << numCols << std::endl;
    std::cout << "  Bytes per matrix:  " << matrixBytes << std::endl;
    std::cout << "  Total allocated:   " << allMatrices.size() << " (" << totalBytes / 1024 << " KB)" << std::endl;
    std::cout << "  In use:            " << numInUse() << " (" << inUseBytes / 1024 << " KB)" << std::endl;
    std::cout << "  Free:              " << numFree() << " (" << freeBytes / 1024 << " KB)" << std::endl;
}

void MatrixPool::release(DoubleMatrix* m) {

    // no validation for speed - caller is responsible for only
    // releasing matrices that came from this pool
    // no zeroing - caller can zero if needed before use
    freeList.push_back(m);
}

void MatrixPool::releaseBatch(DoubleMatrix* const* matrices, size_t count) {

    freeList.reserve(freeList.size() + count);
    for (size_t i = 0; i < count; ++i)
        freeList.push_back(matrices[i]);
}

void MatrixPool::reserve(size_t capacity) {

    if (capacity > allMatrices.size())
        grow(capacity - allMatrices.size());
}

void MatrixPool::shrinkToFit(size_t targetFreeCount) {

    // remove excess free matrices beyond targetFreeCount
    // this can help reduce memory footprint after a spike in usage
    
    if (freeList.size() <= targetFreeCount)
        return;
    
    size_t toRemove = freeList.size() - targetFreeCount;
    
    // collect matrices to delete
    std::vector<DoubleMatrix*> toDelete;
    toDelete.reserve(toRemove);
    for (size_t i = 0; i < toRemove; ++i)
        {
        toDelete.push_back(freeList.back());
        freeList.pop_back();
        }
    
    // remove from allMatrices
    for (DoubleMatrix* m : toDelete)
        {
        auto it = std::find(allMatrices.begin(), allMatrices.end(), m);
        if (it != allMatrices.end())
            {
            // swap with last and pop for O(1) removal
            std::swap(*it, allMatrices.back());
            allMatrices.pop_back();
            }
        delete m;
        }
}

