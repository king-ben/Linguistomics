#ifndef TransitionMatrixMap_hpp
#define TransitionMatrixMap_hpp

#include <cstdint>
#include <cstddef>
#include "Container.hpp"
class MatrixPool;



static constexpr double BRANCH_LENGTH_SCALE = 1.0e9;

inline uint64_t branchLengthToKey(double branchLength) {
    return static_cast<uint64_t>(branchLength * BRANCH_LENGTH_SCALE);
}

inline double keyToBranchLength(uint64_t key) {
    return static_cast<double>(key) / BRANCH_LENGTH_SCALE;
}



class TransitionMatrixMap {

    public:
                            TransitionMatrixMap(void) = delete;
                            TransitionMatrixMap(size_t numStates, size_t initialCapacity, MatrixPool* pool);
                            TransitionMatrixMap(const TransitionMatrixMap&) = delete;
                           ~TransitionMatrixMap(void);
        TransitionMatrixMap& operator=(const TransitionMatrixMap&) = delete;
        
                            // basic accessors
        size_t              size(void) const { return numEntries; }
        size_t              capacity(void) const { return tableCapacity; }
        size_t              getNumStates(void) const { return numStates; }
        bool                empty(void) const { return numEntries == 0; }
        
                            // lookup (read-only, does not mark dirty)
        DoubleMatrix*       find(uint64_t key) const;
        DoubleMatrix*       find(double branchLength) const;
        
                            // get or create entry (for initial population)
        DoubleMatrix*       getOrCreate(uint64_t key);
        DoubleMatrix*       getOrCreate(double branchLength);
        
                            // get matrix for modification - backs up before first modification
        DoubleMatrix*       getForModification(uint64_t key);
        DoubleMatrix*       getForModification(double branchLength);
        
                            // mark a specific entry as needing recalculation (backs up first)
        void                markForUpdate(uint64_t key);
        void                markForUpdate(double branchLength);
        
                            // mark ALL entries for update (full recalculation)
        void                markAllForUpdate(void);
        
                            // MCMC support
        void                keep(void);      // Accept proposal: clear dirty flags
        void                restore(void);   // Reject proposal: restore from backup
        
                            // iteration over all entries (O(n))
        uint64_t            getEntryKey(size_t entryIndex) const;
        DoubleMatrix*       getEntryValue(size_t entryIndex) const;
        double              getEntryBranchLength(size_t entryIndex) const;
        
                            // iteration over dirty entries only
        size_t              getDirtyCount(void) const { return dirtyCount; }
        size_t              getDirtyEntryIndex(size_t i) const { return dirtyList[i]; }
        double              getDirtyBranchLength(size_t i) const;
        DoubleMatrix*       getDirtyMatrix(size_t i) const;
        
                            // structure modification
        void                clear(void);
        bool                erase(uint64_t key);
        bool                erase(double branchLength);
        void                reserve(size_t capacity);
        
        void                print(void) const;
        
    private:
        void                addToEntryList(size_t slot);
        void                backupIfNeeded(size_t entryIndex);
        void                growTable(void);
        size_t              hash(uint64_t key) const;
        void                removeFromEntryList(size_t slot);
        
                            // hash table
        uint64_t*           keys;
        DoubleMatrix**      values;
        DoubleMatrix**      backupValues;   // parallel array for MCMC backup
        bool*               occupied;
        size_t              tableCapacity;
        size_t              numEntries;
        
                            // dense entry list for O(n) iteration
        size_t*             entrySlots;
        size_t*             slotToEntry;
        
                            // external matrix pool (not owned)
        MatrixPool*         matrixPool;
        size_t              numStates;
        
                            // dirty tracking
        bool*               dirty;          // dirty[entryIndex] = true if modified
        size_t*             dirtyList;      // indices of dirty entries
        size_t              dirtyCount;
};

#endif
