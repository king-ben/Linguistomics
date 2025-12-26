#include <cstring>
#include <iomanip>
#include <iostream>
#include "TransitionMatrixManager.hpp"
#include "TransitionMatrixMap.hpp"



TransitionMatrixMap::TransitionMatrixMap(size_t numStates, size_t initialCapacity) : 
    keys(nullptr),
    values(nullptr),
    occupied(nullptr),
    tableCapacity(0),
    numEntries(0),
    entrySlots(nullptr),
    slotToEntry(nullptr),
    pool(nullptr),
    backup(nullptr),
    freeList(nullptr),
    poolCapacity(0),
    freeCount(0),
    numStates(numStates),
    dirty(nullptr),
    dirtyList(nullptr),
    dirtyCount(0) {

    // ensure the matrix manager is initialized with the correct number of states
    TransitionMatrixManager& manager = TransitionMatrixManager::getInstance();
    manager.initialize(numStates, initialCapacity * 2);
    
    reserve(initialCapacity);
    reservePool(initialCapacity);
}

TransitionMatrixMap::~TransitionMatrixMap(void) {

    // Return all matrices to the manager instead of deleting them
    TransitionMatrixManager& manager = TransitionMatrixManager::getInstance();
    
    // Return matrices that are in the pool arrays (both used and free)
    for (size_t i = 0; i < poolCapacity; i++)
        {
        if (pool[i] != nullptr)
            manager.returnMatrix(pool[i]);
        if (backup[i] != nullptr)
            manager.returnMatrix(backup[i]);
        }
    
    // Delete the pointer arrays themselves
    delete[] keys;
    delete[] values;
    delete[] occupied;
    delete[] entrySlots;
    delete[] slotToEntry;
    delete[] pool;
    delete[] backup;
    delete[] freeList;
    delete[] dirty;
    delete[] dirtyList;
}

void TransitionMatrixMap::reserve(size_t requestedCapacity) {

    size_t newCapacity = 16;
    while (newCapacity < requestedCapacity)
        newCapacity <<= 1;
    
    if (newCapacity <= tableCapacity)
        return;
    
    uint64_t* oldKeys = keys;
    DoubleMatrix** oldValues = values;
    bool* oldOccupied = occupied;
    size_t* oldEntrySlots = entrySlots;
    bool* oldDirty = dirty;
    size_t* oldDirtyList = dirtyList;
    size_t oldNumEntries = numEntries;
    size_t oldDirtyCount = dirtyCount;
    
    keys = new uint64_t[newCapacity];
    values = new DoubleMatrix*[newCapacity];
    occupied = new bool[newCapacity];
    entrySlots = new size_t[newCapacity];
    slotToEntry = new size_t[newCapacity];
    dirty = new bool[newCapacity];
    dirtyList = new size_t[newCapacity];
    tableCapacity = newCapacity;
    
    std::memset(occupied, 0, newCapacity * sizeof(bool));
    std::memset(dirty, 0, newCapacity * sizeof(bool));
    for (size_t i = 0; i < newCapacity; i++)
        slotToEntry[i] = SIZE_MAX;
    numEntries = 0;
    dirtyCount = 0;
    
    if (oldKeys)
        {
        // Build mapping from old entry index to new entry index
        size_t* oldToNew = new size_t[oldNumEntries];
        
        for (size_t i = 0; i < oldNumEntries; i++)
            {
            size_t oldSlot = oldEntrySlots[i];
            uint64_t key = oldKeys[oldSlot];
            DoubleMatrix* value = oldValues[oldSlot];
            
            size_t idx = hash(key);
            while (occupied[idx])
                idx = (idx + 1) & (tableCapacity - 1);
            
            keys[idx] = key;
            values[idx] = value;
            occupied[idx] = true;
            
            oldToNew[i] = numEntries;
            addToEntryList(idx);
            }
        
        // Rebuild dirty list with new indices
        for (size_t i = 0; i < oldDirtyCount; i++)
            {
            size_t oldEntryIdx = oldDirtyList[i];
            size_t newEntryIdx = oldToNew[oldEntryIdx];
            dirty[newEntryIdx] = true;
            dirtyList[dirtyCount++] = newEntryIdx;
            }
        
        delete[] oldToNew;
        delete[] oldKeys;
        delete[] oldValues;
        delete[] oldOccupied;
        delete[] oldEntrySlots;
        delete[] oldDirty;
        delete[] oldDirtyList;
        }
}

void TransitionMatrixMap::reservePool(size_t capacity) {

    if (capacity <= poolCapacity)
        return;
    growPool(capacity);
}

void TransitionMatrixMap::growPool(size_t newCapacity) {

    TransitionMatrixManager& manager = TransitionMatrixManager::getInstance();
    
    // allocate new pointer arrays
    DoubleMatrix** newPool = new DoubleMatrix*[newCapacity];
    DoubleMatrix** newBackup = new DoubleMatrix*[newCapacity];
    size_t* newFreeList = new size_t[newCapacity];
    
    // initialize all pointers to nullptr
    std::memset(newPool, 0, newCapacity * sizeof(DoubleMatrix*));
    std::memset(newBackup, 0, newCapacity * sizeof(DoubleMatrix*));
    
    // copy existing pointers (matrices stay the same, just copy the pointers)
    if (pool)
        {
        for (size_t i = 0; i < poolCapacity; i++)
            {
            newPool[i] = pool[i];
            newBackup[i] = backup[i];
            }
        std::memcpy(newFreeList, freeList, freeCount * sizeof(size_t));
        }
    
    // allocate new matrices from the manager for the new slots
    for (size_t i = poolCapacity; i < newCapacity; i++)
        {
        newPool[i] = manager.getMatrix();
        newBackup[i] = manager.getMatrix();
        newFreeList[freeCount++] = i;
        }
    
    // delete old pointer arrays (not the matrices themselves!)
    delete [] pool;
    delete [] backup;
    delete [] freeList;
    
    pool = newPool;
    backup = newBackup;
    freeList = newFreeList;
    poolCapacity = newCapacity;
}

DoubleMatrix* TransitionMatrixMap::allocateFromPool(void) {

    if (freeCount == 0)
        growPool(poolCapacity + poolCapacity / 2 + 16);
    
    size_t idx = freeList[--freeCount];
    pool[idx]->setZero();
    return pool[idx];
}

void TransitionMatrixMap::returnToPool(DoubleMatrix* matrix) {

    size_t idx = matrixIndex(matrix);
    freeList[freeCount++] = idx;
}

size_t TransitionMatrixMap::matrixIndex(DoubleMatrix* matrix) const {

    // search for the matrix pointer in the pool array
    for (size_t i = 0; i < poolCapacity; i++)
        {
        if (pool[i] == matrix)
            return i;
        }
    return SIZE_MAX;  // not found (error condition)
}

void TransitionMatrixMap::addToEntryList(size_t slot) {

    entrySlots[numEntries] = slot;
    slotToEntry[slot] = numEntries;
    ++numEntries;
}

void TransitionMatrixMap::removeFromEntryList(size_t slot) {

    size_t entryIdx = slotToEntry[slot];
    
    // if this entry was dirty, remove from dirty list
    if (dirty[entryIdx])
        {
        for (size_t i = 0; i < dirtyCount; i++)
            {
            if (dirtyList[i] == entryIdx)
                {
                dirtyList[i] = dirtyList[--dirtyCount];
                break;
                }
            }
        dirty[entryIdx] = false;
        }
    
    // swap with last entry if not already last
    if (entryIdx != numEntries - 1)
        {
        size_t lastSlot = entrySlots[numEntries - 1];
        entrySlots[entryIdx] = lastSlot;
        slotToEntry[lastSlot] = entryIdx;
        
        // update dirty tracking for moved entry
        if (dirty[numEntries - 1])
            {
            for (size_t i = 0; i < dirtyCount; i++)
                {
                if (dirtyList[i] == numEntries - 1)
                    {
                    dirtyList[i] = entryIdx;
                    break;
                    }
                }
            dirty[entryIdx] = true;
            dirty[numEntries - 1] = false;
            }
        }
    
    slotToEntry[slot] = SIZE_MAX;
    --numEntries;
}

void TransitionMatrixMap::clear(void) {

    for (size_t i = 0; i < numEntries; i++)
        {
        size_t slot = entrySlots[i];
        returnToPool(values[slot]);
        occupied[slot] = false;
        slotToEntry[slot] = SIZE_MAX;
        dirty[i] = false;
        }
    numEntries = 0;
    dirtyCount = 0;
}

size_t TransitionMatrixMap::hash(uint64_t key) const {

    size_t h = 14695981039346656037ULL;
    h ^= key;
    h *= 1099511628211ULL;
    h ^= (key >> 32);
    h *= 1099511628211ULL;
    return h & (tableCapacity - 1);
}

DoubleMatrix* TransitionMatrixMap::find(uint64_t key) const {

    if (tableCapacity == 0)
        return nullptr;
    
    size_t idx = hash(key);
    size_t startIdx = idx;
    
    do {
        if (!occupied[idx])
            return nullptr;
        if (keys[idx] == key)
            return values[idx];
        idx = (idx + 1) & (tableCapacity - 1);
    } while (idx != startIdx);
    
    return nullptr;
}

DoubleMatrix* TransitionMatrixMap::find(double branchLength) const {

    return find(branchLengthToKey(branchLength));
}

DoubleMatrix* TransitionMatrixMap::getOrCreate(uint64_t key) {

    // Check if already exists
    if (tableCapacity > 0)
        {
        size_t idx = hash(key);
        size_t startIdx = idx;
        
        do {
            if (!occupied[idx])
                break;
            if (keys[idx] == key)
                return values[idx];
            idx = (idx + 1) & (tableCapacity - 1);
        } while (idx != startIdx);
        }
    
    // Need to create new entry
    if (numEntries >= tableCapacity / 2)
        growTable();
    
    size_t idx = hash(key);
    while (occupied[idx])
        idx = (idx + 1) & (tableCapacity - 1);
    
    keys[idx] = key;
    values[idx] = allocateFromPool();
    occupied[idx] = true;
    addToEntryList(idx);
    
    return values[idx];
}

DoubleMatrix* TransitionMatrixMap::getOrCreate(double branchLength) {

    return getOrCreate(branchLengthToKey(branchLength));
}

void TransitionMatrixMap::backupIfNeeded(size_t entryIndex) {

    if (dirty[entryIndex])
        return;
    
    size_t slot = entrySlots[entryIndex];
    size_t poolIdx = matrixIndex(values[slot]);
    backup[poolIdx]->copy(*pool[poolIdx]);
    
    dirty[entryIndex] = true;
    dirtyList[dirtyCount++] = entryIndex;
}

DoubleMatrix* TransitionMatrixMap::getForModification(uint64_t key) {

    if (tableCapacity == 0)
        return nullptr;
    
    size_t idx = hash(key);
    size_t startIdx = idx;
    
    do {
        if (!occupied[idx])
            return nullptr;
        if (keys[idx] == key)
            {
            size_t entryIdx = slotToEntry[idx];
            backupIfNeeded(entryIdx);
            return values[idx];
            }
        idx = (idx + 1) & (tableCapacity - 1);
    } while (idx != startIdx);
    
    return nullptr;
}

DoubleMatrix* TransitionMatrixMap::getForModification(double branchLength) {

    return getForModification(branchLengthToKey(branchLength));
}

void TransitionMatrixMap::markForUpdate(uint64_t key) {

    if (tableCapacity == 0)
        return;
    
    size_t idx = hash(key);
    size_t startIdx = idx;
    
    do {
        if (!occupied[idx])
            return;
        if (keys[idx] == key)
            {
            backupIfNeeded(slotToEntry[idx]);
            return;
            }
        idx = (idx + 1) & (tableCapacity - 1);
    } while (idx != startIdx);
}

void TransitionMatrixMap::markForUpdate(double branchLength) {

    markForUpdate(branchLengthToKey(branchLength));
}

void TransitionMatrixMap::markAllForUpdate(void) {

    for (size_t i = 0; i < numEntries; i++)
        backupIfNeeded(i);
}

void TransitionMatrixMap::keep(void) {

    // Accept: clear dirty flags (backup becomes stale, that's fine)
    for (size_t i = 0; i < dirtyCount; i++)
        dirty[dirtyList[i]] = false;
    dirtyCount = 0;
}

void TransitionMatrixMap::restore(void) {

    // Reject: restore dirty entries from backup
    for (size_t i = 0; i < dirtyCount; i++)
        {
        size_t entryIdx = dirtyList[i];
        size_t slot = entrySlots[entryIdx];
        size_t poolIdx = matrixIndex(values[slot]);
        pool[poolIdx]->copy(*backup[poolIdx]);
        dirty[entryIdx] = false;
        }
    dirtyCount = 0;
}

double TransitionMatrixMap::getDirtyBranchLength(size_t i) const {

    size_t entryIdx = dirtyList[i];
    return keyToBranchLength(keys[entrySlots[entryIdx]]);
}

DoubleMatrix* TransitionMatrixMap::getDirtyMatrix(size_t i) const {

    size_t entryIdx = dirtyList[i];
    return values[entrySlots[entryIdx]];
}

bool TransitionMatrixMap::erase(uint64_t key) {

    if (tableCapacity == 0)
        return false;
    
    size_t idx = hash(key);
    size_t startIdx = idx;
    
    while (occupied[idx])
        {
        if (keys[idx] == key)
            {
            returnToPool(values[idx]);
            occupied[idx] = false;
            removeFromEntryList(idx);
            
            // Rehash subsequent entries
            idx = (idx + 1) & (tableCapacity - 1);
            while (occupied[idx])
                {
                uint64_t rehashKey = keys[idx];
                DoubleMatrix* rehashValue = values[idx];
                
                occupied[idx] = false;
                removeFromEntryList(idx);
                
                size_t newIdx = hash(rehashKey);
                while (occupied[newIdx])
                    newIdx = (newIdx + 1) & (tableCapacity - 1);
                
                keys[newIdx] = rehashKey;
                values[newIdx] = rehashValue;
                occupied[newIdx] = true;
                addToEntryList(newIdx);
                
                idx = (idx + 1) & (tableCapacity - 1);
                }
            
            return true;
            }
        idx = (idx + 1) & (tableCapacity - 1);
        if (idx == startIdx)
            break;
        }
    
    return false;
}

bool TransitionMatrixMap::erase(double branchLength) {

    return erase(branchLengthToKey(branchLength));
}

uint64_t TransitionMatrixMap::getEntryKey(size_t entryIndex) const {

    return keys[entrySlots[entryIndex]];
}

DoubleMatrix* TransitionMatrixMap::getEntryValue(size_t entryIndex) const {

    return values[entrySlots[entryIndex]];
}

double TransitionMatrixMap::getEntryBranchLength(size_t entryIndex) const {

    return keyToBranchLength(keys[entrySlots[entryIndex]]);
}

void TransitionMatrixMap::growTable(void) {

    size_t newCapacity = (tableCapacity == 0) ? 16 : tableCapacity * 2;
    reserve(newCapacity);
}

void TransitionMatrixMap::shrinkPoolIfNeeded(size_t targetUtilization) {

    TransitionMatrixManager& manager = TransitionMatrixManager::getInstance();
    
    // Calculate target pool size based on current usage
    size_t usedCount = poolCapacity - freeCount;
    size_t targetCapacity = usedCount * targetUtilization;
    
    // Minimum pool size to avoid thrashing
    if (targetCapacity < 32)
        targetCapacity = 32;
    
    // Only shrink if we can reclaim at least 25% of the pool
    if (poolCapacity <= targetCapacity || poolCapacity < targetCapacity + targetCapacity / 4)
        return;
    
    // Build a mapping from old pool indices to new pool indices
    // Only the indices that are currently in use (not in freeList) will be kept
    size_t* oldToNew = new size_t[poolCapacity];
    bool* isUsed = new bool[poolCapacity];
    
    std::memset(isUsed, true, poolCapacity * sizeof(bool));
    
    // Mark free slots as unused
    for (size_t i = 0; i < freeCount; i++)
        isUsed[freeList[i]] = false;
    
    // Return unused matrices to the manager
    for (size_t i = 0; i < poolCapacity; i++)
        {
        if (!isUsed[i])
            {
            manager.returnMatrix(pool[i]);
            manager.returnMatrix(backup[i]);
            pool[i] = nullptr;
            backup[i] = nullptr;
            }
        }
    
    // Build the mapping - only used slots get new indices
    size_t newIdx = 0;
    for (size_t i = 0; i < poolCapacity; i++)
        {
        if (isUsed[i])
            oldToNew[i] = newIdx++;
        else
            oldToNew[i] = SIZE_MAX;
        }
    
    size_t newCapacity = newIdx * targetUtilization;
    if (newCapacity < 32)
        newCapacity = 32;
    
    // Allocate new pointer arrays
    DoubleMatrix** newPool = new DoubleMatrix*[newCapacity];
    DoubleMatrix** newBackup = new DoubleMatrix*[newCapacity];
    size_t* newFreeList = new size_t[newCapacity];
    
    // Initialize to nullptr
    std::memset(newPool, 0, newCapacity * sizeof(DoubleMatrix*));
    std::memset(newBackup, 0, newCapacity * sizeof(DoubleMatrix*));
    
    // Move used matrix pointers to new arrays (compacting them)
    for (size_t i = 0; i < poolCapacity; i++)
        {
        if (isUsed[i])
            {
            newPool[oldToNew[i]] = pool[i];
            newBackup[oldToNew[i]] = backup[i];
            }
        }
    
    // Allocate new matrices from manager for the extra slots
    for (size_t i = newIdx; i < newCapacity; i++)
        {
        newPool[i] = manager.getMatrix();
        newBackup[i] = manager.getMatrix();
        }
    
    // Update all values[] pointers to point to new pool locations
    for (size_t i = 0; i < numEntries; i++)
        {
        size_t slot = entrySlots[i];
        if (values[slot] != nullptr)
            {
            // Find the old pool index by searching
            size_t oldPoolIdx = SIZE_MAX;
            for (size_t j = 0; j < poolCapacity; j++)
                {
                if (pool[j] == values[slot])
                    {
                    oldPoolIdx = j;
                    break;
                    }
                }
            if (oldPoolIdx != SIZE_MAX)
                values[slot] = newPool[oldToNew[oldPoolIdx]];
            }
        }
    
    // Build new free list with the extra slots
    size_t newFreeCount = 0;
    for (size_t i = newIdx; i < newCapacity; i++)
        newFreeList[newFreeCount++] = i;
    
    // Clean up
    delete[] oldToNew;
    delete[] isUsed;
    delete[] pool;
    delete[] backup;
    delete[] freeList;
    
    pool = newPool;
    backup = newBackup;
    freeList = newFreeList;
    poolCapacity = newCapacity;
    freeCount = newFreeCount;
}

void TransitionMatrixMap::print(void) const {

    std::cout << "TransitionMatrixMap:" << std::endl;
    std::cout << "  Entries: " << numEntries << ", capacity " << tableCapacity << std::endl;
    std::cout << "  Pool: " << (poolCapacity - freeCount) << " used, " << poolCapacity << " total" << std::endl;
    std::cout << "  Dirty: " << dirtyCount << " entries" << std::endl;
    
    if (numEntries == 0)
        return;
    
    std::cout << "  Contents:" << std::endl;
    for (size_t i = 0; i < numEntries; i++)
        {
        double bl = keyToBranchLength(keys[entrySlots[i]]);
        std::cout << "    [" << i << "] bl=" << std::fixed << std::setprecision(6) << bl
                  << (dirty[i] ? " (dirty)" : "") << std::endl;
        }
}
