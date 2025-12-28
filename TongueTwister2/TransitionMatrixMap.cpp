#include <cstring>
#include <iomanip>
#include <iostream>
#include "MatrixPool.hpp"
#include "TransitionMatrixMap.hpp"



TransitionMatrixMap::TransitionMatrixMap(size_t ns, size_t initialCapacity, MatrixPool* pool) : 
    keys(nullptr),
    values(nullptr),
    backupValues(nullptr),
    occupied(nullptr),
    tableCapacity(0),
    numEntries(0),
    entrySlots(nullptr),
    slotToEntry(nullptr),
    matrixPool(pool),
    numStates(ns),
    dirty(nullptr),
    dirtyList(nullptr),
    dirtyCount(0) {

    reserve(initialCapacity);
}

TransitionMatrixMap::~TransitionMatrixMap(void) {

    // release all matrices back to the pool
    for (size_t i = 0; i < numEntries; i++)
        {
        size_t slot = entrySlots[i];
        if (values[slot] != nullptr)
            matrixPool->release(values[slot]);
        if (backupValues[slot] != nullptr)
            matrixPool->release(backupValues[slot]);
        }
    
    // delete the pointer arrays
    delete[] keys;
    delete[] values;
    delete[] backupValues;
    delete[] occupied;
    delete[] entrySlots;
    delete[] slotToEntry;
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
    DoubleMatrix** oldBackupValues = backupValues;
    bool* oldOccupied = occupied;
    size_t* oldEntrySlots = entrySlots;
    bool* oldDirty = dirty;
    size_t* oldDirtyList = dirtyList;
    size_t oldNumEntries = numEntries;
    size_t oldDirtyCount = dirtyCount;
    
    keys = new uint64_t[newCapacity];
    values = new DoubleMatrix*[newCapacity];
    backupValues = new DoubleMatrix*[newCapacity];
    occupied = new bool[newCapacity];
    entrySlots = new size_t[newCapacity];
    slotToEntry = new size_t[newCapacity];
    dirty = new bool[newCapacity];
    dirtyList = new size_t[newCapacity];
    tableCapacity = newCapacity;
    
    std::memset(occupied, 0, newCapacity * sizeof(bool));
    std::memset(dirty, 0, newCapacity * sizeof(bool));
    std::memset(values, 0, newCapacity * sizeof(DoubleMatrix*));
    std::memset(backupValues, 0, newCapacity * sizeof(DoubleMatrix*));
    for (size_t i = 0; i < newCapacity; i++)
        slotToEntry[i] = SIZE_MAX;
    numEntries = 0;
    dirtyCount = 0;
    
    if (oldKeys)
        {
        // build mapping from old entry index to new entry index
        size_t* oldToNew = new size_t[oldNumEntries];
        
        for (size_t i = 0; i < oldNumEntries; i++)
            {
            size_t oldSlot = oldEntrySlots[i];
            uint64_t key = oldKeys[oldSlot];
            DoubleMatrix* value = oldValues[oldSlot];
            DoubleMatrix* backupValue = oldBackupValues[oldSlot];
            
            size_t idx = hash(key);
            while (occupied[idx])
                idx = (idx + 1) & (tableCapacity - 1);
            
            keys[idx] = key;
            values[idx] = value;
            backupValues[idx] = backupValue;
            occupied[idx] = true;
            
            oldToNew[i] = numEntries;
            addToEntryList(idx);
            }
        
        // rebuild dirty list with new indices
        for (size_t i = 0; i < oldDirtyCount; i++)
            {
            size_t oldEntryIdx = oldDirtyList[i];
            size_t newEntryIdx = oldToNew[oldEntryIdx];
            dirty[newEntryIdx] = true;
            dirtyList[dirtyCount++] = newEntryIdx;
            }
        
        delete [] oldToNew;
        delete [] oldKeys;
        delete [] oldValues;
        delete [] oldBackupValues;
        delete [] oldOccupied;
        delete [] oldEntrySlots;
        delete [] oldDirty;
        delete [] oldDirtyList;
        }
}

void TransitionMatrixMap::addToEntryList(size_t slot) {

    entrySlots[numEntries] = slot;
    slotToEntry[slot] = numEntries;
    numEntries++;
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
        
        // release matrices back to pool
        matrixPool->release(values[slot]);
        if (backupValues[slot] != nullptr)
            {
            matrixPool->release(backupValues[slot]);
            backupValues[slot] = nullptr;
            }
        values[slot] = nullptr;
        
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

    // check if already exists
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
    
    // need to create new entry
    if (numEntries >= tableCapacity / 2)
        growTable();
    
    size_t idx = hash(key);
    while (occupied[idx])
        idx = (idx + 1) & (tableCapacity - 1);
    
    // acquire matrix from pool and zero it
    DoubleMatrix* matrix = matrixPool->acquire();
    matrix->setZero();
    
    keys[idx] = key;
    values[idx] = matrix;
    backupValues[idx] = nullptr;  // backup allocated lazily when needed
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
    
    // allocate backup matrix if not already allocated
    if (backupValues[slot] == nullptr)
        backupValues[slot] = matrixPool->acquire();
    
    // copy current values to backup
    backupValues[slot]->copy(*values[slot]);
    
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

    // accept: clear dirty flags (backup becomes stale, that's fine)
    for (size_t i = 0; i < dirtyCount; i++)
        dirty[dirtyList[i]] = false;
    dirtyCount = 0;
}

void TransitionMatrixMap::restore(void) {

    // reject: restore dirty entries from backup
    for (size_t i = 0; i < dirtyCount; i++)
        {
        size_t entryIdx = dirtyList[i];
        size_t slot = entrySlots[entryIdx];
        
        // copy backup back to working matrix
        values[slot]->copy(*backupValues[slot]);
        
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
            // release matrices back to pool
            matrixPool->release(values[idx]);
            if (backupValues[idx] != nullptr)
                {
                matrixPool->release(backupValues[idx]);
                backupValues[idx] = nullptr;
                }
            values[idx] = nullptr;
            
            occupied[idx] = false;
            removeFromEntryList(idx);
            
            // Rehash subsequent entries
            idx = (idx + 1) & (tableCapacity - 1);
            while (occupied[idx])
                {
                uint64_t rehashKey = keys[idx];
                DoubleMatrix* rehashValue = values[idx];
                DoubleMatrix* rehashBackup = backupValues[idx];
                
                occupied[idx] = false;
                removeFromEntryList(idx);
                
                size_t newIdx = hash(rehashKey);
                while (occupied[newIdx])
                    newIdx = (newIdx + 1) & (tableCapacity - 1);
                
                keys[newIdx] = rehashKey;
                values[newIdx] = rehashValue;
                backupValues[newIdx] = rehashBackup;
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

void TransitionMatrixMap::print(void) const {

    std::cout << "TransitionMatrixMap:" << std::endl;
    std::cout << "  Entries: " << numEntries << ", capacity " << tableCapacity << std::endl;
    std::cout << "  Pool: " << matrixPool->numInUse() << " in use, " 
              << matrixPool->numFree() << " free, "
              << matrixPool->totalAllocated() << " total" << std::endl;
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
