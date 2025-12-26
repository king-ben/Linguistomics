#ifndef BitMask_H
#define BitMask_H

#include <set>
#include <string>
#include "BitSet.hpp"



namespace BitMask {

    bool        isSet(unsigned x, size_t n);
    bool        isSet(uint64_t x, size_t n);
    size_t      numSet(const unsigned& bs, size_t n);
    void        setBit(unsigned& bs, size_t n);
    void        setBit(uint64_t& bs, size_t n);
    std::string unsignedString(unsigned bs, size_t n);
    std::string unsignedString(std::set<unsigned>& bs, size_t n);
    std::string unsignedString(uint64_t bs, size_t n);
    unsigned    unsignedBitSet(BitSet& b);
    unsigned    unsignedBitSet(size_t n);
    unsigned    unsignedBitSet(std::set<unsigned>& b, size_t n);
}

#endif
