#include "BitMask.hpp"
#include "Msg.hpp"



bool BitMask::isSet(unsigned x, size_t n) {

    bool isOn = (x & (1u << n)) != 0;
    return isOn;
}

bool BitMask::isSet(uint64_t x, size_t n) {

    // use 1ULL for 64-bit shift. Using 1u (32-bit) with n >= 32 is undefined behavior.
    bool isOn = (x & (1ULL << n)) != 0;
    return isOn;
}

size_t BitMask::numSet(const unsigned& bs, size_t n) {

    size_t nSet = 0;
    for (size_t i=0; i<n; i++)
        {
        if (isSet(bs, i) == true)
            nSet++;
        }
    return nSet;
}

void BitMask::setBit(unsigned& bs, size_t n) {

    // use >= instead of > since bit positions are 0-indexed
    if (n >= sizeof(unsigned)*8)
        Msg::error("Out of bounds bit set");
    unsigned part = 1u << n;
    bs |= part;
}

void BitMask::setBit(uint64_t& bs, size_t n) {

    // use >= instead of > since bit positions are 0-indexed
    if (n >= sizeof(uint64_t)*8)
        Msg::error("Out of bounds bit set");
    // use 1ULL for 64-bit shift. Using 1u (32-bit) with n >= 32 is undefined behavior.
    uint64_t part = 1ULL << n;
    bs |= part;
}

std::string BitMask::unsignedString(unsigned bs, size_t n) {

    std::string str = "";
    for (size_t i=0; i<n; i++)
        {
        if (bs & (1u << i)) 
            str += "1";
        else
            str += "0";
        }
    return str;
}

std::string BitMask::unsignedString(std::set<unsigned>& bs, size_t n) {

    std::string str = "";
    for (unsigned x : bs)
        {
        for (size_t i=0; i<n; i++)
            {
            if (x & (1u << i)) 
                str += "1";
            else
                str += "0";
            }
        }
    return str;
}

std::string BitMask::unsignedString(uint64_t bs, size_t n) {

    std::string str = "";
    for (size_t i=0; i<n; i++)
        {
        // use 1ULL for 64-bit shift
        if (bs & (1ULL << i)) 
            str += "1";
        else
            str += "0";
        }
    return str;
}

unsigned BitMask::unsignedBitSet(BitSet& b) {

    if (b.size() > sizeof(unsigned)*8)
        Msg::error("Bitset overflow");

    unsigned x = 0;
    for (size_t i=0, n=b.size(); i<n; i++)
        {
        if (b.isSet(i))  // only set the bit if it's set in the BitSet
            {
            unsigned part = 1u << i;
            x |= part;
            }
        }
    return x;
}

unsigned BitMask::unsignedBitSet(size_t n) {

    if (n >= sizeof(unsigned)*8)
        Msg::error("Bitset overflow");

    unsigned x = 1u << n;
    return x;
}

unsigned BitMask::unsignedBitSet(std::set<unsigned>& b, size_t n) {

    if (b.size() * n > sizeof(unsigned)*8)
        Msg::error("Bitset overflow");
        
    unsigned x = 0;
    size_t k = 0;
    for (unsigned elem : b)
        {
        for (size_t i=0; i<n; i++, k++)
            {
            unsigned b1 = 1u << i;
            if ( (b1 & elem) != 0)
                {
                unsigned b2 = 1u << k;
                x |= b2;
                }
            }
        }
    return x;
}
