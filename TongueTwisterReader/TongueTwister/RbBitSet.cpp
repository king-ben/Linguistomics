#include <iostream>
#include <limits>
#include <sstream>
#include "Msg.hpp"
#include "RbBitSet.h"

int RbBitSet::bitsPerUint = sizeof(unsigned) * 8;
unsigned RbBitSet::maxUInt = std::numeric_limits<unsigned>::max() ;
unsigned RbBitSet::highBit = (unsigned(1) << (RbBitSet::bitsPerUint-1));



/** RbBitSet constructor */
RbBitSet::RbBitSet(void) {

    numElements = 0;
    numUints = 0;
    mask = 0;
    v.clear();
}

/** RbBitSet copy constructor
* @param b reference to the RbBitSet to be copied */
RbBitSet::RbBitSet(const RbBitSet& b) {

    clone(b);
}

/** RbBitSet constructor
* @param n the number of elements in the bitset */
RbBitSet::RbBitSet(size_t n) {

    // calculate size of vector and allocate memory for bitfield
    numElements = n;
    numUints = (numElements / RbBitSet::bitsPerUint);
    if (numElements % RbBitSet::bitsPerUint != 0)
        numUints++;
    v.resize(numUints);
    for (size_t i=0; i<numUints; i++)
        v[i] = 0;
    
    // initialize the mask
    if (numElements % RbBitSet::bitsPerUint != 0)
        {
        mask = (RbBitSet::highBit >> (numElements % RbBitSet::bitsPerUint));
        mask--;
        mask ^= RbBitSet::maxUInt;
        }
    else
        {
        mask = RbBitSet::maxUInt;
        }
}

RbBitSet::RbBitSet(std::vector<bool>& bs) {

    numElements = bs.size();
    numUints = (numElements / RbBitSet::bitsPerUint);
    if (numElements % RbBitSet::bitsPerUint != 0)
        numUints++;
    v.resize(numUints);
    for (size_t i=0; i<numUints; i++)
        v[i] = 0;
    
    // initialize the mask
    if (numElements % RbBitSet::bitsPerUint != 0)
        {
        mask = (RbBitSet::highBit >> (numElements % RbBitSet::bitsPerUint));
        mask--;
        mask ^= RbBitSet::maxUInt;
        }
    else
        {
        mask = RbBitSet::maxUInt;
        }

    // set all of the bits
    for (size_t i=0; i<numElements; i++)
        {
        if (bs[i] == true)
            set(i);
        else
            unset(i);
        }
}

RbBitSet::RbBitSet(size_t n, bool def) {

    // calculate size of vector and allocate memory for bitfield
    numElements = n;
    numUints = (numElements / RbBitSet::bitsPerUint);
    if (numElements % RbBitSet::bitsPerUint != 0)
        numUints++;
    v.resize(numUints);
    for (size_t i=0; i<numUints; i++)
        v[i] = 0;
    
    // initialize the mask
    if (numElements % RbBitSet::bitsPerUint != 0)
        {
        mask = (RbBitSet::highBit >> (numElements % RbBitSet::bitsPerUint));
        mask--;
        mask ^= RbBitSet::maxUInt;
        }
    else
        {
        mask = RbBitSet::maxUInt;
        }

    // set all of the bits
    for (size_t i=0; i<numElements; i++)
        {
        if (def == true)
            set(i);
        else
            unset(i);
        }
}

/** RbBitSet subscript operator
* @param idx the bit to return */
bool RbBitSet::operator[](size_t idx) {

    if (idx >= numElements)
        Msg::error("Subscript index out-of-range");
    return (( v[idx/RbBitSet::bitsPerUint] & (RbBitSet::highBit >> (idx % RbBitSet::bitsPerUint))) != 0);
}

/** RbBitSet subscript operator
* @param idx the bit to return */
bool RbBitSet::operator[](size_t idx) const {

    if (idx >= numElements)
        Msg::error("Subscript index out-of-range");
    return (( v[idx/RbBitSet::bitsPerUint] & (RbBitSet::highBit >> (idx % RbBitSet::bitsPerUint))) != 0);
}

/** RbBitSet assignment operator
* @param b the bitset to equate this one to */
RbBitSet& RbBitSet::operator=(const RbBitSet& b) {

    if (this != &b)
        clone(b);
    return *this;
}

/** RbBitSet equality operator
* @param b the bitset to compare to this one */
bool RbBitSet::operator==(const RbBitSet& b) const {

    if (numElements != b.numElements)
        return false;
    for (int i=0; i<numUints; i++)
        if (v[i] != b.v[i])
            return false;
    return true;
}

/** RbBitSet not operator
* @param b the bitset whose bits are to be flipped */
RbBitSet operator!(const RbBitSet& b) {

    RbBitSet newB(b);
    newB.flip();
    return newB;
}

/** RbBitSet inequality operator
* @param lhs the left hand side bitset
* @param rhs the right hand side bitset */
bool operator!=(const RbBitSet& lhs, const RbBitSet& rhs) {

    return !(lhs == rhs);
}

/** RbBitSet less than or equal to operator
* @param b the bitset to compare to this one */
bool operator<=(const RbBitSet& a, const RbBitSet& b) {

    if (a > b)
        return false;
    else
        return true;
}

/** RbBitSet less than operator
* @param lhs the left hand side bitset
* @param rhs the right hand side bitset */
bool operator<(const RbBitSet& lhs, const RbBitSet& rhs) {

    if (lhs > rhs || lhs == rhs)
        return false;
    else
        return true;
}

/** RbBitSet greater than or equal to operator
* @param lhs the left hand side bitset
* @param rhs the right hand side bitset */
bool operator>=(const RbBitSet& lhs, const RbBitSet& rhs) {

    if (lhs > rhs || lhs == rhs)
        return true;
    else
        return false;
    
}

/** RbBitSet greater than operator
* @param b the right hand side bitset */
bool RbBitSet::operator>(const RbBitSet& b) const {

    size_t x = (numUints < b.numUints) ? numUints : b.numUints;
    for (size_t i=0; i<x; i++)
        {
        if (v[i] > b.v[i])
            return true;
        else if (v[i] < b.v[i])
            return false;
        }
    return false;
}

/** RbBitSet and operator
* @param a the right hand side bitset */
RbBitSet RbBitSet::operator&(const RbBitSet& a) const {

    if (a.numElements != numElements )
        return RbBitSet();
    else
        {
        RbBitSet b(numElements);
        for (size_t i=0; i<numUints; i++)
            b.v[i] = v[i] & a.v[i];
        return b;
        }
}

/** RbBitSet or operator
* @param a the right hand side bitset */
RbBitSet RbBitSet::operator|(const RbBitSet& a) const {

    if (a.numElements != numElements )
        return RbBitSet();
    else
        {
        RbBitSet b(numElements);
        for (size_t i=0; i<numUints; i++)
            b.v[i] = v[i] | a.v[i];
        return b;
        }
}

/** RbBitSet xor operator
* @param a the right hand side bitset */
RbBitSet RbBitSet::operator^(const RbBitSet& a) const {

    if (a.numElements != numElements )
        return RbBitSet();
    else
        {
        RbBitSet b(numElements);
        for (size_t i=0; i<numUints; i++)
            b.v[i] = v[i] ^ a.v[i];
        return b;
        }
}

/** Unary not */
RbBitSet& RbBitSet::operator~() {

    flip();
    return *this;
}

/** RbBitSet stream operator
* @param s the stream
* @param b the right hand side bitset */
std::ostream& operator<<(std::ostream& s, const RbBitSet& b) {

    s << "[" << b.bitString() << "]";
    return s;
}

/** RbBitSet and equals  operator
* @param lhs the left hand side bitset
* @param rhs the right hand side bitset */
RbBitSet& operator&=(RbBitSet& lhs, const RbBitSet& rhs) {

    if (lhs.size() == rhs.size())
        {
        for (size_t i=0; i<lhs.getNumUnints(); i++)
            {
            unsigned lvi = lhs.getUintElement(i);
            unsigned rvi = rhs.getUintElement(i);
            lhs.setUintElement(i, lvi &= rvi);
            }
        }
    return lhs;
}

/** RbBitSet or equals  operator
* @param lhs the left hand side bitset
* @param rhs the right hand side bitset */
RbBitSet& operator|=(RbBitSet& lhs, const RbBitSet& rhs) {

    if (lhs.size() == rhs.size())
        {
        for (size_t i=0; i<lhs.getNumUnints(); i++)
            {
            unsigned lvi = lhs.getUintElement(i);
            unsigned rvi = rhs.getUintElement(i);
            lhs.setUintElement(i, lvi |= rvi);
            }
        }
    return lhs;
}

/** RbBitSet xor equals  operator
* @param lhs the left hand side bitset
* @param rhs the right hand side bitset */
RbBitSet& operator^=(RbBitSet& lhs, const RbBitSet& rhs) {

    if (lhs.size() == rhs.size())
        {
        for (size_t i=0; i<lhs.getNumUnints(); i++)
            {
            unsigned lvi = lhs.getUintElement(i);
            unsigned rvi = rhs.getUintElement(i);
            lhs.setUintElement(i, lvi ^= rvi);
            }
        }
    return lhs;
}

/** Return a string with the bitset */
std::string RbBitSet::bitString(void) const {

    std::string str = "";
    for (size_t i=0; i<numElements; i++)
        {
        bool tf = (*this)[i];
        if (tf == true)
            str += "1";
        else
            str += "0";
        }
    return str;
}

/** Clone function
* @param b the RbBitSet to be copied */
void RbBitSet::clone(const RbBitSet& b) {

    numElements = b.numElements;
    numUints    = b.numUints;
    mask        = b.mask;
    v           = b.v;
}

/** Erase all bits */
void RbBitSet::clear(void) {

    numElements = 0;
    numUints = 0;
    mask = 0;
    v.clear();
}

/** Is the bitset of size zero */
bool RbBitSet::empty(void) const {

    return v.empty();
}

/** Flip the bit at the index
* @param idx the index of the bit to flip */
void RbBitSet::flip(size_t idx) {

    v[idx/RbBitSet::bitsPerUint] ^= ((RbBitSet::highBit >> (idx % RbBitSet::bitsPerUint)));
}

/** Flip all bits */
void RbBitSet::flip(void) {

    for (size_t i=0; i<numUints; i++)
        v[i] = v[i] ^ RbBitSet::maxUInt;
    v[numUints-1] &= mask;
}

/** The position of the first on bit
* @return The index of the first on bit */
size_t RbBitSet::getFirstSetBit(void) const {
    
    size_t i, j;
    for (i=0; i<numUints; i++)
        {
        if (v[i] != 0)
            break;
        }
    if (i == numUints)
        return (size_t)-1;
    for (j=0; j<RbBitSet::bitsPerUint; j++)
        {
        if ((v[i] & (RbBitSet::highBit >> j)) != 0)
            break;
        }
    return (i * RbBitSet::bitsPerUint + j);

}

/** The number of on bits
* @return The number of on bits */
size_t RbBitSet::getNumberSetBits(void) const {
    
    size_t count = 0;
    for (size_t i=0; i<numElements; i++)
        {
        if ( (v[i/RbBitSet::bitsPerUint] & (RbBitSet::highBit >> (i % RbBitSet::bitsPerUint))) != 0 )
            count++;
        }
    return count;
}

/** Insert a bit
* @param pos the position in the bitset to insert
* @param tf the sign of the bit */
void RbBitSet::insert(size_t pos, bool tf) {

    // check position of insertion
    if (pos > numElements)
        {
        std::cout << "Cannot insert bit at position " << pos << " because there are only " << numElements << " bits" << std::endl;
        return;
        }
    // resize, if necessary
    size_t n = numElements + 1;
    size_t nu = (n / RbBitSet::bitsPerUint);
    if (n % RbBitSet::bitsPerUint != 0)
        nu++;
    if (nu != numUints)
        v.push_back(0);
    
    numElements = n;
    numUints = nu;
    
    // initialize the mask
    if (numElements % RbBitSet::bitsPerUint != 0)
        {
        mask = (RbBitSet::highBit >> (numElements % RbBitSet::bitsPerUint));
        mask--;
        mask ^= RbBitSet::maxUInt;
        }
    else
        {
        mask = RbBitSet::maxUInt;
        }

    // set the bits
    for (size_t i=numElements-1; i>pos; i--)
        {
        if ( (*this)[i-1] == true )
            set(i);
        else
            unset(i);
        }
    if (tf == true)
        set(pos);
    else
        unset(pos);
}

/** The number of on bits
* @param idx the index of the bit
* @return  whether the bit is on */
bool RbBitSet::isSet(size_t idx) const {

    return ((v[idx/RbBitSet::bitsPerUint] & (RbBitSet::highBit >> (idx % RbBitSet::bitsPerUint))) != 0);
}

/** Print in a long format the bitset */
void RbBitSet::print(void) {

    std::string str = "";
    str += "RbBitSet\n";
    str += "   Size     = " + std::to_string(numElements) + "\n";
    str += "   Capacity = " + std::to_string(numUints*RbBitSet::bitsPerUint) + "\n   ";
    str += bitString();
    std::cout << str << std::endl;
}

/** Add a bit to the end of the set
* @param tf the value of the bit to add */
void RbBitSet::push_back(bool tf) {

    size_t n = numElements + 1;
    size_t nu = (n / RbBitSet::bitsPerUint);
    if (n % RbBitSet::bitsPerUint != 0)
        nu++;
    if (nu != numUints)
        v.push_back(0);
    
    numElements = n;
    numUints = nu;
    
    // initialize the mask
    if (numElements % RbBitSet::bitsPerUint != 0)
        {
        mask = (RbBitSet::highBit >> (numElements % RbBitSet::bitsPerUint));
        mask--;
        mask ^= RbBitSet::maxUInt;
        }
    else
        {
        mask = RbBitSet::maxUInt;
        }
    if (tf == true)
        set(numElements-1);
    else
        unset(numElements-1);
}

/** Add a bit to the front of the set
* @param tf the value of the bit to add */
void RbBitSet::push_front(bool tf) {

    size_t n = numElements + 1;
    size_t nu = (n / RbBitSet::bitsPerUint);
    if (n % RbBitSet::bitsPerUint != 0)
        nu++;
    if (nu != numUints)
        v.push_back(0);
    
    numElements = n;
    numUints = nu;
    
    // initialize the mask
    if (numElements % RbBitSet::bitsPerUint != 0)
        {
        mask = (RbBitSet::highBit >> (numElements % RbBitSet::bitsPerUint));
        mask--;
        mask ^= RbBitSet::maxUInt;
        }
    else
        {
        mask = RbBitSet::maxUInt;
        }

    for (size_t i=numElements-1; i>=1; i--)
        {
        if ( (*this)[i-1] == true )
            set(i);
        else
            unset(i);
        }
    if (tf == true)
        set(0);
    else
        unset(0);
}

void RbBitSet::resize(size_t n) {

    // size the vector
    numElements = n;
    numUints = (numElements / RbBitSet::bitsPerUint);
    if (numElements % RbBitSet::bitsPerUint != 0)
        numUints++;
    v.resize(numUints);
    for (size_t i=0; i<numUints; i++)
        v[i] = 0;
    
    // initialize the mask
    if (numElements % RbBitSet::bitsPerUint != 0)
        {
        mask = (RbBitSet::highBit >> (numElements % RbBitSet::bitsPerUint));
        mask--;
        mask ^= RbBitSet::maxUInt;
        }
    else
        {
        mask = RbBitSet::maxUInt;
        }
}

void RbBitSet::resize(size_t n, bool tf) {

    resize(n); // this sets all bits to false
    if (tf == true)
        {
        for (size_t i=0; i<numElements; i++)
            set(i);
        }
}

/** Set all bist to on  */
void RbBitSet::set(void) {

    for (size_t i=0; i<numElements; i++)
        v[i/RbBitSet::bitsPerUint] |= (RbBitSet::highBit >> (i % RbBitSet::bitsPerUint));
}

/** Set a bit to on
* @param idx the index of the bit to turn on */
void RbBitSet::set(size_t idx) {

    if (idx >= numElements)
        Msg::error("Subscript index out-of-range");
    v[idx/RbBitSet::bitsPerUint] |= (RbBitSet::highBit >> (idx % RbBitSet::bitsPerUint));
}

/** Set all bits to 0  */
void RbBitSet::unset(void) {

    for (size_t i=0; i<numElements; i++)
        v[i/RbBitSet::bitsPerUint] &= ((RbBitSet::highBit >> (i % RbBitSet::bitsPerUint))^RbBitSet::maxUInt);
}

/** Set a bit to 0
* @param idx the index of the bit to change */
void RbBitSet::unset(size_t idx) {

    if (idx >= numElements)
        Msg::error("Subscript index out-of-range");
    v[idx/RbBitSet::bitsPerUint] &= ((RbBitSet::highBit >> (idx % RbBitSet::bitsPerUint))^RbBitSet::maxUInt);
}

/** Check the functionality of the class */
void RbBitSet::test(void) {

    // first check that edge cases are OK
    RbBitSet x1(10);
    x1.set(3);
    std::cout << x1 << " should be: (0001000000)" << std::endl;
    ~x1;
    std::cout << x1 << " should be: (1110111111)" << std::endl;
    RbBitSet x2(32);
    x2.set(31);
    std::cout << x2 << " should be: (00000000000000000000000000000001)" << std::endl;
    RbBitSet x3(33);
    x3.set(32);
    std::cout << x3 << " should be: (000000000000000000000000000000001)" << std::endl;
    
    // check operators
    RbBitSet x4(10);
    x4.set(4);
    std::cout << (x1 == x4) << " should be 0" << std::endl;
    RbBitSet x5(33);
    x5.set(32);
    std::cout << (x3 == x5) << " should be 1" << std::endl;
    std::cout << (x1 > x4) << " should be 1" << std::endl;
    std::cout << (x1 | x4) << " should be (0001100000)" << std::endl;
    std::cout << (x1 & x4) << " should be (0000000000)" << std::endl;
    std::cout << (x1 ^ x4) << " should be (0001100000)" << std::endl;
    std::cout << x1[3] << x2[31] << x3[32] << x4[4] << x5[32] << " should be 11111" << std::endl;
    std::cout << x1[0] << x2[14] << x3[30] << x4[8] << x5[29] << " should be 00000" << std::endl;
    
    RbBitSet x6(1);
    for (int i=0; i<100; i++)
        {
        x6.push_back(true);
        std::cout << x6 << std::endl;
        }

    RbBitSet x7(1);
    for (int i=0; i<100; i++)
        {
        x7.push_front(true);
        std::cout << x7 << std::endl;
        }
        
    RbBitSet x8(32);
    x8.insert(1, true);
    std::cout << x8 << std::endl;
    x8.print();

    RbBitSet x9(32);
    x9.insert(33, true);
    std::cout << x9 << std::endl;
    x9.print();
}
