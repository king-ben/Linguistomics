/*!
*
* \brief RbBitSet class
*
* The RbBitSet class implements a dynamic array of bits. The bits are densely packed
* into an array of unsigned integers. This array is a STL vector of unsigned ints. The
* bitset can be dynamically resized by pushing bits to the front or back or bby arbitrarily
* inserting a bit. Bits can also be deleted.
*/

#ifndef RbBitSet_H
#define RbBitSet_H

#include <string>
#include <vector>

class RbBitSet {

    public:
                                RbBitSet(void);                                         //!< Default constructor (blank bitset)
                                RbBitSet(const RbBitSet& b);                            //!< Copy constructor
                                RbBitSet(size_t n);                                     //!< Constructor setting up an empty bitset of size n
                                RbBitSet(size_t n, bool def);                           //!< Constructor setting up an array of size n with all bits set to def
                                RbBitSet(std::vector<bool>& bs);
        bool                    operator[](size_t idx);                                 //!< Indexing of bits
        bool                    operator[](size_t idx) const;                           //!< Indexing of bits  (for const bitset)
        RbBitSet&               operator=(const RbBitSet& b);                           //!< Assignment operator (deep copy)
        bool                    operator==(const RbBitSet& b) const;                    //!< Equality operator
        bool                    operator>(const RbBitSet& b) const;                     //!< Larger than operator
        RbBitSet                operator&(const RbBitSet& b) const;                     //!< Operator &
        RbBitSet                operator|(const RbBitSet& b) const;                     //!< Operator |
        RbBitSet                operator^(const RbBitSet& b) const;                     //!< Operator ^
        RbBitSet&               operator~();                                            //!< Unary not
        friend std::ostream&    operator<<(std::ostream& s, const RbBitSet& b);         //!< Operator <<
        std::string             bitString(void) const;                                  //!< Get a string representation of the bitset
        void                    clear(void);                                            //!< Remove all elements, setting size to zero
        size_t                  getFirstSetBit(void) const;                             //!< Return the index of the first bit that is on (1)
        bool                    empty(void) const;                                      //!< Whether the bitset is empty (is of size zero)
        void                    flip(size_t idx);                                       //!< Flip bit at position idx
        void                    flip(void);                                             //!< Flip all of the bits from 0 to 1 and from 1 to 0
        size_t                  getNumberSetBits(void) const;                           //!< Return the number of bits that are on
        size_t                  getNumUnints(void) const { return numUints; }           //!< Return the number of unsigned ints used to store the bits in the set
        unsigned                getUintElement(size_t idx) const { return v[idx]; }     //!< Return the unsigned ints used to store the bits in the set at location idx
        void                    insert(size_t pos, bool tf);                            //!< Insert a bit at position pos
        bool                    isSet(size_t idx) const;                                //!< Is bit at position idx set (on)?
        void                    print(void);                                            //!< Print a long description of the bitset
        void                    push_back(bool tf);                                     //!< Add a bit to the end of the vector of bits
        void                    push_front(bool tf);                                    //!< Add a bit to the front of the vector of bits
        void                    resize(size_t n);                                       //!< Resize the vector , setting all bits to false
        void                    resize(size_t n, bool tf);                              //!< Resize the vector and set them all to tf
        void                    set(void);                                              //!< set all bits to on (1)
        void                    set(size_t idx);                                        //!< set bit at position idx to on (1)
        void                    setUintElement(size_t idx, unsigned x) { v[idx] = x; }  //!< Set the unsigned ints used to store the bits in the set at location idx
        size_t                  size(void) const { return numElements; }                //!< Return the number of bits in the set
        void                    test(void);                                             //!< Run a series of tests on the class
        void                    unset(void);                                            //!< set all bits to off (0)
        void                    unset(size_t idx);                                      //!< set bit at position idx to off (0)

    private:
        void                    clone(const RbBitSet& b);
        static int              bitsPerUint;                                            //!< The number of bits per unsigned integer
        static unsigned         maxUInt;                                                //!< The maximum value of an unsigned integer
        static unsigned         highBit;                                                //!< An unsigned integer with only the highest bit on
        size_t                  numElements;                                            //!< The number of bits
        size_t                  numUints;                                               //!< The number of unsigned integers needed to store the bits
        unsigned                mask;                                                   //!< A mask for the last unsigned integer in the array
        std::vector<unsigned>   v;                                                      //!< The array of Uints that store the bits
};

        RbBitSet                operator!(const RbBitSet& b);                           //!< Operator !
        RbBitSet&               operator&=(RbBitSet& lhs, const RbBitSet& rhs);         //!< Operator &=
        RbBitSet&               operator|=(RbBitSet& lhs, const RbBitSet& rhs);         //!< Operator |=
        RbBitSet&               operator^=(RbBitSet& lhs, const RbBitSet& rhs);         //!< Operator ^=
        bool                    operator!=(const RbBitSet& lhs, const RbBitSet& rhs);   //!< Inequality operator
        bool                    operator< (const RbBitSet& lhs, const RbBitSet& rhs);   //!< Operator <
        bool                    operator<=(const RbBitSet& lhs, const RbBitSet& rhs);   //!< Operator <=
        bool                    operator>=(const RbBitSet& lhs, const RbBitSet& rhs);   //!< Operator >=

struct CompBitSet {

    bool operator()(const RbBitSet* n1, const RbBitSet* n2) const {
//        is n1 bigger than n2 true or false
        
        if (n1->size() == n2->size())
            {
            for (size_t i=0; i<n1->size(); i++)
                {
                bool v1 = n1->isSet(i);
                bool v2 = n2->isSet(i);
                if (v1 != v2)
                    {
                    if (v1 == true && v2 == false)
                        return true;
                    else
                        return false;
                    }
                }
            }
        else if (n1->size() > n2->size())
            return true;
        else
            return false;
        return false;
        }
};

#endif
