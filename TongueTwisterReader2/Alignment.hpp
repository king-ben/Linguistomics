#ifndef Alignment_hpp
#define Alignment_hpp

#include <cstdint>
#include <functional>
#include <map>
#include <string>
#include <vector>
#include "json.hpp"
#include "Sequence.hpp"



class Alignment {

    public:
                                        Alignment(void) = delete;
                                        Alignment(nlohmann::json j);
        Sequence&                       operator[](size_t i) { return matrix[i]; }
        const Sequence&                 operator[](size_t i) const { return matrix[i]; }
        bool                            operator==(const Alignment& aln) const;
        int                             getCharCode(int tIdx, int cIdx) const { return matrix[tIdx][cIdx]; }
        int                             getNumChar(void) const { return numChar; }
        int                             getNumTaxa(void) const { return numTaxa; }
        std::vector<double>             getGapSpectrum(void) const;
        int                             lengthOfLongestName(void) const;
        std::map<double,double>         gapInfo(void) const;
        void                            print(void) const;
        void                            print(std::string h) const;
        size_t                          size(void) const { return matrix.size(); }
        nlohmann::json                  toFile(std::ostream& file) const;
    
    private:
        static int                      max_int8;
        int                             numTaxa;
        int                             numChar;
        std::vector<Sequence>           matrix;
};

// Hash functor for Alignment - works with std::unordered_map/set
struct HashAlignment {

    size_t operator()(const Alignment& aln) const {
        
        size_t hash = 0;
        for (size_t i = 0; i < aln.size(); i++)
            {
            const Sequence& seq = aln[i];
            for (size_t j = 0; j < seq.size(); j++)
                {
                // boost::hash_combine approach
                hash ^= std::hash<int8_t>{}(seq[j]) + 0x9e3779b9 + (hash << 6) + (hash >> 2);
                }
            }
        return hash;
    }

    size_t operator()(const Alignment* aln) const {
        
        if (aln == nullptr)
            return 0;
        return (*this)(*aln);
    }
};

// Equality functor for Alignment pointers
struct EqualAlignment {

    bool operator()(const Alignment* a1, const Alignment* a2) const {
        
        if (a1 == a2)
            return true;
        if (a1 == nullptr || a2 == nullptr)
            return false;
        return *a1 == *a2;
    }
};

// Less-than comparator for Alignment - works with std::map/set
struct CompAlignment {

    bool operator()(const Alignment& a1, const Alignment& a2) const {
        
        if (a1.size() != a2.size())
            return a1.size() > a2.size();
        
        if (a1.size() == 0)
            return false;
        
        if (a1[0].size() != a2[0].size())
            return a1[0].size() > a2[0].size();
        
        for (size_t i = 0; i < a1.size(); i++)
            {
            for (size_t j = 0; j < a1[i].size(); j++)
                {
                if (a1[i][j] != a2[i][j])
                    return a1[i][j] > a2[i][j];
                }
            }
        return false;
    }

    bool operator()(const Alignment* a1, const Alignment* a2) const {

        if (a1 == nullptr && a2 == nullptr)
            return false;
        if (a1 == nullptr)
            return false;
        if (a2 == nullptr)
            return true;
        return (*this)(*a1, *a2);
    }
};

#endif
