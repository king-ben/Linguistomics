#ifndef Partition_hpp
#define Partition_hpp

#include <set>
#include "json.hpp"
class Subset;


class Partition {

    public:
                            Partition(void) = delete;
                            Partition(nlohmann::json js);
                           ~Partition(void);
        void                addSubset(std::string s, std::vector<int> v);
        Subset*             findSubsetIndexed(int x);
        Subset*             findSubsetWithValue(int x);
        int                 indexOfSubsetWithValue(int x);
        int                 numSubsets(void) { return (int)subsets.size(); }
        void                print(void);
        int                 maxValue(void);
    
    private:
        std::set<Subset*>   subsets;
};

#endif
