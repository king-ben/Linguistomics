#ifndef Subset_hpp
#define Subset_hpp

#include <set>
#include <string>
#include <vector>



class Subset {

    public:
                        Subset(void) = delete;
                        Subset(std::string s, std::vector<int> v);
        bool            containsValue(int x);
        int             getIndex(void) { return index; }
        std::string     getLabel(void) { return label; }
        int             getNumElements(void) { return (int)vals.size(); }
        std::set<int>&  getValues(void) { return vals; }
        void            setIndex(int x) { index = x; }
    
    private:
        std::string     label;
        int             index;
        std::set<int>   vals;
};

#endif
