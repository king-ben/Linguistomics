#include "Subset.hpp"



Subset::Subset(std::string s, std::vector<int> v) {

    label = s;
    for (int i=0; i<v.size(); i++)
        vals.insert(v[i]);
}

bool Subset::containsValue(int x) {

    std::set<int>::iterator it = vals.find(x);
    if (it == vals.end())
        return false;
    return true;
}
