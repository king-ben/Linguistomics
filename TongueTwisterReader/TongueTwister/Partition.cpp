#include <iostream>
#include "Msg.hpp"
#include "Partition.hpp"
#include "Subset.hpp"



Partition::Partition(nlohmann::json js) {

    for (int i=0; i<js.size(); i++)
        {
        std::string pName = js[i]["Name"];
        std::vector<int> pVals = js[i]["Set"];
        addSubset(pName, pVals);
        }
        
    int rgfIdx = 0;
    for (int i=0; i<maxValue(); i++)
        {
        Subset* s = findSubsetWithValue(i);
        if (s != NULL)
            {
            if (s->getIndex() == -1)
                s->setIndex(++rgfIdx);
            }
        }
        
    //print();
}

Partition::~Partition(void) {

    for (Subset* s : subsets)
        delete s;
}

void Partition::addSubset(std::string s, std::vector<int> v) {

    Subset* ss = new Subset(s, v);
    ss->setIndex(-1);
    subsets.insert(ss);
}

Subset* Partition::findSubsetIndexed(int x) {

    for (Subset* s : subsets)
        {
        if (s->getIndex() == x)
            return s;
        }
    return NULL;
}

Subset* Partition::findSubsetWithValue(int x) {

    for (Subset* s : subsets)
        {
        if (s->containsValue(x) == true)
            return s;
        }
    return NULL;
}

int Partition::indexOfSubsetWithValue(int x) {

    Subset* s = findSubsetWithValue(x);
    if (s == NULL)
        Msg::error("Could not find subset with value " + std::to_string(x));
    return s->getIndex();
}

void Partition::print(void) {

    int len = 0;
    for (Subset* s : subsets)
        {
        std::string n = s->getLabel();
        if (n.length() > len)
            len = (int)n.length();
        }
        
    for (int i=1; i<=subsets.size(); i++)
        {
        Subset* s = findSubsetIndexed(i);
        std::set<int> vals = s->getValues();
        std::cout << s->getIndex() << " " << s->getLabel();
        for (int j=0; j<len-s->getLabel().length(); j++)
            std::cout << " ";
        std::cout << " -- ";
        for (int j : vals)
            std::cout << j << " ";
        std::cout << std::endl;
        }
}

int Partition::maxValue(void) {

    int m = -1;
    for (Subset* s : subsets)
        {
        std::set<int> vals = s->getValues();
        for (int i : vals)
            {
            if (i > m)
                m = i;
            }
        }

    return m;
}
