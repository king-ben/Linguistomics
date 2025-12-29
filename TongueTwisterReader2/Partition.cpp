#include <iomanip>
#include <iostream>
#include "Msg.hpp"
#include "Partition.hpp"
#include "Subset.hpp"



Partition::Partition(const Partition& p) {
            
    for (size_t n=0; n<p.subsets.size(); n++)
        this->subsets.push_back(new Subset(*p.subsets[n]));
}
                                
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

bool Partition::operator==(const Partition& rhs) const {

    if (this->subsets.size() != rhs.subsets.size())
        return false;
        
    for (Subset* lhsSubset : this->subsets)
        {
        int lhsIndex = lhsSubset->getIndex();
        std::string lhsLabel = lhsSubset->getLabel();
        std::set<int>& lhsVals = lhsSubset->getValues();
        
        Subset* rhsSubset = rhs.findSubsetIndexed(lhsIndex);
        if (rhsSubset == NULL)
            {
            return false;
            }
        else
            {
            int rhsIndex = rhsSubset->getIndex();
            std::string rhsLabel = rhsSubset->getLabel();
            std::set<int>& rhsVals = rhsSubset->getValues();
            
            if (lhsIndex != rhsIndex)
                return false;
            if (lhsLabel != rhsLabel)
                return false;
            if (lhsVals != rhsVals)
                return false;
            }
        }
    return true;
}

bool Partition::isEqualTo(const Partition& rhs) const {

    return (*this == rhs);
}

void Partition::addSubset(std::string s, std::vector<int> v) {

    Subset* ss = new Subset(s, v);
    ss->setIndex((int)subsets.size()+1);
    subsets.push_back(ss);
}

Subset* Partition::findSubsetIndexed(int x) {

    for (Subset* s : subsets)
        {
        if (s->getIndex() == x)
            return s;
        }
    return NULL;
}

Subset* Partition::findSubsetIndexed(int x) const {

    for (Subset* s : subsets)
        {
        if (s->getIndex() == x)
            return s;
        }
    return NULL;
}

Subset* Partition::findSubsetLabeled(std::string sName) {

    for (Subset* s : subsets)
        {
        if (s->getLabel() == sName)
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

int Partition::getNumElements(void) {

    int n = 0;
    for (Subset* s : subsets)
        n += s->getNumElements();
    return n;
}

int Partition::indexOfSubsetWithValue(int x) {

    Subset* s = findSubsetWithValue(x);
    if (s == NULL)
        Msg::error("Could not find subset with value " + std::to_string(x));
    return s->getIndex();
}

size_t Partition::longestLabelLength(void) {

    size_t len = 0;
    for (Subset* s : subsets)
        {
        std::string n = s->getLabel();
        if (n.length() > len)
            len = (int)n.length();
        }
    return len;
}

void Partition::print(void) {

    size_t len = longestLabelLength();
        
    std::cout << "   * Natural Class Partition:" << std::endl;
    for (int i=1; i<=subsets.size(); i++)
        {
        std::cout << "      ";
        Subset* s = findSubsetIndexed(i);
        std::set<int>& vals = s->getValues();
        std::cout << std::setw(2) << s->getIndex() << " " << s->getLabel();
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
        std::set<int>& vals = s->getValues();
        for (int i : vals)
            {
            if (i > m)
                m = i;
            }
        }

    return m;
}

nlohmann::json Partition::toJson(void) {

    nlohmann::json j = nlohmann::json::object();
    
    nlohmann::json jSubsets = nlohmann::json::array();
    for (int i=1; i<=subsets.size(); i++)
        {
        Subset* s = findSubsetIndexed(i);
        std::set<int>& vals = s->getValues();
        
        nlohmann::json jIdx = nlohmann::json::array();
        for (int v : vals)
            jIdx.push_back(v);

        nlohmann::json jPart = nlohmann::json::object();
        jPart["index"] = s->getIndex();
        jPart["set"] = jIdx;
        jPart["name"] = s->getLabel();
        
        jSubsets.push_back(jPart);
        }
    j["partition"] = jSubsets;
    
    return j;
}

nlohmann::json Partition::toFile(std::map<int,double>& partFreqs, std::ostream& file) {

    file << "\nEquilibriumFrequencies = [\n  ";
    nlohmann::json j = nlohmann::json::object();
    
    nlohmann::json jSubsets = nlohmann::json::array();
    for (int i=1; i<=subsets.size(); i++)
        {
        Subset* s = findSubsetIndexed(i);
        std::set<int>& vals = s->getValues();
        
        std::map<int,double>::iterator it = partFreqs.find(i);
        if (it == partFreqs.end())
            Msg::error("Could not find subset " + std::to_string(i));
        double subsetFreq = it->second;
        
        nlohmann::json jIdx = nlohmann::json::array();
        for (int v : vals)
            jIdx.push_back(v);

        nlohmann::json jPart = nlohmann::json::object();
        jPart["index"] = s->getIndex();
        jPart["set"] = jIdx;
        jPart["name"] = s->getLabel();
        jPart["freq"] = subsetFreq;
        
        jSubsets.push_back(jPart);
    
        if (i > 1)
            file << ", ";
        file << subsetFreq;

        }
    j["partition"] = jSubsets;
    
    file << "\n];\n\n";

    return j;

}
