#include <iomanip>
#include <iostream>
#include <vector>
#include "Alignment.hpp"
#include "AlignmentDistribution.hpp"
#include "Msg.hpp"
#include "Partition.hpp"
#include "RandomVariable.hpp"



bool cmp(std::pair<Alignment*, int>& a, std::pair<Alignment*, int>& b) {

    return a.second > b.second;
}



AlignmentDistribution::AlignmentDistribution(RandomVariable* r, std::vector<Alignment*>& alns, std::string n) : rv(r), name(n) {

    for (size_t i=0; i<alns.size(); i++)
        addAlignmentCopy(alns[i]);
}

AlignmentDistribution::~AlignmentDistribution(void) {

    for (auto [key,val] : samples)
        delete key;
}

AlignmentDistribution& AlignmentDistribution::operator+=(const AlignmentDistribution& rhs) {

    // check that the partition information is the same in each
    //if ( this->partition->isEqualTo(*(rhs.partition)) == false )
    // combine the information
    for (std::map<Alignment*,int>::const_iterator it = rhs.samples.begin(); it != rhs.samples.end(); it++)
        {
        Alignment* rhsAln = it->first;
        this->addAlignmentCopy(rhsAln);
        }
    
    return *this;
}

void AlignmentDistribution::addAlignment(Alignment* aln) {

    bool foundAlignmentInMap = false;
    for (std::map<Alignment*,int>::iterator it = samples.begin(); it != samples.end(); it++)
        {
        if (*it->first == *aln)
            {
            it->second++;
            foundAlignmentInMap = true;
            delete aln;
            break;
            }
        }
        
    if (foundAlignmentInMap == false)
        {
        samples.insert( std::make_pair(aln,1) );
        }
}

void AlignmentDistribution::addAlignmentCopy(Alignment* aln) {

    bool foundAlignmentInMap = false;
    for (std::map<Alignment*,int>::iterator it = samples.begin(); it != samples.end(); it++)
        {
        if (*it->first == *aln)
            {
            it->second++;
            foundAlignmentInMap = true;
            break;
            }
        }
        
    if (foundAlignmentInMap == false)
        {
        Alignment* newAlignment = new Alignment(*aln);
        samples.insert( std::make_pair(newAlignment,1) );
        }
}

int AlignmentDistribution::numSamples(void) {

    int n = 0;
    for (auto& it : samples)
        n += it.second;
    return n;
}

int AlignmentDistribution::ciSize(void) {

    std::vector<std::pair<Alignment*, int> > v;
    for (auto& it : samples)
        v.push_back(it);
  
    sort(v.begin(), v.end(), cmp);

    int n = numSamples();
    double cumulativeProb = 0.0;
    int i = 0;
    for (auto& it : v)
        {
        double prob = (double)it.second / n;
        cumulativeProb += prob;
        if (cumulativeProb < 0.95)
            {
            i++;
            }
        else
            {
            double lowerVal = cumulativeProb - prob;
            double u = rv->uniformRv();
            if (u < (0.95-lowerVal)/(cumulativeProb-lowerVal))
                {
                i++;
                }
            break;
           
            }
        }
    return i;
}

Alignment* AlignmentDistribution::getMapAlignment(void) {

    std::vector<std::pair<Alignment*, int> > v;
    for (auto& it : samples)
        v.push_back(it);
    sort(v.begin(), v.end(), cmp);
    return v[0].first;
}

void AlignmentDistribution::print(void) {

    print(true, nullptr);
}

void AlignmentDistribution::print(bool showAlignments) {

    print(showAlignments, nullptr);
}

void AlignmentDistribution::print(bool showAlignments, Partition* part) {

    std::vector<std::pair<Alignment*, int> > v;
    for (auto& it : samples)
        v.push_back(it);

    int n = numSamples();
  
    sort(v.begin(), v.end(), cmp);
    int credibleSetSize = 0;
    double cumulativeProb = 0.0;
    for (auto& it : v)
        {
        double prob = (double)it.second / n;
        cumulativeProb += prob;
        if (cumulativeProb < 0.95)
            {
            credibleSetSize++;
            }
        else
            {
            double lowerVal = cumulativeProb - prob;
            double u = rv->uniformRv();
            if (u < (0.95-lowerVal)/(cumulativeProb-lowerVal))
                credibleSetSize++;
            break;
            }
        }
    if (credibleSetSize == 0)
        credibleSetSize = 1;
        
    std::cout << std::endl << "   Cognate: " << name << std::endl;
    std::cout << "      Number of alignments in 95% credible set: " << credibleSetSize << std::endl;

    if (showAlignments == false)
        return;

    cumulativeProb = 0.0;
    int i = 0;
    for (auto& it : v)
        {
        double prob = (double)it.second / n;
        cumulativeProb += prob;
        if (cumulativeProb < 0.95)
            {
            //it.first->print( std::to_string(++i) + " Probability: " + std::to_string(prob) + " (" + std::to_string(cumulativeProb) + ")" );
            std::cout << ++i << " -- Probability: " << prob << " (" << cumulativeProb << ")" << std::endl;
            print(it.first, part);
            }
        else
            {
            double lowerVal = cumulativeProb - prob;
            double u = rv->uniformRv();
            if (u < (0.95-lowerVal)/(cumulativeProb-lowerVal))
                {
                //it.first->print( std::to_string(++i) + " Probability: " + std::to_string(prob) + " (" + std::to_string(cumulativeProb) + ")" );
                std::cout << ++i << " -- Probability: " << prob << " (" << cumulativeProb << ")" << std::endl;
                print(it.first, part);
                }
            break;
            }
        }
}

void AlignmentDistribution::print(Alignment* aln, Partition* part) {

    int len = aln->lengthOfLongestName();
    for (int i=0; i<aln->getNumTaxa(); i++)
        {
        std::cout << "      ";
        std::string tName = (*aln)[i].getName();
        std::cout << tName << " ";
        for (int j=0; j<len-tName.length(); j++)
            std::cout << " ";
        for (int j=0; j<aln->getNumChar(); j++)
            {
            int8_t charCode = aln->getCharCode(i, j);
            int subsetIndex = static_cast<int>(charCode);
            if (subsetIndex != -1 && part != nullptr)
                subsetIndex = part->indexOfSubsetWithValue(charCode);
            if (charCode == -1)
                std::cout << " - ";
            else
                std::cout << std::setw(2) << subsetIndex << " ";
            }
        std::cout << std::endl;
        }
}

