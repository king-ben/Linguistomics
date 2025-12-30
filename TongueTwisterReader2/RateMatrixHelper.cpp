#include <iomanip>
#include <iostream>
#include <set>
#include <map>
#include "Msg.hpp"
#include "Partition.hpp"
#include "RateMatrixHelper.hpp"
#include "Subset.hpp"



RateMatrixHelper::RateMatrixHelper(size_t ns, Partition* part) : 
    numStates(ns),
    numGroups(0),
    numRates(0),
    stateToGroup(nullptr),
    groupPairToRate(nullptr),
    rateIndexMap(nullptr) {

    initialize(part);
}

RateMatrixHelper::~RateMatrixHelper(void) {

    delete [] stateToGroup;
    delete [] groupPairToRate;
    delete [] rateIndexMap;
}

void RateMatrixHelper::initialize(Partition* part) {

    numGroups = part->numSubsets();
    
    // allocate lookup arrays
    stateToGroup = new int[numStates];
    groupPairToRate = new int[numGroups * numGroups];
    rateIndexMap = new int[numStates * numStates];
    
    // initialize to -1 (invalid)
    for (size_t i=0; i<numStates; i++)
        stateToGroup[i] = -1;
    for (size_t i=0; i < numGroups * numGroups; i++)
        groupPairToRate[i] = -1;
    for (size_t i=0; i<numStates * numStates; i++)
        rateIndexMap[i] = -1;
    
    // build stateToGroup mapping and collect group names
    for (size_t g=0; g<numGroups; g++)
        {
        Subset* ss = part->findSubsetIndexed(static_cast<int>(g + 1));
        groupNames.push_back(ss->getLabel());
        
        const std::set<int>& groupSet = ss->getValues();
        for (int state : groupSet)
            {
            if (state < 0 || static_cast<size_t>(state) >= numStates)
                Msg::error("Found state that is not sensible in state groups: " + std::to_string(state));
            if (stateToGroup[state] != -1)
                Msg::error("State " + std::to_string(state) + " appears in multiple groups");
            stateToGroup[state] = static_cast<int>(g);
            }
        }
    for (size_t i=0; i<numStates; i++)
        {
        if (stateToGroup[i] == -1)
            Msg::error("Could not find state " + std::to_string(i) + " in state groups");
        }
    
    // Build groupPairToRate mapping
    // First pass: discover unique group pairs and assign rate indices
    std::map<std::pair<int,int>, int> tempGroupPairs;
    for (size_t i=0; i<numStates; i++)
        {
        for (size_t j=i+1; j<numStates; j++)
            {
            int gi = stateToGroup[i];
            int gj = stateToGroup[j];
            
            // canonical order: smaller group index first
            int g1 = (gi < gj) ? gi : gj;
            int g2 = (gi < gj) ? gj : gi;
            
            std::pair<int,int> key = std::make_pair(g1, g2);
            if (tempGroupPairs.find(key) == tempGroupPairs.end())
                tempGroupPairs[key] = 0;  // placeholder
            }
        }
    
    // assign sequential rate indices
    numRates = 0;
    for (auto& kv : tempGroupPairs)
        {
        kv.second = numRates++;
        int g1 = kv.first.first;
        int g2 = kv.first.second;
        // store in both directions for symmetric lookup
        groupPairToRate[g1 * numGroups + g2] = kv.second;
        groupPairToRate[g2 * numGroups + g1] = kv.second;
        }
    
    // build the full rateIndexMap for fast lookup during updateRateMatrix
    for (size_t i=0; i<numStates; i++)
        {
        int gi = stateToGroup[i];
        for (size_t j=0; j<numStates; j++)
            {
            if (i == j)
                {
                rateIndexMap[i * numStates + j] = -1;  // diagonal
                }
            else
                {
                int gj = stateToGroup[j];
                rateIndexMap[i * numStates + j] = groupPairToRate[gi * numGroups + gj];
                }
            }
        }
}

std::vector<std::string> RateMatrixHelper::getLabels(void) const {

    std::vector<std::string> labels(numRates);
    
    // reconstruct labels from groupPairToRate
    for (size_t g1=0; g1<numGroups; g1++)
        {
        for (size_t g2=g1+1; g2<numGroups; g2++)
            {
            int rateIdx = groupPairToRate[g1 * numGroups + g2];
            if (rateIdx >= 0 && rateIdx < numRates)
                {
                labels[rateIdx] = std::to_string(g1 + 1) + "-" + std::to_string(g2 + 1);
                }
            }
        }
    
    // handle within-group rates if they exist
    for (size_t g=0; g<numGroups; g++)
        {
        int rateIdx = groupPairToRate[g * numGroups + g];
        if (rateIdx >= 0 && rateIdx < numRates)
            {
            labels[rateIdx] = std::to_string(g + 1) + "-" + std::to_string(g + 1);
            }
        }
    
    return labels;
}

void RateMatrixHelper::print(void) const {

    std::cout << "   * Custom rate matrix with states grouped into " << numGroups << " sets:" << std::endl;
    
    for (size_t g=0; g<numGroups; g++)
        {
        std::cout << "     Group " << (g + 1) << " \"" << groupNames[g] << "\" = ( ";
        for (size_t s=0; s<numStates; s++)
            {
            if (stateToGroup[s] == static_cast<int>(g))
                std::cout << s << " ";
            }
        std::cout << ")" << std::endl;
        }
    
    std::cout << "     The model has " << numRates << " rate parameters" << std::endl;
}

void RateMatrixHelper::printMap(void) const {

    for (size_t i=0; i<numStates; i++)
        {
        for (size_t j=0; j<numStates; j++)
            {
            int idx = rateIndexMap[i * numStates + j];
            if (idx < 0)
                std::cout << " - ";
            else
                std::cout << std::setw(2) << idx << " ";
            }
        std::cout << std::endl;
        }
}
