#include <iomanip>
#include <iostream>
#include "Msg.hpp"
#include "Partition.hpp"
#include "RateMatrixHelper.hpp"
#include "Subset.hpp"
#include "UserSettings.hpp"



RateMatrixHelper::RateMatrixHelper(void) {
    numStates = 0;
    numGroups = 0;
    isInitialized = false;
    m = NULL;
}

RateMatrixHelper::~RateMatrixHelper(void) {

    if (m != NULL)
        {
        delete [] m[0];
        delete [] m;
        }
}

std::vector<std::string> RateMatrixHelper::getLabels(void) {

    std::vector<std::string> labels(groupIndices.size());
    for (std::map<GroupPair, int>::iterator it = groupIndices.begin(); it != groupIndices.end(); it++)
        {
        int idx = it->second;
        int i = it->first.group1;
        int j = it->first.group2;
        std::string str = std::to_string(i+1);
        str += "-";
        str += std::to_string(j+1);
        labels[idx] = str;
        }
    return labels;
}

int RateMatrixHelper::groupIdForState(int stateIdx) {

    for (int i=0; i<numGroups; i++)
        {
        std::set<int>::iterator it = stateGroupings[i].find(stateIdx);
        if (it != stateGroupings[i].end())
            return i;
        }
    return -1;
}

void RateMatrixHelper::initialize(int ns, Partition* part) {

    if (isInitialized == true)
        return;
                
    // allocate the matrix holding the transition types
    numStates = ns;
    m = new int*[numStates];
    m[0] = new int[numStates*numStates];
    for (int i=1; i<numStates; i++)
        m[i] = m[i-1] + numStates;
    for (int i=0; i<numStates; i++)
        for (int j=0; j<numStates; j++)
            m[i][j] = 0;
        
    // get the state groupings from the json object
//    numGroups = (int)js.size();
//    for (int i=0; i<numGroups; i++)
//        {
//        std::string groupName = js[i]["Name"];
//        std::vector<int> group = js[i]["Set"];
//        std::set<int> groupSet;
//        for (int j=0; j<group.size(); j++)
//            groupSet.insert(group[j]);
//        stateGroupings.push_back(groupSet);
//        stateGroupingsNames.push_back(groupName);
//        }
    numGroups = part->numSubsets();
    for (int i=1; i<=numGroups; i++)
        {
        Subset* ss = part->findSubsetIndexed(i);
        std::string groupName = ss->getLabel();
        std::set<int> groupSet = ss->getValues();
        stateGroupings.push_back(groupSet);
        stateGroupingsNames.push_back(groupName);
        }
        
    // check that all of the states are found in the stateGroupings object
    for (int i=0; i<numStates; i++)
        {
        if (groupIdForState(i) == -1)
            Msg::error("Could not find state " + std::to_string(i) + " in state groups");
        }
        
    // check that there are no out-of-bound states in the state groupings
    for (int i=0; i<numGroups; i++)
        {
        std::set<int>& group = getGroup(i);
        for (int j : group)
            {
            if (j < 0 || j >= numStates)
                Msg::error("Found state that is not sensible in state groups");
            }
        }
        
    // find the change types for this model
    for (int i=0; i<numStates; i++)
        {
        for (int j=i+1; j<numStates; j++)
            {
            int groupI = groupIdForState(i);
            int groupJ = groupIdForState(j);
//            int groupI = part->indexOfSubsetWithValue(i);
//            int groupJ = part->indexOfSubsetWithValue(j);
            GroupPair key(groupI, groupJ);
            std::map<GroupPair,int,CompGroupPair>::iterator it = groupIndices.find(key);
            if (it == groupIndices.end())
                groupIndices.insert( std::make_pair(GroupPair(groupI, groupJ), 0) );
            }
        }
    int idx = 0;
    for (std::map<GroupPair, int>::iterator it = groupIndices.begin(); it != groupIndices.end(); it++)
        {
        it->second = idx++;
        //std::cout << it->first.group1 << "-" << it->first.group2 << " -- " << it->second << std::endl;
        }
        
    // set up the map
    for (int i=0; i<numStates; i++)
        {
        for (int j=i+1; j<numStates; j++)
            {
            int groupI = groupIdForState(i);
            int groupJ = groupIdForState(j);
//            int groupI = part->indexOfSubsetWithValue(i);
//            int groupJ = part->indexOfSubsetWithValue(j);
            GroupPair key(groupI, groupJ);
            std::map<GroupPair,int,CompGroupPair>::iterator it = groupIndices.find(key);
            if (it == groupIndices.end())
                Msg::error("Could not find group");
            m[i][j] = it->second;
            m[j][i] = it->second;
            }
        }
    
    isInitialized = true;
    
#   if 0
    for (int i=0; i<numStates; i++)
        {
        for (int j=0; j<numStates; j++)
            {
            if (i == j)
                std::cout << " -";
            else
                std::cout << std::setw(2) << m[i][j];
            }
        std::cout << std::endl;
        }
#   endif
}

void RateMatrixHelper::print(void) {

    std::cout << "   * Custom rate matrix with states grouped into " << numGroups << " sets:" << std::endl;
    for (int i=0; i<stateGroupings.size(); i++)
        {
        std::cout << "     Group " << i+1 << " \"" << stateGroupingsNames[i] << "\" = ( ";
        for (int j : stateGroupings[i])
            std::cout << j << " ";
        std::cout << ")" << std::endl;
        }
    std::cout << "     The model has " << groupIndices.size() << " rate parameters" << std::endl;
}

void RateMatrixHelper::printMap(void) {

    for (int i=0; i<numStates; i++)
        {
        for (int j=0; j<numStates; j++)
            {
            std::cout << std::setw(2) << m[i][j] << " ";
            }
        std::cout << std::endl;
        }
}

GroupPair::GroupPair(int i, int j) {

     if (i < j)
        {
        group1 = i;
        group2 = j;
        }
    else if (i > j)
        {
        group1 = j;
        group2 = i;
        }
    else
        {
        group1 = i;
        group2 = j;
        }
}
