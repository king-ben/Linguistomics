#ifndef RateMatrixHelper_hpp
#define RateMatrixHelper_hpp

#include <map>
#include <set>
#include <string>
#include <vector>

struct GroupPair {

                                    GroupPair(int i, int j);
        int                         group1;
        int                         group2;
};

struct CompGroupPair {

    bool operator()(const GroupPair& p1, const GroupPair& p2) const {
        
        if (p1.group1 < p2.group1)
            return true;
        else if (p1.group1 == p2.group1)
            {
            if (p1.group2 < p2.group2)
                return true;
            }
        return false;
        }
};

class Partition;

class RateMatrixHelper {

    public:
                                                RateMatrixHelper(void);
                                               ~RateMatrixHelper(void);
        std::set<int>&                          getGroup(int idx) { return stateGroupings[idx]; }
        std::map<GroupPair,int, CompGroupPair>& getGroupIndices(void) { return groupIndices; }
        std::vector<std::string>                getLabels(void);
        int**                                   getMap(void) { return m; }
        int                                     getNumRates(void) { return (int)groupIndices.size(); }
        int                                     getNumStates(void) { return numStates; }
        int                                     groupIdForState(int stateIdx);
        void                                    initialize(int ns, Partition* part);
        void                                    print(void);
        void                                    printMap(void);
        
    private:
        int**                                   m;
        int                                     numStates;
        bool                                    isInitialized;
        std::vector<std::set<int> >             stateGroupings;
        std::vector<std::string>                stateGroupingsNames;
        int                                     numGroups;
        std::map<GroupPair,int, CompGroupPair>  groupIndices;
};

#endif
