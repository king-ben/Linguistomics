#ifndef RateMatrixHelper_hpp
#define RateMatrixHelper_hpp

#include <map>
#include <set>
#include <string>
#include <vector>
#include "Container.hpp"

struct GroupPair {

                                    GroupPair(int i, int j);
        size_t                      group1;
        size_t                      group2;
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
                                                RateMatrixHelper(void) = delete;
                                                RateMatrixHelper(size_t ns, Partition* part);
                                               ~RateMatrixHelper(void);
        std::set<int>&                          getGroup(size_t idx) { return stateGroupings[idx]; }
        std::map<GroupPair,int, CompGroupPair>& getGroupIndices(void) { return groupIndices; }
        std::vector<std::string>                getLabels(void);
        IntMatrix&                              getMap(void) { return m; }
        int                                     getNumRates(void) { return static_cast<int>(groupIndices.size()); }
        size_t                                  getNumStates(void) { return numStates; }
        int                                     groupIdForState(size_t stateIdx);
        void                                    print(void);
        void                                    printMap(void);
        
    private:
        void                                    initialize(Partition* part);

                                                // accessed during rate matrix construction
        size_t                                  numStates;
        size_t                                  numGroups;
        
                                                // matrix accessed for lookups
        IntMatrix                               m;
        
                                                // containers (keep together, all are 24-56 bytes each)
        std::vector<std::set<int>>              stateGroupings;
        std::vector<std::string>                stateGroupingsNames;
        std::map<GroupPair,int, CompGroupPair>  groupIndices;
};

#endif
