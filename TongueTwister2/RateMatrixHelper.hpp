#ifndef RateMatrixHelper_hpp
#define RateMatrixHelper_hpp

#include <string>
#include <vector>
#include "Container.hpp"
class Partition;



class RateMatrixHelper {

    public:
                                    RateMatrixHelper(void) = delete;
                                    RateMatrixHelper(size_t ns, Partition* part);
                                   ~RateMatrixHelper(void);
        
                                    // which rate parameter governs the i->j transition
        int                         getRateIndex(size_t i, size_t j) const {
                                        return rateIndexMap[i * numStates + j];
                                    }
        
                                    // which group does state belong to
        int                         getGroupForState(size_t state) const {
                                        return stateToGroup[state];
                                    }
        
                                    // rate index for group pair
        int                         getRateIndexForGroups(size_t g1, size_t g2) const {
                                        return groupPairToRate[g1 * numGroups + g2];
                                    }
        
                                    // accessors
        size_t                      getNumStates(void) const { return numStates; }
        size_t                      getNumGroups(void) const { return numGroups; }
        int                         getNumRates(void) const { return numRates; }
        const int*                  getRateIndexMapPtr(void) const { return rateIndexMap; }
        
                                    // For compatibility - returns pointer to rate index for state pair
                                    // Usage: rateIdx = helper->getRateIndex(i, j);
        std::vector<std::string>    getLabels(void) const;
        void                        print(void) const;
        void                        printMap(void) const;
        
    private:
        void                        initialize(Partition* part);
        size_t                      numStates;
        size_t                      numGroups;
        int                         numRates;
        int*                        stateToGroup;       // [numStates] -> group index
        int*                        groupPairToRate;    // [numGroups * numGroups] -> rate index
        int*                        rateIndexMap;       // [numStates * numStates] -> rate index
        std::vector<std::string>    groupNames;         // for labels/printing only
};

#endif
