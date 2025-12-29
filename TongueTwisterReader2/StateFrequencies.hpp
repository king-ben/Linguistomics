#ifndef StateFrequencies_hpp
#define StateFrequencies_hpp

#include "Simplex.hpp"
class McmcSummary;
class Partition;



class StateFrequencies : public Simplex {

    public:
                                            StateFrequencies(void) = delete;
                                            StateFrequencies(McmcSummary& samples, double bf);
        size_t                              getFirstIndex(size_t idx) { return stateValue[idx]; }
        size_t                              getSecondIndex(size_t idx) { return 0; }
        size_t                              getNumStates(void) { return dimension; }
        void                                print(void);
        void                                print(Partition* part);
        
    private:
        std::vector<int>                    stateValue;
};

#endif
