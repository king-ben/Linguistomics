#ifndef States_hpp
#define States_hpp

#include <iostream>
class Partition;



class States {

    public:
                                    States(void);
                                   ~States(void);
        size_t                      getNumStates(void) { return numStates; }
        Partition*                  getPartition(void) { return statePartitions; }
        void                        initialzieStates(void);
        void                        print(void);
    
    private:
                                    States(States&);
        States&                     operator=(const States&);
        Partition*                  statePartitions;
        size_t                      numStates;
        std::vector<std::string>    stateNames;
        std::vector<std::string>    segmentIPA;
};

#endif
