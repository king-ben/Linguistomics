#ifndef TransitionProbabilities_hpp
#define TransitionProbabilities_hpp

#include <map>
#include <string>
#include <vector>
#include "Container.hpp"
#include "RbBitSet.h"
#include "UserSettings.hpp"
#include "Threads.hpp"
class Alignment;
class Model;
class thread_pool;
class Tree;


struct TransitionProbabilitiesInfo {

    int             numMatrices;
    int             numStates;
    DoubleMatrix**  probs[2];
};



class TransitionProbabilities {

    public:
                                                        TransitionProbabilities(void);
                                                       ~TransitionProbabilities(void);
        bool                                            areTransitionProbabilitiesValid(double tolerance);
        void                                            flipActive(void);
        int                                             getNumNodes(void) { return numNodes; }
        int                                             getNumStates(void) { return numStates; }
        std::vector<double>&                            getStationaryFrequencies(void) { return stationaryFreqs[activeProbs]; }
        void                                            getStationaryFrequencies(std::vector<double>& f);
        DoubleMatrix**                                  getTransitionProbabilities(RbBitSet& bs);
        DoubleMatrix&                                   getTransitionProbabilities(RbBitSet& bs, int nodeIdx);
        void                                            initialize(Model* m, ThreadPool* pool, std::vector<Alignment*>& alns, int nn, int ns, int sm);
        void                                            print(void);
        void                                            setNeedsUpdate(bool tf) { needsUpdate = tf; }
        void                                            setTransitionProbabilities(void);
    
    private:
                                                        TransitionProbabilities(const TransitionProbabilities& tp) = delete;
        void                                            setTransitionProbabilitiesJc69(void);
        void                                            setTransitionProbabilitiesUsingPadeMethod(void);
        bool                                            isInitialized;
        Model*                                          modelPtr;
        ThreadPool*                                     threadPool;
        int                                             numNodes;
        int                                             numStates;
        int                                             activeProbs;
        std::map<RbBitSet,TransitionProbabilitiesInfo>  transProbs;
        std::vector<double>                             stationaryFreqs[2];
        bool                                            needsUpdate;
        int                                             substitutionModel;
        int                                             numRateCategories;
};

#endif
