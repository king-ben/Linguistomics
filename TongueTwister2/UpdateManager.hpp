#ifndef UpdateManager_hpp
#define UpdateManager_hpp

#include <vector>
class Model;
class Parameter;
class RateMatrix;
class TransitionProbabilityManager;
class Update;



class UpdateManager {

    public:
                                        UpdateManager(void) = delete;
                                        UpdateManager(Model* m, RandomVariable* r);
                                       ~UpdateManager(void);
        void                            accept(Update* u);
        void                            print(void);
        Update*                         randomlyChooseUpdate(void);
        void                            reject(Update* u);
        void                            summary(void);
        void                            updateDependants(Update* u);
    
    private:
        void                            buildAliasTable(void);
        void                            setProposalProbabilities(void);
        
        Model*                          model;
        RandomVariable*                 rng;
        RateMatrix*                     rateMatrix;
        TransitionProbabilityManager*   tiProbs;
        
        std::vector<Update*>            updates;
        std::vector<Update*>            alignmentUpdates;
        std::vector<Update*>            otherUpdates;
        
                                        // proposal probability management (decoupled from Update objects)
        std::vector<double>             proposalProbabilities;
        
                                        // Walker's alias method tables for O(1) selection
        std::vector<double>             aliasProbability;
        std::vector<size_t>             aliasIndex;
};

#endif
