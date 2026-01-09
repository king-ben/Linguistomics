#ifndef UpdateManager_hpp
#define UpdateManager_hpp

#include <unordered_map>
#include <vector>
#include "UpdateStatistics.hpp"
class Model;
class Parameter;
class RateMatrix;
class RandomVariable;
class TransitionProbabilities;
class Update;
class UpdateAlignment;



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
        
                                        // accessed on every MCMC cycle
        RandomVariable*                 rng;
        
                                        // Walker's alias method tables for O(1) selection
                                        // accessed on every randomlyChooseUpdate() call
        std::vector<double>             aliasProbability;
        std::vector<size_t>             aliasIndex;
        
                                        // update vectors (accessed after selection)
        std::vector<Update*>            updates;
        
                                        // statistics vectors (accessed on accept/reject)
        std::vector<double>             proposalProbabilities;
        UpdateStatistics                updateInfo;
        
                                        // lookup map (accessed for index lookup)
        std::unordered_map<Update*, size_t> updateIndex;
        
                                        // specialized update vectors
        std::vector<UpdateAlignment*>   alignmentUpdates;
        std::vector<Update*>            otherUpdates;
        
                                        // model references (rarely dereferenced after init)
        Model*                          model;
        RateMatrix*                     rateMatrix;
        TransitionProbabilities*        tiProbs;
};

#endif
