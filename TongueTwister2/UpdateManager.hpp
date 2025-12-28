#ifndef UpdateManager_hpp
#define UpdateManager_hpp

#include <unordered_map>
#include <vector>
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
        
        Model*                          model;
        RandomVariable*                 rng;
        RateMatrix*                     rateMatrix;
        TransitionProbabilities*        tiProbs;
        
        std::vector<Update*>            updates;
        std::vector<UpdateAlignment*>   alignmentUpdates;
        std::vector<Update*>            otherUpdates;
        
                                        // map from Update pointer to index (for O(1) lookup)
        std::unordered_map<Update*, size_t> updateIndex;
        
                                        // proposal probability management (decoupled from Update objects)
        std::vector<double>             proposalProbabilities;
        
                                        // acceptance statistics (decoupled from Update objects)
        std::vector<int>                numTries;
        std::vector<int>                numAcceptances;
        
                                        // Walker's alias method tables for O(1) selection
        std::vector<double>             aliasProbability;
        std::vector<size_t>             aliasIndex;
};

#endif
