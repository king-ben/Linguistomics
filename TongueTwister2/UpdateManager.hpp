#ifndef UpdateManager_hpp
#define UpdateManager_hpp

#include <unordered_map>
#include <vector>
#include "UpdateStatistics.hpp"
#include "UserSettings.hpp"
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
        void                                accept(Update* u);
        void                                markCognatesDirty(Update* u);
        void                                print(void);
        Update*                             randomlyChooseUpdate(void);
        void                                reject(Update* u);
        void                                summary(void);
        void                                tune(void);
        void                                updateDependants(Update* u);
        void                                zeroOut(void);
    
    private:
        void                                buildAliasTable(void);
        void                                setProposalProbabilities(void);
        Model*                              model;
        RateMatrix*                         rateMatrix;
        TransitionProbabilities*            tiProbs;
        RandomVariable*                     rng;
        UpdateStatistics                    updateInfo;
        std::vector<double>                 aliasProbability;
        std::vector<double>                 proposalProbabilities;
        std::vector<size_t>                 aliasIndex;
        std::vector<Update*>                updates;
        std::vector<UpdateAlignment*>       alignmentUpdates;
        std::vector<Update*>                otherUpdates;
        std::unordered_map<Update*, size_t> updateIndex;
        RateMatrix* rateMatrixForUpdate(Update* u);
};

#endif
