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
        void                            setProposalProbabilities(void);
        Model*                          model;
        RandomVariable*                 rng;
        RateMatrix*                     rateMatrix;
        TransitionProbabilityManager*   tiProbs;
        std::vector<Update*>            updates;
        std::vector<Update*>            alignmentUpdates;
        std::vector<Update*>            otherUpdates;
};

#endif
