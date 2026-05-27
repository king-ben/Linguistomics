#ifndef Model_hpp
#define Model_hpp

#include <set>
#include <vector>
#include "UserSettings.hpp"
#include "FamilyData.hpp"
#include "SubModel.hpp"
class LikelihoodCalculator;
class Parameter;
class Partition;
class RandomVariable;
class RateMatrix;
class States;
class ThreadPool;
class TransitionProbabilities;



class Model {

    public:
                                                Model(void) = delete;
                                                Model(RandomVariable* r, SubstitutionModel mt, ThreadPool* p);
                                               ~Model(void);
        template<typename T> T*                 findParameter(void);
        template<typename T> bool               isType(Parameter* parm);
        const std::vector<LikelihoodCalculator*>& getCalculators(void) const { return calculators; }
        size_t                                  getNumCalculators(void) const { return calculators.size(); }
        size_t                                  getNumStates(void) { return numStates; }
        const std::vector<Parameter*>&          getParameters(void) { return parameters; }
        RateMatrix*                             getRateMatrix(void) { return rateMatrix; }
        TransitionProbabilities*                getTiProbs(void) { return tiProbs; }
        double                                  lnLikelihood(void);
        double                                  lnPrior(void);
        void                                    markCognateDirty(size_t cognateIdx);                // mark a specific cognate as needing recalculation
        void                                    markAllCognatesDirty(void);                         // mark all cognates as needing recalculation
        size_t                                  getCognateIndex(LikelihoodCalculator* calc) const;  // get the index of a cognate by its LikelihoodCalculator pointer
        void                                    initializeLikelihoodCache(void);                    // initialize the likelihood cache
        void                                    keepLikelihoodCache(void);                          // MCMC state management for likelihood cache
        void                                    restoreLikelihoodCache(void);
        bool                                    hasAnyCognateDirty(void) const;                     // check if any cognates are dirty
        size_t                                  getNumDirtyCognates(void) const;                    // get count of dirty cognates (for diagnostics)
        void                                    registerParameter(Parameter* p);
        void                                    registerCalculator(LikelihoodCalculator* c);
        States*                                 getStates(void) const { return states; }
        int                                     getNumFamilies(void) const { return (int)subModels.size(); }
        SubModel*                               getSubModel(int i) const { return subModels[i]; }

    private:
        bool                                    alignmentsAreConsistent(std::vector<class ParameterAlignment*>& alns);
        bool                                    initializeAlignments(void);
        bool                                    initializeCalculators(void);
        bool                                    initializeIndelRates(void);
        bool                                    initializeStates(void);
        bool                                    initializeSubstitutionModel(void);
        bool                                    initializeTree(void);
        bool                                    initializeSingleFamily(void);
        void                                    initializeMultiFamily(void);
        void                                    plotAscii(const std::map<int,int>& data);
        RandomVariable*                         rng;
        SubstitutionModel                       modelType;
        ThreadPool*                             pool;
        States*                                 states;
        RateMatrix*                             rateMatrix;
        TransitionProbabilities*                tiProbs;
        size_t                                  numStates;
        std::vector<Parameter*>                 parameters;
        std::vector<LikelihoodCalculator*>      calculators;
        std::set<unsigned>                      uniqueTaxonCombinations;
        
                                                // current state: cached log-likelihoods (the accepted/known-good values)
        std::vector<double>                     cachedLnL;          // one entry per cognate
        std::vector<bool>                       dirtyFlags;         // which cognates need recomputation
        double                                  cachedTotalLnL;     // sum of cachedLnL
        
                                                /* Proposed state: new values computed during lnLikelihood()
                                                   These only become the cached values if the proposal is accepted */
        std::vector<double>                     proposedLnL;        // proposed new values
        double                                  proposedTotalLnL;   // proposed new total
        bool                                    cacheInitialized;   // flag to track if cache is initialized
        std::vector<size_t>                     dirtyIndices;       // temporary storage for dirty calculator 
        std::vector<LikelihoodCalculator*>      dirtyCalcs;         // indices (avoid repeated allocation)
        std::vector<FamilyData*>                familyData;
        std::vector<SubModel*>                  subModels;
};



template<typename T>
T* Model::findParameter(void) {

    for (size_t i=0; i<parameters.size(); i++)
        {
        T* derivedPtr = dynamic_cast<T*>(parameters[i]);
        if (derivedPtr)
            return derivedPtr;
        }
    return nullptr;
}

template<typename T>
bool Model::isType(Parameter* parm) {

    T* derivedPtr = dynamic_cast<T*>(parm);
    if (derivedPtr)
        return true;
    return false;
}

#endif
