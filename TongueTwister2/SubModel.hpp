#pragma once
#include <set>
#include <vector>
#include "UserSettings.hpp"

class FamilyData;
class LikelihoodCalculator;
class Model;
class ParameterAlignment;
class ParameterExchangeabilities;
class ParameterFrequencies;
class ParameterIndelRates;
class ParameterTree;
class RandomVariable;
class RateMatrix;
class ThreadPool;
class TransitionProbabilities;

class SubModel {

    public:
                        SubModel(Model* parentModel,
                                 RandomVariable* rng,
                                 FamilyData* familyData,
                                 int familyIndex,
                                 ThreadPool* pool);
                       ~SubModel(void);

        // initialisation (called in constructor order)
        bool            initializeAlignments(ParameterIndelRates* sharedIndelRates);
        bool            initializeTree(void);
        bool            initializeSubstitutionModel(
                            ParameterFrequencies* sharedFreqs,
                            ParameterExchangeabilities* sharedExch,
                            RateMatrix* sharedRateMatrix,
                            SubstitutionModel modelType);
        bool            initializeCalculators(
                            ParameterIndelRates* indelRates,
                            ParameterFrequencies* freqs);

        // accessors
        ParameterTree*              getTree(void)           const { return tree; }
        TransitionProbabilities*    getTiProbs(void)        const { return tiProbs; }
        ParameterIndelRates*        getIndelRates(void)     const { return indelRates; }
        ParameterFrequencies*       getFrequencies(void)    const { return frequencies; }
        ParameterExchangeabilities* getExchangeabilities(void) const { return exchangeabilities; }
        RateMatrix*                 getRateMatrix(void)     const { return rateMatrix; }
        
        const std::vector<ParameterAlignment*>&  getAlignments(void)   const { return alignments; }
        const std::vector<LikelihoodCalculator*>& getCalculators(void) const { return calculators; }
        const std::set<unsigned>&   getUniqueTaxonCombinations(void)   const { return uniqueTaxonCombinations; }
        int                         getFamilyIndex(void)               const { return familyIndex; }

    private:
        // back-pointer to parent model (for parameter registration)
        Model*                          parentModel;
        RandomVariable*                 rng;
        FamilyData*                     familyData;
        ThreadPool*                     pool;
        int                             familyIndex;
        size_t                          numStates;

        // owned parameters (only when not shared)
        ParameterTree*                  tree;
        ParameterIndelRates*            indelRates;
        ParameterFrequencies*           frequencies;
        ParameterExchangeabilities*     exchangeabilities;
        RateMatrix*                     rateMatrix;
        TransitionProbabilities*        tiProbs;
        
        std::vector<ParameterAlignment*>    alignments;
        std::vector<LikelihoodCalculator*>  calculators;
        std::set<unsigned>                  uniqueTaxonCombinations;
        
        bool                            ownsIndelRates;
        bool                            ownsFrequencies;
        bool                            ownsExchangeabilities;
        bool                            ownsRateMatrix;
};