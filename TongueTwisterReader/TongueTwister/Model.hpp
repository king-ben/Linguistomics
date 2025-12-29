#ifndef Model_hpp
#define Model_hpp

#include <map>
#include <set>
#include <string>
#include <vector>
#include "json.hpp"
#include "Parameter.hpp"
#include "RbBitSet.h"
#include "Threads.hpp"
class Alignment;
class LikelihoodCalculator;
class ParameterAlignment;
class ParameterEquilibirumFrequencies;
class ParameterExchangabilityRates;
class ParameterTree;
class Partition;
class RandomVariable;
class RateMatrix;
class RateMatrixHelper;
class TransitionProbabilities;
class Tree;
class UserSettings;


class WordLnLikeTask: public ThreadTask {

    public:
                              WordLnLikeTask(void);
        void                  Init(LikelihoodCalculator* calculator, double* threadLnL, double* wordLnL);
        virtual void          Run(MathCache& cache);

    private:
        LikelihoodCalculator* Calculator;
        double*               ThreadLnL;
        double*               WordLnL;
};


class Model {

    public:
        static const int                        PrintNumberSize = 12;

                                                Model(void) = delete;
                                       explicit Model(RandomVariable* r, ThreadPool& threadPool);
                                               ~Model(void);
        void                                    accept(void);
        void                                    fillParameterValues(double* x, int n);
        void                                    flipActiveLikelihood(void);
        void                                    flipActiveLikelihood(int idx);
        ParameterAlignment*                     getAlignment(int idx);
        std::vector<ParameterAlignment*>        getAlignments(void);
        std::vector<std::string>&               getCanonicalTaxonList(void) { return canonicalTaxonList; }
        double                                  getDeletionRate(void);
        std::vector<double>&                    getEquilibriumFrequencies(void);
        std::vector<double>&                    getExchangabilityRates(void);
        std::vector<double>&                    getIndelGammaRates(void);
        std::vector<double>&                    getRatesGammaRates(void);
        int                                     getIndex(void) { return index; }
        double                                  getInsertionRate(void);
        std::pair<Parameter*,int>&              getLastUpdate(void);
        int                                     getNumAlignments(void);
        std::string                             getParameterHeader(void);
        std::string                             getParameterString(void);
        ParameterTree*                          getParameterTree(void);
        int                                     getNumParameterValues(void);
        RateMatrix*                             getRateMatrix(void) { return rateMatrix; }
        std::string                             getStateSetsJsonString(void);
        TransitionProbabilities*                getTransitionProbabilities(void) { return transitionProbabilities; }
        Tree*                                   getTree(void);
        Tree*                                   getTree(RbBitSet& mask);
        Tree*                                   getTree(const RbBitSet& mask);
        std::string                             getUpdatedParameterName(void);
        double                                  lnLikelihood(void);
        double                                  lnPriorProbability(void);
        void                                    reject(void);
        void                                    setIndex(int x) { index = x; }
        void                                    setUpdateLikelihood(void);
        void                                    setUpdateLikelihood(int idx);
        double                                  update(int iter);
            
    private:
        WordLnLikeTask*                         getTaskList(size_t count);
        std::vector<Alignment*>                 initializeAlignments(nlohmann::json& j);
        void                                    initializeParameters(std::vector<Alignment*>& wordAlignments, nlohmann::json& j);
        void                                    initializeParameterString(void);
        void                                    initializeStateSets(nlohmann::json& j);
        void                                    initializeTransitionProbabilities(std::vector<Alignment*>& wordAlignments);
        nlohmann::json                          parseJsonFile(void);
        void                                    wordLnLike(int i);

        RandomVariable*                         rv;
        double*                                 threadLnL;
        std::vector<bool>                       updateLikelihood;
        std::vector<int>                        activeLikelihood;
        std::vector<double>                     wordLnLikelihoods[2];
        std::vector<Parameter*>                 parameters;
        std::vector<ParameterAlignment*>        wordParameterAlignments;
        std::vector<LikelihoodCalculator*>      wordLikelihoodCalculators;
        Partition*                              partitionInfo;
        std::map<std::string,std::set<int> >    stateSets;
        std::vector<std::string>                canonicalTaxonList;
        WordLnLikeTask*                         taskList;
        ThreadPool                              &threadPool;
        RateMatrix*                             rateMatrix;
        RateMatrixHelper*                       rateMatrixHelper;
        TransitionProbabilities*                transitionProbabilities;
        size_t                                  taskMax;
        int                                     updatedParameterIdx;
        int                                     substitutionModel;
        int                                     index;
        ParameterTree*                          treeParameter;
        ParameterExchangabilityRates*           exchangeabilityParameter;
        ParameterEquilibirumFrequencies*        equilibirumFreqParameter;
        char*                                   parameterString;
};

#endif
