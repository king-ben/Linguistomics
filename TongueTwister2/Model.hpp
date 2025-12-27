#ifndef Model_hpp
#define Model_hpp

#include <set>
#include <string>
#include <vector>
#include "UserSettings.hpp"
class Alignment;
class LikelihoodCalculator;
class Parameter;
class ParameterAlignment;
class RandomVariable;
class RateMatrix;
class States;
class ThreadPool;
class TransitionProbabilityManager;



class Model {

    public:
                                            Model(void) = delete;
                                            Model(RandomVariable* r, SubstitutionModel mt, ThreadPool* p);
                                           ~Model(void);
        template <typename T>
        T*                                  findParameter(void);
        size_t                              getNumStates(void) { return numStates; }
        const std::vector<Parameter*>&      getParameters(void) const { return parameters; }
        template <typename T>
        const std::vector<T*>               getParametersOfType(void) const;
        RateMatrix*                         getRateMatrix(void) { return rateMatrix; }
        TransitionProbabilityManager*       getTiProbs(void) { return tiProbs; }
        double                              lnLikelihood(void);
        double                              lnPrior(void);
    
    private:
        bool                                alignmentsAreConsistent(std::vector<ParameterAlignment*>& alignments);
        size_t                              calculateMaximumAlignmentLength(nlohmann::json j);
        bool                                initializeAlignments(void);
        bool                                initializeCalculators(void);
        bool                                initializeIndelRates(void);
        bool                                initializeStates(void);
        bool                                initializeSubstitutionModel(void);
        bool                                initializeTree(void);
        template <typename T>
        bool                                isType(const Parameter* p) { return dynamic_cast<const T*>(p) != nullptr; }
        void                                plotAscii(const std::map<int,int>& data);

        RandomVariable*                     rng;
        SubstitutionModel                   modelType;
        ThreadPool*                         pool;
        States*                             states;
        RateMatrix*                         rateMatrix;
        TransitionProbabilityManager*       tiProbs;
        size_t                              numStates;
        size_t                              gapCode;
        std::vector<Parameter*>             parameters;
        std::vector<LikelihoodCalculator*>  calculators;
        std::set<unsigned>                  uniqueTaxonCombinations;
};

template <typename T>
T* Model::findParameter(void) {

    for (auto* p : parameters) 
        {
        if (auto* q = dynamic_cast<T*>(p))
            return q;
        }
    return nullptr;
}

template <typename T>
const std::vector<T*> Model::getParametersOfType(void) const {

    std::vector<T*> vec;
    for (Parameter* p : parameters)
        {
        if (dynamic_cast<T*>(p) != nullptr)
            vec.push_back(dynamic_cast<T*>(p));
        }
    return vec;
}

#endif
