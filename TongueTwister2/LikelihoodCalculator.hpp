#ifndef LikelihoodCalculator_hpp
#define LikelihoodCalculator_hpp

#include "MathCache.hpp"
#include "Threads.hpp"
class Alignment;
class Node;
class ParameterAlignment;
class ParameterFrequencies;
class ParameterIndelRates;
class ParameterTree;
class TransitionProbabilityManager;
class Tree;



struct TKF91Probabilities {

    double      insertionRate;
    double      deletionRate;
    double      immortalProbability;
    double*     beta;
    double*     birthProbability;
    double*     extinctionProbability;
    double*     homologousProbability;
    double*     nonHomologousProbability;
};



struct TKF91Combinatorics {

    int*        nodeHomology;            // which homology class each node belongs to
    int*        numHomologousEmissions;  // count of homologous emissions in subtree
};



class LikelihoodCalculator : public LikelihoodTask {

    public:
                                        LikelihoodCalculator(void) = delete;
                                        LikelihoodCalculator(TransitionProbabilityManager* tpc, ParameterAlignment* a, ParameterTree* t, ParameterIndelRates* r, ParameterFrequencies* f);
                                       ~LikelihoodCalculator(void);
                            
    protected:
        double                          computeLnLikelihood(void) override;
    
    private:
        void                            allocateTKF91Probabilities(size_t nn);
        void                            allocateTKF91Combinatorics(size_t nn);
        void                            allocateConditionalProbs(void);
        void                            cacheTransitionMatrices(void);
        double                          computeRootFI(int* signature, size_t site);
        void                            computeFForInternalNode(Node* p);
        void                            computeFForLeafNode(Node* p, int* signature, size_t site);
        void                            freeTKF91Probabilities(void);
        void                            freeTKF91Combinatorics(void);
        void                            freeConditionalProbs(void);
        int                             getLeafIndex(Node* n);
        void                            initializeLeafIndices(void);
        void                            setBirthDeathProbabilities(void);
        void                            setStationaryFrequencies(void);
        TransitionProbabilityManager*   tiProbs;
        ParameterAlignment*             myAlignment;
        ParameterTree*                  myTree;
        ParameterIndelRates*            myIndelRates;
        ParameterFrequencies*           myFrequencies;
        unsigned                        taxonMask;
        Tree*                           tree;
        Alignment*                      alignment;
        size_t                          numTaxa;
        size_t                          numNodes;
        size_t                          numStates;
        size_t                          numSites;
        TKF91Probabilities              tkf91Probs;
        TKF91Combinatorics              tkf91Combos;
        std::vector<int>                leafIndexMap;              // maps node index to sequence index
        double**                        fH;                        // fH[nodeIdx][state] for homologous contributions (inherited from ancestor)
        double**                        fI;                        // fI[nodeIdx][state] for state 0..numStates-1, fI[nodeIdx][numStates] for gap
        double*                         equilibriumFrequencies;
        int*                            signature;
        
        // cached transition probability matrix pointers (indexed by node index)
        // avoids repeated map lookups in inner loops
        double**                        cachedTiMatrices;          // raw pointers to matrix data for each node
};

#endif
