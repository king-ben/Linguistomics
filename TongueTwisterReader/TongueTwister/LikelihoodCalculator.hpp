#ifndef LikelihoodCalculator_hpp
#define LikelihoodCalculator_hpp

#include <set>
#include <string>
#include <vector>
#include "IntVector.hpp"
#include "RbBitSet.h"
#include "TransitionProbabilities.hpp"
class IndelMatrix;
class Model;
class Node;
class ParameterAlignment;
//class TransitionProbabilities;
class Tree;
typedef std::map<IntVector*, double, CompIntVector> PartialProbabilitiesLookup;

struct IndelProbabilities {

    double      insertionRate;
    double      deletionRate;
    double      immortalProbability;
    double*     beta;
    double*     birthProbability;
    double*     extinctionProbability;
    double*     homologousProbability;
    double*     nonHomologousProbability;
};

struct IndelCombinatorics {

    int         maximumSequenceLength;
    int*        state;
    int*        possibleVectorIndices;
    int*        nodeHomology;
    int*        numHomologousEmissions;
    int*        numHomologousEmissionsForClass;
};



class LikelihoodCalculator {

    public:
                                        LikelihoodCalculator(void) = delete;
                                        LikelihoodCalculator(ParameterAlignment* a, Model* m);
                                       ~LikelihoodCalculator(void);
        std::string&                    alignmentName(void);
        double                          lnLikelihood(MathCache& cache);
    
    private:
        void                            allcoateIndelCombinatorics(int nn, int maxSeqLen);
        void                            allocateIndelProbabilities(int nn);
        void                            clearPpTable(void);
        void                            drainPool(void);
        void                            freeIndelCombinatorics(void);
        void                            freeIndelProbabilities(void);
        IntVector*                      getVector(void);
        IntVector*                      getVector(IntVector& vec);
        int                             getNumAllocated(void) { return (int) allocated.size(); }
        int                             getNumInPool(void) { return (int)pool.size(); }
        void                            initialize(void);
        double                          partialProbability(IntVector* signature, IntVector* pos, int site);
        void                            printTable(void);
        void                            returnToPool(IntVector* n);
        void                            setBirthDeathProbabilities(void);
        double                          prune(IntVector* signature, IntVector* pos, std::vector<Node*>& dpSequence, int site);
        void                            pruneBranch(Node* p, int site);


        RbBitSet                        taxonMask;
        ParameterAlignment*             data;
        Model*                          modelPtr;
        
        Tree*                           tree;
        std::vector<Node*>              des;       // instantiated once on construction and modified via reference
        int                             numTaxa;
        int                             numNodes;
        
        PartialProbabilitiesLookup      partialProbabilities;
        std::vector<IntVector*>         pool;
        std::set<IntVector*>            allocated;
        
        IndelMatrix*                    alignment;  // allocated once with a fixed size larger than an alignment should get
        std::vector<std::vector<int> >  sequences;  // instantiated once on construction of object, never modified
        int                             unalignableRegionSize;
        const int                       maxUnalignableDimension  = 10, maxUnalignableDimension1 = maxUnalignableDimension + 1;
        
        DoubleMatrix**                  transitionProbabilities;
        std::vector<double>             stateEquilibriumFrequencies; // resized once on construction of object and modified via reference
        int                             numStates, numStates1;
        
        enum                            StateLabels { freeToUse, possible, edgeUsed, used };
        IndelProbabilities              indelProbs;
        IndelCombinatorics              indelCombos;
        double**                        fH;    // fH and fI can be replaced by Shawn's container class later
        double**                        fI;
};

#endif
