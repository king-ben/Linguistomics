#ifndef AlignmentProposal_hpp
#define AlignmentProposal_hpp

#include <map>
#include <set>
#include <vector>
#include "Container.hpp"
#include "IntVector.hpp"
#include "RbBitSet.h"
class AlnMatrix;
class IndelMatrix;
class Model;
class ParameterAlignment;
class RandomVariable;
class Tree;


/* MCMC proposal that modifies a multiple sequence alignment by selecting a
   contiguous window of columns, stripping the gaps, and re-aligning the
   ungapped sequences using a stochastic progressive alignment guided by
   the current phylogeny.
   
   This is the same proposal mechanism used in 
   
   Lunter G, Miklós I, Drummond A, Jensen JL, Hein J. 2005. Bayesian coestimation 
      of phylogeny and sequence alignment. BMC Bioinformatics, 1;6:83. 
            
   The proposal follows a profile-profile alignment strategy:
     1. Select a region of the alignment (geometric length distribution).
     2. Extract ungapped "profiles" for each leaf from the selected columns.
     3. Progressively align profiles bottom-up along the tree using a 
        Needleman–Wunsch-style DP scored by transition probabilities,
        with a stochastic (Boltzmann-weighted) traceback.
     4. Splice the newly aligned region back into the full alignment.

   The Hastings ratio accounts for:
     - The probability of the forward stochastic traceback path.
     - The probability of the reverse path (deterministic traceback that
       follows the gap pattern of the old alignment).
     - The number of valid column orderings (countPaths) in each direction.
     - The geometric extension probability when the region length differs. */
class AlignmentProposal {

    public:
                                            AlignmentProposal(void) = delete;
                                            AlignmentProposal(ParameterAlignment* a, Tree* t, RandomVariable* r, Model* m, double expnt, double gp);
                                           ~AlignmentProposal(void);
        double                              propose(AlnMatrix* newAlignment, AlnMatrix* oldAlignment, double iP);
                
    private:
        void                                assembleNewAlignment(AlnMatrix* newAlignment, AlnMatrix* oldAlignment, int pos, int oldLen, int newLen);
        void                                buildChildIndices(std::vector<int>& parents, std::vector<int>& lftChildren, std::vector<int>& rhtChildren);
        void                                buildScoringMatrix(DoubleMatrix** tiProbs, std::vector<double>& freqs, int lftChild, int rhtChild);
        double                              calculateReverseProposal(AlnMatrix* aln, int pos, int len, DoubleMatrix** tiProbs, std::vector<double>& freqs, std::vector<int>& lftChildren, std::vector<int>& rhtChildren);
        void                                cleanTable(std::map<IntVector*,int,CompIntVector>& m);
        double                              computeProfileMatchScore(int lftIdx, int rhtIdx, int i, int j);
        void                                computeTracebackProbs(int lftIdx, int rhtIdx, int ii, int jj, double& probGapRight, double& probGapLeft, double& probMatch);
        int                                 countPaths(AlnMatrix* inputAlignment, int startCol, int endCol);
        void                                debugPrint(void);
        double                              deterministicTracebackLogProb(int nodeIdx, int lftIdx, int rhtIdx);
        void                                drainPool(void);
        void                                fillDpTable(int lftIdx, int rhtIdx);
        void                                freeProfile(int*** x, int n);
        void                                getIndelMatrix(AlnMatrix* inputAlignment, IndelMatrix* indelMat);
        int                                 getNumAllocated(void) { return (int) allocated.size(); }
        int                                 getNumInPool(void) { return (int)pool.size(); }
        IntVector*                          getVector(void);
        IntVector*                          getVector(IntVector& vec);
        void                                initialize(void);
        void                                initializeTipProfiles(AlnMatrix* aln, int pos, int len);
        void                                initializeTipProfilesWithGaps(AlnMatrix* aln, int pos, int len);
        void                                mergeDescendantProfiles(std::vector<int>& lftChildren, std::vector<int>& rhtChildren, AlnMatrix* aln, int len);
        void                                mergeProfilesAfterTraceback(int nodeIdx, int lftIdx, int rhtIdx);
        void                                print(std::string s, std::vector<int>& x);
        void                                print(std::string s, std::vector<std::vector<int> >& x);
        void                                print(std::string s, int*** x, int a, int b, int c);
        void                                removeAllGapColumns(AlnMatrix* aln, int len);
        void                                resetProfileDimensions(void);
        void                                returnToPool(IntVector* n);
        void                                selectRegionToRealign(int numSites, double extProb, int& pos, int& len);
        void                                sortChildrenByLeafOrder(std::vector<int>& lftChildren, std::vector<int>& rhtChildren, int nTaxa);
        double                              stochasticTraceback(int nodeIdx, int lftIdx, int rhtIdx);
        RandomVariable*                     rv;
        ParameterAlignment*                 alignmentParm;
        Model*                              model;
        RbBitSet                            taxonMask;            // bitmask identifying which taxa this alignment covers
        Tree*                               tree;
        int***                              profile;              // profile[nodeIdx][taxonInProfile][site]
        IndelMatrix*                        alignment;            // binary indel matrix (columns × taxa)
        AlnMatrix*                          profileNumber;        // maps profile row -> original taxon index
        int*                                xProfile;             // number of columns in each node's profile
        int*                                yProfile;             // number of taxa in each node's profile
        int                                 numStates;            // alphabet size (e.g. 4 for DNA)
        double**                            scoring;              // log-score matrix: scoring[i][j] = log Σ_k P(i|k)P(k|j)/π_j
        IntMatrix*                          tempProfile;          // scratch space for building merged profiles
        double**                            dp;                   // Needleman–Wunsch DP table
        double                              gap;                  // gap penalty for DP (must be < 0)
        double                              basis;                // Boltzmann temperature for stochastic traceback
        int                                 numTaxa;
        int                                 numNodes;             // total nodes (taxa + internals, including root)
        int                                 gapCode;              // integer code for the gap character
        int                                 maxLength;            // maximum number of alignment columns
        int                                 maxUnalignDimension;  // cap on compatible columns before giving up
        static int                          bigUnalignableRegion; // counter for regions that exceed maxUnalignDimension
        std::vector<int>                    possibles;            // compatible column indices for countPaths
        std::vector<int>                    state;                // column state labels for countPaths
        std::vector<IntVector*>             pool;
        std::set<IntVector*>                allocated;
        enum                                StateLabels { freeToUse, possible, edgeUsed, used };
};
#endif
