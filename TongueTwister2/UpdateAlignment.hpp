#ifndef UpdateAlignment_hpp
#define UpdateAlignment_hpp

#include "Update.hpp"
class Model;
class ParameterAlignment;
class ParameterFrequencies;
class ParameterTree;
class TransitionProbabilities;


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
class UpdateAlignment : public Update {

    public:
                                        UpdateAlignment(void) = delete;
                                        UpdateAlignment(Model* m, RandomVariable* r, ParameterAlignment* p);
                                       ~UpdateAlignment(void);
        size_t                          getUpdateIdx(void) { return 0; }
        uint64_t                        getUpdateId(void) { return updateInfo[0].updateHash; }
        UpdateType                      getUpdateType(void) { return updateInfo[0].updateType;}
        double                          getTuningParameter(void) { return updateInfo[0].tuningParameter; }
        std::string                     getUpdateName(void);
        std::string                     parameterType(void) { return "ParameterAlignment"; }
        void                            setDependants(void);
        double                          update(void);
        double                          update(double power);
        double                          updateFromPrior(void);
        
                                        /* Realign the entire alignment (extensionProb = 1.0).
                                           Intended for use by UpdateTopology to jointly update
                                           the tree topology and all alignments together.
                                           Returns the log proposal ratio for the alignment changes.
                                           Does NOT call setDependants() -- the caller manages that. */
        double                          realignFull(void);
        
                                        // accessor for the alignment parameter (needed by UpdateTopology)
        ParameterAlignment*             getAlignmentParameter(void) { return myParm; }
    
    private:
    
                                        // scoring and DP helpers
        void                            buildScoringMatrix(Node* lftChild, Node* rhtChild);
        double                          computeProfileMatchScore(int lftIdx, int rhtIdx, int i, int j);
        void                            fillDpTable(int lftIdx, int rhtIdx);
        void                            computeTracebackProbs(int lftIdx, int rhtIdx, int ii, int jj, double& probGapRight, double& probGapLeft, double& probMatch);
        
                                        // traceback methods
        double                          stochasticTraceback(int idx, int lftIdx, int rhtIdx);
        double                          deterministicTracebackLogProb(int idx, int lftIdx, int rhtIdx);
        void                            mergeProfilesAfterTraceback(int idx, int lftIdx, int rhtIdx);
    
                                        // profile initialization
        void                            initializeTipProfiles(const std::vector<Node*>& postOrder, Alignment* aln, int pos, int len);
        void                            initializeTipProfilesWithGaps(const std::vector<Node*>& postOrder, Alignment* aln, int pos, int len);
        void                            mergeDescendantProfiles(const std::vector<Node*>& postOrder, int len);
        void                            removeAllGapColumns(int len);
        
                                        // alignment assembly
        void                            assembleNewAlignment(Alignment* curAlignment, int rootIdx, int pos, int oldLen, int newLen, int numSegments);
                                        
                                        // Hastings ratio components
        double                          calculateReverseProposal(Alignment* aln, int pos, int len);
        int                             countPaths(Alignment* aln, int startCol, int endCol);
        void                            getIndelMatrix(Alignment* aln, int startCol, int len);
        
                                        // tree and proposal functions
        void                            initializeTreeStructure(void);
        double                          propose(void);
        double                          propose(double extProb, bool alignmentsMustBeDifferent);
        void                            selectRegionToRealign(int numSegments, double extProb, int& pos, int& len);
        void                            resetProfileDimensions(void);

        static constexpr double         basis = 1.5;         // Boltzmann temperature for stochastic traceback
        static constexpr double         gapPenalty = -5.0;   // affine gap open/extend penalty for DP
        double                          extensionProb;       // probability of extending region by one column
        
        ParameterAlignment*             myParm;
        ParameterFrequencies*           freqsParm;
        ParameterTree*                  treeParm;
        TransitionProbabilities*        tiProbs;
        Tree*                           tree;
        int***                          profile;             // profile[nodeIdx][taxonInProfile][site]
        int*                            xProfile;            // number of (non-gap) columns in each node's profile
        int*                            yProfile;            // number of taxa represented in each node's profile
        double**                        dp;                  // Needleman–Wunsch DP table
        double**                        scoring;             // log-score matrix: scoring[i][j] = log Σ_k P(i|k)P(k|j)/π_j
        int**                           indelMatrix;         // binary matrix: 1 if taxon j has a character at column i
        int*                            profileNumber;       // maps profile row -> original taxon index
        int**                           tempProfile;         // scratch space for building merged profiles
        Node**                          lftDescendants;      // left child of each internal node
        Node**                          rhtDescendants;      // right child of each internal node
        int*                            sorter;              // smallest leaf index below each node (for ordering)
        int*                            possibles;           // compatible column indices for countPaths
        int*                            state;               // column state labels for countPaths
        int*                            pathPos;             // current position per taxon in countPaths
        int*                            pathNewPos;          // candidate next position per taxon
        int*                            pathMask;            // running column-occupancy mask
        int*                            pathFinalPos;        // terminal position for each taxon
        int                             numTaxa;
        int                             numNodes;
        int                             numStates;           // alphabet size (e.g. 4 for DNA)
        int                             gapCode;             // integer code for the gap character (= numStates)
        int                             maxLength;           // maximum number of alignment columns
        int                             maxUnalignDimension; // cap on compatible columns before giving up
        unsigned                        taxonMask;           // bitmask identifying which taxa this alignment covers
        std::vector<int>                pathKey;             // reusable key for the dpTable map
        std::map<std::vector<int>, int> dpTable;             // memoization table for countPaths
                        
        enum StateLabels { freeToUse, possible, edgeUsed, used };
};

#endif 
