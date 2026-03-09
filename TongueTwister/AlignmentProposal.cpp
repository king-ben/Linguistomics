#include <cmath>
#include <iomanip>
#include <iostream>
#include "AlnMatrix.hpp"
#include "AlignmentProposal.hpp"
#include "IndelMatrix.hpp"
#include "IntVector.hpp"
#include "Model.hpp"
#include "Msg.hpp"
#include "ParameterAlignment.hpp"
#include "RandomVariable.hpp"
#include "TransitionProbabilities.hpp"
#include "Tree.hpp"

int AlignmentProposal::bigUnalignableRegion = 0;



AlignmentProposal::AlignmentProposal(ParameterAlignment* a, Tree* t, RandomVariable* r, Model* m, double expnt, double gp) {

    if (expnt < 1.0)
        Msg::error("Exponenet parameter must be greater than 1");
    if (gp >= 0.0)
        Msg::error("Gap penalty parameter must be less than 0");
        
    // store pointers to important objects
    alignmentParm = a;
    tree = t;
    rv = r;
    model = m;
    
    // initialize instance variables
    numTaxa = t->getNumTaxa();
    numNodes = (int)tree->getAncestorIndices().size();
    taxonMask = alignmentParm->getTaxonMask();
    gap = gp;
    basis = expnt;
    gapCode = alignmentParm->getGapCode();
    maxLength = 20;
    maxUnalignDimension = 17;
    
    // alignment:     IndelMatrix  (numSites × numTaxa) — binary indel pattern
    // profileNumber: AlnMatrix    (numNodes × numSites) — maps profile rows to taxon indices
    // tempProfile:   IntMatrix    (numNodes × maxLength) — scratch for merged profiles
    alignment = new IndelMatrix(numTaxa, alignmentParm->lengthOfLongestSequence()*5);
    profileNumber = new AlnMatrix(numNodes, alignmentParm->lengthOfLongestSequence()*5);
    tempProfile = new IntMatrix(numNodes, maxLength);
        
    // allocate the profile array: profile[nodeIndex][taxonInProfile][site]
    // each node gets a numNodes × maxLength block (worst case: root has numTaxa rows)
    profile = new int**[numNodes];
    for (int i=0; i<numNodes; i++)
        {
        profile[i] = new int*[numNodes];
        profile[i][0] = new int[numNodes*maxLength];
        for (int j=1; j<numNodes; j++)
            profile[i][j] = profile[i][j-1] + maxLength;
        for (int j = 0; j < numNodes; j++)
            memset(profile[i][j], 0, maxLength * sizeof(int));
        }
                
    // pre-size working vectors for countPaths (sized once, never reallocated)
    possibles.resize(maxUnalignDimension);
    state.resize(3*maxUnalignDimension);
    numStates = alignmentParm->getNumStates();

    // allocate Needleman–Wunsch DP table (maxLength × maxLength)
    dp = new double*[maxLength];
    dp[0] = new double[maxLength*maxLength];
    for (int i=1; i<maxLength; i++)
        dp[i] = dp[i-1] + maxLength;
    for (int i=0; i<maxLength; i++)
        for (int j=0; j<maxLength; j++)
            dp[i][j] = 0.0;
        
    // allocate scoring matrix (numStates × numStates)
    scoring = new double*[numStates];
    scoring[0] = new double[numStates*numStates];
    for (int i=1; i<numStates; i++)
        scoring[i] = scoring[i-1] + numStates;
    for (int i=0; i<numStates; i++)
        for (int j=0; j<numStates; j++)
            scoring[i][j] = 0.0;
            
    // allocate profile dimension arrays
    xProfile = new int[numNodes];
    yProfile = new int[numNodes];
    for (int i=0; i<numNodes; i++)
        {
        xProfile[i] = 0;
        yProfile[i] = 0;
        }
}

AlignmentProposal::~AlignmentProposal(void) {

    freeProfile(profile, numNodes);
    
    delete [] scoring[0];
    delete [] scoring;
    delete [] xProfile;
    delete [] yProfile;
    
    delete alignment;
    delete profileNumber;
    delete tempProfile;
}

/* Assemble the new alignment from three pieces: the prefix (columns before
   the realigned region), the newly aligned region from the root profile,
   and the suffix (columns after the old region). The root is always
   numNodes-1 in this index-based scheme. */
void AlignmentProposal::assembleNewAlignment(AlnMatrix* newAlignment, AlnMatrix* oldAlignment, int pos, int oldLen, int newLen) {

    int rootIdx = numNodes - 1;
    newAlignment->setNumSites(oldAlignment->getNumSites() + newLen - oldLen);
    
    // copy prefix (columns before the realigned region)
    for (int i=0; i<pos; i++)
        for (int j=0; j<numTaxa; j++)
            (*newAlignment)(j,i) = (*oldAlignment)(j,i);
    
    // copy newly aligned columns from the root profile,
    // using profileNumber to map profile rows back to original taxon indices
    for (int i=pos; i<pos+newLen; i++)
        for (int j=0; j<numTaxa; j++)
            (*newAlignment)((*profileNumber)(rootIdx,j), i) = profile[rootIdx][j][i-pos];
    
    // copy suffix (columns after the old region)
    for (int i=pos+oldLen; i<oldAlignment->getNumSites(); i++)
        for (int j=0; j<numTaxa; j++)
            (*newAlignment)(j, i+newLen-oldLen) = (*oldAlignment)(j,i);
}

/* Build left-child and right-child index arrays from the parent index
   array. For each non-root node i, parents[i] gives its parent. The
   first child encountered becomes the left child, the second becomes
   the right child. This relies on nodes being indexed such that
   children always have lower indices than their parents. */
void AlignmentProposal::buildChildIndices(std::vector<int>& parents, std::vector<int>& lftChildren, std::vector<int>& rhtChildren) {

    lftChildren.assign(numNodes, 0);
    rhtChildren.assign(numNodes, 0);
    for (int i=0; i<numNodes-1; i++)
        {
        if (lftChildren[parents[i]] == 0)
            lftChildren[parents[i]] = i;
        else
            rhtChildren[parents[i]] = i;
        }
}

/* Build the profile-profile scoring matrix from the current transition
   probability matrices. The score for aligning character i against
   character j is:
   
       scoring[i][j] = log( Σ_k P(i|k,t_L) · P(k|j,t_R) / π_j )
   
   where t_L and t_R are the branch lengths to the left and right children,
   and π_j is the stationary frequency. The tiProbs array is indexed by
   node index: tiProbs[lftChild] gives the matrix for the left branch. */
void AlignmentProposal::buildScoringMatrix(DoubleMatrix** tiProbs, std::vector<double>& freqs, int lftChild, int rhtChild) {

    for (int i=0; i<numStates; i++)
        {
        for (int j=0; j<numStates; j++)
            {
            scoring[i][j] = 0.0;
            for (int l=0; l<numStates; l++)
                scoring[i][j] += (*tiProbs[lftChild])(i,l) * (*tiProbs[rhtChild])(l,j);
            }
        }
    for (int i=0; i<numStates; i++)
        for (int j=0; j<numStates; j++)
            scoring[i][j] = log(scoring[i][j] / freqs[j]);
}

/* Calculate the log-probability of the reverse alignment proposal. Given
   the new alignment (after the forward proposal), this function asks:
   "what is the probability of recovering the old alignment from the new
   one?" This probability enters the Hastings ratio to ensure detailed
   balance.

   The calculation reconstructs the profile structure from the supplied
   alignment, strips all-gap columns, fills the DP table at each internal
   node, then performs a deterministic traceback (following the known gap
   pattern) to accumulate the log-probability of each traceback move. */
double AlignmentProposal::calculateReverseProposal(AlnMatrix* aln, int pos, int len, DoubleMatrix** tiProbs, std::vector<double>& freqs, std::vector<int>& lftChildren, std::vector<int>& rhtChildren) {

    double reverseLnProb = 0.0;
    
    // step 1: build tip profiles from the alignment (keeping gaps)
    initializeTipProfilesWithGaps(aln, pos, len);
    
    // step 2: merge descendant profiles at internal nodes
    mergeDescendantProfiles(lftChildren, rhtChildren, aln, len);
    
    // step 3: remove all-gap columns
    removeAllGapColumns(aln, len);
    
    // step 4: at each internal node, fill the DP table and compute the
    //         log-probability of the deterministic traceback path
    for (int k=aln->getNumTaxa(); k<numNodes; k++)
        {
        int lft = lftChildren[k];
        int rht = rhtChildren[k];
        
        buildScoringMatrix(tiProbs, freqs, lft, rht);
        fillDpTable(lft, rht);
        reverseLnProb += deterministicTracebackLogProb(k, lft, rht);
        }
    
    return reverseLnProb;
}

void AlignmentProposal::cleanTable(std::map<IntVector*, int, CompIntVector>& m) {

    for (std::map<IntVector*,int,CompIntVector>::iterator it = m.begin(); it != m.end(); it++)
        returnToPool(it->first);
    m.clear();
}

/* Compute the average pairwise alignment score between column i of the
   left profile and column j of the right profile. Only non-gap pairs
   contribute to the average; if all pairs involve at least one gap,
   returns 0.0 (neutral contribution to the DP). This averaging over
   profile rows implements a sum-of-pairs objective for the progressive
   alignment. */
double AlignmentProposal::computeProfileMatchScore(int lftIdx, int rhtIdx, int i, int j) {

    int count = 0;
    double score = 0.0;
    for (int l=0; l<yProfile[lftIdx]; l++)
        {
        for (int m=0; m<yProfile[rhtIdx]; m++)
            {
            if (profile[lftIdx][l][i] != gapCode &&
                profile[rhtIdx][m][j] != gapCode)
                {
                score += scoring[profile[lftIdx][l][i]]
                                [profile[rhtIdx][m][j]];
                count++;
                }
            }
        }
    if (count > 0)
        return score / count;
    return 0.0;
}

/* Compute the Boltzmann-weighted probabilities of the three possible
   traceback moves at position (ii, jj) in the DP table:
     - probGapRight: advance left profile only  (gap in right descendant)
     - probGapLeft:  advance right profile only  (gap in left descendant)
     - probMatch:    advance both profiles       (match/mismatch)
   
   Each move's score is exponentiated as basis^score (computed via
   exp(score * log(basis))), then normalized to form a probability
   distribution. Moves that are not available (ii==0 or jj==0) get
   probability 0. */
void AlignmentProposal::computeTracebackProbs(int lftIdx, int rhtIdx, int ii, int jj, double& probGapRight, double& probGapLeft, double& probMatch) {

    double logBasis = log(basis);
    
    double scoreGapRight = 0.0, scoreGapLeft = 0.0, scoreMatch = 0.0;
    
    if (ii > 0)
        scoreGapRight = dp[ii-1][jj] + gap;
    if (jj > 0)
        scoreGapLeft = dp[ii][jj-1] + gap;
    if (ii > 0 && jj > 0)
        scoreMatch = dp[ii-1][jj-1] + computeProfileMatchScore(lftIdx, rhtIdx, ii-1, jj-1);
    
    // Boltzmann weighting: convert DP scores to unnormalized probabilities
    double sum = 0.0;
    if (ii > 0)
        sum += exp(scoreGapRight * logBasis);
    if (jj > 0)
        sum += exp(scoreGapLeft * logBasis);
    if (ii > 0 && jj > 0)
        sum += exp(scoreMatch * logBasis);
    
    // normalize
    probGapRight = (ii > 0)           ? exp(scoreGapRight * logBasis) / sum : 0.0;
    probGapLeft  = (jj > 0)           ? exp(scoreGapLeft  * logBasis) / sum : 0.0;
    probMatch    = (ii > 0 && jj > 0) ? exp(scoreMatch    * logBasis) / sum : 0.0;
}

/* Count the number of distinct column orderings that yield the same
   alignment pattern in the region [startCol, endCol]. Two orderings
   are equivalent if they produce the same per-taxon character sequence
   (i.e., columns with non-overlapping taxon sets can be freely permuted).

   This is needed for the Hastings ratio because the proposal generates
   one specific column ordering, but multiple orderings yield the same
   alignment. The ratio of forward and reverse path counts enters the
   acceptance probability.

   The algorithm uses a DP table keyed by per-taxon position vectors
   (IntVector objects managed via the pool). It enumerates "compatible sets"
   of columns — columns whose non-gap taxa don't overlap — and counts
   combinations. Returns -1 if the number of compatible columns exceeds
   maxUnalignDimension (too complex to count). */
int AlignmentProposal::countPaths(AlnMatrix* inputAlignment, int startCol, int endCol) {
		
    int len = endCol - startCol + 1;

    // build the binary indel matrix for this alignment
    getIndelMatrix(inputAlignment, alignment);

    // seed the DP with an initial position vector (all zeros) and count = 1
    std::map<IntVector*,int,CompIntVector> dpTable;
    IntVector* pos = getVector();
    IntVector* newVec = getVector(*pos);
    dpTable.insert( std::make_pair( newVec, 1) );
		
    // reset working arrays
    for (int i=0; i<maxUnalignDimension; i++)
        possibles[i] = 0;
    state.resize(len);
    for (int i=0; i<len; i++)
        state[i] = 0;
    
    do
        {
        // a column is "compatible" if its non-gap taxa don't overlap with
        // any column already accumulated in the mask
        IntVector* mask = getVector();
        int ptr, numPossible = 0, firstNotUsed = 0;
        for (ptr = firstNotUsed; mask->zeroEntry() && ptr<len; ptr++)
            {
            if (state[ ptr ] != used)
                {
                if (mask->innerProduct( alignment->getRow(ptr), alignment->getNumTaxa() ) == 0)
                    {
                    state[ ptr ] = possible;
                    if (numPossible == maxUnalignDimension)
                        {
                        // too many compatible columns — bail out
                        bigUnalignableRegion++;
                        returnToPool(pos);
                        cleanTable(dpTable);
                        return -1;
                        }
                    possibles[numPossible++] = ptr;
                    }
                mask->add( alignment->getRow(ptr), alignment->getNumTaxa() );
                }
            }
        returnToPool(mask);
        
        // uses binary counting: for each subset, advance the position
        // vector and update the path count DP
        IntVector* newPos = getVector(*pos);
        IntVector* signature = getVector();
        int posPtr;
        bool unusedPos;
        bool foundNonZero;
        do
            {
            // find next combination via binary counting over possibles
            foundNonZero = false;
            for (posPtr = numPossible - 1; posPtr >= 0; --posPtr)
                {
                int curPtr = possibles[ posPtr ];
                if (state[ curPtr ] == possible)
                    {
                    // this digit goes 0->1: mark as edge-used and advance position
                    state[ curPtr ] = edgeUsed;
                    newPos->add( alignment->getRow(curPtr), alignment->getNumTaxa() );
                    signature->addMultiple( alignment->getRow(curPtr), alignment->getNumTaxa(), posPtr+1 );
                    foundNonZero = true;
                    posPtr = 0;  // break after finding the first unset digit
                    }
                else
                    {
                    // this digit goes 1->0: reset and continue (carry)
                    state[ curPtr ] = possible;
                    newPos->subtract( alignment->getRow(curPtr), alignment->getNumTaxa() );
                    signature->addMultiple( alignment->getRow(curPtr), alignment->getNumTaxa(), -posPtr-1 );
                    }
                }
        
            if (foundNonZero)
                {
                // look up the path count at the current position
                std::map<IntVector*, int, CompIntVector>::iterator it = dpTable.find(pos);
                if (it == dpTable.end())
                    Msg::error("Could not find pos in dpTable map.");
                int left = it->second;
                
                // look up the path count at the new position
                int right;
                it = dpTable.find(newPos);
                if (it == dpTable.end())
                    {
                    right = 0;
                    unusedPos = true;
                    }
                else
                    {
                    right = it->second;
                    unusedPos = false;
                    }
                right += left;

                // store the updated count —- use a fresh key for new positions
                if (unusedPos)
                    {
                    IntVector* v = getVector(*newPos);
                    dpTable.insert( std::make_pair(v, right) );
                    }
                else
                    {
                    std::map<IntVector*, int, CompIntVector>::iterator it2 = dpTable.find(newPos);
                    if (it2 == dpTable.end())
                        Msg::error("We should have found newPos in dpTable for the alignment proposal!");
                    it2->second = right;
                    }
                }
        
            } while (foundNonZero);
        
        // walk backwards through columns; the first column still in state
        // "possible" becomes the next "used" column
        --ptr;
        while (ptr >= 0 && state[ptr] != possible)
            {
            if (state[ptr] == used)
                {
                pos->subtract( alignment->getRow(ptr), alignment->getNumTaxa() );
                state[ptr] = freeToUse;
                }
            --ptr;
            }
        
        if (ptr == -1)
            {
            // we've exhausted all orderings; look up the terminal position
            std::map<IntVector*, int, CompIntVector>::iterator it = dpTable.find(newPos);
            if (it == dpTable.end())
                Msg::error("Could not find newPos in dpTable map.");
            int result = it->second;
            returnToPool(pos);
            returnToPool(newPos);
            returnToPool(signature);
            cleanTable(dpTable);
            return result;
            }
        else
            {
            returnToPool(newPos);
            returnToPool(signature);
            }
        
        // mark this column as used and advance the base position
        state[ptr] = used;
        pos->add( alignment->getRow(ptr), alignment->getNumTaxa() );
        if (ptr <= firstNotUsed)
            firstNotUsed++;
        
        } while (true);

    return 0;   // unreachable, but keeps the compiler happy
}

void AlignmentProposal::debugPrint(void) {

    std::cout << "xProfile: ";
    for (int i=0; i<numNodes; i++)
        std::cout << std::setw(3) << xProfile[i];
    std::cout << std::endl;

    std::cout << "yProfile: ";
    for (int i=0; i<numNodes; i++)
        std::cout << std::setw(3) << yProfile[i];
    std::cout << std::endl;

    for (int i=0; i<numNodes; i++)
        {
        for (int j=0; j<numNodes; j++)
            {
            std::cout << "(" << i << "," << j << "): ";
            for (int k=0; k<maxLength; k++)
                std::cout << std::setw(3) << profile[i][j][k];
            std::cout << std::endl;
            }
        }
}

/* Follow the existing gap pattern in the merged profile at nodeIdx to
   determine which traceback move was made at each step, and accumulate
   the log probability of each move. This gives the reverse proposal
   probability for the Hastings ratio.

   The logic at each column is:
     - If the left taxa are all gaps -> "gap in left" (advance right).
     - Else if the right taxa are all gaps -> "gap in right" (advance left).
     - Otherwise -> "match" (both advanced). */
double AlignmentProposal::deterministicTracebackLogProb(int nodeIdx, int lftIdx, int rhtIdx) {

    int i = xProfile[lftIdx];
    int j = xProfile[rhtIdx];
    int n = xProfile[nodeIdx];
    double logProb = 0.0;
    
    while (n > 0)
        {
        double probGapRight, probGapLeft, probMatch;
        computeTracebackProbs(lftIdx, rhtIdx, i, j, probGapRight, probGapLeft, probMatch);
        
        // determine which move was made by inspecting the gap pattern
        bool lftAllGap = true;
        for (int l=0; l<yProfile[lftIdx] && lftAllGap; l++)
            lftAllGap = (profile[nodeIdx][l][n-1] == gapCode);
        
        if (lftAllGap)
            {
            // left taxa all gaps -> this was a gap-in-left (advance right) move
            logProb += log(probGapLeft);
            j--;
            }
        else
            {
            bool rhtAllGap = true;
            for (int l=yProfile[lftIdx]; l<yProfile[lftIdx]+yProfile[rhtIdx] && rhtAllGap; l++)
                rhtAllGap = (profile[nodeIdx][l][n-1] == gapCode);
            
            if (rhtAllGap)
                {
                // right taxa all gaps -> gap-in-right (advance left) move
                logProb += log(probGapRight);
                i--;
                }
            else
                {
                // both sides have non-gaps -> match move
                logProb += log(probMatch);
                i--;
                j--;
                }
            }
        n--;
        }
    
    return logProb;
}

void AlignmentProposal::drainPool(void) {

    for (std::vector<IntVector*>::iterator v=pool.begin(); v != pool.end(); v++)
        {
        allocated.erase(*v);
        delete (*v);
        }
}

/* Fill the Needleman–Wunsch DP table for aligning the profiles of the
   left child (xProfile[lftIdx] columns) against the right child
   (xProfile[rhtIdx] columns). The recurrence is:
   
       dp[i][j] = max { dp[i-1][j]   + gap,              (gap in right)
                        dp[i][j-1]   + gap,              (gap in left)
                        dp[i-1][j-1] + matchScore(i,j) } (match)
   
   The match score is the average over all taxon pairs in the two profiles.
   The borders initialize with cumulative gap penalties. */
void AlignmentProposal::fillDpTable(int lftIdx, int rhtIdx) {

    // borders: cumulative gap penalties
    dp[0][0] = 0.0;
    for (int i=1; i<=xProfile[lftIdx]; i++)
        dp[i][0] = dp[i-1][0] + gap;
    for (int j=1; j<=xProfile[rhtIdx]; j++)
        dp[0][j] = dp[0][j-1] + gap;

    // interior: three-way max
    for (int i=1; i<=xProfile[lftIdx]; i++)
        {
        for (int j=1; j<=xProfile[rhtIdx]; j++)
            {
            double scoreGapRight = dp[i-1][j]   + gap;
            double scoreGapLeft  = dp[i][j-1]   + gap;
            double scoreMatch    = dp[i-1][j-1] + computeProfileMatchScore(lftIdx, rhtIdx, i-1, j-1);

            dp[i][j] = scoreGapRight;
            if (scoreGapLeft > dp[i][j])
                dp[i][j] = scoreGapLeft;
            if (scoreMatch > dp[i][j])
                dp[i][j] = scoreMatch;
            }
        }
}

void AlignmentProposal::freeProfile(int*** x, int n) {

    for (int i=0; i<n; i++)
        {
        delete [] x[i][0];
        delete [] x[i];
        }
    delete [] x;
}

/* Construct the binary indel matrix from the given alignment. The indel
   matrix has dimensions (numSites × numTaxa): entry (i,j) = 1 if taxon j
   has a non-gap character at column i, and 0 otherwise. */
void AlignmentProposal::getIndelMatrix(AlnMatrix* inputAlignment, IndelMatrix* indelMat) {

    int nt = inputAlignment->getNumTaxa();
    int ns = inputAlignment->getNumSites();
    if (indelMat->getNumTaxa() != nt)
        Msg::error("Mismatch in the number of taxa when initializing the indel matrix");
    indelMat->setNumSites(ns);
    for (int i=0; i<ns; i++)
        for (int j=0; j<nt; j++)
            {
            if ((*inputAlignment)(j,i) == gapCode)
                (*indelMat)(i,j) = 0;
            else
                (*indelMat)(i,j) = 1;
            }
}

IntVector* AlignmentProposal::getVector(void) {

    if (pool.empty() == true)
        {
        IntVector* v = new IntVector(numTaxa);
        allocated.insert(v);
        return v;
        }
    
    IntVector* v = pool.back();
    pool.pop_back();
    return v;
}

IntVector* AlignmentProposal::getVector(IntVector& vec) {

    if (pool.empty() == true)
        {
        IntVector* v = new IntVector(numTaxa);
        allocated.insert(v);
        for (int i=0; i<vec.size(); ++i)
            (*v)[i] = vec[i];
        return v;
        }
    IntVector* v = pool.back();
    pool.pop_back();
    for (int i=0; i<vec.size(); ++i)
        (*v)[i] = vec[i];
    return v;
}


/* Reinitialize state that may have changed since the last proposal:
   refresh the tree pointer, zero out all profiles, and verify that
   all IntVectors have been returned to the pool. */
void AlignmentProposal::initialize(void) {

    tree = model->getTree(taxonMask);
    
    for (int i=0; i<numNodes; i++)
        for (int j=0; j<numNodes; j++)
            for (int k=0; k<maxLength; k++)
                profile[i][j][k] = 0;
                
    // sanity check: all IntVectors should have been returned to the pool
    if (allocated.size() != pool.size())
        {
        std::cout << "allocated.size() = " << allocated.size() << std::endl;
        std::cout << "pool.size()      = " << pool.size() << std::endl;
        Msg::error("Expecting all IntVectors to have been returned to the pool");
        }
}

/* Extract ungapped characters from the selected region [pos, pos+len) for
   each leaf (taxa are indexed 0..numTaxa-1 in this index-based scheme).
   Gap characters are skipped so that xProfile[i] counts only the non-gap
   columns. Also initializes profileNumber so we can later reconstruct
   which taxon each profile row belongs to. */
void AlignmentProposal::initializeTipProfiles(AlnMatrix* aln, int pos, int len) {

    for (int i=0; i<aln->getNumTaxa(); i++)
        {
        yProfile[i] = 1;
        for (int j=pos; j<pos+len; j++)
            {
            if ((*aln)(i,j) != gapCode)
                {
                profile[i][0][xProfile[i]] = (*aln)(i,j);
                xProfile[i]++;
                }
            }
        (*profileNumber)(i,0) = i;
        }
}

/* Copy all characters (including gaps) from the selected region [pos, pos+len)
   for each leaf. Unlike initializeTipProfiles, this preserves gap positions
   because the reverse proposal calculation needs the complete aligned columns
   to determine which traceback path was taken. */
void AlignmentProposal::initializeTipProfilesWithGaps(AlnMatrix* aln, int pos, int len) {

    for (int i=0; i<aln->getNumTaxa(); i++)
        {
        yProfile[i] = 1;
        for (int j=pos; j<pos+len; j++)
            profile[i][0][j-pos] = (*aln)(i,j);
        }
}

/* After a stochastic traceback has built the merged profile in tempProfile
   (in reverse column order), this function:
     1. Sets yProfile for the parent to the sum of its children's counts.
     2. Copies profileNumber mappings from both descendants so we can later
        map profile rows back to the original taxon indices.
     3. Reverses the column order and copies from tempProfile into profile. */
void AlignmentProposal::mergeProfilesAfterTraceback(int nodeIdx, int lftIdx, int rhtIdx) {

    yProfile[nodeIdx] = yProfile[lftIdx] + yProfile[rhtIdx];

    // concatenate profile number mappings: left taxa first, then right
    for (int i=0; i<yProfile[lftIdx]; i++)
        (*profileNumber)(nodeIdx,i) = (*profileNumber)(lftIdx,i);
    for (int i=yProfile[lftIdx]; i<yProfile[lftIdx]+yProfile[rhtIdx]; i++)
        (*profileNumber)(nodeIdx,i) = (*profileNumber)(rhtIdx,i-yProfile[lftIdx]);

    // reverse the traceback output and copy into the profile array
    for (int i=0; i<xProfile[nodeIdx]; i++)
        for (int j=0; j<yProfile[nodeIdx]; j++)
            profile[nodeIdx][j][i] = (*tempProfile)(j, xProfile[nodeIdx]-i-1);
}

/* Build internal node profiles by concatenating the profiles of each
   node's two children column by column. Internal nodes are indexed
   from numTaxa to numNodes-1 (bottom-up by construction). */
void AlignmentProposal::mergeDescendantProfiles(std::vector<int>& lftChildren, std::vector<int>& rhtChildren, AlnMatrix* aln, int len) {

    for (int i=aln->getNumTaxa(); i<numNodes; i++)
        {
        int lft = lftChildren[i];
        int rht = rhtChildren[i];
        yProfile[i] = yProfile[lft] + yProfile[rht];
        for (int j=0; j<len; j++)
            {
            for (int k=0; k<yProfile[lft]; k++)
                profile[i][k][j] = profile[lft][k][j];
            for (int k=yProfile[lft]; k<yProfile[lft]+yProfile[rht]; k++)
                profile[i][k][j] = profile[rht][k-yProfile[lft]][j];
            }
        }
}

void AlignmentProposal::print(std::string s, std::vector<int>& x) {

    std::cout << s << std::endl;
    for (int i=0; i<x.size(); i++)
        std::cout << x[i] << " ";
    std::cout << std::endl;
}

void AlignmentProposal::print(std::string s, std::vector<std::vector<int> >& x) {

    std::cout << s << std::endl;
    for (int i=0; i<x.size(); i++)
        {
        for (int j=0; j<x[i].size(); j++)
            std::cout << x[i][j] << " ";
        std::cout << std::endl;
        }
}

void AlignmentProposal::print(std::string s, int*** x, int a, int b, int c) {

    std::cout << s << std::endl;
    for (int i=0; i<a; i++)
        {
        std::cout << i << ":" << std::endl;
        for (int j=0; j<b; j++)
            {
            for (int k=0; k<c; k++)
                std::cout << x[i][j][k] << " ";
            std::cout << std::endl;
            }
        }
}

/* The main alignment proposal:

   1. Initialize tree structure and transition probabilities.
   2. Select a contiguous region of geometric length.
   3. Extract ungapped tip profiles from the selected columns.
   4. Progressively align profiles bottom-up through the tree using a
      stochastic traceback (Boltzmann-weighted sampling from the DP table).
   5. Assemble the new alignment from prefix + new region + suffix.
   6. Retry if the alignment didn't change.
   7. Compute the Hastings ratio from forward/reverse proposal probabilities
      and the forward/reverse path counts. */
double AlignmentProposal::propose(AlnMatrix* newAlignment, AlnMatrix* oldAlignment, double extensionProb) {
    
    if (extensionProb <= 0.0 || extensionProb > 1.0)
        Msg::error("Extension parameter must be in the range (0,1]");
        
    // refresh tree pointer, zero profiles, check IntVector pool integrity
    initialize();
    
    // --- build tree traversal structure ---
    std::vector<int> parents = tree->getAncestorIndices();
    std::vector<int> lftChildrenIndices, rhtChildrenIndices;
    buildChildIndices(parents, lftChildrenIndices, rhtChildrenIndices);

    // get transition probabilities and stationary frequencies
    TransitionProbabilities* tip = model->getTransitionProbabilities();
    DoubleMatrix** tiProbs = tip->getTransitionProbabilities(taxonMask);
    std::vector<double>& stationaryFrequencies = tip->getStationaryFrequencies();

    // sort children so the left child always has the smaller minimum leaf index
    sortChildrenByLeafOrder(lftChildrenIndices, rhtChildrenIndices, newAlignment->getNumTaxa());

    // retry until alignment changes
    bool foundDifferentAlignment = false;
    int len = 0, len2 = 0, pos = 0;
    double proposalLoglikelihood = 0.0;
    
    do
        {
        // choose region to realign
        selectRegionToRealign(newAlignment->getNumSites(), extensionProb, pos, len);

        // reset profile dimensions for this attempt
        profileNumber->setNumSites(newAlignment->getNumTaxa());
        resetProfileDimensions();

        // initialize leaf profiles (ungapped characters only)
        initializeTipProfiles(newAlignment, pos, len);
        
        // clear scoring matrix before the forward pass
        for (int i=0; i<numStates; i++)
            for (int j=0; j<numStates; j++)
                scoring[i][j] = 0.0;

        proposalLoglikelihood = 0.0;
        
        // internal nodes are indexed from numTaxa to numNodes-1
        for (int k=newAlignment->getNumTaxa(); k<numNodes; k++)
            {
            int lftChild = lftChildrenIndices[k];
            int rhtChild = rhtChildrenIndices[k];
            
            // build the scoring matrix from transition probabilities
            buildScoringMatrix(tiProbs, stationaryFrequencies, lftChild, rhtChild);

            // fill the Needleman–Wunsch DP table
            fillDpTable(lftChild, rhtChild);

            // stochastic traceback: sample a new alignment of the two profiles
            // and accumulate the negative log proposal probability
            proposalLoglikelihood += stochasticTraceback(k, lftChild, rhtChild);
            
            // fix up the merged profile: set sizes, copy profile numbers, reverse
            mergeProfilesAfterTraceback(k, lftChild, rhtChild);
            }
        
        // the root profile now contains the newly aligned region
        len2 = xProfile[numNodes - 1];

        // assemble the new alignment from prefix + new region + suffix
        assembleNewAlignment(newAlignment, oldAlignment, pos, len, len2);

        // check if the alignment actually changed
        foundDifferentAlignment = !(*newAlignment == *oldAlignment);
        
        } while (foundDifferentAlignment == false);
        
    // how many column orderings produce each alignment?
    int iCF = countPaths(newAlignment, pos, pos+len2-1);
    if (iCF == -1)
        {
        freeProfile(profile, numNodes);
        Msg::error("Something is seriously wrong with the alignment");
        proposalLoglikelihood -= exp(100.0);
        iCF = 1;
        }
    proposalLoglikelihood -= log((double)iCF);

    iCF = countPaths(oldAlignment, pos, pos+len-1);
    if (iCF == -1)
        {
        freeProfile(profile, numNodes);
        Msg::error("Something is seriously wrong with the alignment");
        proposalLoglikelihood -= exp(100.0);
        iCF = 1;
        }
    proposalLoglikelihood += log((double)iCF);

    // reconstruct profiles from the new alignment and compute the probability
    // of recovering the old alignment via deterministic traceback
    proposalLoglikelihood += calculateReverseProposal(newAlignment, pos, len, tiProbs, stationaryFrequencies, lftChildrenIndices, rhtChildrenIndices);

    // when extensionProb < 1, account for the probability of drawing regions
    // of different lengths; when extensionProb = 1, this term vanishes
    return proposalLoglikelihood + (len - len2) * log(1.0 - extensionProb);
}

/* Remove columns in which every taxon has a gap character. After merging
   descendant profiles, some columns may be all-gaps at certain nodes.
   These carry no alignment information and must be stripped before the DP.
   
   After this function, xProfile[k] gives the number of retained (non-all-gap)
   columns for each node k. */
void AlignmentProposal::removeAllGapColumns(AlnMatrix*, int len) {

    for (int i=0; i<numNodes; i++)
        {
        // copy current profile to scratch space
        for (int j=0; j<len; j++)
            for (int k=0; k<yProfile[i]; k++)
                (*tempProfile)(k,j) = profile[i][k][j];
        
        // compact: keep only columns with at least one non-gap
        xProfile[i] = 0;
        for (int j=0; j<len; j++)
            {
            bool hasNonGap = false;
            for (int k=0; k<yProfile[i]; k++)
                {
                if ((*tempProfile)(k,j) != gapCode)
                    {
                    hasNonGap = true;
                    break;
                    }
                }
            if (hasNonGap)
                {
                for (int k=0; k<yProfile[i]; k++)
                    profile[i][k][xProfile[i]] = (*tempProfile)(k,j);
                xProfile[i]++;
                }
            }
        }
}

void AlignmentProposal::resetProfileDimensions(void) {

    for (int i=0; i<numNodes; i++)
        {
        xProfile[i] = 0;
        yProfile[i] = 0;
        }
}

void AlignmentProposal::returnToPool(IntVector* v) {

    v->clean();
    pool.push_back(v);
}

/* Choose a contiguous region of the alignment to realign. The region
   length follows a geometric distribution: starting at 1, each extension
   has probability extProb. The starting position is drawn uniformly. */
void AlignmentProposal::selectRegionToRealign(int numSites, double extProb, int& pos, int& len) {

    for (len=1; len<numSites && rv->uniformRv() < extProb; len++)
        ;
    pos = rv->uniformRvInt(numSites - len + 1);
}

/* Sort children so that the left child always has the smaller minimum
   leaf index. This ensures the profile row ordering is deterministic
   and consistent: the root's profile will list taxa in index order.
   If the left child has a larger sorter value, the children are swapped. */
void AlignmentProposal::sortChildrenByLeafOrder(std::vector<int>& lftChildren, std::vector<int>& rhtChildren, int nTaxa) {

    std::vector<int> sorter(numNodes);
    for (int i=0; i<nTaxa; i++)
        sorter[i] = i;
    for (int i=nTaxa; i<numNodes; i++)
        {
        if (sorter[lftChildren[i]] > sorter[rhtChildren[i]])
            sorter[i] = sorter[rhtChildren[i]];
        else
            sorter[i] = sorter[lftChildren[i]];
        
        // swap if left child has a larger minimum leaf index
        if (sorter[lftChildren[i]] > sorter[rhtChildren[i]])
            {
            int temp = lftChildren[i];
            lftChildren[i] = rhtChildren[i];
            rhtChildren[i] = temp;
            }
        }
}

/* Stochastically trace back through the DP table to produce a new
   alignment of the left and right child profiles. At each step, the
   three possible moves are weighted by their Boltzmann probabilities
   and one is sampled. The merged profile is built column-by-column
   in tempProfile (in reverse order — fixed by mergeProfilesAfterTraceback).

   Returns the negative log proposal probability of the sampled path,
   which is accumulated into the Hastings ratio. */
double AlignmentProposal::stochasticTraceback(int nodeIdx, int lftIdx, int rhtIdx) {

    int iI = xProfile[lftIdx];
    int jJ = xProfile[rhtIdx];
    xProfile[nodeIdx] = 0;
    double lnProb = 0.0;
    
    while (iI > 0 || jJ > 0)
        {
        double probGapRight, probGapLeft, probMatch;
        computeTracebackProbs(lftIdx, rhtIdx, iI, jJ, probGapRight, probGapLeft, probMatch);
        
        double u = rv->uniformRv();
        if (u < probGapRight)
            {
            // gap in right: left characters advance, right gets gaps
            // the probability is in the denominator, so we subtract
            lnProb -= log(probGapRight);
            for (int l=0; l<yProfile[lftIdx]; l++)
                (*tempProfile)(l,xProfile[nodeIdx]) = profile[lftIdx][l][iI-1];
            for (int l=yProfile[lftIdx]; l<yProfile[lftIdx]+yProfile[rhtIdx]; l++)
                (*tempProfile)(l,xProfile[nodeIdx]) = gapCode;
            xProfile[nodeIdx]++;
            iI--;
            }
        else if (u < probGapRight + probGapLeft)
            {
            // gap in left: right characters advance, left gets gaps
            lnProb -= log(probGapLeft);
            for (int l=0; l<yProfile[lftIdx]; l++)
                (*tempProfile)(l,xProfile[nodeIdx]) = gapCode;
            for (int l=yProfile[lftIdx]; l<yProfile[lftIdx]+yProfile[rhtIdx]; l++)
                (*tempProfile)(l,xProfile[nodeIdx]) = profile[rhtIdx][l-yProfile[lftIdx]][jJ-1];
            xProfile[nodeIdx]++;
            jJ--;
            }
        else
            {
            // match: both profiles advance
            if (probMatch <= 0.0)
                {
                freeProfile(profile, numNodes);
                Msg::error("Problem proposing new alignment");
                }
            lnProb -= log(probMatch);
            for (int l=0; l<yProfile[lftIdx]; l++)
                (*tempProfile)(l,xProfile[nodeIdx]) = profile[lftIdx][l][iI-1];
            for (int l=yProfile[lftIdx]; l<yProfile[lftIdx]+yProfile[rhtIdx]; l++)
                (*tempProfile)(l,xProfile[nodeIdx]) = profile[rhtIdx][l-yProfile[lftIdx]][jJ-1];
            xProfile[nodeIdx]++;
            jJ--;
            iI--;
            }
        }
    
    return lnProb;
}
