#include "Alignment.hpp"
#include "Model.hpp"
#include "Msg.hpp"
#include "Node.hpp"
#include "ParameterAlignment.hpp"
#include "ParameterFrequencies.hpp"
#include "ParameterTree.hpp"
#include "RandomVariable.hpp"
#include "TransitionProbabilities.hpp"
#include "Tree.hpp"
#include "UpdateAlignment.hpp"



/* Allocates all working memory for the progressive alignment proposal.
   Arrays are allocated as contiguous blocks where possible to improve
   cache locality. The profile array is three-dimensional:
     profile[nodeIndex][taxonInProfile][site]
   where "taxonInProfile" ranges over the taxa below that node. For a
   leaf this is just one row; for the root it has numTaxa rows. */
UpdateAlignment::UpdateAlignment(Model* m, RandomVariable* r, ParameterAlignment* p) : 
    Update(m, r), myParm(p) {
    
    // register this update's metadata for checkpoint matching
    updateInfo.resize(1);
    updateInfo[0].updateIdx = 0;
    updateInfo[0].updateName = getUpdateName();
    updateInfo[0].updateHash = hashUpdateName(getUpdateName());
    updateInfo[0].updateType = probability;
    updateInfo[0].tuningParameter = 0.5;
        
    extensionProb = updateInfo[0].tuningParameter;

    // locate model components (freqsParm may be NULL for equal-frequency models)
    freqsParm = model->findParameter<ParameterFrequencies>();

    treeParm = model->findParameter<ParameterTree>();
    if (treeParm == nullptr)
        Msg::error("Could not find tree parameter for the alignment update");

    tiProbs = model->getTiProbs();
    
    // get taxon mask and tree
    taxonMask = myParm->getTaxonMask();
    tree = treeParm->getTree(taxonMask);
    numTaxa = static_cast<int>(tree->getNumTaxa());
    numNodes = static_cast<int>(tree->getNumNodes());
    
    // get alignment information
    numStates = static_cast<int>(myParm->getNumStates());
    gapCode = numStates;
    
    // set maximum dimensions
    Alignment* aln = myParm->getAlignment(0);
    maxLength = static_cast<int>(aln->getMaximumNumberOfSegments());
    maxUnalignDimension = 17;
    
    // allocate profile array: profile[nodeIndex][taxonInProfile][site]
    // each node gets a numNodes × maxLength block (worst case: root has numTaxa rows)
    profile = new int**[numNodes];
    for (int i=0; i<numNodes; i++)
        {
        profile[i] = new int*[numNodes];
        profile[i][0] = new int[numNodes * maxLength];
        for (int j=1; j<numNodes; j++)
            profile[i][j] = profile[i][j-1] + maxLength;
        for (int j=0; j<numNodes; j++)
            memset(profile[i][j], 0, maxLength * sizeof(int));
        }
    
    // xProfile[i] = number of columns in node i's profile
    // yProfile[i] = number of taxa represented in node i's profile
    xProfile = new int[numNodes];
    yProfile = new int[numNodes];
    for (int i=0; i<numNodes; i++)
        {
        xProfile[i] = 0;
        yProfile[i] = 0;
        }
    
    // allocate Needleman–Wunsch DP table (maxLength × maxLength)
    dp = new double*[maxLength];
    dp[0] = new double[maxLength * maxLength];
    for (int i=1; i<maxLength; i++)
        dp[i] = dp[i-1] + maxLength;
    for (int i=0; i<maxLength; i++)
        for (int j=0; j<maxLength; j++)
            dp[i][j] = 0.0;
    
    // allocate scoring matrix (numStates × numStates) 
    // stores log-transformed pairwise substitution scores
    scoring = new double*[numStates];
    scoring[0] = new double[numStates * numStates];
    for (int i=1; i<numStates; i++)
        scoring[i] = scoring[i-1] + numStates;
    for (int i=0; i<numStates; i++)
        for (int j=0; j<numStates; j++)
            scoring[i][j] = 0.0;
    
    // allocate helper arrays as contiguous blocks 
    
    // indelMatrix[column][taxon]: binary indicator of character presence
    indelMatrix = new int*[maxLength];
    indelMatrix[0] = new int[maxLength * numTaxa];
    for (int i=1; i<maxLength; i++)
        indelMatrix[i] = indelMatrix[i-1] + numTaxa;
    for (int i=0; i<maxLength; i++)
        for (int j=0; j<numTaxa; j++)
            indelMatrix[i][j] = 0;
    
    // profileNumber: maps (nodeIdx * numNodes + profileRow) -> original taxon index
    profileNumber = new int[numNodes * numNodes];
    for (int i=0; i<numNodes*numNodes; i++)
        profileNumber[i] = 0;
    
    // tempProfile: scratch space for building merged profiles during traceback
    tempProfile = new int*[numNodes];
    tempProfile[0] = new int[numNodes * maxLength];
    for (int i=1; i<numNodes; i++)
        tempProfile[i] = tempProfile[i-1] + maxLength;
    for (int i=0; i<numNodes; i++)
        for (int j=0; j<maxLength; j++)
            tempProfile[i][j] = 0;
    
    // descendant pointers and sorter for deterministic left-right ordering
    lftDescendants = new Node*[numNodes];
    rhtDescendants = new Node*[numNodes];
    sorter = new int[numNodes];
    for (int i=0; i<numNodes; i++)
        {
        lftDescendants[i] = nullptr;
        rhtDescendants[i] = nullptr;
        sorter[i] = 0;
        }
    
    // allocate arrays for countPaths
    possibles = new int[maxUnalignDimension];
    for (int i=0; i<maxUnalignDimension; i++)
        possibles[i] = 0;
    
    state = new int[maxLength];
    for (int i=0; i<maxLength; i++)
        state[i] = 0;
    
    pathPos = new int[numTaxa];
    pathNewPos = new int[numTaxa];
    pathMask = new int[numTaxa];
    pathFinalPos = new int[numTaxa];
    pathKey.resize(numTaxa);
    for (int i=0; i<numTaxa; i++)
        {
        pathPos[i] = 0;
        pathNewPos[i] = 0;
        pathMask[i] = 0;
        pathFinalPos[i] = 0;
        pathKey[i] = 0;
        }
}

UpdateAlignment::~UpdateAlignment(void) {

    for (int i=0; i<numNodes; i++)
        {
        delete [] profile[i][0];
        delete [] profile[i];
        }
    delete [] profile;
    
    delete [] xProfile;
    delete [] yProfile;
    
    delete [] dp[0];
    delete [] dp;
    
    delete [] scoring[0];
    delete [] scoring;
    
    delete [] indelMatrix[0];
    delete [] indelMatrix;
    
    delete [] profileNumber;
    
    delete [] tempProfile[0];
    delete [] tempProfile;
    
    delete [] lftDescendants;
    delete [] rhtDescendants;
    delete [] sorter;
    
    delete [] possibles;
    delete [] state;
    delete [] pathPos;
    delete [] pathNewPos;
    delete [] pathMask;
    delete [] pathFinalPos;
}

/* Build the profile-profile scoring matrix from the current transition
   probability matrices. The score for aligning character i against
   character j is:
   
       scoring[i][j] = log( Σ_k P(i|k,t_L) · P(k|j,t_R) / π_j )
   
   where t_L and t_R are the branch lengths to the left and right
   descendants, and π_j is the stationary frequency. This is the
   log-likelihood of observing characters i and j at the tips of a
   two-branch subtree, normalized by the equilibrium frequency. */
void UpdateAlignment::buildScoringMatrix(Node* lftDescendant, Node* rhtDescendant) {

    // get transition probability matrices for the two branches
    DoubleMatrix& tiL = tiProbs->getTransitionProbability(lftDescendant->getBranchLength());
    DoubleMatrix& tiR = tiProbs->getTransitionProbability(rhtDescendant->getBranchLength());
    
    // compute the raw substitution scores: Σ_k P(i|k) · P(k|j)
    for (int i=0; i<numStates; i++)
        {
        for (int j=0; j<numStates; j++)
            {
            scoring[i][j] = 0.0;
            for (int k=0; k<numStates; k++)
                scoring[i][j] += tiL(i,k) * tiR(k,j);
            }
        }
    
    // divide by stationary frequencies and take the log
    if (freqsParm != nullptr)
        {
        std::vector<double>& freqs = freqsParm->getFrequencies();
        for (int i=0; i<numStates; i++)
            for (int j=0; j<numStates; j++)
                scoring[i][j] = log(scoring[i][j] / freqs[j]);
        }
    else
        {
        // equal frequencies: π_j = 1/numStates for all j
        double eqFreq = 1.0 / numStates;
        for (int i=0; i<numStates; i++)
            for (int j=0; j<numStates; j++)
                scoring[i][j] = log(scoring[i][j] / eqFreq);
        }
}

/* Compute the average pairwise alignment score between column i of the
   left profile and column j of the right profile. Only non-gap pairs
   contribute to the average; if all pairs involve at least one gap,
   returns 0.0 (neutral contribution to the DP). This averaging over
   profile rows implements a sum-of-pairs objective for the progressive
   alignment. */
double UpdateAlignment::computeProfileMatchScore(int lftIdx, int rhtIdx, int i, int j) {

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

/* Fill the Needleman–Wunsch DP table for aligning the profiles of the
   left child (xProfile[lftIdx] columns) against the right child
   (xProfile[rhtIdx] columns). The recurrence is:
   
       dp[i][j] = max { dp[i-1][j]   + gapPenalty,       (gap in right)
                        dp[i][j-1]   + gapPenalty,       (gap in left)
                        dp[i-1][j-1] + matchScore(i,j) } (match)
   
   The match score is the average over all taxon pairs in the two profiles
   (computed by computeProfileMatchScore). The borders initialize with
   cumulative gap penalties. */
void UpdateAlignment::fillDpTable(int lftIdx, int rhtIdx) {

    // borders: cumulative gap penalties
    dp[0][0] = 0.0;
    for (int i=1; i<=xProfile[lftIdx]; i++)
        dp[i][0] = dp[i-1][0] + gapPenalty;
    for (int j=1; j<=xProfile[rhtIdx]; j++)
        dp[0][j] = dp[0][j-1] + gapPenalty;
    
    // interior: three-way max
    for (int i=1; i<=xProfile[lftIdx]; i++)
        {
        for (int j=1; j<=xProfile[rhtIdx]; j++)
            {
            double scoreGapRight = dp[i-1][j]   + gapPenalty;
            double scoreGapLeft  = dp[i][j-1]   + gapPenalty;
            double scoreMatch    = dp[i-1][j-1] + computeProfileMatchScore(lftIdx, rhtIdx, i-1, j-1);
            
            dp[i][j] = scoreGapRight;
            if (scoreGapLeft > dp[i][j])
                dp[i][j] = scoreGapLeft;
            if (scoreMatch > dp[i][j])
                dp[i][j] = scoreMatch;
            }
        }
}

/* Compute the Boltzmann-weighted probabilities of the three possible
   traceback moves at position (ii, jj) in the DP table:
     - probGapRight: advance left profile only   (gap in right descendant)
     - probGapLeft:  advance right profile only  (gap in left descendant)
     - probMatch:    advance both profiles       (match/mismatch)
   
   Each move's score is exponentiated as basis^score (computed via
   exp(score * log(basis))), then normalized to form a probability
   distribution. Moves that are not available (ii==0 or jj==0) get
   probability 0. */
void UpdateAlignment::computeTracebackProbs(int lftIdx, int rhtIdx, int ii, int jj, double& probGapRight, double& probGapLeft, double& probMatch) {

    static const double logBasis = log(basis);
    
    double scoreGapRight = 0.0, scoreGapLeft = 0.0, scoreMatch = 0.0;
    
    if (ii > 0)
        scoreGapRight = dp[ii-1][jj] + gapPenalty;
    if (jj > 0)
        scoreGapLeft = dp[ii][jj-1] + gapPenalty;
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

/* Stochastically trace back through the DP table to produce a new
   alignment of the left and right child profiles. At each step, the
   three possible moves are weighted by their Boltzmann probabilities
   and one is sampled. The merged profile is built column-by-column
   in tempProfile (in reverse order —- fixed by mergeProfilesAfterTraceback).

   Returns the negative log proposal probability of the sampled path,
   which is accumulated into the Hastings ratio. */
double UpdateAlignment::stochasticTraceback(int idx, int lftIdx, int rhtIdx) {

    int ii = xProfile[lftIdx];
    int jj = xProfile[rhtIdx];
    xProfile[idx] = 0;
    double lnProb = 0.0;
    
    while (ii > 0 || jj > 0)
        {
        double probGapRight, probGapLeft, probMatch;
        computeTracebackProbs(lftIdx, rhtIdx, ii, jj, probGapRight, probGapLeft, probMatch);
        
        double u = rng->uniformRv();
        if (u < probGapRight)
            {
            // gap in right: left characters advance, right gets gaps
            lnProb -= log(probGapRight);
            for (int l=0; l<yProfile[lftIdx]; l++)
                tempProfile[l][xProfile[idx]] = profile[lftIdx][l][ii-1];
            for (int l=yProfile[lftIdx]; l<yProfile[lftIdx]+yProfile[rhtIdx]; l++)
                tempProfile[l][xProfile[idx]] = gapCode;
            xProfile[idx]++;
            ii--;
            }
        else if (u < probGapRight + probGapLeft)
            {
            // gap in left: right characters advance, left gets gaps
            lnProb -= log(probGapLeft);
            for (int l=0; l<yProfile[lftIdx]; l++)
                tempProfile[l][xProfile[idx]] = gapCode;
            for (int l=yProfile[lftIdx]; l<yProfile[lftIdx]+yProfile[rhtIdx]; l++)
                tempProfile[l][xProfile[idx]] = profile[rhtIdx][l-yProfile[lftIdx]][jj-1];
            xProfile[idx]++;
            jj--;
            }
        else
            {
            // match: both profiles advance
            if (probMatch <= 0.0)
                Msg::error("Problem in alignment proposal: probMatch <= 0");
            lnProb -= log(probMatch);
            for (int l=0; l<yProfile[lftIdx]; l++)
                tempProfile[l][xProfile[idx]] = profile[lftIdx][l][ii-1];
            for (int l=yProfile[lftIdx]; l<yProfile[lftIdx]+yProfile[rhtIdx]; l++)
                tempProfile[l][xProfile[idx]] = profile[rhtIdx][l-yProfile[lftIdx]][jj-1];
            xProfile[idx]++;
            jj--;
            ii--;
            }
        }
    
    return lnProb;
}

/* After a stochastic traceback has built the merged profile in tempProfile
   (in reverse column order), this function:
     1. Sets yProfile for the parent to the sum of its children's counts.
     2. Copies profileNumber mappings from both descendants so we can later
        map profile rows back to the original taxon indices.
     3. Reverses the column order and copies from tempProfile into profile. */
void UpdateAlignment::mergeProfilesAfterTraceback(int idx, int lftIdx, int rhtIdx) {

    yProfile[idx] = yProfile[lftIdx] + yProfile[rhtIdx];
    
    // concatenate profile number mappings: left taxa first, then right
    for (int i=0; i<yProfile[lftIdx]; i++)
        profileNumber[idx*numNodes + i] = profileNumber[lftIdx*numNodes + i];
    for (int i=yProfile[lftIdx]; i<yProfile[lftIdx]+yProfile[rhtIdx]; i++)
        profileNumber[idx*numNodes + i] = profileNumber[rhtIdx*numNodes + i-yProfile[lftIdx]];
    
    // reverse the traceback output and copy into the profile array
    for (int i=0; i<xProfile[idx]; i++)
        for (int j=0; j<yProfile[idx]; j++)
            profile[idx][j][i] = tempProfile[j][xProfile[idx]-i-1];
}

/* Follow the existing gap pattern in the merged profile at node idx to
   determine which traceback move was made at each step, and accumulate
   the log probability of each move. This gives the reverse proposal
   probability for the Hastings ratio.

   The logic at each column is:
     - If the left taxa are all gaps -> the move was "gap in left" (advance right).
     - Else if the right taxa are all gaps -> "gap in right" (advance left).
     - Otherwise -> "match" (both advanced). */
double UpdateAlignment::deterministicTracebackLogProb(int idx, int lftIdx, int rhtIdx) {

    int ii = xProfile[lftIdx];
    int jj = xProfile[rhtIdx];
    int n = xProfile[idx];
    double logProb = 0.0;
    
    while (n > 0)
        {
        double probGapRight, probGapLeft, probMatch;
        computeTracebackProbs(lftIdx, rhtIdx, ii, jj, probGapRight, probGapLeft, probMatch);
        
        // determine which move was made by inspecting the gap pattern
        bool lftAllGap = true;
        for (int l=0; l<yProfile[lftIdx]; l++)
            {
            if (profile[idx][l][n-1] != gapCode)
                {
                lftAllGap = false;
                break;
                }
            }
        
        if (lftAllGap)
            {
            // left taxa all gaps -> this was a gap-in-left (advance right) move
            logProb += log(probGapLeft);
            jj--;
            }
        else
            {
            bool rhtAllGap = true;
            for (int l=yProfile[lftIdx]; l<yProfile[lftIdx]+yProfile[rhtIdx]; l++)
                {
                if (profile[idx][l][n-1] != gapCode)
                    {
                    rhtAllGap = false;
                    break;
                    }
                }
            
            if (rhtAllGap)
                {
                // right taxa all gaps -> gap-in-right (advance left) move
                logProb += log(probGapRight);
                ii--;
                }
            else
                {
                // both sides have non-gaps -> match move
                logProb += log(probMatch);
                ii--;
                jj--;
                }
            }
        n--;
        }
    
    return logProb;
}

/* Extract ungapped characters from the selected region [pos, pos+len) for
   each leaf node. Gap characters are skipped so that xProfile[leafIdx]
   counts only the non-gap columns. Also initializes profileNumber so
   we can later reconstruct which taxon each profile row belongs to. */
void UpdateAlignment::initializeTipProfiles(const std::vector<Node*>& postOrder, Alignment* aln, int pos, int len) {

    for (Node* p : postOrder)
        {
        if (p->getIsLeaf() == true)
            {
            int idx = p->getIndex();
            yProfile[idx] = 1;
            for (int j=pos; j<pos+len; j++)
                {
                int c = (*aln)(idx, j);
                if (c != gapCode)
                    {
                    profile[idx][0][xProfile[idx]] = c;
                    xProfile[idx]++;
                    }
                }
            profileNumber[idx * numNodes] = idx;
            }
        }
}

/* Copy all characters (including gaps) from the selected region [pos, pos+len)
   for each leaf node. Unlike initializeTipProfiles, this preserves gap
   positions because the reverse proposal calculation needs the complete
   aligned columns to determine which traceback path was taken. */
void UpdateAlignment::initializeTipProfilesWithGaps(const std::vector<Node*>& postOrder, Alignment* aln, int pos, int len) {

    for (Node* p : postOrder)
        {
        if (p->getIsLeaf() == true)
            {
            int idx = p->getIndex();
            yProfile[idx] = 1;
            for (int j=pos; j<pos+len; j++)
                profile[idx][0][j-pos] = (*aln)(idx, j);
            }
        }
}

/* Build internal node profiles by concatenating the profiles of each
   node's two descendants column by column. This reconstructs the complete
   profile structure from the existing alignment — needed for the reverse
   proposal to determine what traceback path produced the current state. */
void UpdateAlignment::mergeDescendantProfiles(const std::vector<Node*>& postOrder, int len) {

    for (Node* p : postOrder)
        {
        if (p->getIsLeaf() == true)
            continue;
        
        int idx = p->getIndex();
        Node* lft = lftDescendants[idx];
        Node* rht = rhtDescendants[idx];
        if (lft == nullptr || rht == nullptr)
            continue;
        
        int lftIdx = lft->getIndex();
        int rhtIdx = rht->getIndex();
        
        yProfile[idx] = yProfile[lftIdx] + yProfile[rhtIdx];
        
        for (int j=0; j<len; j++)
            {
            for (int l=0; l<yProfile[lftIdx]; l++)
                profile[idx][l][j] = profile[lftIdx][l][j];
            for (int l=yProfile[lftIdx]; l<yProfile[lftIdx]+yProfile[rhtIdx]; l++)
                profile[idx][l][j] = profile[rhtIdx][l-yProfile[lftIdx]][j];
            }
        }
}

/* Remove columns in which every taxon has a gap character. After merging
   descendant profiles, some columns may be all-gaps at certain nodes
   (because the descendants' gap patterns don't overlap). These columns
   carry no alignment information and must be stripped before the DP.
   
   After this function, xProfile[k] gives the number of retained (non-all-gap)
   columns for each node k. */
void UpdateAlignment::removeAllGapColumns(int len) {

    for (int k=0; k<numNodes; k++)
        {
        // copy current profile to scratch space
        for (int j=0; j<len; j++)
            for (int l=0; l<yProfile[k]; l++)
                tempProfile[l][j] = profile[k][l][j];
        
        // compact: keep only columns with at least one non-gap
        xProfile[k] = 0;
        for (int j=0; j<len; j++)
            {
            bool hasNonGap = false;
            for (int l=0; l<yProfile[k]; l++)
                {
                if (tempProfile[l][j] != gapCode)
                    {
                    hasNonGap = true;
                    break;
                    }
                }
            if (hasNonGap)
                {
                for (int l=0; l<yProfile[k]; l++)
                    profile[k][l][xProfile[k]] = tempProfile[l][j];
                xProfile[k]++;
                }
            }
        }
}

/* Splice the newly aligned region back into the full alignment. The region
   [pos, pos+oldLen) in the current alignment is replaced with the newLen
   columns from the root profile. If the region changed length, the suffix
   (columns after the old region) is shifted left or right accordingly. */
void UpdateAlignment::assembleNewAlignment(Alignment* curAlignment, int rootIdx, int pos, int oldLen, int newLen, int numSegments) {

    int newNumSegments = numSegments + newLen - oldLen;
    
    // shift suffix columns to accommodate the new region length
    if (newLen < oldLen)
        {
        // region shrank: shift suffix left (forward direction to avoid overwriting)
        for (int i=pos+oldLen; i<numSegments; i++)
            for (int j=0; j<numTaxa; j++)
                (*curAlignment)(j, i+newLen-oldLen) = (*curAlignment)(j, i);
        }
    else if (newLen > oldLen)
        {
        // region grew: shift suffix right (backward direction to avoid overwriting)
        for (int i=numSegments-1; i>=pos+oldLen; i--)
            for (int j=0; j<numTaxa; j++)
                (*curAlignment)(j, i+newLen-oldLen) = (*curAlignment)(j, i);
        }
    
    curAlignment->setNumSegments(newNumSegments);
    
    // copy newly aligned columns from the root profile into the alignment,
    // using profileNumber to map profile rows back to original taxon indices
    for (int i=pos; i<pos+newLen; i++)
        {
        for (int j=0; j<numTaxa; j++)
            {
            int taxonIdx = profileNumber[rootIdx*numNodes + j];
            (*curAlignment)(taxonIdx, i) = profile[rootIdx][j][i-pos];
            }
        }
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
double UpdateAlignment::calculateReverseProposal(Alignment* aln, int pos, int len) {

    const std::vector<Node*>& postOrder = tree->getPostOrder();
    
    double reverseLnProb = 0.0;
    
    // step 1: build tip profiles from the alignment (keeping gaps)
    initializeTipProfilesWithGaps(postOrder, aln, pos, len);
    
    // step 2: merge descendant profiles at internal nodes
    mergeDescendantProfiles(postOrder, len);
    
    // step 3: remove all-gap columns (they don't participate in the DP)
    removeAllGapColumns(len);
    
    // step 4: at each internal node, fill the DP table and compute the
    //         log-probability of the deterministic traceback path
    for (Node* p : postOrder)
        {
        if (p->getIsLeaf() == true)
            continue;
        
        int idx = p->getIndex();
        Node* lft = lftDescendants[idx];
        Node* rht = rhtDescendants[idx];
        if (lft == nullptr || rht == nullptr)
            continue;
        
        int lftIdx = lft->getIndex();
        int rhtIdx = rht->getIndex();
        
        buildScoringMatrix(lft, rht);
        fillDpTable(lftIdx, rhtIdx);
        reverseLnProb += deterministicTracebackLogProb(idx, lftIdx, rhtIdx);
        }
    
    return reverseLnProb;
}

/* Construct the binary indel matrix for columns [startCol, startCol+len).
   indelMatrix[i][j] = 1 if taxon j has a non-gap character at column
   (startCol+i), and 0 otherwise. This matrix is used by countPaths to
   determine how many distinct column orderings produce the same alignment. */
void UpdateAlignment::getIndelMatrix(Alignment* aln, int startCol, int len) {

    for (int i=0; i<len; i++)
        {
        for (int j=0; j<numTaxa; j++)
            {
            if ((*aln)(j, startCol+i) == gapCode)
                indelMatrix[i][j] = 0;
            else
                indelMatrix[i][j] = 1;
            }
        }
}

/* Count the number of distinct column orderings that yield the same
   alignment pattern in the region [startCol, endCol]. Two orderings
   are equivalent if they produce the same per-taxon character sequence
   (i.e., columns with non-overlapping taxon sets can be freely permuted).

   This is needed for the Hastings ratio because the proposal generates
   one specific column ordering, but multiple orderings yield the same
   alignment. The ratio of forward and reverse path counts enters the
   acceptance probability.

   The algorithm enumerates "compatible sets" of columns -- columns whose
   non-gap taxa don't overlap -- and counts combinations using a DP table
   keyed by the per-taxon position vector. Returns -1 if the number of
   compatible columns exceeds maxUnalignDimension (too complex to count). */
int UpdateAlignment::countPaths(Alignment* aln, int startCol, int endCol) {

    int len = endCol - startCol + 1;
    if (len <= 0)
        return 1;
    
    getIndelMatrix(aln, startCol, len);
    
    // clear the persistent dpTable for this invocation
    dpTable.clear();
    
    // reset pre-allocated working arrays
    for (int i=0; i<numTaxa; i++)
        {
        pathPos[i] = 0;
        pathNewPos[i] = 0;
        pathMask[i] = 0;
        pathFinalPos[i] = 0;
        }
    
    for (int i=0; i<maxUnalignDimension; i++)
        possibles[i] = 0;
    for (int i=0; i<len; i++)
        state[i] = freeToUse;
    
    // seed the DP with the initial position (all zeros)
    for (int i=0; i<numTaxa; i++)
        pathKey[i] = pathPos[i];
    dpTable[pathKey] = 1;
    
    int firstNotUsed = 0;
    
    while (true)
        {
        // find compatible columns from the current position:
        // a column is "compatible" if its non-gap taxa don't overlap with
        // any column already accumulated in pathMask
        for (int j=0; j<numTaxa; j++)
            pathMask[j] = 0;
        int numPossible = 0;
        int ptr = 0;
        
        for (ptr=firstNotUsed; ptr<len; ptr++)
            {
            if (state[ptr] != used)
                {
                bool compatible = true;
                for (int j=0; j<numTaxa; j++)
                    {
                    if (pathMask[j] != 0 && indelMatrix[ptr][j] != 0)
                        {
                        compatible = false;
                        break;
                        }
                    }
                
                if (compatible)
                    {
                    state[ptr] = possible;
                    if (numPossible >= maxUnalignDimension)
                        return -1;     // too many compatible columns —- bail out
                    possibles[numPossible++] = ptr;
                    }
                
                // accumulate occupancy mask regardless of compatibility
                for (int j=0; j<numTaxa; j++)
                    pathMask[j] += indelMatrix[ptr][j];
                
                // once every taxon is covered, no more columns can be compatible
                bool allNonzero = true;
                for (int j=0; j<numTaxa; j++)
                    {
                    if (pathMask[j] == 0)
                        {
                        allNonzero = false;
                        break;
                        }
                    }
                if (allNonzero)
                    break;
                }
            }
        
        // enumerate all non-empty subsets of compatible columns:
        // for each subset, update the path count DP
        for (int combo=1; combo<(1<<numPossible); combo++)
            {
            // compute the new position by advancing through the selected columns
            for (int j=0; j<numTaxa; j++)
                pathNewPos[j] = pathPos[j];
            for (int p=0; p<numPossible; p++)
                {
                if (combo & (1<<p))
                    {
                    int col = possibles[p];
                    for (int j=0; j<numTaxa; j++)
                        pathNewPos[j] += indelMatrix[col][j];
                    }
                }
            
            // look up the count for the current position
            for (int i=0; i<numTaxa; i++)
                pathKey[i] = pathPos[i];
            std::map<std::vector<int>, int>::iterator it = dpTable.find(pathKey);
            int left = (it != dpTable.end()) ? it->second : 0;
            
            // look up the count for the new position (may not exist yet)
            for (int i=0; i<numTaxa; i++)
                pathKey[i] = pathNewPos[i];
            it = dpTable.find(pathKey);
            int right = (it != dpTable.end()) ? it->second : 0;
            
            // accumulate: the new position inherits paths from the current position
            dpTable[pathKey] = right + left;
            }
        
        // advance to the next position:
        // walk backwards through the compatible columns; the first column
        // still in state "possible" becomes the next "used" column
        --ptr;
        while (ptr >= 0 && state[ptr] != possible)
            {
            if (state[ptr] == used)
                {
                for (int j=0; j<numTaxa; j++)
                    pathPos[j] -= indelMatrix[ptr][j];
                state[ptr] = freeToUse;
                }
            --ptr;
            }
        
        if (ptr == -1)
            {
            // we've exhausted all orderings; look up the terminal position
            for (int j=0; j<numTaxa; j++)
                pathFinalPos[j] = 0;
            for (int i=0; i<len; i++)
                for (int j=0; j<numTaxa; j++)
                    pathFinalPos[j] += indelMatrix[i][j];
            
            for (int i=0; i<numTaxa; i++)
                pathKey[i] = pathFinalPos[i];
            std::map<std::vector<int>, int>::iterator it = dpTable.find(pathKey);
            if (it != dpTable.end())
                return it->second;
            return 1;
            }
        
        // mark this column as used and advance the position
        state[ptr] = used;
        for (int j=0; j<numTaxa; j++)
            pathPos[j] += indelMatrix[ptr][j];
        if (ptr <= firstNotUsed)
            firstNotUsed++;
        }
    
    return 1;   // unreachable, but keeps the compiler happy
}

std::string UpdateAlignment::getUpdateName(void) {

    return "Alignment Update: " + myParm->getName();
}

/* Initialize (or reinitialize) the tree structure for this proposal.
   This must be called at the start of each propose() in case the tree
   topology has changed since the last call. The function:
     1. Gets a fresh tree pointer for this alignment's taxon set.
     2. Zeros out all profiles.
     3. Builds left/right descendant pointers from the tree's post-order.
     4. Sorts descendants so that the left child always has the smaller
        minimum-leaf-index —- this ensures a deterministic profile ordering
        that is consistent across proposals. */
void UpdateAlignment::initializeTreeStructure(void) {

    // get fresh tree pointer (in case topology changed)
    tree = treeParm->getTree(taxonMask);
    
    // reset all profiles to zero
    for (int i=0; i<numNodes; i++)
        for (int j=0; j<numNodes; j++)
            for (int k=0; k<maxLength; k++)
                profile[i][j][k] = 0;
    
    // build descendant pointers from the tree's post-order traversal
    for (int i=0; i<numNodes; i++)
        {
        lftDescendants[i] = nullptr;
        rhtDescendants[i] = nullptr;
        }
    
    const std::vector<Node*>& postOrder = tree->getPostOrder();
    for (Node* p : postOrder)
        {
        if (p->getIsLeaf() == false)
            {
            int idx = p->getIndex();
            p->orderDescendantsByOffset();
            if (p->getDescendant(0) != nullptr)
                lftDescendants[idx] = p->getDescendant(0);
            if (p->getDescendant(1) != nullptr)
                rhtDescendants[idx] = p->getDescendant(1);
            }
        }
    
    // Sort descendants so the left child always has the smaller minimum
    // leaf index. This makes the profile row ordering deterministic
    for (int i=0; i<numTaxa; i++)
        sorter[i] = i;
    for (Node* p : postOrder)
        {
        if (p->getIsLeaf() == false)
            {
            int idx = p->getIndex();
            Node* lft = lftDescendants[idx];
            Node* rht = rhtDescendants[idx];
            if (lft != nullptr && rht != nullptr)
                {
                int lftSort = sorter[lft->getIndex()];
                int rhtSort = sorter[rht->getIndex()];
                if (lftSort > rhtSort)
                    {
                    sorter[idx] = rhtSort;
                    // swap descendants so left has the smaller sorter value
                    lftDescendants[idx] = rht;
                    rhtDescendants[idx] = lft;
                    }
                else
                    {
                    sorter[idx] = lftSort;
                    }
                }
            }
        }
}

/* Choose a contiguous region of the alignment to realign. The region
   length follows a geometric distribution: starting at 1, each extension
   has probability extProb. The starting position is drawn uniformly from
   all valid positions for the chosen length. */
void UpdateAlignment::selectRegionToRealign(int numSegments, double extProb, int& pos, int& len) {

    len = 1;
    while (len < numSegments && rng->uniformRv() < extProb)
        len++;
    pos = static_cast<int>(rng->uniformRv() * (numSegments - len + 1));
}

/* Zero out the xProfile and yProfile arrays for all nodes.
   Called at the start of each alignment attempt within propose(). */
void UpdateAlignment::resetProfileDimensions(void) {

    for (int i=0; i<numNodes; i++)
        {
        xProfile[i] = 0;
        yProfile[i] = 0;
        }
}

/* Default propose: uses the member extensionProb and requires that the
   proposed alignment differ from the current one. */
double UpdateAlignment::propose(void) {

    return propose(extensionProb, true);
}

/* The main alignment proposal. This is the heart of the MCMC move:
   
   1. Select a contiguous region of geometric length.
   2. Extract ungapped tip profiles from the selected columns.
   3. Progressively align profiles bottom-up through the tree using a
      stochastic traceback (Boltzmann-weighted sampling from the DP table).
   4. Splice the new aligned region into the alignment.
   5. Compute the Hastings ratio from forward/reverse proposal probabilities
      and the forward/reverse path counts.
   
   Parameters:
     extProb — probability of extending the region by one column (0 < p ≤ 1).
               When extProb = 1.0, the entire alignment is realigned.
     alignmentsMustBeDifferent — if true, the proposal retries until the new
               alignment differs from the old one (avoids wasted MCMC steps). */
double UpdateAlignment::propose(double extProb, bool alignmentsMustBeDifferent) {

    initializeTreeStructure();
    
    // get current and backup alignments
    Alignment* curAlignment = myParm->getAlignment(0);
    Alignment* oldAlignment = myParm->getAlignment(1);
    int numSegments = static_cast<int>(curAlignment->getNumSegments());
    
    tree = treeParm->getTree(taxonMask);
    const std::vector<Node*>& postOrder = tree->getPostOrder();
    Node* rootNode = tree->getRoot();
    int rootIdx = rootNode->getIndex();
    
    bool foundDifferentAlignment = false;
    int len = 0, len2 = 0, pos = 0;
    double proposalLnProb = 0.0;
    
    do {
        // restore from backup if we're retrying after getting the same alignment
        if (foundDifferentAlignment == false && len > 0)
            (*curAlignment) = (*oldAlignment);
        
        // choose region to realign
        selectRegionToRealign(numSegments, extProb, pos, len);
        
        // reset profile dimensions for this attempt
        resetProfileDimensions();
        
        // initialize leaf profiles (ungapped characters only)
        initializeTipProfiles(postOrder, curAlignment, pos, len);
        
        proposalLnProb = 0.0;
        
        // progressive alignment: process internal nodes in post order sequence
        for (Node* p : postOrder)
            {
            if (p->getIsLeaf() == true)
                continue;
            
            int idx = p->getIndex();
            Node* lftDescendant = lftDescendants[idx];
            Node* rhtDescendant = rhtDescendants[idx];
            if (lftDescendant == nullptr || rhtDescendant == nullptr)
                continue;
            
            int lftIdx = lftDescendant->getIndex();
            int rhtIdx = rhtDescendant->getIndex();
            
            // build the scoring matrix from transition probabilities
            buildScoringMatrix(lftDescendant, rhtDescendant);
            
            // fill the Needleman–Wunsch DP table
            fillDpTable(lftIdx, rhtIdx);
            
            // stochastic traceback: sample a new alignment of the two profiles
            // and accumulate the negative log proposal probability
            proposalLnProb += stochasticTraceback(idx, lftIdx, rhtIdx);
            
            // fix up the merged profile: set sizes, copy profile numbers, reverse
            mergeProfilesAfterTraceback(idx, lftIdx, rhtIdx);
            }
        
        // the root profile now contains the newly aligned region
        len2 = xProfile[rootIdx];
        
        // splice the new aligned region into the full alignment
        assembleNewAlignment(curAlignment, rootIdx, pos, len, len2, numSegments);
        
        // check if the alignment actually changed
        if (*curAlignment == *oldAlignment)
            foundDifferentAlignment = false;
        else
            foundDifferentAlignment = true;
            
        } while (foundDifferentAlignment == false && alignmentsMustBeDifferent == true);
    
    // Hastings ratio
    
    // path count correction: how many column orderings produce each alignment?
    int cf1 = countPaths(curAlignment, pos, pos+len2-1);
    if (cf1 == -1)
        {
        // too complex to count, so penalize heavily so this proposal is rejected
        proposalLnProb -= 1e10;
        cf1 = 1;
        }
    proposalLnProb -= log((double)cf1);
    
    int cf2 = countPaths(oldAlignment, pos, pos+len-1);
    if (cf2 == -1)
        {
        // too complex to count, so reward heavily so this proposal is accepted
        proposalLnProb += 1e10;
        cf2 = 1;
        }
    proposalLnProb += log((double)cf2);
    
    // geometric length correction: when extProb < 1, account for the
    // probability of drawing regions of different lengths. When extProb = 1
    // both forward and reverse always select the full alignment, so this is 0.
    if (extProb < 1.0)
        proposalLnProb += (len - len2) * log(1.0 - extProb);
    
    // reverse proposal probability: probability of recovering the old
    // alignment from the new one via deterministic traceback
    proposalLnProb += calculateReverseProposal(curAlignment, pos, len);
    
    return proposalLnProb;
}

void UpdateAlignment::setDependants(void) {

    clearDependencyFlags();

    updatedParameter      = myParm;
    rateMatrixNeedsUpdate = false;
}

/* Standard update: set dependency flags and propose a new alignment. */
double UpdateAlignment::update(void) {
    
    setDependants();
    return propose();
}

/* Stepping-stone update: with probability determined by the Hill function,
   sample from the prior instead of the posterior. */
double UpdateAlignment::update(double power) {

    if (rng->uniformRv() < priorSampleProb(power))
        return updateFromPrior();
    return update();
}

double UpdateAlignment::updateFromPrior(void) {

    setDependants();
    return propose();
}

/* Realign the entire alignment by using extensionProb = 1.0, which
   forces len = numSegments and pos = 0 so every column is realigned.
   Does NOT call setDependants() because the caller (UpdateTopology)
   manages dependencies when updating multiple parameters jointly.
   Also does not require the alignment to differ (alignmentsMustBeDifferent = false)
   because in a joint topology+alignment move, a same-alignment outcome
   is still informative. */
double UpdateAlignment::realignFull(void) {
    
    return propose(1.0, false);
}
