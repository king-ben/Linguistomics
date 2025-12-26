#include "Alignment.hpp"
#include "Model.hpp"
#include "Msg.hpp"
#include "Node.hpp"
#include "ParameterAlignment.hpp"
#include "ParameterFrequencies.hpp"
#include "ParameterTree.hpp"
#include "RandomVariable.hpp"
#include "TransitionProbabilityManager.hpp"
#include "Tree.hpp"
#include "UpdateAlignment.hpp"



UpdateAlignment::UpdateAlignment(Model* m, RandomVariable* r, ParameterAlignment* p, double exponent, double gapPen, double extProb) : 
    Update(m, r), myParm(p) {

    if (exponent < 1.0)
        Msg::error("Exponent parameter must be greater than 1");
    if (gapPen >= 0.0)
        Msg::error("Gap penalty parameter must be less than 0");
    if (extProb <= 0.0 || extProb > 1.0)
        Msg::error("Extension probability must be in range (0,1]");

    freqsParm = model->findParameter<ParameterFrequencies>(); // OK if this is NULL

    treeParm = model->findParameter<ParameterTree>();
    if (treeParm == nullptr)
        Msg::error("Could not find tree parameter for the alignment update");

    tiProbs = model->getTiProbs();

    basis = exponent;
    gapPenalty = gapPen;
    extensionProb = extProb;
    
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
    maxLength = static_cast<int>(aln->getMaximumNumberOfSites());
    maxUnalignDimension = 17;
    
    // allocate profile array: profile[nodeIndex][taxonInProfile][site]
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
    
    xProfile = new int[numNodes];
    yProfile = new int[numNodes];
    for (int i=0; i<numNodes; i++)
        {
        xProfile[i] = 0;
        yProfile[i] = 0;
        }
    
    // allocate dp table
    dp = new double*[maxLength];
    dp[0] = new double[maxLength * maxLength];
    for (int i=1; i<maxLength; i++)
        dp[i] = dp[i-1] + maxLength;
    for (int i=0; i<maxLength; i++)
        for (int j=0; j<maxLength; j++)
            dp[i][j] = 0.0;
    
    // allocate scoring matrix
    scoring = new double*[numStates];
    scoring[0] = new double[numStates * numStates];
    for (int i=1; i<numStates; i++)
        scoring[i] = scoring[i-1] + numStates;
    for (int i=0; i<numStates; i++)
        for (int j=0; j<numStates; j++)
            scoring[i][j] = 0.0;
    
    // allocate helper structures as contiguous blocks
    indelMatrix = new int*[maxLength];
    indelMatrix[0] = new int[maxLength * numTaxa];
    for (int i=1; i<maxLength; i++)
        indelMatrix[i] = indelMatrix[i-1] + numTaxa;
    for (int i=0; i<maxLength; i++)
        for (int j=0; j<numTaxa; j++)
            indelMatrix[i][j] = 0;
    
    profileNumber = new int[numNodes * numNodes];
    for (int i=0; i<numNodes*numNodes; i++)
        profileNumber[i] = 0;
    
    tempProfile = new int*[numNodes];
    tempProfile[0] = new int[numNodes * maxLength];
    for (int i=1; i<numNodes; i++)
        tempProfile[i] = tempProfile[i-1] + maxLength;
    for (int i=0; i<numNodes; i++)
        for (int j=0; j<maxLength; j++)
            tempProfile[i][j] = 0;
    
    lftChildren = new Node*[numNodes];
    rhtChildren = new Node*[numNodes];
    sorter = new int[numNodes];
    for (int i=0; i<numNodes; i++)
        {
        lftChildren[i] = nullptr;
        rhtChildren[i] = nullptr;
        sorter[i] = 0;
        }
    
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
    
    delete [] lftChildren;
    delete [] rhtChildren;
    delete [] sorter;
    
    delete [] possibles;
    delete [] state;
    delete [] pathPos;
    delete [] pathNewPos;
    delete [] pathMask;
    delete [] pathFinalPos;
}

void UpdateAlignment::buildScoringMatrix(Node* lftChild, Node* rhtChild) {

    // get transition probabilities for both branches
    DoubleMatrix& tiL = tiProbs->getTransitionProbability(lftChild->getBranchLength());
    DoubleMatrix& tiR = tiProbs->getTransitionProbability(rhtChild->getBranchLength());
    
    // build scoring matrix: score[i][j] = log(sum_k P(i|k)*P(k|j) / pi_j)
    for (int i=0; i<numStates; i++)
        {
        for (int j=0; j<numStates; j++)
            {
            scoring[i][j] = 0.0;
            for (int k=0; k<numStates; k++)
                scoring[i][j] += tiL(i,k) * tiR(k,j);
            }
        }
    
    // divide by stationary frequencies
    if (freqsParm != nullptr)
        {
        std::vector<double>& freqs = freqsParm->getFrequencies();
        for (int i=0; i<numStates; i++)
            for (int j=0; j<numStates; j++)
                scoring[i][j] = log(scoring[i][j] / freqs[j]);
        }
    else
        {
        // equal frequencies
        double eqFreq = 1.0 / numStates;
        for (int i=0; i<numStates; i++)
            for (int j=0; j<numStates; j++)
                scoring[i][j] = log(scoring[i][j] / eqFreq);
        }
}

double UpdateAlignment::calculateReverseProposal(Alignment* aln, int pos, int len) {

    const std::vector<Node*>& postOrder = tree->getPostOrder();
    
    double reverseLnProb = 0.0;
    
    // build profiles from alignment for tips
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
    
    // build internal node profiles by merging children
    for (Node* p : postOrder)
        {
        if (p->getIsLeaf() == true)
            continue;
        
        int idx = p->getIndex();
        Node* lftChild = lftChildren[idx];
        Node* rhtChild = rhtChildren[idx];
        if (lftChild == nullptr || rhtChild == nullptr)
            continue;
        
        int lftIdx = lftChild->getIndex();
        int rhtIdx = rhtChild->getIndex();
        
        yProfile[idx] = yProfile[lftIdx] + yProfile[rhtIdx];
        
        for (int j=0; j<len; j++)
            {
            for (int l=0; l<yProfile[lftIdx]; l++)
                profile[idx][l][j] = profile[lftIdx][l][j];
            for (int l=yProfile[lftIdx]; l<yProfile[lftIdx]+yProfile[rhtIdx]; l++)
                profile[idx][l][j] = profile[rhtIdx][l-yProfile[lftIdx]][j];
            }
        }
    
    // remove all-gap columns
    for (int k=0; k<numNodes; k++)
        {
        for (int j=0; j<len; j++)
            for (int l=0; l<yProfile[k]; l++)
                tempProfile[l][j] = profile[k][l][j];
        
        xProfile[k] = 0;
        for (int j=0; j<len; j++)
            {
            bool notAllGaps = false;
            for (int l=0; l<yProfile[k]; l++)
                {
                if (tempProfile[l][j] != gapCode)
                    {
                    notAllGaps = true;
                    break;
                    }
                }
            if (notAllGaps)
                {
                for (int l=0; l<yProfile[k]; l++)
                    profile[k][l][xProfile[k]] = tempProfile[l][j];
                xProfile[k]++;
                }
            }
        }
    
    // calculate reverse traceback probabilities
    double logBasis = log(basis);
    
    for (Node* p : postOrder)
        {
        if (p->getIsLeaf() == true)
            continue;
        
        int idx = p->getIndex();
        Node* lftChild = lftChildren[idx];
        Node* rhtChild = rhtChildren[idx];
        if (lftChild == nullptr || rhtChild == nullptr)
            continue;
        
        int lftIdx = lftChild->getIndex();
        int rhtIdx = rhtChild->getIndex();
        
        // build scoring matrix
        buildScoringMatrix(lftChild, rhtChild);
        
        // fill dp table
        dp[0][0] = 0.0;
        for (int i=1; i<=xProfile[lftIdx]; i++)
            dp[i][0] = dp[i-1][0] + gapPenalty;
        for (int j=1; j<=xProfile[rhtIdx]; j++)
            dp[0][j] = dp[0][j-1] + gapPenalty;
        
        for (int i=1; i<=xProfile[lftIdx]; i++)
            {
            for (int j=1; j<=xProfile[rhtIdx]; j++)
                {
                double score1 = dp[i-1][j] + gapPenalty;
                double score2 = dp[i][j-1] + gapPenalty;
                
                int count = 0;
                double score3 = 0.0;
                for (int l=0; l<yProfile[lftIdx]; l++)
                    {
                    for (int m=0; m<yProfile[rhtIdx]; m++)
                        {
                        if (profile[lftIdx][l][i-1] != gapCode && 
                            profile[rhtIdx][m][j-1] != gapCode)
                            {
                            score3 += scoring[profile[lftIdx][l][i-1]]
                                             [profile[rhtIdx][m][j-1]];
                            count++;
                            }
                        }
                    }
                if (count > 0)
                    score3 = dp[i-1][j-1] + score3/count;
                else
                    score3 = dp[i-1][j-1];
                
                dp[i][j] = score1;
                if (score2 > dp[i][j])
                    dp[i][j] = score2;
                if (score3 > dp[i][j])
                    dp[i][j] = score3;
                }
            }
        
        // deterministic traceback following gap pattern
        int ii = xProfile[lftIdx];
        int jj = xProfile[rhtIdx];
        int n = xProfile[idx];
        
        while (n > 0)
            {
            double score1 = 0.0, score2 = 0.0, score3 = 0.0;
            
            if (ii > 0)
                score1 = dp[ii-1][jj] + gapPenalty;
            if (jj > 0)
                score2 = dp[ii][jj-1] + gapPenalty;
            if (ii > 0 && jj > 0)
                {
                int count = 0;
                for (int l=0; l<yProfile[lftIdx]; l++)
                    {
                    for (int m=0; m<yProfile[rhtIdx]; m++)
                        {
                        if (profile[lftIdx][l][ii-1] != gapCode && 
                            profile[rhtIdx][m][jj-1] != gapCode)
                            {
                            score3 += scoring[profile[lftIdx][l][ii-1]]
                                             [profile[rhtIdx][m][jj-1]];
                            count++;
                            }
                        }
                    }
                if (count > 0)
                    score3 = dp[ii-1][jj-1] + score3/count;
                else
                    score3 = dp[ii-1][jj-1];
                }
            
            // compute probabilities
            double sum = 0.0;
            if (ii > 0)
                sum += exp(score1 * logBasis);
            if (jj > 0)
                sum += exp(score2 * logBasis);
            if (ii > 0 && jj > 0)
                sum += exp(score3 * logBasis);
            
            double prob1 = (ii > 0) ? exp(score1 * logBasis) / sum : 0.0;
            double prob2 = (jj > 0) ? exp(score2 * logBasis) / sum : 0.0;
            double prob3 = (ii > 0 && jj > 0) ? exp(score3 * logBasis) / sum : 0.0;
            
            // determine move from gap pattern
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
                reverseLnProb += log(prob2);
                jj--;
                n--;
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
                    reverseLnProb += log(prob1);
                    ii--;
                    n--;
                    }
                else
                    {
                    reverseLnProb += log(prob3);
                    ii--;
                    jj--;
                    n--;
                    }
                }
            }
        }
    
    return reverseLnProb;
}

int UpdateAlignment::countPaths(Alignment* aln, int startCol, int endCol) {

    int len = endCol - startCol + 1;
    if (len <= 0)
        return 1;
    
    getIndelMatrix(aln, startCol, len);
    
    // clear the persistent dpTable
    dpTable.clear();
    
    // reset pre-allocated arrays
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
    
    // insert initial position
    for (int i=0; i<numTaxa; i++)
        pathKey[i] = pathPos[i];
    dpTable[pathKey] = 1;
    
    int firstNotUsed = 0;
    
    while (true)
        {
        // find compatible columns from current position
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
                        return -1;
                    possibles[numPossible++] = ptr;
                    }
                
                for (int j=0; j<numTaxa; j++)
                    pathMask[j] += indelMatrix[ptr][j];
                
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
        
        // enumerate combinations
        for (int combo=1; combo<(1<<numPossible); combo++)
            {
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
            
            // lookup pathPos
            for (int i=0; i<numTaxa; i++)
                pathKey[i] = pathPos[i];
            std::map<std::vector<int>, int>::iterator it = dpTable.find(pathKey);
            int left = (it != dpTable.end()) ? it->second : 0;
            
            // lookup pathNewPos
            for (int i=0; i<numTaxa; i++)
                pathKey[i] = pathNewPos[i];
            it = dpTable.find(pathKey);
            int right = (it != dpTable.end()) ? it->second : 0;
            
            // insert/update pathNewPos
            dpTable[pathKey] = right + left;
            }
        
        // find next position
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
            for (int j=0; j<numTaxa; j++)
                pathFinalPos[j] = 0;
            for (int i=0; i<len; i++)
                for (int j=0; j<numTaxa; j++)
                    pathFinalPos[j] += indelMatrix[i][j];
            
            // lookup pathFinalPos
            for (int i=0; i<numTaxa; i++)
                pathKey[i] = pathFinalPos[i];
            std::map<std::vector<int>, int>::iterator it = dpTable.find(pathKey);
            if (it != dpTable.end())
                return it->second;
            return 1;
            }
        
        state[ptr] = used;
        for (int j=0; j<numTaxa; j++)
            pathPos[j] += indelMatrix[ptr][j];
        if (ptr <= firstNotUsed)
            firstNotUsed++;
        }
    
    return 1;
}

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

std::string UpdateAlignment::getUpdateName(void) {

    return myParm->getName() + " Update";
}

void UpdateAlignment::initializeTreeStructure(void) {

    // get fresh tree pointer (in case it changed)
    tree = treeParm->getTree(taxonMask);
    
    // reset profiles
    for (int i=0; i<numNodes; i++)
        for (int j=0; j<numNodes; j++)
            for (int k=0; k<maxLength; k++)
                profile[i][j][k] = 0;
    
    // build child pointers from tree structure
    for (int i=0; i<numNodes; i++)
        {
        lftChildren[i] = nullptr;
        rhtChildren[i] = nullptr;
        }
    
    const std::vector<Node*>& postOrder = tree->getPostOrder();
    for (Node* p : postOrder)
        {
        if (p->getIsLeaf() == false)
            {
            int idx = p->getIndex();
            std::set<Node*,NodeComparator>& des = p->getDescendants().getNodes();
            auto it = des.begin();
            if (it != des.end())
                {
                lftChildren[idx] = *it;
                ++it;
                }
            if (it != des.end())
                rhtChildren[idx] = *it;
            }
        }
    
    // sort children to ensure profile order matches taxon order
    for (int i=0; i<numTaxa; i++)
        sorter[i] = i;
    for (Node* p : postOrder)
        {
        if (p->getIsLeaf() == false)
            {
            int idx = p->getIndex();
            Node* lft = lftChildren[idx];
            Node* rht = rhtChildren[idx];
            if (lft != nullptr && rht != nullptr)
                {
                int lftSort = sorter[lft->getIndex()];
                int rhtSort = sorter[rht->getIndex()];
                if (lftSort > rhtSort)
                    {
                    sorter[idx] = rhtSort;
                    // swap children so left has smaller sorter value
                    lftChildren[idx] = rht;
                    rhtChildren[idx] = lft;
                    }
                else
                    {
                    sorter[idx] = lftSort;
                    }
                }
            }
        }
}

void UpdateAlignment::notifyDependants(void) {

}

double UpdateAlignment::propose(void) {

    initializeTreeStructure();
    
    // get alignments
    Alignment* curAlignment = myParm->getAlignment(0);
    Alignment* oldAlignment = myParm->getAlignment(1);
    int numSites = static_cast<int>(curAlignment->getNumSites());
    
    const std::vector<Node*>& postOrder = tree->getPostOrder();
    Node* rootNode = tree->getRoot();
    int rootIdx = rootNode->getIndex();
    
    bool foundDifferentAlignment = false;
    int len = 0, len2 = 0, pos = 0;
    double proposalLnProb = 0.0;
    
    do
        {
        // restore alignment from backup if retrying
        if (foundDifferentAlignment == false && len > 0)
            (*curAlignment) = (*oldAlignment);
        
        // choose region to realign using geometric distribution
        len = 1;
        while (len < numSites && rng->uniformRv() < extensionProb)
            len++;
        //pos = rv->uniformRvInt(numSites - len + 1);
        pos = static_cast<int>(rng->uniformRv() * (numSites - len + 1));
        
        // reset profile dimensions
        for (int i=0; i<numNodes; i++)
            {
            xProfile[i] = 0;
            yProfile[i] = 0;
            }
        
        // initialize tip profiles (extract non-gap characters from selected region)
        for (Node* p : postOrder)
            {
            if (p->getIsLeaf() == true)
                {
                int idx = p->getIndex();
                yProfile[idx] = 1;
                for (int j=pos; j<pos+len; j++)
                    {
                    int c = (*curAlignment)(idx, j);
                    if (c != gapCode)
                        {
                        profile[idx][0][xProfile[idx]] = c;
                        xProfile[idx]++;
                        }
                    }
                profileNumber[idx * numNodes] = idx;
                }
            }
        
        proposalLnProb = 0.0;
        
        // progressive alignment: process internal nodes in post-order
        for (Node* p : postOrder)
            {
            if (p->getIsLeaf() == true)
                continue;
            
            int idx = p->getIndex();
            Node* lftChild = lftChildren[idx];
            Node* rhtChild = rhtChildren[idx];
            if (lftChild == nullptr || rhtChild == nullptr)
                continue;
            
            int lftIdx = lftChild->getIndex();
            int rhtIdx = rhtChild->getIndex();
            
            // build scoring matrix from transition probabilities
            buildScoringMatrix(lftChild, rhtChild);
            
            // initialize dp table borders
            dp[0][0] = 0.0;
            for (int i=1; i<=xProfile[lftIdx]; i++)
                dp[i][0] = dp[i-1][0] + gapPenalty;
            for (int j=1; j<=xProfile[rhtIdx]; j++)
                dp[0][j] = dp[0][j-1] + gapPenalty;
            
            // fill dp table
            for (int i=1; i<=xProfile[lftIdx]; i++)
                {
                for (int j=1; j<=xProfile[rhtIdx]; j++)
                    {
                    double score1 = dp[i-1][j] + gapPenalty;
                    double score2 = dp[i][j-1] + gapPenalty;
                    
                    int count = 0;
                    double score3 = 0.0;
                    for (int l=0; l<yProfile[lftIdx]; l++)
                        {
                        for (int m=0; m<yProfile[rhtIdx]; m++)
                            {
                            if (profile[lftIdx][l][i-1] != gapCode && 
                                profile[rhtIdx][m][j-1] != gapCode)
                                {
                                score3 += scoring[profile[lftIdx][l][i-1]]
                                                 [profile[rhtIdx][m][j-1]];
                                count++;
                                }
                            }
                        }
                    if (count > 0)
                        score3 = dp[i-1][j-1] + score3/count;
                    else
                        score3 = dp[i-1][j-1];
                    
                    dp[i][j] = score1;
                    if (score2 > dp[i][j])
                        dp[i][j] = score2;
                    if (score3 > dp[i][j])
                        dp[i][j] = score3;
                    }
                }
            
            // stochastic traceback
            int ii = xProfile[lftIdx];
            int jj = xProfile[rhtIdx];
            double logBasis = log(basis);
            
            while (ii > 0 || jj > 0)
                {
                double score1 = 0.0, score2 = 0.0, score3 = 0.0;
                
                if (ii > 0)
                    score1 = dp[ii-1][jj] + gapPenalty;
                if (jj > 0)
                    score2 = dp[ii][jj-1] + gapPenalty;
                if (ii > 0 && jj > 0)
                    {
                    int count = 0;
                    for (int l=0; l<yProfile[lftIdx]; l++)
                        {
                        for (int m=0; m<yProfile[rhtIdx]; m++)
                            {
                            if (profile[lftIdx][l][ii-1] != gapCode && 
                                profile[rhtIdx][m][jj-1] != gapCode)
                                {
                                score3 += scoring[profile[lftIdx][l][ii-1]]
                                                 [profile[rhtIdx][m][jj-1]];
                                count++;
                                }
                            }
                        }
                    if (count > 0)
                        score3 = dp[ii-1][jj-1] + score3/count;
                    else
                        score3 = dp[ii-1][jj-1];
                    }
                
                // convert scores to probabilities
                double sum = 0.0;
                if (ii > 0)
                    sum += exp(score1 * logBasis);
                if (jj > 0)
                    sum += exp(score2 * logBasis);
                if (ii > 0 && jj > 0)
                    sum += exp(score3 * logBasis);
                
                double prob1 = (ii > 0) ? exp(score1 * logBasis) / sum : 0.0;
                double prob2 = (jj > 0) ? exp(score2 * logBasis) / sum : 0.0;
                double prob3 = (ii > 0 && jj > 0) ? exp(score3 * logBasis) / sum : 0.0;
                
                // stochastic choice
                double u = rng->uniformRv();
                if (u < prob1)
                    {
                    proposalLnProb -= log(prob1);
                    for (int l=0; l<yProfile[lftIdx]; l++)
                        tempProfile[l][xProfile[idx]] = profile[lftIdx][l][ii-1];
                    for (int l=yProfile[lftIdx]; l<yProfile[lftIdx]+yProfile[rhtIdx]; l++)
                        tempProfile[l][xProfile[idx]] = gapCode;
                    xProfile[idx]++;
                    ii--;
                    }
                else if (u < prob1 + prob2)
                    {
                    proposalLnProb -= log(prob2);
                    for (int l=0; l<yProfile[lftIdx]; l++)
                        tempProfile[l][xProfile[idx]] = gapCode;
                    for (int l=yProfile[lftIdx]; l<yProfile[lftIdx]+yProfile[rhtIdx]; l++)
                        tempProfile[l][xProfile[idx]] = profile[rhtIdx][l-yProfile[lftIdx]][jj-1];
                    xProfile[idx]++;
                    jj--;
                    }
                else
                    {
                    if (prob3 <= 0.0)
                        Msg::error("Problem in alignment proposal: prob3 <= 0");
                    proposalLnProb -= log(prob3);
                    for (int l=0; l<yProfile[lftIdx]; l++)
                        tempProfile[l][xProfile[idx]] = profile[lftIdx][l][ii-1];
                    for (int l=yProfile[lftIdx]; l<yProfile[lftIdx]+yProfile[rhtIdx]; l++)
                        tempProfile[l][xProfile[idx]] = profile[rhtIdx][l-yProfile[lftIdx]][jj-1];
                    xProfile[idx]++;
                    jj--;
                    ii--;
                    }
                }
            
            yProfile[idx] = yProfile[lftIdx] + yProfile[rhtIdx];
            
            // update profile numbers
            for (int i=0; i<yProfile[lftIdx]; i++)
                profileNumber[idx*numNodes + i] = profileNumber[lftIdx*numNodes + i];
            for (int i=yProfile[lftIdx]; i<yProfile[lftIdx]+yProfile[rhtIdx]; i++)
                profileNumber[idx*numNodes + i] = profileNumber[rhtIdx*numNodes + i-yProfile[lftIdx]];
            
            // reverse profile and copy from tempProfile
            for (int i=0; i<xProfile[idx]; i++)
                for (int j=0; j<yProfile[idx]; j++)
                    profile[idx][j][i] = tempProfile[j][xProfile[idx]-i-1];
            }
        
        len2 = xProfile[rootIdx];
        
        // build new alignment
        int newNumSites = numSites + len2 - len;
        
        // shift suffix if needed
        if (len2 < len)
            {
            for (int i=pos+len; i<numSites; i++)
                for (int j=0; j<numTaxa; j++)
                    (*curAlignment)(j, i+len2-len) = (*curAlignment)(j, i);
            }
        else if (len2 > len)
            {
            for (int i=numSites-1; i>=pos+len; i--)
                for (int j=0; j<numTaxa; j++)
                    (*curAlignment)(j, i+len2-len) = (*curAlignment)(j, i);
            }
        
        curAlignment->setNumSites(newNumSites);
        
        // copy new aligned region from root profile
        for (int i=pos; i<pos+len2; i++)
            {
            for (int j=0; j<numTaxa; j++)
                {
                int taxonIdx = profileNumber[rootIdx*numNodes + j];
                (*curAlignment)(taxonIdx, i) = profile[rootIdx][j][i-pos];
                }
            }
        
        // check if alignment changed
        if (*curAlignment == *oldAlignment)
            foundDifferentAlignment = false;
        else
            foundDifferentAlignment = true;
            
        } while (foundDifferentAlignment == false);
    
    // calculate Hastings ratio components
    int cf1 = countPaths(curAlignment, pos, pos+len2-1);
    if (cf1 == -1)
        {
        proposalLnProb -= 1e10;
        cf1 = 1;
        }
    proposalLnProb -= log((double)cf1);
    
    int cf2 = countPaths(oldAlignment, pos, pos+len-1);
    if (cf2 == -1)
        {
        proposalLnProb += 1e10;
        cf2 = 1;
        }
    proposalLnProb += log((double)cf2);
    
    // length change contribution
    proposalLnProb += (len - len2) * log(1.0 - extensionProb);
    
    // reverse proposal probability
    proposalLnProb += calculateReverseProposal(curAlignment, pos, len);
    
    return proposalLnProb;
}

void UpdateAlignment::setDependants(void) {

    clearDependencyFlags();

    updatedParameter      = myParm;
    rateMatrixNeedsUpdate = false;
}

double UpdateAlignment::update(void) {

    numTries++;
    
    setDependants();
    return propose();
}

double UpdateAlignment::update(double power) {

    if (rng->uniformRv() < priorSampleProb(power))
        return updateFromPrior();
    return update();
}

double UpdateAlignment::updateFromPrior(void) {

    setDependants();
    return propose();
}

