#include <cstring>
#include "Alignment.hpp"
#include "LikelihoodCalculator.hpp"
#include "Msg.hpp"
#include "Node.hpp"
#include "ParameterAlignment.hpp"
#include "ParameterFrequencies.hpp"
#include "ParameterIndelRates.hpp"
#include "ParameterTree.hpp"
#include "TransitionProbabilitiesGpu.hpp"
#include "Tree.hpp"



LikelihoodCalculator::LikelihoodCalculator(TransitionProbabilities* tpc, ParameterAlignment* a, ParameterTree* t, ParameterIndelRates* r, ParameterFrequencies* f) : 
    tiProbs(tpc), myAlignment(a), myTree(t), myIndelRates(r), myFrequencies(f) {

    // set the alignment for this calculator
    alignment = myAlignment->getAlignment();
    
    // and the taxon mask
    taxonMask = myAlignment->getTaxonMask();

    // and the number of states 
    numStates = myAlignment->getNumStates();
    
    // allocate the stationary frequencies
    equilibriumFrequencies = new double[numStates];

    // tree information
    tree = myTree->getTree(taxonMask);
    if (tree == nullptr)
        Msg::error("Could not find tree when initializing TKF91 calculator");
    numTaxa = tree->getNumTaxa();
    numNodes = 2 * numTaxa - 1;
    if (tree->isBinary() == false)
        Msg::error("Expecting a rooted and binary tree");

    // determine number of taxa, sites, and states from alignment
    numSites = alignment->getNumSites();
    if (alignment->getNumTaxa() != numTaxa)
        {
        tree->print();
        alignment->print();
        Msg::error("Inconsistent number of taxa in the tree and the alignment");
        }
            
    // store TKF91 rates
    tkf91Probs.insertionRate = myIndelRates->getInsertionRate();
    tkf91Probs.deletionRate = myIndelRates->getDeletionRate();
    
    // allocate probability arrays
    allocateTKF91Probabilities(numNodes);
    allocateTKF91Combinatorics(numNodes);
    allocateConditionalProbs();
    
    // allocate cached transition matrix pointer array
    cachedTiMatrices = new double*[numNodes];
    for (size_t i = 0; i < numNodes; i++)
        cachedTiMatrices[i] = nullptr;
    
    // allocate signature array (reused across all likelihood calculations)
    signature = new int[numTaxa];
    
    // initialize leaf index mapping
    initializeLeafIndices();
}

LikelihoodCalculator::~LikelihoodCalculator(void) {

    freeTKF91Probabilities();
    freeTKF91Combinatorics();
    freeConditionalProbs();
    delete [] cachedTiMatrices;
    delete [] signature;
    delete [] equilibriumFrequencies;
}

void LikelihoodCalculator::allocateTKF91Probabilities(size_t nn) {

    tkf91Probs.beta = new double[nn];
    tkf91Probs.birthProbability = new double[nn];
    tkf91Probs.extinctionProbability = new double[nn];
    tkf91Probs.homologousProbability = new double[nn];
    tkf91Probs.nonHomologousProbability = new double[nn];
    for (size_t i=0; i<nn; i++)
        {
        tkf91Probs.beta[i] = 0.0;
        tkf91Probs.birthProbability[i] = 0.0;
        tkf91Probs.extinctionProbability[i] = 0.0;
        tkf91Probs.homologousProbability[i] = 0.0;
        tkf91Probs.nonHomologousProbability[i] = 0.0;
        }
}

void LikelihoodCalculator::allocateTKF91Combinatorics(size_t nn) {

    tkf91Combos.nodeHomology = new int[nn];
    tkf91Combos.numHomologousEmissions = new int[nn];
    for (size_t i=0; i<nn; i++)
        {
        tkf91Combos.nodeHomology[i] = 0;
        tkf91Combos.numHomologousEmissions[i] = 0;
        }
}

void LikelihoodCalculator::allocateConditionalProbs(void) {

    // fH[nodeIdx][state] for state 0..numStates-1 (homologous contributions)
    // fI[nodeIdx][state] for state 0..numStates-1, plus fI[nodeIdx][numStates] for gap
    size_t numStates1 = numStates + 1;
    
    fH = new double*[numNodes];
    fH[0] = new double[numNodes * numStates];
    for (size_t i=1; i<numNodes; i++)
        fH[i] = fH[i-1] + numStates;
    
    fI = new double*[numNodes];
    fI[0] = new double[numNodes * numStates1];
    for (size_t i=1; i<numNodes; i++)
        fI[i] = fI[i-1] + numStates1;
}

void LikelihoodCalculator::cacheTransitionMatrices(void) {

    // cache raw pointers to transition probability matrices for all non-root nodes
    // this avoids repeated map/hash lookups in the inner site loop
    const std::vector<Node*>& postOrder = tree->getPostOrder();
    for (Node* p : postOrder)
        {
        if (p->getAncestor() != nullptr)  // skip root (no branch)
            {
            const int pIdx = p->getIndex();
            DoubleMatrix& tp = tiProbs->getTransitionProbability(p->getBranchLength());
            cachedTiMatrices[pIdx] = tp.begin();
            }
        }
}

void LikelihoodCalculator::computeFForInternalNode(Node* p) {

    // get descendant nodes
    Node* lftDescendant = p->getDescendant(0);
    Node* rhtDescendant = p->getDescendant(1);
    if (lftDescendant == nullptr || rhtDescendant == nullptr)
        Msg::error("Did not find descendants of node indexed " + std::to_string(p->getIndex()));
    const int pIdx = p->getIndex();
    const int lftChildIdx = lftDescendant->getIndex();
    const int rhtChildIdx = rhtDescendant->getIndex();
    
    // get homology class for this node and children
    const int nodeHomologyI = tkf91Combos.nodeHomology[pIdx];
    const int nodeHomologyLft = tkf91Combos.nodeHomology[lftChildIdx];
    const int nodeHomologyRht = tkf91Combos.nodeHomology[rhtChildIdx];
    
    // use cached transition probability matrix pointers (avoids map lookup)
    double* const tpLftData = cachedTiMatrices[lftChildIdx];
    double* const tpRhtData = cachedTiMatrices[rhtChildIdx];
    const size_t n = numStates;
    
    // pointers to descendant conditional likelihoods
    double* const fHLeft  = fH[lftChildIdx];
    double* const fHRight = fH[rhtChildIdx];
    double* const fILeft  = fI[lftChildIdx];
    double* const fIRight = fI[rhtChildIdx];
    
    // pointers to parent conditional likelihoods
    double* const fHIdx = fH[pIdx];
    double* const fIIdx = fI[pIdx];

    if ((nodeHomologyI != 0) && (nodeHomologyI == nodeHomologyLft) && (nodeHomologyI == nodeHomologyRht))
        {
        // Case 1: Homologous nucleotide travels down both edges
        const double probLeft  = tkf91Probs.homologousProbability[lftChildIdx];
        const double probRight = tkf91Probs.homologousProbability[rhtChildIdx];
        
        for (size_t i = 0; i < n; i++)
            {
            const double* tpLftRow = tpLftData + i * n;
            const double* tpRhtRow = tpRhtData + i * n;
            double lft = 0.0, rht = 0.0;
            for (size_t j = 0; j < n; j++)
                {
                lft += fHLeft[j]  * tpLftRow[j];
                rht += fHRight[j] * tpRhtRow[j];
                }
            fHIdx[i] = (lft * probLeft) * (rht * probRight);
            }
        // fI values remain zero for this case
        }
    else if (nodeHomologyI != 0)
        {
        // Case 2: Homologous nucleotide travels down one specific edge
        int homologousChildIdx, inhomologousChildIdx;
        double* tpHomData;
        double* tpInhomData;
        double* fH_hom;
        double* fH_inhom;
        double* fI_inhom;
        
        if (nodeHomologyI == nodeHomologyLft)
            {
            homologousChildIdx   = lftChildIdx;
            inhomologousChildIdx = rhtChildIdx;
            tpHomData   = tpLftData;
            tpInhomData = tpRhtData;
            fH_hom   = fHLeft;
            fH_inhom = fHRight;
            fI_inhom = fIRight;
            }
        else
            {
            homologousChildIdx   = rhtChildIdx;
            inhomologousChildIdx = lftChildIdx;
            tpHomData   = tpRhtData;
            tpInhomData = tpLftData;
            fH_hom   = fHRight;
            fH_inhom = fHLeft;
            fI_inhom = fILeft;
            }
        
        const double birthProb            = tkf91Probs.birthProbability[inhomologousChildIdx];
        const double extinctionProb       = tkf91Probs.extinctionProbability[inhomologousChildIdx];
        const double homologousProb_hom   = tkf91Probs.homologousProbability[homologousChildIdx];
        const double homologousProb_inhom = tkf91Probs.homologousProbability[inhomologousChildIdx];
        const double nonHomologousProb    = tkf91Probs.nonHomologousProbability[inhomologousChildIdx];
        const double rhtStartProb         = extinctionProb * fI_inhom[n];
        const double factor1              = nonHomologousProb - extinctionProb * birthProb;
        
        // precompute frequency-weighted sum (constant across i)
        double freqSum = 0.0;
        for (size_t j = 0; j < n; j++)
            freqSum += (fH_inhom[j] + fI_inhom[j]) * equilibriumFrequencies[j];
        const double freqTerm = factor1 * freqSum;
        
        for (size_t i = 0; i < n; i++)
            {
            const double* tpHomRow = tpHomData + i * n;
            const double* tpInhomRow = tpInhomData + i * n;
            double lft = 0.0;
            double rht = 0.0;
            for (size_t j = 0; j < n; j++)
                {
                lft += fH_hom[j] * tpHomRow[j];
                rht += fI_inhom[j] * tpInhomRow[j];
                }
            fHIdx[i] = (lft * homologousProb_hom) * (rhtStartProb + freqTerm + rht * homologousProb_inhom);
            }
        // fI values remain zero for this case
        }
    else
        {
        // Case 3: No homology family's spanning tree intersects either child edge
        const double extinctionProbabilityLeft  = tkf91Probs.extinctionProbability[lftChildIdx];
        const double extinctionProbabilityRight = tkf91Probs.extinctionProbability[rhtChildIdx];
        const double homologousProbabilityLeft  = tkf91Probs.homologousProbability[lftChildIdx];
        const double homologousProbabilityRight = tkf91Probs.homologousProbability[rhtChildIdx];
        const double birthProbabilityLeft  = tkf91Probs.birthProbability[lftChildIdx];
        const double birthProbabilityRight = tkf91Probs.birthProbability[rhtChildIdx];
        const double nonHomologousProbabilityLeft  = tkf91Probs.nonHomologousProbability[lftChildIdx];
        const double nonHomologousProbabilityRight = tkf91Probs.nonHomologousProbability[rhtChildIdx];
        
        const double lfactor = nonHomologousProbabilityLeft  - extinctionProbabilityLeft  * birthProbabilityLeft;
        const double rfactor = nonHomologousProbabilityRight - extinctionProbabilityRight * birthProbabilityRight;
        
        const double rht1Base = extinctionProbabilityLeft  * fILeft[n];
        const double rht2Base = extinctionProbabilityRight * fIRight[n];
        
        // precompute frequency-weighted sums (constant across i)
        double freqSumLeft = 0.0, freqSumRight = 0.0;
        for (size_t j = 0; j < n; j++)
            {
            freqSumLeft  += (fHLeft[j]  + fILeft[j])  * equilibriumFrequencies[j];
            freqSumRight += (fHRight[j] + fIRight[j]) * equilibriumFrequencies[j];
            }
        const double freqTermLeft  = lfactor * freqSumLeft;
        const double freqTermRight = rfactor * freqSumRight;
        
        for (size_t i = 0; i < n; i++)
            {
            const double* tpLftRow = tpLftData + i * n;
            const double* tpRhtRow = tpRhtData + i * n;
            
            double lft1 = 0.0, lft2 = 0.0;
            double rht1 = 0.0, rht2 = 0.0;
            
            for (size_t j = 0; j < n; j++)
                {
                const double tpLftIJ = tpLftRow[j];
                const double tpRhtIJ = tpRhtRow[j];
                
                lft1 += fHLeft[j]  * tpLftIJ;
                lft2 += fHRight[j] * tpRhtIJ;
                rht1 += fILeft[j]  * tpLftIJ;
                rht2 += fIRight[j] * tpRhtIJ;
                }
            
            lft1 *= homologousProbabilityLeft;
            lft2 *= homologousProbabilityRight;
            rht1 = rht1Base + freqTermLeft  + rht1 * homologousProbabilityLeft;
            rht2 = rht2Base + freqTermRight + rht2 * homologousProbabilityRight;
            
            fHIdx[i] = lft1 * rht2 + lft2 * rht1;
            fIIdx[i] = rht1 * rht2;
            }
        
        // Compute fI[p][numStates] (gap state)
        const double left  = fILeft[n]  - birthProbabilityLeft  * freqSumLeft;
        const double right = fIRight[n] - birthProbabilityRight * freqSumRight;
        fIIdx[n] = left * right;
        }
}

void LikelihoodCalculator::computeFForLeafNode(Node* p, int* sig, size_t site) {

    const int pIdx = p->getIndex();
    const int leafIdx = getLeafIndex(p);

    // initialize all values to zero
    double* const fHPtr = fH[pIdx];
    double* const fIPtr = fI[pIdx];
    std::memset(fHPtr, 0, numStates * sizeof(double));
    std::memset(fIPtr, 0, (numStates + 1) * sizeof(double));
    
    // set homology tracking for this leaf
    tkf91Combos.nodeHomology[pIdx] = sig[leafIdx];
    tkf91Combos.numHomologousEmissions[pIdx] = (sig[leafIdx] != 0) ? 1 : 0;
    
    // set values based on signature (gap vs character)
    if (sig[leafIdx] == 0)
        {
        // gap at this position
        fIPtr[numStates] = 1.0;
        }
    else
        {
        // character at this position
        const size_t charState = (*alignment)(leafIdx, site);
        if (charState < numStates)
            {
            fHPtr[charState] = 1.0;
            }
        else
            {
            // this shouldn't happen if signature is correct
            Msg::warning("Character state >= numStates at leaf " + std::to_string(pIdx));
            fIPtr[numStates] = 1.0;
            }
        }
}

double LikelihoodCalculator::computeLnLikelihood(void) {

    // get the number of columns in the alignment
    numSites = alignment->getNumSites();
    if (numSites == 0)
        return 0.0;

    // set TKF91 parameters
    setBirthDeathProbabilities();
    
    // set the stationary frequencies
    setStationaryFrequencies();

    // get the tree
    tree = myTree->getTree(taxonMask);
    if (tree == nullptr)
        Msg::error("Could not find tree when initializing TKF91 calculator");
    
    // cache transition probability matrix pointers (avoids repeated lookups per site)
    cacheTransitionMatrices();
    
    // compute the null emission factor (all gaps)
    // (use member signature array, zero it out for null signature)
    std::memset(signature, 0, numTaxa * sizeof(int));
    double nullFactor = computeRootFI(signature, -1);
    
    // initialize probability with immortal link contribution
    // P(0) = prod_{n in T} (1 - B_n) / G^0(r, -)
    double prob = log(tkf91Probs.immortalProbability) - log(nullFactor);
    
#   if defined(DEBUG_TKF91)
    std::cout << "immortalProbability = " << tkf91Probs.immortalProbability << std::endl;
    std::cout << "Null emission factor = " << nullFactor << std::endl;
    std::cout << "Initial P(0) = " << prob << std::endl;
#   endif
    
    // process each column of the alignment
    // P(K) = P(K-v) * (-G^v(r,-)) / G^0(r,-)
    for (size_t site=0; site<numSites; site++)
        {
        // build signature for this column (1 if character, 0 if gap)
        for (size_t i = 0; i < numTaxa; i++)
            {
            const size_t charState = (*alignment)(i, site);
            signature[i] = (charState == numStates) ? 0 : 1;
            }
        
        // compute emission factor for this column
        const double Gv = computeRootFI(signature, site);
        
        // update probability: P(K) = P(K-v) * (-G^v) / G^0
        const double factor = -Gv / nullFactor;
        prob += log(factor);
        
#       if defined(DEBUG_TKF91)
        std::cout << "Site " << site << ": signature=(";
        for (size_t i = 0; i < numTaxa; i++)
            std::cout << signature[i] << (i < numTaxa-1 ? "," : "");
        std::cout << std::scientific << std::setprecision(10);
        std::cout << "), Gv=" << Gv << ", factor=" << factor << ", prob=" << prob << std::endl;
#       endif
        }
        
    return prob;
}

double LikelihoodCalculator::computeRootFI(int* sig, size_t site) {

    // initialize all fH, fI, and homology tracking values
    const size_t numStates1 = numStates + 1;
    std::memset(fH[0], 0, numNodes * numStates * sizeof(double));
    std::memset(fI[0], 0, numNodes * numStates1 * sizeof(double));
    std::memset(tkf91Combos.nodeHomology, 0, numNodes * sizeof(int));
    std::memset(tkf91Combos.numHomologousEmissions, 0, numNodes * sizeof(int));

    // count total emissions for the homology class (class 1 = has character)
    int numHomologousEmissionsForClass1 = 0;
    for (size_t i = 0; i < numTaxa; i++)
        {
        if (sig[i] != 0)
            numHomologousEmissionsForClass1++;
        }
    
    // get post-order traversal (tips to root)
    const std::vector<Node*>& postOrder = tree->getPostOrder();
    
    // first pass: initialize leaves
    for (Node* p : postOrder)
        {
        if (p->getIsLeaf() == true)
            computeFForLeafNode(p, sig, site);
        }
    
    // second pass: propagate homology information from descendants to ancestors
    for (Node* p : postOrder)
        {
        if (p->getAncestor() != nullptr)
            {
            const int pIdx = p->getIndex();
            const int pAncIdx = p->getAncestor()->getIndex();
            
            const int numHomologousEmission = tkf91Combos.numHomologousEmissions[pIdx];
            const int nodeHomologyI = tkf91Combos.nodeHomology[pIdx];
            
            if (nodeHomologyI == 0)
                {
                // gap node - do nothing
                }
            else if (numHomologousEmission == numHomologousEmissionsForClass1)
                {
                // this node is the MRCA of all homologous taxa - do not propagate
                }
            else
                {
                // not yet MRCA, propagate to ancestor
                const int nodeHomologyAnc = tkf91Combos.nodeHomology[pAncIdx];
                if (nodeHomologyAnc == 0 || nodeHomologyAnc == nodeHomologyI)
                    {
                    tkf91Combos.nodeHomology[pAncIdx] = nodeHomologyI;
                    tkf91Combos.numHomologousEmissions[pAncIdx] += numHomologousEmission;
                    }
                }
            }
        }
    
    // third pass: compute conditional probabilities
    for (Node* p : postOrder)
        {
        if (p->getIsLeaf() == false)
            computeFForInternalNode(p);
        }
    
    // compute result at root with correction
    const int rootIdx = tree->getRoot()->getIndex();
    const double* const fIRoot = fI[rootIdx];
    const double* const fHRoot = fH[rootIdx];
    const double birthProbRoot = tkf91Probs.birthProbability[rootIdx];
    
    double res = fIRoot[numStates];
    double freqSum = 0.0;
    for (size_t i=0; i<numStates; i++)
        freqSum += (fIRoot[i] + fHRoot[i]) * equilibriumFrequencies[i];
    res -= birthProbRoot * freqSum;
    
    return res;
}

void LikelihoodCalculator::freeTKF91Probabilities(void) {

    delete [] tkf91Probs.beta;
    delete [] tkf91Probs.birthProbability;
    delete [] tkf91Probs.extinctionProbability;
    delete [] tkf91Probs.homologousProbability;
    delete [] tkf91Probs.nonHomologousProbability;
}

void LikelihoodCalculator::freeTKF91Combinatorics(void) {

    delete [] tkf91Combos.nodeHomology;
    delete [] tkf91Combos.numHomologousEmissions;
}

void LikelihoodCalculator::freeConditionalProbs(void) {

    delete [] fH[0];
    delete [] fH;
    delete [] fI[0];
    delete [] fI;
}

int LikelihoodCalculator::getLeafIndex(Node* n) {

    return leafIndexMap[n->getIndex()];
}

void LikelihoodCalculator::initializeLeafIndices(void) {

    leafIndexMap.resize(2*numTaxa-1, -1);
    
    // leaf nodes have indices 0, 1, 2, ... which directly corresponds to sequence order
    const std::vector<Node*>& postOrder = tree->getPostOrder();
    for (Node* p : postOrder)
        {
        if (p->getIsLeaf() == true)
            {
            const int pIdx = p->getIndex();
            leafIndexMap[pIdx] = pIdx;
            }
        }
#   if defined(DEBUG_TKF91)
    for (size_t i = 0; i < leafIndexMap.size(); i++)
        std::cout << i << " -> " << leafIndexMap[i] << std::endl;
#   endif
}

void LikelihoodCalculator::setBirthDeathProbabilities(void) {

    tkf91Probs.insertionRate = myIndelRates->getInsertionRate();
    tkf91Probs.deletionRate = myIndelRates->getDeletionRate();

    const double lambda = tkf91Probs.insertionRate;
    const double mu     = tkf91Probs.deletionRate;
    
    tkf91Probs.immortalProbability = 1.0;
    const std::vector<Node*>& postOrder = tree->getPostOrder();
    for (Node* p : postOrder)
        {
        const int pIdx = p->getIndex();
        const double v = p->getBranchLength();
        if (p->getAncestor() == nullptr)
            {
            // root node
            tkf91Probs.beta[pIdx] = 1.0 / mu;
            tkf91Probs.homologousProbability[pIdx] = 0.0;
            }
        else
            {
            // internal node or tip
            const double expPart = std::exp((lambda - mu) * v);
            tkf91Probs.beta[pIdx] = (1.0 - expPart) / (mu - lambda * expPart);
            tkf91Probs.homologousProbability[pIdx] = std::exp(-mu * v) * (1.0 - lambda * tkf91Probs.beta[pIdx]);
            }
        tkf91Probs.birthProbability[pIdx]         = lambda * tkf91Probs.beta[pIdx];
        tkf91Probs.extinctionProbability[pIdx]    = mu * tkf91Probs.beta[pIdx];
        tkf91Probs.nonHomologousProbability[pIdx] = (1.0 - mu * tkf91Probs.beta[pIdx]) * 
                                                    (1.0 - tkf91Probs.birthProbability[pIdx]) - 
                                                    tkf91Probs.homologousProbability[pIdx];
        tkf91Probs.immortalProbability *= (1.0 - tkf91Probs.birthProbability[pIdx]);
        }
    
#   if defined(DEBUG_TKF91)
    std::cout << "TKF91 event probabilities:" << std::endl;
    std::cout << std::fixed << std::setprecision(6);
    for (size_t i = 0; i < numNodes; i++)
        {
        std::cout << "  Node " << i << ": B=" << tkf91Probs.birthProbability[i]
                  << " E=" << tkf91Probs.extinctionProbability[i]
                  << " H=" << tkf91Probs.homologousProbability[i]
                  << " N=" << tkf91Probs.nonHomologousProbability[i] << std::endl;
        }
    std::cout << "  Immortal prob = " << tkf91Probs.immortalProbability << std::endl;
#   endif
}

void LikelihoodCalculator::setStationaryFrequencies(void) {

    if (myFrequencies == nullptr)
        {
        const double x = 1.0 / numStates;
        for (double* p = equilibriumFrequencies, *end = equilibriumFrequencies + numStates; p < end; p++)
            (*p) = x;
        }
    else
        {
        std::vector<double>& x = myFrequencies->getFrequencies();
        for (double* p = equilibriumFrequencies, *q = x.data(), *end = equilibriumFrequencies + numStates; p < end; p++, q++)
            (*p) = (*q);
        }
}
