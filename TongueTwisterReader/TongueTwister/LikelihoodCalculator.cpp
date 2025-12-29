#include <cmath>
#include <iomanip>
#include <iostream>
#include "IndelMatrix.hpp"
#include "IntVector.hpp"
#include "LikelihoodCalculator.hpp"
#include "Model.hpp"
#include "Msg.hpp"
#include "Node.hpp"
#include "ParameterAlignment.hpp"
#include "TransitionProbabilities.hpp"
#include "Tree.hpp"

#undef DEBUG_LIKELIHOOD



LikelihoodCalculator::LikelihoodCalculator(ParameterAlignment* a, Model* m) {

    // memory addresses of important objects
    data = a;
    modelPtr = m;
    sequences = data->getRawSequenceMatrix();

    numStates = data->getNumStates();
    numStates1 = numStates + 1;
    taxonMask = data->getTaxonMask();
    numTaxa = (int)taxonMask.getNumberSetBits();
    numNodes = 2 * numTaxa - 1;
    stateEquilibriumFrequencies.resize(numStates);
    
    allocateIndelProbabilities(numNodes);
    allcoateIndelCombinatorics(numNodes, data->lengthOfLongestSequence()*5);
    alignment = new IndelMatrix(numTaxa, data->lengthOfLongestSequence()*5); // numSites X numTaxa

    fI = new double* [numNodes];
    fI[0] = new double[numNodes * numStates1];
    fH = new double*[numNodes];
    fH[0] = new double[numNodes * numStates];
    for (int i=1; i<numNodes; i++)
        {
        fH[i] = fH[i-1] + numStates;
        fI[i] = fI[i-1] + numStates1;
        }
}

LikelihoodCalculator::~LikelihoodCalculator(void) {

    drainPool();
    freeIndelCombinatorics();
    freeIndelProbabilities();
    delete [] fH[0];
    delete [] fH;
    delete [] fI[0];
    delete [] fI;
}

std::string& LikelihoodCalculator::alignmentName(void) {

    return data->parmName;
}

void LikelihoodCalculator::allcoateIndelCombinatorics(int nn, int maxSeqLen) {

    indelCombos.maximumSequenceLength = maxSeqLen;
    indelCombos.state = new int[maxSeqLen];
    indelCombos.possibleVectorIndices = new int[maxUnalignableDimension];
    indelCombos.nodeHomology = new int[nn];
    indelCombos.numHomologousEmissions = new int[nn];
    indelCombos.numHomologousEmissionsForClass = new int[maxUnalignableDimension1];
}

void LikelihoodCalculator::allocateIndelProbabilities(int nn) {

    indelProbs.beta = new double[nn];
    indelProbs.birthProbability = new double[nn];
    indelProbs.extinctionProbability = new double[nn];
    indelProbs.homologousProbability = new double[nn];
    indelProbs.nonHomologousProbability = new double[nn];
    for (int i=0; i<nn; i++)
        {
        indelProbs.beta[i] = 0.0;
        indelProbs.birthProbability[i] = 0.0;
        indelProbs.extinctionProbability[i] = 0.0;
        indelProbs.homologousProbability[i] = 0.0;
        indelProbs.nonHomologousProbability[i] = 0.0;
        }
}
        
void LikelihoodCalculator::clearPpTable(void) {

    for (PartialProbabilitiesLookup::iterator it = partialProbabilities.begin(); it != partialProbabilities.end(); it++)
        returnToPool(it->first);
    partialProbabilities.clear();
}


void LikelihoodCalculator::drainPool(void) {

    for (std::vector<IntVector*>::iterator v=pool.begin(); v != pool.end(); v++)
        {
        allocated.erase(*v);
        delete (*v);
        }
}

void LikelihoodCalculator::freeIndelCombinatorics(void) {

    delete [] indelCombos.state;
    delete [] indelCombos.possibleVectorIndices;
    delete [] indelCombos.nodeHomology;
    delete [] indelCombos.numHomologousEmissions;
    delete [] indelCombos.numHomologousEmissionsForClass;
}

void LikelihoodCalculator::freeIndelProbabilities(void) {

    delete [] indelProbs.beta;
    delete [] indelProbs.birthProbability;
    delete [] indelProbs.extinctionProbability;
    delete [] indelProbs.homologousProbability;
    delete [] indelProbs.nonHomologousProbability;
}

IntVector* LikelihoodCalculator::getVector(void) {

    if (pool.empty() == true)
        {
        /* If the vector pool is empty, we allocate a new vector and return it. We
           do not need to add it to the vector pool. */
        IntVector* v = new IntVector(numTaxa);
        allocated.insert(v);
#       if defined (TRACK_ALLOCS)
        numAllocs++;
#       endif
        return v;
        }
    
    // Return a vector from the vector pool, remembering to remove it from the pool.
    IntVector* v = pool.back();
    pool.pop_back();
    return v;
}

IntVector* LikelihoodCalculator::getVector(IntVector& vec) {

    if (pool.empty() == true)
        {
        IntVector* v = new IntVector(numTaxa);
        allocated.insert(v);
        for (int i=0; i<vec.size(); ++i)
            (*v)[i] = vec[i];
#       if defined (TRACK_ALLOCS)
        numAllocs++;
#       endif
        return v;
        }
    IntVector* v = pool.back();
    pool.pop_back();
    for (int i=0; i<vec.size(); ++i)
        (*v)[i] = vec[i];
    return v;
}

void LikelihoodCalculator::initialize(void) {

    //alignment = data->getIndelMatrix();
    data->getIndelMatrix(alignment);
    unalignableRegionSize = 0;
    tree = modelPtr->getTree(taxonMask);

    RbBitSet bs = data->getTaxonMask();
    transitionProbabilities = modelPtr->getTransitionProbabilities()->getTransitionProbabilities(bs);
    std::vector<Node*>& downPassSequence = tree->getDownPassSequence();
    for (int n=0; n<downPassSequence.size(); n++)
        {
        Node* p = downPassSequence[n];
        p->setTransitionProbability( transitionProbabilities[p->getIndex()] );
        }
    modelPtr->getTransitionProbabilities()->getStationaryFrequencies(stateEquilibriumFrequencies);
    setBirthDeathProbabilities();
}

double LikelihoodCalculator::lnLikelihood(MathCache&) {

    //std::cout << "Calculating likelihood" << std::endl;
    initialize();
    
    // null emissions probability
	IntVector* pos = getVector();
    //std::cout << "pos: " << *pos << std::endl;
    double nullEmissionFactor = partialProbability(pos, pos, -1);
    double f = indelProbs.immortalProbability / nullEmissionFactor;
    partialProbabilities.insert( std::make_pair(pos, f) );

    // array of possible vector indices, used in inner loop
    for (int i=0; i<maxUnalignableDimension; i++)
        indelCombos.possibleVectorIndices[i] = 0;
    int firstNotUsed = 0;
    //int len = (int)alignment.size();
    int len = alignment->getNumSites();
    int currentAlignmentColumn = 0;
    if (len > indelCombos.maximumSequenceLength)
        Msg::error("Alignment unexpectedly exceeded 3X the length of the longest sequence");
    for (int i=0; i<len; i++)
        indelCombos.state[i] = 0;
    do
        {
        // Find all possible vectors from current position, pos
        int ptr;
        int numPossible = 0;
        IntVector* mask = getVector();
	    for (ptr=firstNotUsed; mask->zeroEntry() && ptr<len; ptr++)
            {
            if (indelCombos.state[ ptr ] != used)
                {
                if (mask->innerProduct( alignment->getRow(ptr), numTaxa ) == 0)
                    {
                    indelCombos.state[ ptr ] = possible;
                    if (numPossible == maxUnalignableDimension)
                        {
                        // This gets too hairy - bail out.
                        unalignableRegionSize++;
                        Msg::warning("We bailed out because things became too hairy: numPossible = " + std::to_string(numPossible));
                        return -std::numeric_limits<double>::infinity();
                        }
                    indelCombos.possibleVectorIndices[numPossible++] = ptr;
                    }
                //mask->add( alignment[ptr], numTaxa );
                mask->add( alignment->getRow(ptr), numTaxa );
                }
            }
        returnToPool(mask);

        //std::cout << "ptr = " << ptr << " " << *pos << std::endl;
        // Loop over all combinations of possible vectors, which define edges from
        // pos to another possible position, by ordinary binary counting.
        IntVector* newPos = getVector(*pos);
        IntVector* signature = getVector();
	    bool unusedPos, foundNonZero;
        do
            {
            // Find next combination
            foundNonZero = false;
            int site = 0;
			for (int posPtr=numPossible-1; posPtr>=0; --posPtr)
                {
			    int curPtr = indelCombos.possibleVectorIndices[ posPtr ];
                //std::cout << "curPtr = " << curPtr << std::endl;
                site = curPtr;
                currentAlignmentColumn = curPtr;
                if (indelCombos.state[ curPtr ] == possible)
                    {
                    indelCombos.state[ curPtr ] = edgeUsed;
					newPos->add( alignment->getRow(curPtr), numTaxa );
                    // Compute signature vector
					signature->addMultiple( alignment->getRow(curPtr), numTaxa, posPtr + 1 );
                    // Signal: non-zero combination found, and stop
                    foundNonZero = true;
                    posPtr = 0;
                    }
                else
                    {
                    // It was edgeUsed (i.e., digit == 1), so reset digit and continue
					indelCombos.state[ curPtr ] = possible;
					newPos->subtract( alignment->getRow(curPtr), numTaxa );
					signature->addMultiple( alignment->getRow(curPtr), numTaxa, -posPtr - 1 );
                    }
                }

            if (foundNonZero)
                {
                PartialProbabilitiesLookup::iterator it = partialProbabilities.find(pos);
                if (it == partialProbabilities.end())
                    Msg::error("Could not find pos vector in dpTable map.");
                double lft = it->second;
                it = partialProbabilities.find(newPos);
                double rht = 0.0;
                if (it == partialProbabilities.end())
                    {
                    unusedPos = true;
                    rht = 0.0;
                    }
                else
                    {
                    rht = it->second;
                    unusedPos = false;
                    }
                    
                //std::cout << site << " -- " << *signature << " " << *pos << std::endl;
                double transFac = (-partialProbability(signature, pos, site)) / nullEmissionFactor;
                lft *= transFac;
                rht += lft;

                // If we are storing a value at a previously unused position, make sure we use a fresh key object
                if (unusedPos)
                    {
                    IntVector* newPosCopy = getVector(*newPos);
                    partialProbabilities.insert( std::make_pair(newPosCopy, rht) );
                    }
                else
                    {
                    PartialProbabilitiesLookup::iterator it2 = partialProbabilities.find(newPos);
                    if (it2 == partialProbabilities.end())
                        Msg::error("Should have found newPos in table. What happened?");
                    it2->second = rht;
                    //delete rht; // need to remember to delete rht if it's not being inserted back into the table
                    }
                }

            } while (foundNonZero);

    returnToPool(signature);

    // Now find next entry in partial probability table.  Use farthest unused vector
    --ptr;
    while (ptr >= 0 && indelCombos.state[ptr] != possible)
        {
        // Undo any possible used vector that we encounter
        if (indelCombos.state[ptr] == used)
            {
            pos->subtract( alignment->getRow(ptr), numTaxa );
            indelCombos.state[ptr] = freeToUse;
            }
        --ptr;
        }
    
    if (ptr == -1)
        {
        // No more unused vectors, so we also fell through the edge loop above,
        // hence newPos contains the final position
        PartialProbabilitiesLookup::iterator it = partialProbabilities.find(newPos);
        if (it == partialProbabilities.end())
            {
            printTable();
            Msg::error("Could not find newPos in partialProbabilities");
            }
        double lnL = log(it->second);
        if (isinf(lnL) == true)
            {
            data->print();
            printTable();
            }
        clearPpTable();
        returnToPool(newPos);
        if (pool.size() - allocated.size() != 0)
            Msg::warning("Some IntVectors were not returned to the pool");
        return lnL;
        }

    returnToPool(newPos);

    // Now use this farthest-out possible vector
    indelCombos.state[ptr] = used;
    pos->add( alignment->getRow(ptr), numTaxa );
    if (ptr <= firstNotUsed)
        firstNotUsed++;
            
    } while (true);
}

double LikelihoodCalculator::partialProbability(IntVector* signature, IntVector* pos, int site) {

    // make certain that the conditional probabilities for all
    // nodes are zero before we start
    auto fHi = fH;
    auto fIi = fI;

    for (int i=0; i<numNodes; i++)
        {
        memset(*fHi++, 0, numStates*sizeof(double) );
        memset(*fIi++, 0, numStates1*sizeof(double) );
        indelCombos.nodeHomology[i] = 0;
        indelCombos.numHomologousEmissions[i] = 0;
        }

    for (int i=0; i<maxUnalignableDimension1; i++)
        indelCombos.numHomologousEmissionsForClass[i] = 0;

    for (int i=0; i<numTaxa; i++)
        {
        auto sigi = (*signature)[i];

        indelCombos.numHomologousEmissionsForClass[sigi]++;
        indelCombos.nodeHomology[i] = sigi;
        indelCombos.numHomologousEmissions[i] = sigi != 0;
        }
        
    // Loop over all nodes except root, and find out which nodes need carry nucleotides of what
    // homology class.
    bool isHomologyClashing = false;
    std::vector<Node*>& dpSequence = tree->getDownPassSequence();
    for (int n=0; n<dpSequence.size(); n++)
        {
        Node* p = dpSequence[n];
        if (p->getAncestor() != NULL)
            {
            int pIdx = p->getIndex();
            int pAncIdx = p->getAncestor()->getIndex();

            auto numHomologousEmission = indelCombos.numHomologousEmissions[pIdx];
            auto nodeHomologyI = indelCombos.nodeHomology[pIdx];

            if (numHomologousEmission == indelCombos.numHomologousEmissionsForClass[nodeHomologyI] || nodeHomologyI == 0)
                {
                // This node is the MRCA of the homologous nucleotides iHom[i], or a gap --- do nothing
                }
            else
                {
                auto nodeHomologyAnc = indelCombos.nodeHomology[pAncIdx];
                if (nodeHomologyAnc == 0 || nodeHomologyAnc == nodeHomologyI)
                    {
                    // This node is not yet MRCA, so node above must be homologous
                    indelCombos.nodeHomology[pAncIdx] = nodeHomologyI;
                    // If iHom[i]==0, iHomNum[i]==0.
                    indelCombos.numHomologousEmissions[pAncIdx] += numHomologousEmission;
                    }
                else
                    {
                    // Clashing homology; signal
                    isHomologyClashing = true;
                    }
                }
            }
        }
    if (isHomologyClashing)
        return 0.0;
    return prune(signature, pos, dpSequence, site);
}

void LikelihoodCalculator::pruneBranch(Node* p, int /*site*/) {

    int pIdx = p->getIndex();
    auto fIIdx = fI[pIdx];
    auto fHIdx = fH[pIdx];
    auto nodeHomologyI = indelCombos.nodeHomology[pIdx];

    // calculate conditional probabilities for interior nodes
    p->getDescendantsVector(des);
    if (des.size() != 2)
        Msg::error("Expecting two descendants");
    Node* pLft = des[0];
    Node* pRht = des[1];
    int lftChildIdx = pLft->getIndex();
    int rhtChildIdx = pRht->getIndex();
    DoubleMatrix& tpLft = *(pLft->getTransitionProbability());
    DoubleMatrix& tpRht = *(pRht->getTransitionProbability());

    if ((nodeHomologyI != 0) && (nodeHomologyI == indelCombos.nodeHomology[lftChildIdx]) && (nodeHomologyI == indelCombos.nodeHomology[rhtChildIdx]))
        {
        // Case 1: One homology family spanning tree intersects both child edges, i.e.
        //         a 'homologous' nucleotide should travel down both edges.

        auto fHleft  = fH[lftChildIdx];
        auto fHright = fH[rhtChildIdx];
        auto probLeft  = indelProbs.homologousProbability[lftChildIdx];
        auto probRight = indelProbs.homologousProbability[rhtChildIdx];

        auto tpLeft = tpLft.begin();
        auto tpRight = tpRht.begin();

        for (int i=0; i<numStates; i++)
            {
            auto lj = fHleft;
            auto rj = fHright;

            double lft = 0.0, rht = 0.0;
            for (int j=0; j<numStates; j++)
                {
                lft += *lj++ * probLeft  * *tpLeft++;
                rht += *rj++ * probRight * *tpRight++;
                }
            fHIdx[i] = lft * rht;
            }
        }
    else if (nodeHomologyI != 0)
        {
        // Case 2: One homology family Spanning tree intersects one of the child edges, i.e.
        //        'homologous' nucleotide must travel down a specific edge.
        int homologousChildIdx, inhomologousChildIdx;
        DoubleMatrix* tpHomologous   = NULL;
        DoubleMatrix* tpInhomologous = NULL;
        if (nodeHomologyI == indelCombos.nodeHomology[lftChildIdx])
            {
            homologousChildIdx   = lftChildIdx;
            inhomologousChildIdx = rhtChildIdx;
            tpHomologous   = pLft->getTransitionProbability();
            tpInhomologous = pRht->getTransitionProbability();
            }
        else
            {
            homologousChildIdx   = rhtChildIdx;
            inhomologousChildIdx = lftChildIdx;
            tpHomologous   = pRht->getTransitionProbability();
            tpInhomologous = pLft->getTransitionProbability();
            }

        double birthProb            = indelProbs.birthProbability[inhomologousChildIdx];
        double extinctionProb       = indelProbs.extinctionProbability[inhomologousChildIdx];
        double homologousProb_hIdx  = indelProbs.homologousProbability[homologousChildIdx];
        double homologousProb_ihIdx = indelProbs.homologousProbability[inhomologousChildIdx];
        double nonHomologousProb    = indelProbs.nonHomologousProbability[inhomologousChildIdx];
        double rhtStartProb         = extinctionProb * fI[inhomologousChildIdx][numStates];
        double* tiHomologous        = &(*tpHomologous)(0,0); // we assume row-major order
        double* tiNonHomologous     = &(*tpInhomologous)(0,0);
        double* freqsStart          = &stateEquilibriumFrequencies[0];
        double* freqsEnd            = freqsStart + numStates;
        double* clH_hIdxStart       = fH[homologousChildIdx];
        double* clH_ihIdxStart      = fH[inhomologousChildIdx];
        double* clI_ihIdxStart      = fI[inhomologousChildIdx];
        double* clH_start           = fHIdx;
        double* clH_end             = clH_start + numStates;
        
        for (double* clH = clH_start; clH != clH_end; clH++)
            {
            double* clH_hIdx  = clH_hIdxStart;
            double* clH_ihIdx = clH_ihIdxStart;
            double* clI_ihIdx = clI_ihIdxStart;
            double lft = 0.0;
            double rht = rhtStartProb;
            for (double* freq = freqsStart; freq != freqsEnd; freq++)
                {
                lft += (*clH_hIdx) * homologousProb_hIdx * (*tiHomologous);
                rht += ((*clH_ihIdx) + (*clI_ihIdx)) * (nonHomologousProb - extinctionProb * birthProb) * (*freq) +
                       (*clI_ihIdx) * homologousProb_ihIdx * (*tiNonHomologous);
                       
                tiHomologous++;
                tiNonHomologous++;
                clH_hIdx++;
                clH_ihIdx++;
                clI_ihIdx++;
                }
            *clH = lft * rht;
            }
            
        }
    else
        {
        // Case 3: No homology family's spanning tree intersect either child edge, i.e.
        //         a 'homologous' nucleotide may travel down either edge and pop out at a leaf,
        //         or may travel down both edges but has to die in at least one subtree, or
        //         an 'inhomologous' nucleotide may do whatever it likes here.
        auto extinctionProbabilityLeft  = indelProbs.extinctionProbability[lftChildIdx];
        auto extinctionProbabilityRight = indelProbs.extinctionProbability[rhtChildIdx];
        auto fILeft  = fI[lftChildIdx];
        auto fIRight = fI[rhtChildIdx];
        auto fHLeft  = fH[lftChildIdx];
        auto fHRight = fH[rhtChildIdx];
        auto homologousProbabilityLeft  = indelProbs.homologousProbability[lftChildIdx];
        auto homologousProbabilityRight = indelProbs.homologousProbability[rhtChildIdx];
        auto birthProbabilityLeft  = indelProbs.birthProbability[lftChildIdx];
        auto birthProbabilityRight = indelProbs.birthProbability[rhtChildIdx];
        auto nonHomologousProbabilityLeft  = indelProbs.nonHomologousProbability[lftChildIdx];
        auto nonHomologousProbabilityRight = indelProbs.nonHomologousProbability[rhtChildIdx];
        auto rfactor = nonHomologousProbabilityRight - extinctionProbabilityRight * birthProbabilityRight;
        auto lfactor = nonHomologousProbabilityLeft  - extinctionProbabilityLeft  * birthProbabilityLeft;

        for (int i = 0; i < numStates; i++)
            {
            double lft1 = 0.0, lft2 = 0.0;
            double rht1 = extinctionProbabilityLeft * fILeft[numStates];
            double rht2 = extinctionProbabilityRight * fIRight[numStates];
            auto fhleft  = fHLeft;
            auto fhright = fHRight;
            auto fileft  = fILeft;
            auto firight = fIRight;

            for (int j = 0; j < numStates; j++)
                {
                lft1 += *fhleft  * homologousProbabilityLeft  * tpLft(i, j);
                lft2 += *fhright * homologousProbabilityRight * tpRht(i, j);
                auto stateEquilibriumFrequency = stateEquilibriumFrequencies[j];
                rht1 += (*fhleft  + *fileft)  * lfactor * stateEquilibriumFrequency + *fileft  * homologousProbabilityLeft  * tpLft(i, j);
                rht2 += (*fhright + *firight) * rfactor * stateEquilibriumFrequency + *firight * homologousProbabilityRight * tpRht(i, j);
                ++fhleft;
                ++fhright;
                ++firight;
                }

            fHIdx[i] = lft1 * rht2 + lft2 * rht1;
            fIIdx[i] = rht1 * rht2;
            }

        auto left  = fILeft[numStates];
        auto right = fIRight[numStates];
        auto fIL = fILeft;
        auto fHL = fHLeft;
        auto fIR = fIRight;
        auto fHR = fHRight;

        for (int j=0; j<numStates; j++)
            {
            auto stateEquilibriumFrequency = stateEquilibriumFrequencies[j];
            left  -= birthProbabilityLeft  * (*fIL + *fHL) * stateEquilibriumFrequency;
            right -= birthProbabilityRight * (*fIR + *fHR) * stateEquilibriumFrequency;
            ++fIL;
            ++fHL;
            ++fIR;
            ++fHR;
            }

        fIIdx[numStates] = left * right;
        }
}

double LikelihoodCalculator::prune(IntVector* signature, IntVector* pos, std::vector<Node*>& dpSequence, int site) {
 
    // pruning algorithm
    auto nsize = dpSequence.size();
    for (int n=0; n<nsize; n++)
        {
        auto p = dpSequence[n];
        if (p->getIsLeaf() == true)
            {
            int pIdx = p->getIndex();
            auto fIIdx = fI[pIdx];
            auto fHIdx = fH[pIdx];
            
            //std::cout << "(*pos)[pIdx] = " << (*pos)[pIdx] << std::endl;

            // initialize conditional probabilties for the tip
            if ( (*signature)[pIdx] == 0 )
                fIIdx[numStates] = 1.0;  // gap
            else
                fHIdx[ sequences[pIdx][(*pos)[pIdx]] ] = 1.0; // nucleotide
            }
        else
            pruneBranch(p, site);
        }

    // average probability over states at the root
    int rootIdx = tree->getRoot()->getIndex();
    auto fip = fI[rootIdx];
    double res = fip[numStates];
    auto bpr = indelProbs.birthProbability[rootIdx];
    auto fhr = fH[rootIdx];

    for (int i=0; i<numStates; i++) 
        {
        res -= (*fip + *fhr) * bpr * stateEquilibriumFrequencies[i];
        ++fip;
        ++fhr;
        }
    return res;
}

void LikelihoodCalculator::printTable(void) {

    std::cout << "Table" << std::endl;
    std::cout << std::fixed << std::setprecision(10) << std::scientific;
    int i = 0;
    for (PartialProbabilitiesLookup::iterator it = partialProbabilities.begin(); it != partialProbabilities.end(); it++)
        {
        std::cout << i++ << " -- " << *(it->first) << " -- " << it->second << std::endl;
        }
    Msg::error("Found infinite likelihood");
}

void LikelihoodCalculator::returnToPool(IntVector* v) {

    v->clean();
    pool.push_back(v);
}

void LikelihoodCalculator::setBirthDeathProbabilities(void) {

    indelProbs.insertionRate = modelPtr->getInsertionRate();
    indelProbs.deletionRate  = modelPtr->getDeletionRate();

    indelProbs.immortalProbability= 1.0;
    std::vector<Node*>& dpSequence = tree->getDownPassSequence();
    for (int n=0; n<dpSequence.size(); n++)
        {
        Node* p = dpSequence[n];
        int pIdx = p->getIndex();
        double v = p->getBranchLength();
        if (p->getAncestor() == NULL)
            {
            // root
            indelProbs.beta[pIdx] = 1.0 / indelProbs.deletionRate;
            indelProbs.homologousProbability[pIdx] = 0.0;
            }
        else
            {
            // internal node or tip
            indelProbs.beta[pIdx] = exp( (indelProbs.insertionRate - indelProbs.deletionRate) * v );
            indelProbs.beta[pIdx] = (1.0 - indelProbs.beta[pIdx]) / (indelProbs.deletionRate - indelProbs.insertionRate * indelProbs.beta[pIdx]);
            indelProbs.homologousProbability[pIdx] = exp( -indelProbs.deletionRate * v) * (1.0 - indelProbs.insertionRate * indelProbs.beta[pIdx] );
            }
        indelProbs.birthProbability[pIdx]         = indelProbs.insertionRate * indelProbs.beta[pIdx];
        indelProbs.extinctionProbability[pIdx]    = indelProbs.deletionRate * indelProbs.beta[pIdx];
        indelProbs.nonHomologousProbability[pIdx] = (1.0 - indelProbs.deletionRate * indelProbs.beta[pIdx]) * (1.0 - indelProbs.birthProbability[pIdx]) - indelProbs.homologousProbability[pIdx];
        indelProbs.immortalProbability *= (1.0 - indelProbs.birthProbability[pIdx]);
        }
        
#   if defined(DEBUG_LIKELIHOOD)
    std::cout << "TKF91 event probabilties" << std::endl;
    std::cout << std::fixed << std::setprecision(10);
    for (int i=0; i<numNodes; i++)
        {
        std::cout << i << " -- " << birthProbability[i] << " " << extinctionProbability[i] << " " << homologousProbability[i] << " " << nonHomologousProbability[i] << std::endl;
        }
#   endif
}
