#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include "EigenSystem.hpp"
#include "IntVector.hpp"
#include "Msg.hpp"
#include "LikelihoodTkf.hpp"
#include "Model.hpp"
#include "Node.hpp"
#include "ParameterAlignment.hpp"
#include "SiteLikelihood.hpp"
#include "TransitionProbabilities.hpp"
#include "Tree.hpp"
#include "UserSettings.hpp"



#undef DEBUG_TKF91
int LikelihoodTkf::unalignableRegionSize = 0;
int LikelihoodTkf::maxUnalignableDimension = 10;

/* This is code from Beast (now defunct in Beast) for calculating the likelihood
   under the TKF91 model. The main changes were made in the hash map called iTable
   in the original code and called dpTable here. I use a map, sorting on the vector
   contents for the key. The original code hashes the vectors by the vector elements,
   so both should be equivalent, here. I also use a high performance arbitrary
   precision double in place of the roll-it-yourself class (BFloat) used in Beast. */



LikelihoodTkf::LikelihoodTkf(ParameterAlignment* a, Model* m) {

    data = a;
    tree = NULL;
    model = m;
        
    siteProbs = data->getSiteProbs();
    fH = siteProbs->getProbsH();
    fI = siteProbs->getProbsI();
    
    init();
}

LikelihoodTkf::~LikelihoodTkf(void) {

    clearDpTable();
}

void LikelihoodTkf::clearDpTable(void) {

    dpTable.clear();
}

void LikelihoodTkf::init(void) {

    // initialize some useful variables
    numStates = data->getNumStates();
    taxonMask = data->getTaxonMask();
    
    // initialize the tree
    initTree();

    // initialize the alignment array from the given alignment.
    initAlignment();

    // initialize the sequences array from the given alignment.
    initSequences();

    // initialize the transitionProbabilities array from the substitution model
    initTransitionProbabilities();

    // initialise TKF91 probabilities
    initTKF91();
    
    //initPrint();
}

void LikelihoodTkf::initAlignment(void) {
            
    alignment = data->getIndelMatrix();
            
#   if 0
    std::cout << "Variable = alignment" << std::endl;
    std::cout << "name = " << data->getName() << std::endl;
    for (int i=0; i<alignment.size(); i++)
        {
        std::cout << std::setw(3) << i << " -- ";
        for (int j=0; j<alignment[i].size(); j++)
            std::cout << alignment[i][j] << " ";
        std::cout << std::endl;
        }
#   endif
}

void LikelihoodTkf::initPrint(void) {

    for (int i=0; i<60; i++)
        std::cout << "-";
    std::cout << std::endl;
    std::cout << "data-.getName() = " << data->getName() << std::endl;
    std::cout << "taxonMask = " << taxonMask.bitString() << std::endl;
    std::cout << "numNodes = " << siteProbs->getNumNodes() << std::endl;
    data->print();
    tree->print();
    std::cout << "alignment" << std::endl;
    data->print();
    for (int i=0; i<alignment.size(); i++)
        {
        std::cout << std::setw(3) << i << " -- ";
        for (int j=0; j<alignment[i].size(); j++)
            std::cout << alignment[i][j] << " ";
        std::cout << std::endl;
        }
    std::cout << "sequences" << std::endl;
    for (int i=0; i<sequences.size(); i++)
        {
        std::cout << std::setw(3) << i << " -- ";
        for (int j=0; j<sequences[i].size(); j++)
            std::cout << std::setw(2) << sequences[i][j] << " ";
        std::cout << std::endl;
        }
    std::cout << "numStates = " << numStates << std::endl;
    std::cout << "numIndelCategories = " << numIndelCategories << std::endl;
}

void LikelihoodTkf::initSequences(void) {
    
    sequences = data->getRawSequenceMatrix();

#   if 0
    std::cout << "Variable = sequences" << std::endl;
    for (int i=0; i<sequences.size(); i++)
        {
        std::cout << std::setw(3) << i << " -- ";
        for (int j=0; j<sequences[i].size(); j++)
            std::cout << std::setw(2) << sequences[i][j] << " ";
        std::cout << std::endl;
        }
#   endif
}

void LikelihoodTkf::initTransitionProbabilities(void) {
        
    TransitionProbabilities& tip = TransitionProbabilities::transitionProbabilties();
    RbBitSet bs = data->getTaxonMask();
    transitionProbabilities = tip.getTransitionProbabilities(bs);
    stateEquilibriumFrequencies = tip.getStationaryFrequencies();
    
    // we have an alignment of word segments for the complete set of canonical taxa
    std::vector<Node*>& downPassSequence = tree->getDownPassSequence();
    for (int n=0; n<downPassSequence.size(); n++)
        {
        Node* p = downPassSequence[n];
        p->setTransitionProbability( transitionProbabilities[p->getIndex()] );
        }
}

void LikelihoodTkf::initTKF91(void) {

    // initialize lambda and mu
    insertionRate = model->getInsertionRate();
    deletionRate  = model->getDeletionRate();
    
    // initialize parameters of TKF91 model
    UserSettings& settings = UserSettings::userSettings();
    numIndelCategories = settings.getNumIndelCategories();

    // size the vectors for the indel portion of the process
    int numNodes = tree->getNumNodes();
    std::vector<double> beta(numNodes);
    birthProbability.resize(numNodes);
    extinctionProbability.resize(numNodes);
    homologousProbability.resize(numNodes);
    nonHomologousProbability.resize(numNodes);
        
    immortalProbability= 1.0;
    std::vector<Node*> dpSequence = tree->getDownPassSequence();
    for (int n=0; n<dpSequence.size(); n++)
        {
        Node* p = dpSequence[n];
        int pIdx = p->getIndex();
        double v = p->getBranchLength();
        if (p->getAncestor() == NULL)
            {
            // root
            beta[pIdx] = 1.0 / deletionRate;
            homologousProbability[pIdx] = 0.0;
            }
        else
            {
            // internal node or tip
            beta[pIdx] = exp( (insertionRate - deletionRate) * v );
            beta[pIdx] = (1.0 - beta[pIdx]) / (deletionRate - insertionRate * beta[pIdx]);
            homologousProbability[pIdx] = exp( -deletionRate * v) * (1.0 - insertionRate * beta[pIdx] );
            }
        birthProbability[pIdx]         = insertionRate * beta[pIdx];
        extinctionProbability[pIdx]    = deletionRate * beta[pIdx];
        nonHomologousProbability[pIdx] = (1.0 - deletionRate * beta[pIdx]) * (1.0 - birthProbability[pIdx]) - homologousProbability[pIdx];
        immortalProbability *= (1.0 - birthProbability[pIdx]);
        }
        
#   if defined(DEBUG_TKF91)
    std::cout << "TKF91 event probabilties" << std::endl;
    std::cout << std::fixed << std::setprecision(10);
    for (int i=0; i<numNodes; i++)
        {
        std::cout << i << " -- " << birthProbability[i] << " " << extinctionProbability[i] << " " << homologousProbability[i] << " " << nonHomologousProbability[i] << std::endl;
        }
#   endif
}

void LikelihoodTkf::initTree(void) {

    tree = model->getTree(taxonMask);
    //tree->print();
}

void LikelihoodTkf::printTable(void) {

    std::cout << "dpTable" << std::endl;
    std::cout << std::fixed << std::setprecision(10) << std::scientific;
    int i = 0;
    for (TkfLookup::iterator it = dpTable.begin(); it != dpTable.end(); it++)
        {
        std::cout << i++ << " -- " << (it->first) << " -- " << it->second << std::endl;
        }
}

void LikelihoodTkf::printVector(std::string header, std::vector<int>& v) {

    if (header != "")
        std::cout << header << std::endl;
    for (int i=0; i<v.size(); i++)
        std::cout << v[i] << " ";
    std::cout << std::endl;
}

void LikelihoodTkf::printVector(std::string header, std::vector<double>& v) {

    std::cout << std::fixed << std::setprecision(10);
    
    if (header != "")
        std::cout << header << std::endl;
    for (int i=0; i<v.size(); i++)
        std::cout << v[i] << " ";
    std::cout << std::endl;
}

void LikelihoodTkf::printVector(std::string header, std::vector< std::vector<double> >& v) {

    std::cout << std::fixed << std::setprecision(10);

    if (header != "")
        std::cout << header << std::endl;
    for (int i=0; i<v.size(); i++)
        {
        for (int j=0; j<v[i].size(); j++)
            std::cout << v[i][j] << " ";
        std::cout << std::endl;
        }
}

void LikelihoodTkf::printVector(std::string header, std::vector< std::vector<int> >& v) {

    if (header != "")
        std::cout << header << std::endl;
    for (int i=0; i<v.size(); i++)
        {
        for (int j=0; j<v[i].size(); j++)
            std::cout << std::setw(2) << v[i][j] << " ";
        std::cout << std::endl;
        }
}

double LikelihoodTkf::tkfLike(void) {
            
    int len = (int)alignment.size();
    int numLeaves = (int)alignment[0].size();
    int firstNotUsed = 0;                                    // First not-'used' alignment vector (for efficiency)
    std::vector<int> state(len);                             // Helper array, to traverse the region in the DP table corresp. to the alignment
	IntVector pos(numLeaves);                                // Current position; sum of all used vectors
    int currentAlignmentColumn = 0;                          // Keep track of which alignment column is being calculated
    
    // Calculate correction factor for null emissions ("wing folding", or linear equation solving.)
    double nullEmissionFactor = treeRecursion( &pos, &pos, -1 );
    double f = immortalProbability / nullEmissionFactor;
//    mpfr::mpreal f = new mpfr::mpreal;
//    mpfr::mpreal  f = immortalProbability / nullEmissionFactor;
    dpTable.insert( std::make_pair(IntVector(pos), f) );
    //printTable();
 
    // Array of possible vector indices, used in inner loop
    std::vector<int> possibleVectorIndices(maxUnalignableDimension, 0);

    do
        {
        // Find all possible vectors from current position, pos
	    IntVector mask(numLeaves);
        int ptr;
        int numPossible = 0;
	    for (ptr=firstNotUsed; mask.zeroEntry() && ptr<len; ptr++)
            {
            if (state[ ptr ] != used)
                {
                if (mask.innerProduct( alignment[ptr] ) == 0)
                    {
                    state[ ptr ] = possible;
                    if (numPossible == maxUnalignableDimension)
                        {
                        // This gets too hairy - bail out.
                        unalignableRegionSize++;
                        Msg::warning("We bailed out because things became too hairy: numPossible = " + std::to_string(numPossible));
                        return -std::numeric_limits<double>::infinity();
                        }
                    possibleVectorIndices[numPossible++] = ptr;
                    }
                mask.add( alignment[ptr] );
                }
            }

        // Loop over all combinations of possible vectors, which define edges from
        // pos to another possible position, by ordinary binary counting.
	    IntVector newPos(pos);
	    IntVector signature(pos.size());
	    bool unusedPos, foundNonZero;
        do
            {
            // Find next combination
            foundNonZero = false;
			for (int posPtr=numPossible-1; posPtr>=0; --posPtr)
                {
			    int curPtr = possibleVectorIndices[ posPtr ];
                currentAlignmentColumn = curPtr;
                if (state[ curPtr ] == possible)
                    {
                    state[ curPtr ] = edgeUsed;
					newPos.add( alignment[ curPtr ] );
                    // Compute signature vector
					signature.addMultiple( alignment[ curPtr ], posPtr+1 );
                    // Signal: non-zero combination found, and stop
                    foundNonZero = true;
                    posPtr = 0;
                    }
                else
                    {
                    // It was edgeUsed (i.e., digit == 1), so reset digit and continue
					state[ curPtr ] = possible;
					newPos.subtract( alignment[ curPtr ] );
					signature.addMultiple( alignment[ curPtr ], -posPtr-1 );
                    }
                }

            if (foundNonZero)
                {
                TkfLookup::iterator it = dpTable.find(pos);
                if (it == dpTable.end())
                    Msg::error("Could not find pos vector in dpTable map.");
                double lft = it->second;
//                mpfr::mpreal lft = it->second;
                it = dpTable.find(newPos);
                double rht = 0.0;
                //mpfr::mpreal* rht = new mpfr::mpreal;
//                mpfr::mpreal rht = 0.0;
                if (it == dpTable.end())
                    {
                    unusedPos = true;
                    rht = 0.0;
                    }
                else
                    {
                    rht = it->second;
                    unusedPos = false;
                    }
                    
//                std::cout << "currentAlignmentColumn = " << currentAlignmentColumn << std::endl;
                double transFac = (-treeRecursion( &signature, &pos, currentAlignmentColumn )) / nullEmissionFactor;
                lft *= transFac;
                rht += lft;

                // If we are storing a value at a previously unused position, make sure we use a fresh key object
                if (unusedPos)
                    {
                    dpTable.insert( std::make_pair(IntVector(newPos), rht) );
                    //printTable();
                    }
                else
                    {
                    TkfLookup::iterator it2 = dpTable.find(newPos);
                    if (it2 == dpTable.end())
                        Msg::error("Should have found newPos in table. What happened?");
                    it2->second = rht;
                    //printTable();
                    //delete rht; // need to remember to delete rht if it's not being inserted back into the table
                    }
                }

            } while (foundNonZero);

    // Now find next entry in DP table.  Use farthest unused vector
    --ptr;
    while (ptr >= 0 && state[ptr] != possible)
        {
        // Undo any possible used vector that we encounter
        if (state[ptr] == used)
            {
            pos.subtract( alignment[ptr] );
            state[ptr] = free;
            }
        --ptr;
        }
    
    if (ptr == -1)
        {
        // No more unused vectors, so we also fell through the edge loop above,
        // hence newPos contains the final position
        TkfLookup::iterator it = dpTable.find(newPos);
        if (it == dpTable.end())
            Msg::error("Could not find newPos in dpTable");
        
//        mpfr::mpreal res = log(it->second, MPFR_RNDN);
//        double lnL = res.toDouble();
        double lnL = log(it->second);
        
        if (isinf(lnL) == true)
            {
            std::cout << "res = " << it->second << std::endl;
            data->print();
            printTable();
            }

        clearDpTable();
        return lnL;
        }
    
    // Now use this farthest-out possible vector
    state[ptr] = used;
    pos.add( alignment[ptr] );
    if (ptr <= firstNotUsed)
        firstNotUsed++;
            
    } while (true);
}

double LikelihoodTkf::treeRecursion(IntVector* signature, IntVector* pos, int /*siteColumn*/) {

    int numLeaves = (int)signature->size();
    int numNodes = tree->getNumNodes();

    // make certain that the conditional probabilities for all
    // nodes are zero before we start
    siteProbs->zeroOutH();
    siteProbs->zeroOutI();
        
    std::vector<int> nodeHomology(numNodes);                                      // Homology for every node; 0 if node need not be homologous to an emitted nucleotide
    std::vector<int> numHomologousEmissions(numNodes);                            // Number of homologous emissions accounted for by homologous nucleotide at this node
    std::vector<int> numHomologousEmissionsForClass(maxUnalignableDimension + 1); // Number of emissions for each class of homologous nucleotides
    for (int i=0; i<numLeaves; i++)
        {
        numHomologousEmissionsForClass[ (*signature)[ i ] ]++;
        nodeHomology[ i ] = (*signature)[ i ];
        if ( (*signature)[i] == 0)
            numHomologousEmissions[ i ] = 0;
        else
            numHomologousEmissions[ i ] = 1;
        }
        
    // Loop over all nodes except root, and find out which nodes need carry nucleotides of what
    // homology class.
    bool isHomologyClashing = false;
    std::vector<Node*> dpSequence = tree->getDownPassSequence();
    for (int n=0; n<dpSequence.size(); n++)
        {
        Node* p = dpSequence[n];
        if (p->getAncestor() != NULL)
            {
            int pIdx = p->getIndex();
            int pAncIdx = p->getAncestor()->getIndex();
            
            if ((numHomologousEmissions[pIdx] == numHomologousEmissionsForClass[ nodeHomology[pIdx] ]) || (nodeHomology[pIdx] == 0))
                {
                // This node is the MRCA of the homologous nucleotides iHom[i], or a gap - do nothing
                }
            else
                {
                if (nodeHomology[pAncIdx] == 0 || nodeHomology[pAncIdx] == nodeHomology[pIdx])
                    {
                    // This node is not yet MRCA, so node above must be homologous
                    nodeHomology[pAncIdx] = nodeHomology[pIdx];
                    // If iHom[i]==0, iHomNum[i]==0.
                    numHomologousEmissions[pAncIdx] += numHomologousEmissions[pIdx];
                    }
                else
                    {
                    // Clashing homology; signal
                    isHomologyClashing = true;
                    }
                }
            }
        }

    // argh
    if (isHomologyClashing)
        return 0.0;

//    std::cout << *signature << std::endl;
    
    // pruning algorithm
    for (int n=0; n<dpSequence.size(); n++)
        {
        Node* p = dpSequence[n];
        int pIdx = p->getIndex();
        if (p->getIsLeaf() == true)
            {
            // initialize conditional probabilties for the tip
            if ( (*signature)[pIdx] == 0 )
                {
                // gap
                fI[pIdx][numStates] = 1.0;
                }
            else
                {
                // nucleotide
                fH[pIdx][ sequences[pIdx][(*pos)[pIdx]] ] = 1.0;
                }
            }
        else
            {
            // calculate conditional probabilities for interior nodes
            std::vector<Node*> des = p->getDescendantsVector();
            if (des.size() != 2)
                Msg::error("Expecting two descendants");
            Node* pLft = des[0];
            Node* pRht = des[1];
            int lftChildIdx = pLft->getIndex();
            int rhtChildIdx = pRht->getIndex();
            StateMatrix_t* tpLft = pLft->getTransitionProbability();
            StateMatrix_t* tpRht = pRht->getTransitionProbability();
            
            if ( (nodeHomology[pIdx] != 0) && (nodeHomology[pIdx] == nodeHomology[lftChildIdx]) && (nodeHomology[pIdx] == nodeHomology[rhtChildIdx]) )
                {
                // Case 1: One homology family spanning tree intersects both child edges, i.e.
                //         a 'homologous' nucleotide should travel down both edges.
                for (int i=0; i<numStates; i++)
                    {
                    double lft = 0.0;
                    double rht = 0.0;
                    for (int j=0; j<numStates; j++)
                        {
                        lft += fH[lftChildIdx][j] * homologousProbability[lftChildIdx] * (*tpLft)(i,j);
                        rht += fH[rhtChildIdx][j] * homologousProbability[rhtChildIdx] * (*tpRht)(i,j);
                        }
                    fH[pIdx][i] = lft * rht;
                    }
                    // Others: 0.0
                }
            else if (nodeHomology[pIdx] != 0)
                {
                // Case 2: One homology family Spanning tree intersects one of the child edges, i.e.
                //        'homologous' nucleotide must travel down a specific edge.
                int homologousChildIdx, inhomologousChildIdx;
                StateMatrix_t* tpHomologous = NULL;
                StateMatrix_t* tpInhomologous = NULL;
                if (nodeHomology[pIdx] == nodeHomology[lftChildIdx])
                    {
                    homologousChildIdx = lftChildIdx;
                    inhomologousChildIdx = rhtChildIdx;
                    tpHomologous = pLft->getTransitionProbability();
                    tpInhomologous = pRht->getTransitionProbability();
                    }
                else
                    {
                    homologousChildIdx = rhtChildIdx;
                    inhomologousChildIdx = lftChildIdx;
                    tpHomologous = pRht->getTransitionProbability();
                    tpInhomologous = pLft->getTransitionProbability();
                    }

                for (int i=0; i<numStates; i++)
                    {
                    double lft = 0.0;
                    double rht = extinctionProbability[inhomologousChildIdx] * fI[inhomologousChildIdx][numStates];
                    for (int j=0; j<numStates; j++)
                        {
                        lft += fH[homologousChildIdx][j] * homologousProbability[homologousChildIdx] * (*tpHomologous)(i,j);
                        rht +=
                            (fH[inhomologousChildIdx][j] + fI[inhomologousChildIdx][j]) * (nonHomologousProbability[inhomologousChildIdx] - extinctionProbability[inhomologousChildIdx] * birthProbability[inhomologousChildIdx]) * stateEquilibriumFrequencies[j] +
                            fI[inhomologousChildIdx][j] * homologousProbability[inhomologousChildIdx] * (*tpInhomologous)(i,j);
                        }
                    fH[pIdx][i] = lft * rht;
                    }
                // Others: 0.0
                }
            else
                {
                // Case 3: No homology family's spanning tree intersect either child edge, i.e.
                //         a 'homologous' nucleotide may travel down either edge and pop out at a leaf,
                //         or may travel down both edges but has to die in at least one subtree, or
                //         an 'inhomologous' nucleotide may do whatever it likes here.
                for (int i=0; i<numStates; i++)
                    {
                    double lft1 = 0.0, lft2 = 0.0;
                    double rht1 = extinctionProbability[lftChildIdx] * fI[lftChildIdx][numStates];
                    double rht2 = extinctionProbability[rhtChildIdx] * fI[rhtChildIdx][numStates];
                    for (int j=0; j<numStates; j++)
                        {
                        lft1 += fH[lftChildIdx][j] * homologousProbability[lftChildIdx] * (*tpLft)(i,j);
                        lft2 += fH[rhtChildIdx][j] * homologousProbability[rhtChildIdx] * (*tpRht)(i,j);
                        rht1 += (fH[lftChildIdx][j] + fI[lftChildIdx][j]) * (nonHomologousProbability[lftChildIdx] - extinctionProbability[lftChildIdx] * birthProbability[lftChildIdx]) * stateEquilibriumFrequencies[j] + fI[lftChildIdx][j] * homologousProbability[lftChildIdx] * (*tpLft)(i,j);
                        rht2 += (fH[rhtChildIdx][j] + fI[rhtChildIdx][j]) * (nonHomologousProbability[rhtChildIdx] - extinctionProbability[rhtChildIdx] * birthProbability[rhtChildIdx]) * stateEquilibriumFrequencies[j] + fI[rhtChildIdx][j] * homologousProbability[rhtChildIdx] * (*tpRht)(i,j);
                        }
                    fH[pIdx][i] = lft1 * rht2 + lft2 * rht1;
                    fI[pIdx][i] = rht1 * rht2;
                    }
                double lft = fI[lftChildIdx][numStates];
                double rht = fI[rhtChildIdx][numStates];
                for (int j=0; j<numStates; j++)
                    {
                    lft -= birthProbability[lftChildIdx] * (fI[lftChildIdx][j] + fH[lftChildIdx][j]) * stateEquilibriumFrequencies[j];
                    rht -= birthProbability[rhtChildIdx] * (fI[rhtChildIdx][j] + fH[rhtChildIdx][j]) * stateEquilibriumFrequencies[j];
                    }
                fI[pIdx][numStates] = lft * rht;
                }
            
            }
        }

    // average probability over states at the root
    int rootIdx = tree->getRoot()->getIndex();
    double res = fI[rootIdx][numStates];
    for (int i=0; i<numStates; i++)
        res -= (fI[rootIdx][i] + fH[rootIdx][i]) * birthProbability[rootIdx] * stateEquilibriumFrequencies[i];
    return res;
}


