#include <cmath>
#include <iomanip>
#include <iostream>
#include <set>
#include "Alignment.hpp"
#include "Model.hpp"
#include "Msg.hpp"
#include "Node.hpp"
#include "NodeSet.hpp"
#include "ParameterTree.hpp"
#include "Probability.hpp"
#include "RandomVariable.hpp"
#include "TransitionProbabilities.hpp"
#include "Tree.hpp"
#include "UserSettings.hpp"

#undef DEBUG_LOCAL



ParameterTree::ParameterTree(RandomVariable* r, Model* m, std::string treeStr, std::vector<std::string> tNames, std::vector<Alignment*>& wordAlignments, double itl) : Parameter(r, m, "tree") {

    //parmId = TreeParm;

    updateChangesRateMatrix = false;
    lambda = itl;
    isClockConstrained = UserSettings::userSettings().getUseClockConstraint();
    
    if (treeStr != "")
        std::cout << "   * Setting up the tree parameter from information in json file " << std::endl;
    else
        std::cout << "   * Setting up the tree parameter with randomly chosen tree" << std::endl;
        
    if (treeStr != "")
        fullTree.trees[0] = new Tree(treeStr, tNames, lambda, rv, isClockConstrained);
    else
        fullTree.trees[0] = new Tree(tNames, lambda, rv, isClockConstrained);
    fullTree.trees[1] = new Tree(*fullTree.trees[0]);
    
    // initialize the subtrees
    initializeSubtrees(wordAlignments);
    
    if (checkSubtreeCompatibility(0.0001) == false)
        Msg::error("Subtrees are not compatible with the canonical tree");
        
    parmStrLen = 0;
}

ParameterTree::~ParameterTree(void) {

    delete fullTree.trees[0];
    delete fullTree.trees[1];
    clearSubtrees();
}

void ParameterTree::accept(void) {

    *(fullTree.trees[1]) = *(fullTree.trees[0]);
    for (std::map<RbBitSet,TreePair>::iterator it = subTrees.begin(); it != subTrees.end(); it++)
        *(it->second.trees[1]) = *(it->second.trees[0]);
}

bool ParameterTree::checkSubtreeCompatibility(double tolerance) {

    std::map<TaxonPair,double,CompTaxonPair> fullTreeDistances = fullTree.trees[0]->pairwiseDistances();
    int n = (int)fullTree.trees[0]->getTaxonNames().size();
    if (fullTreeDistances.size() != n * (n-1) / 2)
        Msg::error("Incomplete pairwise distance map for full tree");
    
    for (std::map<RbBitSet,TreePair>::iterator i = subTrees.begin(); i != subTrees.end(); i++)
        {
        std::map<TaxonPair,double,CompTaxonPair> subTreeDistances = i->second.trees[0]->pairwiseDistances();
        n = (int)i->second.trees[0]->getTaxonNames().size();
        if (subTreeDistances.size() != n * (n-1) / 2)
            Msg::error("Incomplete pairwise distance map for sub tree");
        
        for (std::map<TaxonPair,double,CompTaxonPair>::iterator j = subTreeDistances.begin(); j != subTreeDistances.end(); j++)
            {
            double stD = j->second;
            
            std::map<TaxonPair,double,CompTaxonPair>::iterator it = fullTreeDistances.find(j->first);
            if (it != fullTreeDistances.end())
                {
                double ftD = it->second;
                double diff = fabs(stD - ftD);
                if (diff > tolerance)
                    return false;
                }
            else
                {
                Msg::error("Could not find taxon pair");
                }
            }
        }
    return true;
}

void ParameterTree::clearSubtrees(void) {

    for (std::map<RbBitSet,TreePair>::iterator it = subTrees.begin(); it != subTrees.end(); it++)
        {
        delete it->second.trees[0];
        delete it->second.trees[1];
        }
    subTrees.clear();
}

int ParameterTree::countMaskBits(std::vector<bool>& m) {

    int n = 0;
    for (size_t i=0; i<m.size(); i++)
        {
        if (m[i] == true)
            n++;
        }
    return n;
}

void ParameterTree::fillParameterValues(double* /*x*/, int& /*start*/, int /*maxNumValues*/) {

}

void ParameterTree::nniArea(std::vector<Node*>& backbone, Node*& incidentNode) {

    Tree* t = fullTree.trees[0];
    Node* r = t->getRoot();
    Node* p = NULL;
    bool goodNode = true;
    do
        {
        p = t->randomNode(rv);
        goodNode = true;
        if (p->getIsLeaf() == true || p == r)
            goodNode = false;
        else if (p->getAncestor() == r)
            goodNode = false;
        } while(goodNode == false);
    Node* pAnc = p->getAncestor();
    
    if (pAnc == NULL)
        Msg::error("Problem choosing random backbone");
    std::vector<Node*> pDes = p->getDescendantsVector();
    if (pDes.size() != 2)
        Msg::error("Problem choosing random backbone");
    std::vector<Node*> pAncDes = pAnc->getDescendantsVector();
    if (pAncDes.size() != 2)
        Msg::error("Problem choosing random backbone");

    std::vector<Node*> possibleIncidentNodes(2);
    Node* a = pDes[rv->uniformRvInt(2)];
    if (a == pDes[0])
        possibleIncidentNodes[0] = pDes[1];
    else
        possibleIncidentNodes[0] = pDes[0];
    Node* b = NULL;
    if (rv->uniformRv() < 0.5)
        {
        b = pAnc;
        if (pAncDes[0] == p)
            possibleIncidentNodes[1] = pAncDes[1];
        else
            possibleIncidentNodes[1] = pAncDes[0];
        }
    else
        {
        if (pAncDes[0] == p)
            b = pAncDes[1];
        else
            b = pAncDes[0];
        possibleIncidentNodes[1] = pAnc;
        }
        
    // we always assume this order, from tip to root
    backbone.clear();
    backbone.push_back(a);
    backbone.push_back(p);
    backbone.push_back(b);
    incidentNode = possibleIncidentNodes[rv->uniformRvInt(2)];
}

Tree* ParameterTree::getActiveTree(RbBitSet& mask) {

    if (mask.getNumberSetBits() == mask.size())
        return fullTree.trees[0];
        
    std::map<RbBitSet,TreePair>::iterator it = subTrees.find(mask);
    if (it != subTrees.end())
        return it->second.trees[0];
    return NULL;
}

Tree* ParameterTree::getActiveTree(const RbBitSet& mask) {

    if (mask.getNumberSetBits() == mask.size())
        return fullTree.trees[0];
        
    std::map<RbBitSet,TreePair>::iterator it = subTrees.find(mask);
    if (it != subTrees.end())
        return it->second.trees[0];
    return NULL;
}

std::string ParameterTree::getJsonString(void) {

    std::string str = "\"Tree\": \"" + fullTree.trees[0]->getNewick(20) + "\"";
    return str;
}

std::string ParameterTree::getString(void) {

    std::string str = fullTree.trees[0]->getNewick(5);
    return str;
}

char* ParameterTree::getCString(void) {

    return parmStr;
}

std::string ParameterTree::getUpdateName(int idx) {

    if (idx == 0)
        return "local";
    else if (idx == 1)
        return "branch length";
    else if (idx == 2)
        return "tree length";
    else if (idx == 3)
        return "random tree";
    return "";
}

void ParameterTree::initializeSubtrees(std::vector<Alignment*>& alns) {
        
    // remove all of the current subtrees
    clearSubtrees();
    
    // loop over all of the alignments, finding those that require subtrees because they only have a subset of the
    // full canonical list of taxa
    for (int i=0; i<alns.size(); i++)
        {
        // get the taxon mask and only proceed if the taxon list is incomplete
        std::vector<bool>& m = alns[i]->taxonMask;
        if (countMaskBits(m) == m.size())
            continue;
        RbBitSet mask(m);

        std::map<RbBitSet,TreePair>::iterator it = subTrees.find(mask);
        if (it == subTrees.end())
            {
            // build the subtree
            Tree* t0 = new Tree(*fullTree.trees[0], mask);
            Tree* t1 = new Tree(*t0);
            TreePair pair(t0, t1);
            subTrees.insert( std::make_pair(mask,pair) );
            }
        }
        
#   if 0
    int i = 0;
    for (std::map<RbBitSet,TreePair>::iterator it = subTrees.begin(); it != subTrees.end(); it++)
        {
        std::cout << "Subtree: " << i << " (" << it->first << ")" << std::endl;
        it->second.trees[0]->print();
        i++;
        }
#   endif
}

double ParameterTree::lnPriorProbability(void) {

#   if 1

    Tree* t = fullTree.trees[0];

    // get a vector with the branch lengths
    std::vector<Node*>& dpSeq = fullTree.trees[0]->getDownPassSequence();
    double lnP = 0.0;
    if (isClockConstrained == false)
        {
        // unconstrained prior probability
        for (Node* p : dpSeq)
            {
            if (p->getAncestor() != t->getRoot())
                {
                if (p == t->getRoot())
                    {
                    std::set<Node*,CompNode>& rootDes = p->getDescendants()->getNodes();
                    double len = 0.0;
                    for (Node* n : rootDes)
                        len += n->getBranchLength();
                    lnP += Probability::Exponential::lnPdf(lambda, len);
                    }
                else
                    {
                    lnP += Probability::Exponential::lnPdf(lambda, p->getBranchLength());
                    }
                }
            }
        }
    else
        {
        // clock-constrained prior probability
        for (Node* p : dpSeq)
            {
            if (p->getIsLeaf() == false)
                {
                double d = p->getNodeTime() - p->oldestDescendantTime();
                lnP += Probability::Exponential::lnPdf(lambda, d);
                }
            }
        }
    return lnP;

#   else

    // Prior from:
    // Rannala, B., T. Zhu, and Z. Yang. 2012. Tail paradox, partial identifiability, and
    //    influential priors in Bayesian branch length inference. Molecular Biology and
    //    Evolution 29(1):325-335.
    // with C = 1.0 and alpha = 1.0.

    // joint prior on tree length and branch lengths from
    Tree* t = fullTree.trees[0];
    double alphaT = 1.0;
    double s = (double)t->getNumTaxa();

    // get a vector with the branch lengths
    double treeLength = t->getTreeLength();
    std::vector<Node*>& dpSeq = fullTree.trees[0]->getDownPassSequence();
    std::vector<double> branchLengths;
    for (int i=0; i<dpSeq.size(); i++)
        {
        Node* p = dpSeq[i];
        if (p->getAncestor() != t->getRoot())
            {
            if (p == t->getRoot())
                {
                std::set<Node*,CompNode>& rootDes = p->getDescendants()->getNodes();
                double len = 0.0;
                for (Node* n : rootDes)
                    len += n->getBranchProportion() * treeLength;
                branchLengths.push_back(len);
                }
            else
                {
                branchLengths.push_back( treeLength * p->getBranchProportion() );
                }
            }
        }
        
    double lnP = 0.0;
    lnP += alphaT * log(betaT) - Probability::Helper::lnGamma(alphaT) - betaT * treeLength;
    lnP += (alphaT - 1.0) * log(treeLength);
    lnP += Probability::Helper::lnGamma(2 * s - 4.0 + 1.0);
    lnP += (-2.0 * s + 4.0) * log(treeLength);
            
    return lnP;
#   endif
}

void ParameterTree::print(void) {

    fullTree.trees[0]->print();
}

void ParameterTree::printNewick(void) {

    std::cout << fullTree.trees[0]->getNewick(6) << ";" << std::endl;
    for (std::map<RbBitSet,TreePair>::iterator it = subTrees.begin(); it != subTrees.end(); it++)
        std::cout << it->second.trees[0]->getNewick(6) << ";" << std::endl;
}

void ParameterTree::reject(void) {

    *(fullTree.trees[0]) = *(fullTree.trees[1]);
    for (std::map<RbBitSet,TreePair>::iterator it = subTrees.begin(); it != subTrees.end(); it++)
        *(it->second.trees[0]) = *(it->second.trees[1]);
    modelPtr->flipActiveLikelihood();
}

double ParameterTree::update(int iter) {
    
    UserSettings& settings = UserSettings::userSettings();
    
    // pick a tree parameter to update
    double probNni   = 0.00;
    double probBrlen = 0.75;
    if (iter < settings.getNumMcmcCycles() * 0.25 || settings.getCalculateMarginalLikelihood() == true || isClockConstrained == true)
        {
        probNni = 0.0;
        probBrlen = 0.75;
        }
    double u = rv->uniformRv();
    if (u <= probNni)
        return updateNni();
    else if (u > probNni && u <= probNni + probBrlen)
        return updateBrlen();
    else
        return updateTreeLength();
}

double ParameterTree::updateBrlen(void) {

    lastUpdateType.first = this;
    lastUpdateType.second = 1;
    
    if (isClockConstrained == false)
        {
        // update an unconstrained tree
        
        // tuning parameter for move
        double tuning = log(4.0);
        
        // pick a branch at random (note that the two branches
        // from the root are treated as a single branch)
        Tree* t = getActiveTree();
        std::vector<Node*>& nodes = t->getDownPassSequence();
        Node* nde = NULL;
        do
            {
            nde = nodes[rv->uniformRvInt(nodes.size())];
            if (nde == t->getRoot())
                nde = NULL;
            } while(nde == NULL);
            
        // determine if the branch is one of the
        // two branches incident to the root
        bool isRootBranch = false;
        if (nde->getAncestor() == t->getRoot())
            {
            isRootBranch = true;
            std::set<Node*,CompNode>& des = t->getRoot()->getDescendants()->getNodes();
            if (des.size() != 2)
                Msg::error("Expecting two descendants at root of the tree");
            }
            
        // get the old branch length
        double oldBrlen = 0.0;
        if (isRootBranch == false)
            oldBrlen = nde->getBranchLength();
        else
            {
            std::set<Node*,CompNode>& des = t->getRoot()->getDescendants()->getNodes();
            for (Node* n : des)
                oldBrlen += n->getBranchLength();
            }

        // propose new branch length
        double randomFactor = exp( -tuning*(rv->uniformRv()-0.5) );
        double newBrlen = oldBrlen * randomFactor;
        
        // set the new branch length
        if (isRootBranch == false)
            nde->setBranchLength(newBrlen);
        else
            {
            std::set<Node*,CompNode>& des = t->getRoot()->getDescendants()->getNodes();
            for (Node* n : des)
                n->setBranchLength(newBrlen*0.5);
            }

        // update the subtrees (before updating transition probabilities)
        updateSubtrees();
        
        // update the transition probabilities
        updateChangesTransitionProbabilities = true;
        TransitionProbabilities* tip = modelPtr->getTransitionProbabilities();
        tip->flipActive();
        tip->setNeedsUpdate(true);
        tip->setTransitionProbabilities();

        modelPtr->setUpdateLikelihood();
        modelPtr->flipActiveLikelihood();

        return log(randomFactor);
        }
    else
        {
        // update a clock-constrained tree

        // pick an interior node at random
        Tree* t = getActiveTree();
        //t->print("before");
        std::vector<Node*>& nodes = t->getDownPassSequence();
        Node* nde = NULL;
        do
            {
            nde = nodes[rv->uniformRvInt(nodes.size())];
            } while(nde->getIsLeaf() == true);
            
        // find the time of the oldest descendant of nde
        double oldestDescendantTime = nde->oldestDescendantTime();
        
        // and the time of the ancestor
        double ancestorTime = 0.0;
        if (nde->getAncestor() != NULL)
            ancestorTime = nde->getAncestor()->getNodeTime();
            
        // pick a new node time
        double lnProposalProbability = 0.0;
        if (nde->getAncestor() != NULL)
            {
            double u = rv->uniformRv();
            /*std::cout << "node                 = " << nde->getIndex() << std::endl;
            std::cout << "oldestDescendantTime = " << oldestDescendantTime << std::endl;
            std::cout << "ancestorTime         = " << ancestorTime << std::endl;
            std::cout << "new time             = " << oldestDescendantTime + u * (ancestorTime - oldestDescendantTime) << std::endl;*/
            nde->setNodeTime( oldestDescendantTime + u * (ancestorTime - oldestDescendantTime) );
            }
        else
            {
            // we treat the root differently
            double tuning = log(4.0);
            double oldD = nde->getNodeTime() - oldestDescendantTime;
            double randomFactor = exp( -tuning*(rv->uniformRv()-0.5) );
            double newD = oldD * randomFactor;
            //std::cout << oldD << " -> " << oldestDescendantTime + newD << std::endl;
            nde->setNodeTime( oldestDescendantTime + newD );
            
            lnProposalProbability = log(randomFactor);
            }

        // adjust the branch lengths
        for (Node* n : t->getDownPassSequence())
            {
            if (n != t->getRoot())
                n->setBranchLength( n->getAncestor()->getNodeTime() - n->getNodeTime() );
            else
                n->setBranchLength(0.0);
            }
        //t->print("after");


        // update the subtrees (before updating transition probabilities)
        updateSubtrees();
        
        // update the transition probabilities
        updateChangesTransitionProbabilities = true;
        TransitionProbabilities* tip = modelPtr->getTransitionProbabilities();
        tip->flipActive();
        tip->setNeedsUpdate(true);
        tip->setTransitionProbabilities();

        modelPtr->setUpdateLikelihood();
        modelPtr->flipActiveLikelihood();

        return lnProposalProbability;
        }
    
}

double ParameterTree::updateNni(void) {

    if (isClockConstrained == true)
        Msg::error("Cannot update a clock constrained tree using NNI");
        
    lastUpdateType.first = this;
    lastUpdateType.second = 0;

    double tuning = log(2.0);
    Tree* t = fullTree.trees[0];
#   if defined(DEBUG_LOCAL)
    t->print("BEFORE");
#   endif
        
    // randomly choose area of rearrangement
    std::vector<Node*> backbone;
    Node* incidentNode = NULL;
    nniArea(backbone, incidentNode);
    Node* p = backbone[1];
    Node* pAnc = p->getAncestor();
    Node* a = backbone[0];
    Node* b = backbone[2];
    if (b == pAnc && b->getAncestor() == t->getRoot())
        {
        std::vector<Node*> rDes = t->getRoot()->getDescendantsVector();
        double rBrlen = rDes[0]->getBranchLength() + rDes[1]->getBranchLength();
        rDes[0]->setBranchLength(0.0);
        rDes[1]->setBranchLength(0.0);
        b->setBranchLength(rBrlen);
        }
        
    // get length of backbone and randomly change it
    double mOld = a->getBranchLength() + p->getBranchLength() + b->getBranchLength();
    double mNew = mOld * exp(tuning*(rv->uniformRv()-0.5));
    double factor = mNew / mOld;
    a->setBranchLength( a->getBranchLength() * factor );
    p->setBranchLength( p->getBranchLength() * factor );
    b->setBranchLength( b->getBranchLength() * factor );

    // rearrange the tree
    if (incidentNode == pAnc)
        {
        // the backbone must be rotated around incident
        double newPt = rv->uniformRv() * mNew;
        if (newPt < a->getBranchLength())
            {
            p->removeDescendant(a);
            p->addDescendant(b);
            pAnc->removeDescendant(b);
            pAnc->addDescendant(a);
            a->setAncestor(pAnc);
            b->setAncestor(p);
            b->setBranchLength( b->getBranchLength() + p->getBranchLength() );
            p->setBranchLength(newPt);
            a->setBranchLength( a->getBranchLength() - newPt );
            }
        else
            {
            // no topology change, but adjust branch lengths for pAnc descendants
            newPt -= a->getBranchLength();
            double sum = b->getBranchLength() + p->getBranchLength();
            b->setBranchLength(newPt);
            p->setBranchLength(sum-newPt);
            }
        }
    else
        {
        // detach and reattach incident as it's above backbone
        Node* incidentAnc = incidentNode->getAncestor();
        Node* incidentSis = incidentNode->getSisterNode();
        Node* incidentAncAnc = incidentAnc->getAncestor();
        incidentAncAnc->removeDescendant(incidentAnc);
        incidentAncAnc->addDescendant(incidentSis);
        incidentSis->setAncestor(incidentAncAnc);
        incidentSis->setBranchLength( incidentSis->getBranchLength() + incidentAnc->getBranchLength() );
        incidentAnc->removeDescendant(incidentSis);
        
        double newPt = rv->uniformRv() * mNew;
        if (newPt < a->getBranchLength())
            {
            // reattach along branch a
            Node* aAnc = a->getAncestor();
            aAnc->removeDescendant(a);
            aAnc->addDescendant(incidentAnc);
            incidentAnc->addDescendant(a);
            incidentAnc->setAncestor(aAnc);
            a->setAncestor(incidentAnc);
            a->setBranchLength( a->getBranchLength() - newPt );
            incidentAnc->setBranchLength( newPt );
            }
        else
            {
            Node* x = NULL;
            if (a->getAncestor() == p)
                x = p;
            else
                {
                if (pAnc == b)
                    x = pAnc;
                else
                    x = b;
                }
            Node* xAnc = x->getAncestor();
            xAnc->removeDescendant(x);
            xAnc->addDescendant(incidentAnc);
            incidentAnc->addDescendant(x);
            incidentAnc->setAncestor(xAnc);
            x->setAncestor(incidentAnc);
            double xPt = newPt - a->getBranchLength();
            x->setBranchLength( x->getBranchLength() - xPt );
            incidentAnc->setBranchLength( xPt );
            }
        }
    
    // reinitialize the down pass sequence for the tree, in case the topology has changed
    t->initializeDownPassSequence();
        
    // break the branch at the root evenly in two
    std::vector<Node*> rDes = t->getRoot()->getDescendantsVector();
    double rBrlenSum = rDes[0]->getBranchLength() + rDes[1]->getBranchLength();
    rDes[0]->setBranchLength(rBrlenSum*0.5);
    rDes[1]->setBranchLength(rBrlenSum*0.5);
    
    // update the subtrees (before updating transition probabilities)
    updateSubtrees();

    // update the transition probabilities
    updateChangesTransitionProbabilities = true;
    TransitionProbabilities* tip = modelPtr->getTransitionProbabilities();
    tip->flipActive();
    tip->setNeedsUpdate(true);
    tip->setTransitionProbabilities();

    modelPtr->setUpdateLikelihood();
    modelPtr->flipActiveLikelihood();

#   if defined(DEBUG_LOCAL)
    t->print("AFTER");
    std::cout << "backbone[0]  = " << backbone[0]->getIndex() << std::endl;
    std::cout << "backbone[1]  = " << backbone[1]->getIndex() << std::endl;
    std::cout << "backbone[2]  = " << backbone[2]->getIndex() << std::endl;
    std::cout << "incidentNode = " << incidentNode->getIndex() << std::endl;
#   endif

#   if 0
    if (checkSubtreeCompatibility(0.0001) == false)
        Msg::error("Subtrees are not compatible with the canonical tree");
#   endif

    return 3.0 * log(factor);
}

double ParameterTree::updateSpr(void) {

    if (isClockConstrained == true)
        Msg::error("Cannot update a clock constrained tree using SPR");

    return 0.0;
}

double ParameterTree::updateTreeLength(void) {

    lastUpdateType.first = this;
    lastUpdateType.second = 2;

    // update the tree length
    Tree* t = getActiveTree();
    double tuning = log(2.0);
    double randomFactor = exp(tuning * (rv->uniformRv()-0.5));
    
    std::vector<Node*>& dpSeq = t->getDownPassSequence();
    if (isClockConstrained == false)
        {
        // update all of the branch lengths (unconstrained)
        for (Node* p : dpSeq)
            {
            if (p->getAncestor() != NULL)
                p->setBranchLength( p->getBranchLength() * randomFactor );
            }
        }
    else
        {
        // update all node times, then branch lengths (clock constrained)
        for (Node* p : dpSeq)
            {
            if (p->getIsLeaf() == false)
                {
                double oldT = p->getNodeTime();
                double newT = oldT * randomFactor;
                p->setNodeTime(newT);
                }
            }
        for (Node* p : dpSeq)
            {
            if (p != t->getRoot())
                p->setBranchLength( p->getAncestor()->getNodeTime() - p->getNodeTime() );
            else
                p->setBranchLength(0.0);
            }
        }

    // update the transition probabilities
    updateChangesTransitionProbabilities = true;
    TransitionProbabilities* tip = modelPtr->getTransitionProbabilities();
    tip->flipActive();
    tip->setNeedsUpdate(true);
    tip->setTransitionProbabilities();

    modelPtr->setUpdateLikelihood();
    modelPtr->flipActiveLikelihood();
    
    // calculate and return the log of the proposal probability
    double lnP = 0.0;
    if (isClockConstrained == false)
        lnP = (double)(dpSeq.size()-1) * log(randomFactor);
    else
        lnP = (t->getNumTaxa()-1) * log(randomFactor);
    return lnP;
}

double ParameterTree::updateTopologyFromPrior(void) {

    if (isClockConstrained == true)
        Msg::error("Cannot update a clock constrained topology from prior");

    lastUpdateType.first = this;
    lastUpdateType.second = 3;
    
    Tree* t = getActiveTree();
    std::vector<std::string> tNames;
    std::vector<std::string>& tNamesToCopy = t->getTaxonNames();
    for (int i=0; i<tNamesToCopy.size(); i++)
        tNames.push_back(tNamesToCopy[i]);
        
    // update the topology and calculate proposal probability
    double lnP1 = lnPriorProbability();
    t->buildRandomTree(tNames, lambda, rv, isClockConstrained);
    double lnP2 = lnPriorProbability();

    // update the transition probabilities
    updateChangesTransitionProbabilities = true;
    TransitionProbabilities* tip = modelPtr->getTransitionProbabilities();
    tip->flipActive();
    tip->setNeedsUpdate(true);
    tip->setTransitionProbabilities();
    modelPtr->setUpdateLikelihood();
    modelPtr->flipActiveLikelihood();

    return lnP1 - lnP2;
}

double ParameterTree::updateBranchlengthsFromPrior(void) {

    if (isClockConstrained == true)
        Msg::error("Cannot update a clock constrained tree's branch lengths from prior");

    //lastUpdateType = "random branch lengths";
    lastUpdateType.first = this;
    lastUpdateType.second = 4;

    Tree* t = getActiveTree();
    double lnP1 = lnPriorProbability();
    std::vector<Node*> dpSeq = t->getDownPassSequence();
    for (int i=0; i<dpSeq.size(); i++)
        {
        Node* p = dpSeq[i];
        if (p != t->getRoot())
            p->setBranchLength( Probability::Exponential::rv(rv, lambda) );
        else
            p->setBranchLength(0.0);
        }
    double rootLen = Probability::Exponential::rv(rv, lambda);
    std::set<Node*,CompNode>& rootDes = t->getRoot()->getDescendants()->getNodes();
    for (Node* n : rootDes)
        n->setBranchLength(rootLen * 0.5);

    double lnP2 = lnPriorProbability();
    return lnP1 - lnP2;
}

void ParameterTree::updateSubtrees(void) {

    Tree* t = fullTree.trees[0];
    for (std::map<RbBitSet,TreePair>::iterator it = subTrees.begin(); it != subTrees.end(); it++)
        {
        Tree* st = it->second.trees[0];
        RbBitSet bs = RbBitSet(it->first);
        st->makeSubtree(*t, bs);
        }
}
