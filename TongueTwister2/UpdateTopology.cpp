#include <cmath>
#include <limits>
#include <set>
#include "Model.hpp"
#include "Msg.hpp"
#include "Node.hpp"
#include "ParameterAlignment.hpp"
#include "ParameterTree.hpp"
#include "Probability.hpp"
#include "RandomVariable.hpp"
#include "TransitionProbabilities.hpp"
#include "Tree.hpp"
#include "UpdateAlignment.hpp"
#include "UpdateTopology.hpp"



UpdateTopology::UpdateTopology(Model* m, RandomVariable* r, ParameterTree* p, const std::vector<UpdateAlignment*>& alnVec) : 
    Update(m, r), myParm(p), myAlignmentUpdates(alnVec) {

    tiProbs = model->getTiProbs();
    maxBrlen = myParm->getMaximumBrlen();
    tuning = 2.0 * log(2.0);  // λ₂ parameter from paper
}

void UpdateTopology::applyNni(Node* u, Node* v, Node* a, Node* c) {

    // NNI swap: move a from u to v, move c from v to u
    // Before: a,b are children of u; u,c are children of v
    // After:  c,b are children of u; u,a are children of v
    
    u->removeDescendant(a);
    v->removeDescendant(c);
    u->addDescendant(c);
    v->addDescendant(a);
    a->setAncestor(v);
    c->setAncestor(u);
}

std::vector<Parameter*> UpdateTopology::getAdditionalModifiedParameters(void) {

    // Return the alignment parameters that are modified during topology updates.
    // These need to be kept/restored along with the tree parameter.
    std::vector<Parameter*> params;
    params.reserve(myAlignmentUpdates.size());
    for (UpdateAlignment* aln : myAlignmentUpdates)
        params.push_back(aln->getAlignmentParameter());
    return params;
}

Node* UpdateTopology::randomlyChooseInternalBranch(Tree* t) {

    const std::vector<Node*>& dpseq = t->getPostOrder();
    size_t n = dpseq.size();
    Node* rootNode = t->getRoot();
    Node* p = nullptr;
    while(p == nullptr)
        {
        p = dpseq[static_cast<size_t>(rng->uniformRv() * n)];
        if (p->getIsLeaf() == true || p == rootNode)
            {
            p = nullptr;
            continue;
            }
        if (p == rootNode->getDescendant(0) || p == rootNode->getDescendant(1))
            {
            p = nullptr;
            continue;
            }
        }
    return p;
}

void UpdateTopology::setDependants(void) {

    clearDependencyFlags();
    
    updatedParameter      = myParm;
    rateMatrixNeedsUpdate = false;
    allTiprobsNeedUpdate  = true;   // topology/branch changes require full recomputation
    singleBranchChanged   = false;
}

double UpdateTopology::update(void) {

    // LOCAL algorithm without molecular clock
    // Larget, B. and Simon, D. L. 1999. Markov chain Monte Carlo algorithms for 
    //   the Bayesian analysis of phylogenetic trees. Mol. Biol. Evol. 16:750-759.
    
    // get the tree
    Tree* t = myParm->getTree();
    if (t->getNumTaxa() <= 3)
        Msg::error("Can only update topology on trees with more than three taxa");
            
    // Find valid internal edges
    // An internal branch connects two internal nodes (neither is a leaf)
    // We exclude: root, descendants of root (to keep root structure unchanged)
    Node* u = randomlyChooseInternalBranch(t);
    Node* v = u->getAncestor();
    
    // label u's descendants as a and b (randomly choose which is a)
    Node* a = u->getDescendant(0);
    Node* b = u->getDescendant(1);
    if (rng->uniformRv() < 0.5)
        std::swap(a, b);
    
    // find c: v's other descendant (not u)
    Node* c = nullptr;
    for (size_t i=0; i<v->numDescendants(); i++)
        {
        Node* des = v->getDescendant(i);
        if (des != u)
            {
            c = des;
            break;
            }
        }
    if (c == nullptr)
        Msg::error("Problem in LOCAL update: could not find c");
            
    // Current path distances (path: a -> u -> v -> c)
    // x = dist(a, u) = a's branch length
    // y = dist(a, v) = a's branch length + u's branch length
    // m = dist(a, c) = a's branch length + u's branch length + c's branch length
    double oldALen = a->getBranchLength();
    double oldULen = u->getBranchLength();
    double oldCLen = c->getBranchLength();
    
    double x = oldALen;
    double y = x + oldULen;
    double m = y + oldCLen;
    
    // propose new total path length: m* = m × exp(λ(U₁ - 0.5))
    double mStar = m * exp(tuning * (rng->uniformRv() - 0.5));
    
    // Propose new node positions
    // Choose to move u or v with equal probability
    double xStar, yStar;
    bool moveU = (rng->uniformRv() < 0.5);
    double U2 = rng->uniformRv();
    
    if (moveU)
        {
        // u moves to random position, v's relative position preserved
        xStar = U2 * mStar;
        yStar = y * (mStar / m);
        }
    else
        {
        // v moves to random position, u's relative position preserved
        xStar = x * (mStar / m);
        yStar = U2 * mStar;
        }
    
    // determine if topology changes
    bool topologyChanges = (xStar > yStar);
    
    // calculate new branch lengths
    double newALen, newULen, newCLen;
    
    if (topologyChanges == false)
        {
        // No topology change: x* < y*
        // Path remains a → u → v → c
        newALen = xStar;                  // dist(a, u*)
        newULen = yStar - xStar;          // dist(u*, v*)
        newCLen = mStar - yStar;          // dist(v*, c)
        }
    else
        {
        // Topology changes: x* > y* (u and v swap relative positions)
        // Path becomes a → v* → u* → c
        // After NNI: a attaches to v, c attaches to u
        newALen = yStar;                  // dist(a, v*) - a's new branch length
        newULen = xStar - yStar;          // dist(v*, u*) - u's new branch length
        newCLen = mStar - xStar;          // dist(u*, c) - c's new branch length
        }
    
    // validate proposed branch lengths
    if (newALen <= 0.0 || newULen <= 0.0 || newCLen <= 0.0 ||
        newALen > maxBrlen || newULen > maxBrlen || newCLen > maxBrlen)
        {
        setDependants();
        Msg::error("Invalid branch lengths in LOCAL");
        }
    
    // apply changes
    if (topologyChanges == true)
        {
        // Mark that topology is changing (for proper backup/restore)
        myParm->setTopologyChanged(true);
        
        // Apply NNI to full tree
        applyNni(u, v, a, c);
        
        // Apply NNI to all subtrees that contain both a and c
        myParm->applyNniToSubtrees();
        
        // Rebuild branch mappings since topology changed
        myParm->initializeBranchMappings();
        }
    
    // update branch lengths (propagates to subtrees via mappings)
    myParm->setBranchLength(a, newALen);
    myParm->setBranchLength(u, newULen);
    myParm->setBranchLength(c, newCLen);
    
    // set dependency flags
    setDependants();
    
    // Update transition probabilities NOW, before realigning.
    // The alignment proposal (realignFull) needs transition probabilities to
    // build the scoring matrix. Normally updateDependants() is called after
    // update() returns, but we need them immediately.
    tiProbs->updateAllBranches();
    
    // clear the flag to prevent redundant computation in updateDependants()
    allTiprobsNeedUpdate = false;
    
    // Update all alignments on the new tree topology.
    // Each realignFull() returns the log proposal ratio for that alignment.
    double lnAlignmentProposalRatio = 0.0;
    for (UpdateAlignment* aln : myAlignmentUpdates)
        lnAlignmentProposalRatio += aln->realignFull();
    
    // Hastings ratio = (m*/m)^3 * alignment_proposal_ratio
    double ratio = mStar / m;
    return 3.0 * log(ratio) + lnAlignmentProposalRatio;
}

double UpdateTopology::update(double power) {

    if (rng->uniformRv() < priorSampleProb(power))
        return updateFromPrior();
    return update();
}

double UpdateTopology::updateFromPrior(void) {

    return 0.0;
}
