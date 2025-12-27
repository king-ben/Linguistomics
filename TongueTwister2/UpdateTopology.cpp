#include <cmath>
#include <limits>
#include <set>
#include "Model.hpp"
#include "Node.hpp"
#include "NodeComparator.hpp"
#include "ParameterTree.hpp"
#include "Probability.hpp"
#include "RandomVariable.hpp"
#include "TransitionProbabilityManager.hpp"
#include "Tree.hpp"
#include "UpdateTopology.hpp"



UpdateTopology::UpdateTopology(Model* m, RandomVariable* r, ParameterTree* p, const std::vector<ParameterAlignment*>& alnVec) : 
    Update(m, r), myParm(p), myAlignments(alnVec) {

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
    
    Tree* t = myParm->getTree();
    Node* root = t->getRoot();
    
    // Find valid internal edges
    // An internal edge connects two internal nodes (neither is a leaf)
    // We exclude: root, children of root (to keep root structure unchanged)
    const std::vector<Node*>& dp = t->getPostOrder();
    std::vector<Node*> validNodes;
    
    for (Node* p : dp)
        {
        // p will be 'u' - the lower node of the internal edge
        if (p->getIsLeaf() == true)
            continue;
        if (p == root)
            continue;
        if (p->numDescendants() != 2)
            continue;
            
        Node* anc = p->getAncestor();
        if (anc == nullptr)
            continue;
        if (anc == root)
            continue;  // don't modify edges directly below root
        if (anc->numDescendants() != 2)
            continue;
        
        // Also exclude if v's parent is the root (keeps root structure unchanged)
        Node* ancAnc = anc->getAncestor();
        if (ancAnc == nullptr)
            continue;
        if (ancAnc == root)
            continue;  // v is a child of root, skip to keep root structure intact
            
        validNodes.push_back(p);
        }
    
    if (validNodes.empty())
        {
        // No valid internal edges (tree too small)
        setDependants();
        return -std::numeric_limits<double>::infinity();
        }
    
    // Randomly select an internal edge
    size_t idx = static_cast<size_t>(rng->uniformRv() * validNodes.size());
    if (idx >= validNodes.size())
        idx = validNodes.size() - 1;
    Node* u = validNodes[idx];
    Node* v = u->getAncestor();
    
    // Label u's children as a and b (randomly choose which is a)
    Node* a = u->getDescendants()[0];
    Node* b = u->getDescendants()[1];
    if (rng->uniformRv() < 0.5)
        std::swap(a, b);
    
    // Find c: v's other child (not u)
    Node* c = nullptr;
    std::set<Node*,NodeComparator>& vChildren = v->getDescendants().getNodes();
    for (Node* child : vChildren)
        {
        if (child != u)
            {
            c = child;
            break;
            }
        }
    
    if (c == nullptr)
        {
        setDependants();
        return -std::numeric_limits<double>::infinity();
        }
    
    // Current path distances (path: a → u → v → c)
    // x = dist(a, u) = a's branch length
    // y = dist(a, v) = a's branch length + u's branch length
    // m = dist(a, c) = a's branch length + u's branch length + c's branch length
    double oldALen = a->getBranchLength();
    double oldULen = u->getBranchLength();
    double oldCLen = c->getBranchLength();
    
    double x = oldALen;
    double y = x + oldULen;
    double m = y + oldCLen;
    
    // Propose new total path length: m* = m × exp(λ(U₁ - 0.5))
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
    
    // Determine if topology changes
    bool topologyChanges = (xStar > yStar);
    
    // Calculate new branch lengths
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
    
    // Validate proposed branch lengths
    if (newALen <= 0.0 || newULen <= 0.0 || newCLen <= 0.0 ||
        newALen > maxBrlen || newULen > maxBrlen || newCLen > maxBrlen)
        {
        setDependants();
        return -std::numeric_limits<double>::infinity();
        }
    
    // Apply changes
    if (topologyChanges == true)
        {
        // Mark that topology is changing (for proper backup/restore)
        myParm->setTopologyChanged(true);
        
        // Apply NNI to full tree
        applyNni(u, v, a, c);
        
        // Apply NNI to all subtrees that contain both a and c
        myParm->applyNniToSubtrees(u, v, a, c);
        
        // Rebuild branch mappings since topology changed
        myParm->initializeBranchMappings();
        }
    
    // Update branch lengths (propagates to subtrees via mappings)
    myParm->setBranchLength(a, newALen);
    myParm->setBranchLength(u, newULen);
    myParm->setBranchLength(c, newCLen);
    
    setDependants();
    
    // Hastings ratio = (m*/m)²
    double ratio = mStar / m;
    return 2.0 * log(ratio);
}

double UpdateTopology::update(double power) {

    if (rng->uniformRv() < priorSampleProb(power))
        return updateFromPrior();
    return update();
}

double UpdateTopology::updateFromPrior(void) {

    // For prior sampling, we do a similar LOCAL move but sample 
    // the total path length from the prior (truncated exponential)
    
    Tree* t = myParm->getTree();
    Node* root = t->getRoot();
    double lambda = myParm->getBrlenLambda();
    
    // Find valid internal edges (same as in update())
    const std::vector<Node*>& dp = t->getPostOrder();
    std::vector<Node*> validNodes;
    
    for (Node* p : dp)
        {
        if (p->getIsLeaf() == true)
            continue;
        if (p == root)
            continue;
        if (p->numDescendants() != 2)
            continue;
            
        Node* anc = p->getAncestor();
        if (anc == nullptr)
            continue;
        if (anc == root)
            continue;
        if (anc->numDescendants() != 2)
            continue;
        
        Node* ancAnc = anc->getAncestor();
        if (ancAnc == nullptr)
            continue;
        if (ancAnc == root)
            continue;
            
        validNodes.push_back(p);
        }
    
    if (validNodes.empty())
        {
        setDependants();
        return -std::numeric_limits<double>::infinity();
        }
    
    // Randomly select an internal edge
    size_t idx = static_cast<size_t>(rng->uniformRv() * validNodes.size());
    if (idx >= validNodes.size())
        idx = validNodes.size() - 1;
    Node* u = validNodes[idx];
    Node* v = u->getAncestor();
    
    // Label nodes
    Node* a = u->getDescendants()[0];
    Node* b = u->getDescendants()[1];
    if (rng->uniformRv() < 0.5)
        std::swap(a, b);
    
    Node* c = nullptr;
    std::set<Node*,NodeComparator>& vChildren = v->getDescendants().getNodes();
    for (Node* child : vChildren)
        {
        if (child != u)
            {
            c = child;
            break;
            }
        }
    
    if (c == nullptr)
        {
        setDependants();
        return -std::numeric_limits<double>::infinity();
        }
    
    // Current values
    double oldALen = a->getBranchLength();
    double oldULen = u->getBranchLength();
    double oldCLen = c->getBranchLength();
    double m = oldALen + oldULen + oldCLen;
    
    // Sample three new branch lengths from truncated exponential prior
    double newALen, newULen, newCLen;
    do { newALen = Probability::Exponential::rv(rng, lambda); } while (newALen > maxBrlen);
    do { newULen = Probability::Exponential::rv(rng, lambda); } while (newULen > maxBrlen);
    do { newCLen = Probability::Exponential::rv(rng, lambda); } while (newCLen > maxBrlen);
    
    double mStar = newALen + newULen + newCLen;
    
    // Since we're sampling from prior, no topology change needed
    // (topology changes would require more complex prior ratio calculation)
    
    myParm->setBranchLength(a, newALen);
    myParm->setBranchLength(u, newULen);
    myParm->setBranchLength(c, newCLen);
    
    setDependants();
    
    // Prior ratio for three independent truncated exponentials
    double lnPriorRatio = -lambda * (newALen + newULen + newCLen - oldALen - oldULen - oldCLen);
    
    // Proposal ratio (sampling from prior means proposal ratio = 1/prior ratio for MH)
    return -lnPriorRatio;
}
