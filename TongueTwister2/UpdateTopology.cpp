#include "Model.hpp"
#include "Node.hpp"
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
}

void UpdateTopology::setDependants(void) {

    clearDependencyFlags();
    
    updatedParameter      = myParm;
    rateMatrixNeedsUpdate = false;
    allTiprobsNeedUpdate  = false;
    singleBranchChanged   = true;
    // changedBranchLength is set in update()
}

double UpdateTopology::update(void) {

    Tree* t = myParm->getTree();
    double tuning = log(4.0);
    
    // pick a branch at random (not the root)
    const std::vector<Node*>& dp = t->getPostOrder();
    Node* p = nullptr;
    do {
        p = dp[static_cast<int>(rng->uniformRv() * dp.size())];
    } while (p == t->getRoot() || p->getIsOutgroup() == true);
            
    // propose a valid branch length using rejection
    double oldV = p->getBranchLength();
    double randomFactor = 1.0;
    do {
        randomFactor = exp(tuning * (rng->uniformRv() - 0.5));
    } while (oldV * randomFactor > maxBrlen);
    double newV = oldV * randomFactor;
    
    // update the length of the branch (this will also update
    // the length of the branches of all of the subtrees of this tree)
    myParm->setBranchLength(p, newV);
    
    setDependants();
    changedBranchLength = newV;  // store for UpdateManager
    
    // calculate the log of the Hastings ratio
    double a = std::exp(-tuning * 0.5);
    double b = std::exp(tuning * 0.5);
    double boundaryFactor = log(std::min(oldV*b,maxBrlen)/(oldV*a)) / log(std::min(newV*b,maxBrlen)/(newV*a));
    
    return log(randomFactor) + log(boundaryFactor);
}

double UpdateTopology::update(double power) {

    if (rng->uniformRv() < priorSampleProb(power))
        return updateFromPrior();
    return update();
}

double UpdateTopology::updateFromPrior(void) {

    Tree* t = myParm->getTree();
    
    // pick a branch at random (not the root)
    const std::vector<Node*>& dp = t->getPostOrder();
    Node* p = nullptr;
    do {
        p = dp[static_cast<int>(rng->uniformRv() * dp.size())];
    } while (p == t->getRoot() || p->getIsOutgroup() == true);
            
    // pick a branch length from the truncated exponential prior
    double lambda = myParm->getBrlenLambda();
    double oldV = p->getBranchLength();
    double newV = 0;
    do {
        newV = Probability::Exponential::rv(rng, lambda);
    } while(newV > maxBrlen);
    
    // update the length of the branch (this will also update
    // the length of the branches of all of the subtrees of this tree)
    myParm->setBranchLength(p, newV);
    
    setDependants();
    changedBranchLength = newV;  // store for UpdateManager
    
    return exp(-lambda * (oldV - newV));
}
