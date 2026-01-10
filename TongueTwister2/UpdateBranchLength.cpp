#include "Model.hpp"
#include "Node.hpp"
#include "ParameterTree.hpp"
#include "Probability.hpp"
#include "RandomVariable.hpp"
#include "TransitionProbabilities.hpp"
#include "Tree.hpp"
#include "UpdateBranchLength.hpp"



UpdateBranchLength::UpdateBranchLength(Model* m, RandomVariable* r, ParameterTree* p) : Update(m, r), myParm(p) {

    updateInfo.resize(1);
    updateInfo[0].updateIdx = 0;
    updateInfo[0].updateName = getUpdateName();
    updateInfo[0].updateHash = hashUpdateName(getUpdateName());
    updateInfo[0].updateType = factor;
    updateInfo[0].tuningParameter = log(4.0);

    tiProbs = model->getTiProbs();
    maxBrlen = myParm->getMaximumBrlen();
}

void UpdateBranchLength::setDependants(void) {

    clearDependencyFlags();
    
    updatedParameter      = myParm;
    rateMatrixNeedsUpdate = false;
    allTiprobsNeedUpdate  = false;
    singleBranchChanged   = true;
    // changedBranchLength is set in update()
}

double UpdateBranchLength::update(void) {

    Tree* t = myParm->getTree();
    double tuning = updateInfo[0].tuningParameter;
    
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

double UpdateBranchLength::update(double power) {

    if (rng->uniformRv() < priorSampleProb(power))
        return updateFromPrior();
    return update();
}

double UpdateBranchLength::updateFromPrior(void) {

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
    double newV = 0.0;
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
