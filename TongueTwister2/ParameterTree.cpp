#include <iomanip>
#include <iostream>
#include <set>
#include "BitMask.hpp"
#include "JsonData.hpp"
#include "Msg.hpp"
#include "Node.hpp"
#include "NodeComparator.hpp"
#include "ParameterTree.hpp"
#include "Probability.hpp"
#include "Tree.hpp"



ParameterTree::ParameterTree(Model* m, RandomVariable* r, std::string n) : Parameter(m, r, n) {

    brlenLambda = 10.0;
    lnProbLessMax = log( 1.0 - exp(-brlenLambda * MAX_BRLEN) );
    topologyChanged = false;
}

ParameterTree::~ParameterTree(void) {

    delete fullTree.trees[0];
    delete fullTree.trees[1];
    for (auto [key, val] : subTrees)
        {
        delete val.trees[0];
        delete val.trees[1];
        }
}

void ParameterTree::applyNniToSubtrees(Node* u, Node* v, Node* a, Node* c) {

    // Apply NNI to all subtrees that contain the relevant nodes
    // The NNI swaps: a from u to v, c from v to u
    
    // Get the taxon indices involved (we need to find which subtrees contain both a and c)
    // For leaf nodes, we can use their index; for internal nodes, we need descendants
    
    for (auto& [mask, treePair] : subTrees)
        {
        Tree* subtree = treePair.trees[0];
        
        // Find corresponding nodes in subtree
        Node* subU = findCorrespondingNode(fullTree.trees[0], subtree, u);
        Node* subV = findCorrespondingNode(fullTree.trees[0], subtree, v);
        Node* subA = findCorrespondingNode(fullTree.trees[0], subtree, a);
        Node* subC = findCorrespondingNode(fullTree.trees[0], subtree, c);
        
        // If all four nodes exist in subtree and have the expected relationships,
        // apply the same NNI
        if (subU != nullptr && subV != nullptr && subA != nullptr && subC != nullptr)
            {
            // Verify relationships before applying NNI
            if (subA->getAncestor() == subU && 
                subC->getAncestor() == subV && 
                subU->getAncestor() == subV)
                {
                // Apply NNI: swap a and c
                subU->removeDescendant(subA);
                subV->removeDescendant(subC);
                subU->addDescendant(subC);
                subV->addDescendant(subA);
                subA->setAncestor(subV);
                subC->setAncestor(subU);
                }
            }
        }
}

bool ParameterTree::checkTipToTipDistances(double threshhold) {

    std::map<std::pair<std::string,std::string>,std::set<double>> allDistances;
    for (auto st : subTrees)
        {
        Tree* t = st.second.trees[0];
        dMap distances = t->pairwiseDistances();
        for (auto [key,val] : distances)
            {
            std::map<std::pair<std::string,std::string>,std::set<double>>::iterator it = allDistances.find(key);
            if (it == allDistances.end())
                {
                std::set<double> x;
                x.insert(val);
                allDistances.insert( std::make_pair(key,x) );
                }
            else
                {
                it->second.insert(val);
                }
            }
        }
    
    for (auto [key,val] : allDistances)
        {
        auto first = val.begin();
        auto last  = val.rbegin();
        double diff = fabs(*last - *first);
        if (fabs(diff) > threshhold)
            return false;
        }
    return true;
}

Node* ParameterTree::findCorrespondingNode(Tree* srcTree, Tree* dstTree, Node* srcNode) {

    // Find the node in dstTree that corresponds to srcNode in srcTree
    // For leaf nodes, match by name; for internal nodes, match by descendant leaf sets
    
    if (srcNode == nullptr)
        return nullptr;
    
    if (srcNode->getIsLeaf())
        {
        // Match leaf by name
        return dstTree->getLeafNamed(srcNode->getName());
        }
    else
        {
        // For internal nodes, we need to find the node with the same set of descendant leaves
        // This is complex; for now, use offset-based matching if trees have same structure
        const std::vector<Node*>& dstPostOrder = dstTree->getPostOrder();
        for (Node* dstNode : dstPostOrder)
            {
            if (dstNode->getIsLeaf() == false && dstNode->getOffset() == srcNode->getOffset())
                return dstNode;
            }
        
        // Fallback: try to match by descendant leaves
        // Get leaves under srcNode
        std::set<std::string> srcLeaves;
        std::vector<Node*> stack;
        stack.push_back(srcNode);
        while (!stack.empty())
            {
            Node* n = stack.back();
            stack.pop_back();
            if (n->getIsLeaf())
                srcLeaves.insert(n->getName());
            else
                {
                std::set<Node*,NodeComparator>& des = n->getDescendants().getNodes();
                for (Node* d : des)
                    stack.push_back(d);
                }
            }
        
        // Find internal node in dstTree with same leaves
        for (Node* dstNode : dstPostOrder)
            {
            if (dstNode->getIsLeaf())
                continue;
                
            std::set<std::string> dstLeaves;
            stack.clear();
            stack.push_back(dstNode);
            while (!stack.empty())
                {
                Node* n = stack.back();
                stack.pop_back();
                if (n->getIsLeaf())
                    dstLeaves.insert(n->getName());
                else
                    {
                    std::set<Node*,NodeComparator>& des = n->getDescendants().getNodes();
                    for (Node* d : des)
                        stack.push_back(d);
                    }
                }
            
            if (srcLeaves == dstLeaves)
                return dstNode;
            }
        }
    
    return nullptr;
}

size_t ParameterTree::getNumTaxa(void) {

    return fullTree.trees[0]->getNumTaxa();
}

Tree* ParameterTree::getTree(unsigned& mask) {

    std::unordered_map<unsigned,TreePair>::iterator it = subTrees.find(mask);
    if (it != subTrees.end())
        return it->second.trees[0];
    return nullptr;
}

Tree* ParameterTree::getTree(const unsigned& mask) {

    std::unordered_map<unsigned,TreePair>::iterator it = subTrees.find(mask);
    if (it != subTrees.end())
        return it->second.trees[0];
    return nullptr;
}

void ParameterTree::getTrees(std::vector<Tree*>& allTrees) {

    allTrees.resize(subTrees.size());
    int i = 0;
    for (auto key : subTrees)
        allTrees[i++] = key.second.trees[0];
}

void ParameterTree::initialize(const std::set<unsigned>& uniqueTaxonCombinations) {

    // check that the necessary keys are in the JSON object
    JsonData& jd = JsonData::jsonInstance();
    if (jd.hasKey("Tree") == false)
        Msg::error("Could not find \"Tree\" key in JSON object");
    if (jd.hasKey("Taxa") == false)
        Msg::error("Could not find \"Taxa\" key in JSON object");
        
    // get the prior for the branch length (default value is 10.0)
    if (jd.hasKey("BranchLengthLambda") == true)
        brlenLambda = jd.getValue<double>("BranchLengthLambda");
        
    // read the taxon names
    canonicalTaxonList = jd.getValue<std::vector<std::string>>("Taxa");

    // read the Newick string
    std::string newickString = jd.getValue<std::string>("Tree");
    
    // instantiate the full tree
    fullTree.trees[0] = new Tree(rng, newickString, canonicalTaxonList, brlenLambda, MAX_BRLEN);
    fullTree.trees[1] = new Tree(*fullTree.trees[0]);
    
    // instantiate the subtrees
    for (const unsigned& bs : uniqueTaxonCombinations)
        {
        Tree* st0 = new Tree(*fullTree.trees[0]);
        makeSubtree(*st0, bs);
        Tree* st1 = new Tree(*st0);
        subTrees.insert( std::make_pair(bs, TreePair(st0, st1)) );
        }
        
    if (checkTipToTipDistances(0.001) == false)
        Msg::error("Problem with subtree branch lengths");
        
    initializeBranchMappings();
}

void ParameterTree::initializeBranchMappings(void) {

    branchMappings.clear();
    
    // build mappings using trees[0] (the working copy)
    for (auto& [mask, treePair] : subTrees)
        {
        BranchMapping mapping;
        mapping.buildMapping(fullTree.trees[0], treePair.trees[0]);
        branchMappings[mask] = mapping;
        }
}

void ParameterTree::keep(void) {

    // Called when MCMC proposal is accepted
    // Copy state from trees[0] to trees[1] for all trees
    
    if (topologyChanged)
        {
        // Topology changed: do full clone to backup the new topology
        saveTopology();
        topologyChanged = false;
        }
    else
        {
        // Only branch lengths changed: just copy branch lengths
        const std::vector<Node*>& srcPostOrder = fullTree.trees[0]->getPostOrder();
        const std::vector<Node*>& dstPostOrder = fullTree.trees[1]->getPostOrder();
        for (size_t i=0; i<srcPostOrder.size(); i++)
            {
            Node* srcNode = srcPostOrder[i];
            Node* dstNode = dstPostOrder[i];
            dstNode->setBranchLength(srcNode->getBranchLength());
            }
        
        // copy subtree branch lengths
        for (auto& [mask, treePair] : subTrees)
            {
            std::unordered_map<unsigned, BranchMapping>::iterator it = branchMappings.find(mask);
            if (it != branchMappings.end())
                it->second.copyBranchLengths(treePair.trees[0], treePair.trees[1]);
            }
        }
}

double ParameterTree::lnPriorProbability(void) {
        
    Tree* t = fullTree.trees[0];
    Node* r = t->getRoot();
    
    // root branch probability
    const NodeSet& rootDescendants = r->getDescendants();
    double rootV = rootDescendants[0]->getBranchLength() + rootDescendants[1]->getBranchLength();
    double lnP = Probability::Exponential::lnPdf(brlenLambda, rootV) - lnProbLessMax;;
    
    // now deal with all of the non-root branches
    const std::vector<Node*>& dp = t->getPostOrder();
    for (Node* p : dp)
        {
        if (p != r)
            {
            if (p->getAncestor() != r)
                lnP += Probability::Exponential::lnPdf(brlenLambda, p->getBranchLength()) - lnProbLessMax;
            }
        }
        
    return lnP;
}

void ParameterTree::makeSubtree(Tree& t, const unsigned& taxonMask) {

    // 1. mark all nodes that are part of the subtree
    t.setAllFlags(false);
    std::map<int,int> tipMap;
    for (size_t i=0, k=0; i<canonicalTaxonList.size(); i++)
        {
        if (BitMask::isSet(taxonMask, i) == true)
            {
            Node* p = t.getLeafIndexed(i);
            if (p == nullptr)
                Msg::error("Could not find leaf to mark path to root");
            do
                {
                p->setFlag(true);
                p = p->getAncestor();
                } while (p != nullptr);
            tipMap.insert( std::make_pair(i,k) );
            k++;
            }
        }
        
    // 2. remove nodes that are not part of the subtree
    const std::vector<Node*>& dpseq1 = t.getPostOrder();
    for (Node* p : dpseq1)
        {
        if (p->getFlag() == false)
            {
            p->removeDescendants();
            Node* a = p->getAncestor();
            if (a != nullptr)
                a->removeDescendant(p);
            p->setAncestor(nullptr);
            }
        }
        
    // 3. remove superfluous nodes from the tree
    for (Node* p : dpseq1)
        {
        if (p->getFlag() == true)
            {
            if (p->numDescendants() == 1)
                {
                Node* d = p->getDescendants()[0];
                if (d == nullptr)
                    Msg::error("Could not find first descendant");
                Node* a = p->getAncestor();
                if (a != nullptr)
                    {
                    d->setAncestor(a);
                    a->removeDescendant(p);
                    a->addDescendant(d);
                    d->setBranchLength( d->getBranchLength() + p->getBranchLength() );
                    }
                else
                    {
                    d->setAncestor(nullptr);
                    d->setBranchLength(0.0);
                    if (p == t.getRoot())
                        t.setRoot(d);
                    }
                p->removeDescendants();
                p->setAncestor(nullptr);
                }
            }
        }
        
    // 4. reinitialize down pass sequence
    t.initializeDownPassSequence();
    const std::vector<Node*>& dpseq2 = t.getPostOrder();
            
    // 5. reindex nodes
    int intIdx = static_cast<int>(BitMask::numSet(taxonMask, canonicalTaxonList.size()));
    for (Node* p : dpseq2)
        {
        if (p->getIsLeaf() == true)
            {
            std::map<int,int>::iterator it = tipMap.find(p->getIndex());
            if (it == tipMap.end())
                Msg::error("Could not find tip index in tipMat for tip reindexing");
            p->setIndex(it->second);
            }
        else
            {
            p->setIndex(intIdx++);
            }
        }
        
    // 6. update the number of taxa and number of nodes
    t.numTaxa = 0;
    for (Node* p : dpseq2)
        {
        if (p->getIsLeaf() == true)
            t.numTaxa++;
        }
    if (t.numTaxa != BitMask::numSet(taxonMask, canonicalTaxonList.size()))
        Msg::error("Inconsistency between the number of tips and number of set bits when pruning taxa");
}

void ParameterTree::print(void) {

    std::cout << "   Canonical tree" << std::endl;
    fullTree.trees[0]->print();
    int i = 0;
    for (auto [key,val] : subTrees)
        {
        std::cout << "   Subtree " << ++i << ": " << BitMask::unsignedString(key, getNumTaxa()) << std::endl;
        val.trees[0]->print();
        }
}

void ParameterTree::restore(void) {

    // Called when MCMC proposal is rejected
    // Copy state from trees[1] back to trees[0] for all trees
    
    if (topologyChanged)
        {
        // Topology changed: do full restore from backup
        restoreTopology();
        topologyChanged = false;
        }
    else
        {
        // Only branch lengths changed: just copy branch lengths back
        const std::vector<Node*>& srcPostOrder = fullTree.trees[1]->getPostOrder();
        const std::vector<Node*>& dstPostOrder = fullTree.trees[0]->getPostOrder();
        for (size_t i=0; i<srcPostOrder.size(); i++)
            {
            Node* srcNode = srcPostOrder[i];
            Node* dstNode = dstPostOrder[i];
            dstNode->setBranchLength(srcNode->getBranchLength());
            }
        
        // restore subtree branch lengths
        for (auto& [mask, treePair] : subTrees)
            {
            std::unordered_map<unsigned, BranchMapping>::iterator it = branchMappings.find(mask);
            if (it != branchMappings.end())
                it->second.copyBranchLengths(treePair.trees[1], treePair.trees[0]);
            }
        }
}

void ParameterTree::restoreTopology(void) {

    // Full topology restore from trees[1] to trees[0]
    // Use the clone operation which copies everything including topology
    
    *fullTree.trees[0] = *fullTree.trees[1];
    
    for (auto& [mask, treePair] : subTrees)
        {
        *treePair.trees[0] = *treePair.trees[1];
        }
    
    // Rebuild branch mappings since we restored topology
    initializeBranchMappings();
}

void ParameterTree::saveTopology(void) {

    // Full topology save from trees[0] to trees[1]
    // Use the clone operation which copies everything including topology
    
    *fullTree.trees[1] = *fullTree.trees[0];
    
    for (auto& [mask, treePair] : subTrees)
        {
        *treePair.trees[1] = *treePair.trees[0];
        }
}

void ParameterTree::setBranchLength(Node* p, double v) {

    // only modify trees[0] (the working copy)
    double delta = v - p->getBranchLength();
    
    // update full tree (trees[0])
    p->setBranchLength(v);
    
    // propagate to all subtrees (trees[0] only)
    for (auto& [mask, treePair] : subTrees)
        {
        std::unordered_map<unsigned, BranchMapping>::iterator it = branchMappings.find(mask);
        if (it != branchMappings.end())
            it->second.propagateChange(p, delta, treePair.trees[0]);
        }
}

size_t ParameterTree::taxonIndex(const std::string& tName) {

    for (size_t i=0, n=canonicalTaxonList.size(); i<n; i++)
        {
        if (tName == canonicalTaxonList[i])
            return i;
        }
    return canonicalTaxonList.size();
}
