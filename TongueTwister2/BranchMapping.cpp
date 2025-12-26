#include <cmath>
#include <iostream>
#include <set>
#include "BranchMapping.hpp"
#include "Msg.hpp"
#include "Node.hpp"
#include "NodeComparator.hpp"
#include "Tree.hpp"

// Minimum allowed branch length to prevent numerical issues
static const double MIN_BRANCH_LENGTH = 1.0e-10;


void BranchMapping::buildMapping(Tree* fullTree, Tree* subTree) {

    subToFullBranches.clear();
    fullToSubBranch.clear();
    
    const std::vector<Node*>& subPostOrder = subTree->getPostOrder();
    Node* subRoot = subTree->getRoot();
    
    // for each subtree node (except root), find contributing full tree branches
    for (Node* subNode : subPostOrder)
        {
        if (subNode == subRoot)
            continue;
            
        int subOffset = subNode->getOffset();
        
        // get taxa descended from this subtree node
        std::vector<std::string> subTaxa;
        getDescendantLeaves(subNode, subTaxa);
        
        // get taxa descended from this node's ancestor in subtree
        std::vector<std::string> ancTaxa;
        Node* subAnc = subNode->getAncestor();
        if (subAnc != nullptr)
            getDescendantLeaves(subAnc, ancTaxa);
        
        // find corresponding nodes in full tree
        Node* fullNode = findMRCA(fullTree, subTaxa);
        Node* fullAnc = findMRCA(fullTree, ancTaxa);
        
        if (fullNode == nullptr)
            Msg::error("Could not find MRCA for subtree node in BranchMapping::buildMapping");
        
        if (fullAnc == nullptr)
            {
            // This can happen if ancTaxa is empty, which means subAnc has no descendants
            // This shouldn't happen for a valid tree, but handle gracefully
            Msg::error("Could not find MRCA for subtree ancestor in BranchMapping::buildMapping");
            }
        
        // collect path from fullNode up to (but not including) fullAnc
        std::vector<int> fullOffsets;
        Node* q = fullNode;
        while (q != nullptr && q != fullAnc)
            {
            int fullOffset = q->getOffset();
            fullOffsets.push_back(fullOffset);
            fullToSubBranch[fullOffset] = subOffset;
            q = q->getAncestor();
            }
        
        // Validate: the mapping should not be empty for non-root subtree nodes
        if (fullOffsets.empty())
            {
            Msg::warning("Empty branch mapping for subtree node " + std::to_string(subOffset) + 
                         " - fullNode and fullAnc may be the same");
            }
            
        subToFullBranches[subOffset] = fullOffsets;
        }
}

void BranchMapping::propagateChange(Node* fullTreeNode, double delta, Tree* subTree) {

    std::unordered_map<int,int>::iterator it = fullToSubBranch.find(fullTreeNode->getOffset());
    if (it == fullToSubBranch.end())
        return;
        
    int subOffset = it->second;
    const std::vector<Node*>& subPostOrder = subTree->getPostOrder();
    
    for (Node* p : subPostOrder)
        {
        if (p->getOffset() == subOffset)
            {
            double newLength = p->getBranchLength() + delta;
            
            // Clamp to minimum positive value to prevent numerical issues
            if (newLength < MIN_BRANCH_LENGTH)
                {
                newLength = MIN_BRANCH_LENGTH;
                }
            
            p->setBranchLength(newLength);
            return;
            }
        }
}

void BranchMapping::copyBranchLengths(Tree* src, Tree* dst) {

    // trees are clones with identical topology, so nodes at same offset correspond
    const std::vector<Node*>& srcPostOrder = src->getPostOrder();
    const std::vector<Node*>& dstPostOrder = dst->getPostOrder();
    
    // postOrder vectors have same size and corresponding nodes at same index
    for (size_t i=0; i<srcPostOrder.size(); i++)
        dstPostOrder[i]->setBranchLength(srcPostOrder[i]->getBranchLength());
}

double BranchMapping::recomputeBranchLength(int subTreeNodeOffset, Tree* fullTree) {

    std::unordered_map<int, std::vector<int>>::iterator it = subToFullBranches.find(subTreeNodeOffset);
    if (it == subToFullBranches.end())
        return 0.0;
        
    const std::vector<Node*>& fullPostOrder = fullTree->getPostOrder();
    
    double sum = 0.0;
    for (int fullOffset : it->second)
        {
        for (Node* p : fullPostOrder)
            {
            if (p->getOffset() == fullOffset)
                {
                sum += p->getBranchLength();
                break;
                }
            }
        }
        
    return sum;
}

bool BranchMapping::affectsSubtree(int fullTreeNodeOffset) {

    return fullToSubBranch.count(fullTreeNodeOffset) > 0;
}

void BranchMapping::getDescendantLeaves(Node* p, std::vector<std::string>& leaves) {

    if (p->getIsLeaf() == true)
        {
        leaves.push_back(p->getName());
        }
    else
        {
        std::set<Node*,NodeComparator>& desc = p->getDescendants().getNodes();
        for (Node* d : desc)
            getDescendantLeaves(d, leaves);
        }
}

Node* BranchMapping::findMRCA(Tree* t, std::vector<std::string>& taxa) {

    if (taxa.empty() == true)
        return nullptr;
        
    const std::vector<Node*>& postOrder = t->getPostOrder();
    
    // count how many target taxa descend from each node
    std::unordered_map<Node*, int> taxaCount;
    for (Node* p : postOrder)
        {
        if (p->getIsLeaf() == true)
            {
            int count = 0;
            for (std::string& s : taxa)
                {
                if (s == p->getName())
                    {
                    count = 1;
                    break;
                    }
                }
            taxaCount[p] = count;
            }
        else
            {
            int count = 0;
            std::set<Node*,NodeComparator>& desc = p->getDescendants().getNodes();
            for (Node* d : desc)
                count += taxaCount[d];
            taxaCount[p] = count;
            }
        }
    
    // find deepest node with count == taxa.size()
    int targetCount = static_cast<int>(taxa.size());
    for (Node* p : postOrder)
        {
        if (taxaCount[p] == targetCount)
            return p;
        }
        
    return t->getRoot();
}

void BranchMapping::print(void) {

    std::cout << "Branch Mapping:" << std::endl;
    std::cout << "  Subtree -> Full Tree branches:" << std::endl;
    for (auto& [subOffset, fullOffsets] : subToFullBranches)
        {
        std::cout << "    Subtree node " << subOffset << " <- Full nodes: ";
        for (int fo : fullOffsets)
            std::cout << fo << " ";
        std::cout << std::endl;
        }
}
