#include "McmcSummary.hpp"
#include "Tree.hpp"
#include "TreeSamples.hpp"
#include <iostream>



TreeSamples::TreeSamples(McmcSummary& samples, double burnFraction) {

    const std::vector<Tree*>& treeList = samples.getTrees();
    size_t startIdx = static_cast<size_t>(treeList.size() * burnFraction);
    for (size_t i=startIdx; i<treeList.size(); i++)
        trees.push_back( new Tree(*treeList[i]) );
}

TreeSamples::~TreeSamples(void) {

    for (Tree* t : trees)
        delete t;
}

void TreeSamples::print(void) {

    for (Tree* t : trees)
        t->print();
}
