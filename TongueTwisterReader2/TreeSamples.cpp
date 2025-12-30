#include <algorithm>
#include <iomanip>
#include <iostream>
#include "McmcSummary.hpp"
#include "ParameterStatistics.hpp"
#include "Statistics.hpp"
#include "Tree.hpp"
#include "TreeSamples.hpp"



TreeSamples::TreeSamples(McmcSummary& samples, double burnFraction) : 
    consensusTree(nullptr), 
    partitionsComputed(false) {

    const std::vector<Tree*>& treeList = samples.getTrees();
    size_t startIdx = static_cast<size_t>(treeList.size() * burnFraction);
    
    // copy trees after burn-in
    for (size_t i=startIdx; i<treeList.size(); i++)
        trees.push_back(new Tree(*treeList[i]));
    
    // get translate map from McmcSummary
    translateMap = samples.getTranslateMap();
    
    getConsensusTree();
    consensusTree->print();
}

TreeSamples::~TreeSamples(void) {

    for (Tree* t : trees)
        delete t;
    
    for (auto& kv : partitions)
        delete kv.second;
    
    if (consensusTree != nullptr)
        delete consensusTree;
}

void TreeSamples::computePartitions(void) {

    if (partitionsComputed)
        return;
    
    // gather all partitions from all trees
    for (Tree* t : trees)
        {
        std::map<BitSet,double> treeParts = t->getPartitions();
        for (auto& kv : treeParts)
            {
            BitSet bs = kv.first;
            double brlen = kv.second;
            
            // add to partitions map
            auto it = partitions.find(bs);
            if (it == partitions.end())
                {
                ParameterStatistics* stats = new ParameterStatistics;
                stats->addValue(brlen);
                partitions[bs] = stats;
                }
            else
                {
                it->second->addValue(brlen);
                }
            }
        }
    
    partitionsComputed = true;
}

ConsensusTree* TreeSamples::getConsensusTree(void) {

    if (consensusTree != nullptr)
        return consensusTree;
    
    if (trees.empty())
        return nullptr;
    
    // make sure partitions are computed
    computePartitions();
    
    // build the consensus tree
    int numSamples = static_cast<int>(trees.size());
    consensusTree = new ConsensusTree(partitions, numSamples, translateMap);
    
    return consensusTree;
}

void TreeSamples::print(void) {

    std::cout << "Consensus Tree:" << std::endl;
    consensusTree->print();
    
//    std::cout << "TreeSamples: " << trees.size() << " trees" << std::endl;
//    for (size_t i=0; i<trees.size(); i++)
//        {
//        std::cout << "Tree " << i+1 << ":" << std::endl;
//        trees[i]->print();
//        }
}

void TreeSamples::printPartitionSummary(void) {

    computePartitions();
    
    int numSamples = static_cast<int>(trees.size());
    
    std::cout << "Partition Summary (" << partitions.size() << " unique partitions from " 
              << numSamples << " trees):" << std::endl;
    std::cout << std::string(70, '-') << std::endl;
    
    // Sort partitions by frequency
    std::vector<std::pair<BitSet, ParameterStatistics*>> sortedParts(partitions.begin(), partitions.end());
    std::sort(sortedParts.begin(), sortedParts.end(),
        [](const auto& a, const auto& b) {
            return a.second->size() > b.second->size();
        });
    
    std::cout << std::left << std::setw(40) << "Partition" 
              << std::right << std::setw(10) << "Count"
              << std::setw(10) << "Prob"
              << std::setw(10) << "Mean BL" << std::endl;
    std::cout << std::string(70, '-') << std::endl;
    
    for (auto& kv : sortedParts)
        {
        double prob = static_cast<double>(kv.second->size()) / numSamples;
        double meanBrlen = Statistics::getMeanAndVariance(kv.second->getValues()).first;
        
        // Print BitSet using operator<< if available, otherwise print summary
        std::cout << std::left << std::setw(40) << kv.first
                  << std::right << std::setw(10) << kv.second->size()
                  << std::setw(10) << std::fixed << std::setprecision(3) << prob
                  << std::setw(10) << std::setprecision(5) << meanBrlen << std::endl;
        }
}
