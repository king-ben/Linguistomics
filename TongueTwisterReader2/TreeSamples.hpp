#ifndef TreeSamples_hpp
#define TreeSamples_hpp

#include <map>
#include <memory>
#include <string>
#include <vector>
#include "BitSet.hpp"
class ConsensusTree;
class McmcSummary;
class ParameterStatistics;
class Tree;



class TreeSamples {

    public:
                                                    TreeSamples(void) = delete;
                                                    TreeSamples(McmcSummary& samples, double burnFraction);
                                                   ~TreeSamples(void);
        ConsensusTree*                              getConsensusTree(void);
        size_t                                      getNumTrees(void) { return trees.size(); }
        Tree*                                       getTree(size_t idx) { return trees[idx]; }
        ConsensusTree*                              getConsensusTree(void) { return consensusTree; }
        void                                        print(void);
        void                                        printPartitionSummary(void);
   
    private:
        void                                        computePartitions(void);
        std::vector<Tree*>                          trees;
        std::map<int,std::string>                   translateMap;
        std::map<BitSet,ParameterStatistics*>       partitions;
        ConsensusTree*                              consensusTree;
        bool                                        partitionsComputed;
};

#endif
