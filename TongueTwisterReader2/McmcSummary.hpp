#ifndef McmcSummary_hpp
#define McmcSummary_hpp

#include <map>
#include <string>
#include "FileManager.hpp"
class Alignment;
class ParameterStatistics;
class Partition;
class ThreadPool;
class Tree;



class McmcSummary {

    public:
                                                McmcSummary(void) = delete;
                                                McmcSummary(ThreadPool* tp, std::string dn);
                                               ~McmcSummary(void);
        ParameterStatistics*                    operator[](size_t idx);
        ParameterStatistics*                    operator[](size_t idx) const;
        std::vector<Alignment*>&                getAlignments(size_t idx) { return alignments[idx]; }
        std::string                             getAlignmentName(size_t idx) { return alignmentNames[idx]; }
        size_t                                  getNumStates(void) { return numStates; }
        Partition*                              getStatePartition(void) { return statePartitions; }
        const std::vector<Tree*>&               getTrees(void) const { return trees; }
        std::map<int,std::string>&              getTranslateMap(void) { return translateMap; }
        bool                                    hasFrequencies(void);
        bool                                    hasExchangeabilities(void);
        bool                                    hasPartition(void);
        bool                                    hasTrees(void);
        size_t                                  numParameterStatistics(void) { return stats.size(); }
        size_t                                  numAlignments(void) { return alignments.size(); }
        void                                    printParameterSummary(void);
        void                                    printParameterSummary(double burnFraction);
        size_t                                  size(void) { return stats.size(); }
    
    private:
        std::vector<std::string>                breakString(std::string str);
        bool                                    hasSemicolon(std::string str);
        std::map<int,std::string>               interpretTranslateString(std::vector<std::string> translateTokens);
        std::string                             interpretTreeString(std::string str);
        void                                    readAlnFile(std::string fn, size_t n, size_t total);
        void                                    readAlignmentFiles(void);
        void                                    readConfigurationFile(void);
        void                                    readParameterFile(void);
        void                                    readTreeFile(void);
        ThreadPool*                             pool;
        FileManager                             fileManager;
        std::vector<std::string>                taxa;
        size_t                                  numStates;
        std::vector<ParameterStatistics*>       stats;
        std::vector<std::vector<Alignment*>>    alignments;
        std::vector<std::string>                alignmentNames;
        Partition*                              statePartitions;
        std::vector<Tree*>                      trees;
        std::map<int,std::string>               translateMap;
};

#endif 
