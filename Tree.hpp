#ifndef Tree_hpp
#define Tree_hpp

#include <map>
#include <sstream>
#include <string>
#include <vector>
#include "RbBitSet.h"
class Node;
class ParameterStatistics;


class Tree {

    public:
                                    Tree(void) = delete;
                                    Tree(std::string newickString, std::map<int, std::string> translateMap);
                                    Tree(std::map<RbBitSet,ParameterStatistics*>& partitions, std::map<int,std::string> translateMap);
        std::string                 getNewick(int brlenPrecision);
        std::string                 getNewick(int brlenPrecision, std::map<RbBitSet,ParameterStatistics*>& partitions);
        std::map<RbBitSet,double>   getPartitions(void);
        void                        initializeDownPassSequence(void);
        void                        print(void);
    
    private:
        Node*                       addNode(void);
        void                        listNode(Node* p, int indent);
        std::vector<std::string>    parseNewickString(std::string ns);
        void                        passDown(Node* p);
        void                        writeTree(Node* p, std::stringstream& ss, int brlenPrecision);
        void                        writeData(Node* p, std::stringstream& ss, int brlenPrecision);
        std::vector<Node*>          nodes;
        std::vector<Node*>          downPassSequence;
        Node*                       root;

    public:
        int                         numTaxa;
};

#endif
