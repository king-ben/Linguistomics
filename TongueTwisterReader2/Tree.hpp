#ifndef Tree_hpp
#define Tree_hpp

#include <map>
#include <sstream>
#include <string>
#include <vector>
#include "BitSet.hpp"
class Node;
class ParameterStatistics;


class Tree {

    public:
                                    Tree(void) = delete;
                                    Tree(const Tree& t);
                                    Tree(std::string newickString, std::map<int, std::string> translateMap);
                                    Tree(std::map<BitSet,ParameterStatistics*>& partitions, std::map<int,std::string> translateMap);
        std::string                 getNewick(int brlenPrecision);
        std::string                 getNewick(int brlenPrecision, std::map<BitSet,ParameterStatistics*>& partitions);
        std::map<BitSet,double>     getPartitions(void);
        void                        initializeDownPassSequence(void);
        void                        print(void);
    
    private:
        Node*                       addNode(void);
        void                        clone(const Tree& t);
        void                        deleteAllNodes(void);
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
