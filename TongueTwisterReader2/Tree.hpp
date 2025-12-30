#ifndef Tree_hpp
#define Tree_hpp

#include <map>
#include <sstream>
#include <string>
#include <vector>
#include "BitSet.hpp"
class Node;
class NodeConTree;
class ParameterStatistics;


class Tree {

    public:
                                    Tree(void);
                                    Tree(const Tree& t);
                                    Tree(std::string newickString, std::map<int, std::string> translateMap);
        virtual                    ~Tree(void);
        Tree&                       operator=(const Tree& t);
        std::string                 getNewick(int brlenPrecision);
        std::map<BitSet,double>     getPartitions(void);
        int                         getNumTaxa(void) { return numTaxa; }
        Node*                       getRoot(void) { return root; }
        void                        initializeDownPassSequence(void);
        void                        print(void);
        void                        writeTree(Node* p, std::stringstream& ss, int brlenPrecision);
    
    protected:
        virtual Node*               allocateNode(void);
        Node*                       addNode(void);
        void                        clone(const Tree& t);
        void                        deleteAllNodes(void);
        void                        listNode(Node* p, int indent);
        std::vector<std::string>    parseNewickString(std::string ns);
        void                        passDown(Node* p);
        virtual void                writeData(Node* p, std::stringstream& ss, int brlenPrecision);
        std::vector<Node*>          nodes;
        std::vector<Node*>          downPassSequence;
        Node*                       root;
        int                         numTaxa;
};

class ConsensusTree : public Tree {

    public:
                                    ConsensusTree(void);
                                    ConsensusTree(std::map<BitSet,ParameterStatistics*>& partitions, 
                                                  int numSamples,
                                                  std::map<int,std::string>& translateMap);
        std::string                 getNewick(int brlenPrecision);
    
    protected:
        Node*                       allocateNode(void) override;
        void                        writeData(Node* p, std::stringstream& ss, int brlenPrecision) override;
};

#endif
