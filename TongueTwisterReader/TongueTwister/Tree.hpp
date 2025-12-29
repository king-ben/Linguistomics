#ifndef Tree_H
#define Tree_H

#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <vector>
#include "RbBitSet.h"
#include "TaxonPair.hpp"
class Node;
class RandomVariable;



class Tree {

    public:
                                                    Tree(void) = delete;
                                                    Tree(std::vector<std::string> tNames, double betaT, RandomVariable* rv, bool ic);
                                                    Tree(std::string treeStr, std::vector<std::string> tNames, double betaT, RandomVariable* rv, bool ic);
                                                    Tree(Tree& t);
                                                    Tree(Tree& t, RbBitSet& taxonMask);
                                                   ~Tree(void);
        Tree&                                       operator=(Tree& t);
        void                                        buildRandomTree(std::vector<std::string> tNames, double betaT, RandomVariable* rv, bool ic);
        void                                        debugPrint(std::string h);
        std::string                                 getNewick(int brlenPrecision);
        int                                         getNumNodes(void) { return (int)nodes.size(); }
        int                                         getNumTaxa(void) { return numTaxa; }
        std::vector<int>                            getAncestorIndices(void);
        std::vector<Node*>&                         getDownPassSequence(void) { return downPassSequence; }
        bool                                        getIsClock(void) { return isClock; }
        Node*                                       getLeafIndexed(int idx);
        Node*                                       getLeafNamed(std::string n);
        Node*                                       getRoot(void) { return root; }
        std::vector<std::string>&                   getTaxonNames(void) { return taxonNames; }
        int                                         getTaxonNameIndex(std::string tName);
        double                                      getTreeLength(void);
        void                                        initializeDownPassSequence(void);
        bool                                        isBinary(void);
        bool                                        isRoot(Node* p) { return ((p == root) ? true : false); }
        bool                                        isTaxonPresent(std::string tn);
        void                                        makeSubtree(Tree& t, RbBitSet& taxonMask);
        void                                        newickStream(std::ofstream& strm, int brlenPrecision);
        std::map<TaxonPair,double,CompTaxonPair>    pairwiseDistances(void);
        void                                        print(void);
        void                                        print(std::string header);
        Node*                                       randomNode(RandomVariable* rv);
        void                                        setAllFlags(bool tf);
        void                                        setIsClock(bool tf) { isClock = tf; }
                
    private:
        Node*                                       addNode(void);
        double                                      calculateDistance(std::string t1, std::string t2);
        void                                        clone(Tree& t);
        void                                        deleteAllNodes(void);
        void                                        initializeClockBranchLengths(RandomVariable* rv, double lam);
        bool                                        isRooted(void);
        void                                        listNodes(Node* p, size_t indent);
        void                                        passDown(Node* p);
        std::vector<std::string>                    tokenizeTreeString(std::string ls);
        void                                        writeTree(Node* p, std::stringstream& ss, int brlenPrecision);
        void                                        writeTree(Node* p, std::ofstream& strm, int brlenPrecision);
        std::vector<Node*>                          nodes;
        std::vector<Node*>                          downPassSequence;
        Node*                                       root;
        int                                         numTaxa;
        std::vector<std::string>                    taxonNames;
        bool                                        isClock;
};

#endif
