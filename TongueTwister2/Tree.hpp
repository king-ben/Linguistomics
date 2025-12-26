#ifndef Tree_hpp
#define Tree_hpp

#include <map>
#include <string>
#include <vector>
class Node;
class RandomVariable;
typedef std::map<std::pair<std::string,std::string>,double> dMap;



class Tree {

    friend class ParameterTree;

    public:
                                    Tree(void) = delete;
                                    Tree(Tree& t);
                                    Tree(RandomVariable* rng, std::string newickString, std::vector<std::string>& canonicalTaxonList, double lambda, double maxBrlen, int outgroupIdx = 0);
                                   ~Tree(void);
        Tree&                       operator=(Tree& t);
        const std::vector<Node*>&   getPostOrder(void) const { return postOrder; }
        Node*                       getLeafIndexed(size_t idx);
        Node*                       getLeafNamed(const char* n);
        size_t                      getNumNodes(void) { return nodes.size(); }
        size_t                      getNumTaxa(void) { return numTaxa; }
        Node*                       getRoot(void) { return root; }
        bool                        isBinary(void);
        dMap                        pairwiseDistances(void);
        void                        print(void);
        void                        print(std::string header); 
        void                        setAllFlags(bool tf);
        void                        setRoot(Node* p) { root = p; }
        
    private:
        Node*                       addNode(void);
        void                        clone(Tree& t);
        void                        deleteAllNodes(void);
        void                        initializeDownPassSequence(void);
        void                        listNodes(Node* p, size_t indent);
        void                        passDown(Node* p);
        void                        tokenizeTreeString(std::string& newickString, std::vector<std::string>& tokens);
        size_t                      numTaxa;
        Node*                       root;
        std::vector<Node*>          nodes;
        std::vector<Node*>          postOrder;
};

#endif
