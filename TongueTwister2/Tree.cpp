#include <iomanip>
#include <iostream>
#include <set>
#include "BitMask.hpp"
#include "Msg.hpp"
#include "Node.hpp"
#include "NodeFactory.hpp"
#include "Probability.hpp"
#include "RandomVariable.hpp"
#include "Tree.hpp"



Tree::Tree(Tree& t) {

    clone(t);
}

Tree::Tree(RandomVariable* rng, std::string newickString, std::vector<std::string>& canonicalTaxonList, double lambda, double maxBrlen, int outgroupIdx) {

    // read the tree file
    root = nullptr;
    numTaxa = 0;

    // tokenize the Newick string
    std::vector<std::string> tokens;
    tokenizeTreeString(newickString, tokens);
    
    bool readingBrlen = false;
    Node* p = nullptr;
    int numOutgroups = 0;
    for (size_t i=0; i<tokens.size(); i++)
        {
        std::string token = tokens[i];
        if (token == "(")
            {
            if (p == nullptr)
                {
                Node* newNode = addNode();
                p = newNode;
                root = newNode;
                }
            else
                {
                Node* newNode = addNode();
                p->addDescendant(newNode);
                newNode->setAncestor(p);
                p = newNode;
                }
            readingBrlen = false;
            }
        else if (token == ")" || token == ",")
            {
            if (p->getAncestor() != nullptr)
                p = p->getAncestor();
            readingBrlen = false;
            }
        else if (token == ";")
            {
            readingBrlen = false;
            }
        else if (token == ":")
            {
            readingBrlen = true;
            }
        else
            {
            if (readingBrlen == false)
                {
                int n = 0;
                for (std::string s : canonicalTaxonList)
                    {
                    if (s == token)
                        break;
                    n++;
                    }
                                   
                Node* newNode = addNode();
                p->addDescendant(newNode);
                newNode->setAncestor(p);
                newNode->setName(token.c_str());
                newNode->setIsLeaf(true);
                newNode->setIndex(n);
                if (n == outgroupIdx)
                    {
                    newNode->setIsOutgroup(true);
                    numOutgroups++;
                    }
                p = newNode;
                numTaxa++;
                }
            else
                {
                double brlen = std::stod(token);
                p->setBranchLength(brlen);
                }
            readingBrlen = false;
            }
        }

    // check that we have a single outgroup specified
    if (numOutgroups != 1)
        Msg::error("Expected one outgroup language");
        
    // initialize down pass sequence
    initializeDownPassSequence();

    // reindex interior nodes
    int intIdx = (int)numTaxa;
    for (Node* p : postOrder)
        {
        if (p->getIsLeaf() == false)
            p->setIndex(intIdx++);
        }
        
    // initialize branch lengths from exponential(lambda) distribution
    // first, initialize all branch lengths, except the root node
    for (Node* p : postOrder)
        {
        if (p != root)
            {
            double v = 0.0;
            do {
                v = Probability::Exponential::rv(rng, lambda);
            } while (v > maxBrlen);
            p->setBranchLength(v);
            }
        else
            p->setBranchLength(0.0);
        }
        
    // The lengths of the two branches emminating from the root are not
    // identificable. Only their sum can be estimated. Here, I use the 
    // convention that the length of the branch leading to the outgroup
    // taxon is 0.0 in length. The other branch, then, is identifiable.
    if (root->numDescendants() != 2)
        Msg::error("Expecting two descendants of the root node");
    Node* d0 = root->getDescendant(0);
    Node* d1 = root->getDescendant(1);
    if (d0->getIsOutgroup() == true)
        d0->setBranchLength(0.0);
    else if (d1->getIsOutgroup() == true)
        d1->setBranchLength(0.0);
    else 
        Msg::error("Expected one descendant of the root to be the outgroup");
}

Tree::~Tree(void) {

    deleteAllNodes();
}

Tree& Tree::operator=(Tree& t) {

    if (this != &t)
        clone(t);
    return *this;
}

Node* Tree::addNode(void) {

    NodeFactory& nf = NodeFactory::nodeFactory();
    Node* newNode = nf.getNode();
    newNode->setOffset( (int)nodes.size() );
    nodes.push_back(newNode);
    return newNode;
}

void Tree::clone(Tree& t) {

    // copy some instance variables
    numTaxa = t.numTaxa;
    
    // make certain we have the saame number of nodes in each tree
    if (nodes.size() != t.nodes.size())
        {
        deleteAllNodes();
        for (size_t i=0; i<t.nodes.size(); i++)
            addNode();
        }

    // copy the nodes
    root = nodes[t.root->getOffset()];
    for (size_t i=0; i<t.nodes.size(); i++)
        {
        Node* pLft = nodes[i];
        Node* pRht = t.nodes[i];
        
        pLft->setIndex( pRht->getIndex() );
        pLft->setIsLeaf( pRht->getIsLeaf() );
        pLft->setIsOutgroup( pRht->getIsOutgroup() );
        pLft->setName( pRht->getName() );
        pLft->setBranchLength( pRht->getBranchLength() );
        
        if (pRht->getAncestor() != NULL)
            pLft->setAncestor( nodes[pRht->getAncestor()->getOffset()] );
        else
            pLft->setAncestor(NULL);
            
        pLft->removeDescendants();
        for (size_t j=0; j<pRht->numDescendants(); j++)
            {
            Node* rightDesc = pRht->getDescendant(j);
            pLft->addDescendant( nodes[rightDesc->getOffset()] );
            }
        }
            
    // copy the down pass sequence
    postOrder.clear();
    for (size_t i=0; i<t.postOrder.size(); i++)
        postOrder.push_back( nodes[t.postOrder[i]->getOffset()] );
}

void Tree::deleteAllNodes(void) {

    NodeFactory& nf = NodeFactory::nodeFactory();
    for (Node* p : nodes)
        nf.returnToPool(p);
    nodes.clear();
}

Node* Tree::getLeafIndexed(size_t idx) {

    for (size_t i=0; i<nodes.size(); i++)
        {
        if ((size_t)nodes[i]->getIndex() == idx)
            {
            if (nodes[i]->getIsLeaf() == true)
                return nodes[i];
            }
        }
    return nullptr;
}

Node* Tree::getLeafNamed(const char* n) {

    for (Node* p : nodes)
        {
        if (p->getIsLeaf() == true)
            {
            if (strcmp(n, p->getName()) == 0)
                return p;
            }
        }
    return nullptr;
}

void Tree::initializeDownPassSequence(void) {

    postOrder.clear();
    passDown(root);
    int idx = (int)numTaxa;
    for (Node* p : postOrder)
        {
        if (p->getIsLeaf() == false)
            p->setIndex(idx++);
        }
}

bool Tree::isBinary(void) {

    for (Node* p : postOrder)
        {
        size_t nd = p->numDescendants();
        if (!(nd == 0 || nd == 2))
            return false;
        }
    return true;
}

void Tree::listNodes(Node* p, size_t indent) {

    if (p != NULL)
        {
        p->orderDescendantsByOffset();
        
        std::cout << "  ";
        for (size_t i=0; i<indent; i++)
            std::cout << " ";
        std::cout << p->getIndex() << " ";
        std::cout << "( ";
        if (p->getAncestor() != NULL)
            std::cout << "a_" << p->getAncestor()->getIndex() << " ";
        else
        std::cout << "a_NULL ";
        for (size_t i=0; i<p->numDescendants(); i++)
            {
            Node* n = p->getDescendant(i);
            std::cout << n->getIndex() << " ";
            }
        std::cout << ") " << std::fixed << std::setprecision(5) << p->getBranchLength() << " ";
    
        if (p->getIsLeaf() == true)
            std::cout << " (" << p->getName() << ")";
            
        if (p == root)
            std::cout << " <- Root";
            
        std::cout << std::endl;

        for (size_t i=0; i<p->numDescendants(); i++)
            {
            Node* n = p->getDescendant(i);
            listNodes( n, indent+3 );
            }
        }
}

dMap Tree::pairwiseDistances(void) {

    // get list of taxa
    std::vector<std::string> taxonList;
    for (Node* p : postOrder)
        {
        if (p->getIsLeaf() == true)
            taxonList.push_back(std::string(p->getName()));
        }

    // enter taxon pairs as keys
    dMap distances;
    for (size_t i=0; i<taxonList.size(); i++)
        {
        for (size_t j=i+1; j<taxonList.size(); j++)
            {
            std::pair<std::string,std::string> key;
            if (taxonList[i] <taxonList[j])
                {
                key.first = taxonList[i];
                key.second = taxonList[j];
                }
            else 
                {
                key.first = taxonList[j];
                key.second = taxonList[i];
                }
            distances.insert( std::make_pair(key, 0.0) );
            }
        }

    // calculate the distances
    for (dMap::iterator it = distances.begin(); it != distances.end(); it++)
        {
        Node* p1 = getLeafNamed(it->first.first.c_str());
        Node* p2 = getLeafNamed(it->first.second.c_str());
        if (p1 == nullptr || p2 == nullptr)
            Msg::error("Could not find nodes for tip taxa");
            
        setAllFlags(false);
        Node* q = p1;
        do {
            q->setFlag(true);
            q = q->getAncestor();
            } while(q != nullptr);
        q = p2;
        do {
            if (q->getFlag() == true)
                q->setFlag(false);
            else
                q->setFlag(true);
            q = q->getAncestor();
            } while(q != nullptr);
        
        double sum = 0.0;
        for (Node* p : postOrder)
            {
            if (p->getFlag() == true)
                sum += p->getBranchLength();
            }
        it->second = sum;
        }
        
#   if 0
    std::cout << "Taxa: ";
    for (std::string& s : taxonList)
        std::cout << s << " ";
    std::cout << std::endl;
    for (auto [key,val] : distances)
        std::cout << "   " << key.first << "-" << key.second << ": " << val << std::endl;
#   endif
    
    return distances;
}

void Tree::passDown(Node* p) {

    if (p != NULL)
        {
        for (size_t i=0; i<p->numDescendants(); i++)
            passDown( p->getDescendant(i) );
        postOrder.push_back(p);
        }
}

void Tree::print(void) {

    std::cout << "     Num. Tips  = " << numTaxa << std::endl;
    std::cout << "     Num. Nodes = " << getNumNodes() << std::endl;
    listNodes(root, 3);
    std::cout << "     DPS: ";
    for (size_t i=0; i<postOrder.size(); i++)
        std::cout << postOrder[i]->getIndex() << " ";
    std::cout << std::endl << std::endl;
}

void Tree::print(std::string header) {

    std::cout << header << std::endl;
    print();
}

void Tree::setAllFlags(bool tf) {

    for (Node* p : nodes)
        p->setFlag(tf);
}

void Tree::tokenizeTreeString(std::string& newickString, std::vector<std::string>& tokens) {
    
    std::string longToken = "";
    for (size_t i=0; i<newickString.size(); i++)
        {
        char x = newickString.at(i);
        if (x ==')' || x == '(' || x == ',' || x == ';' || x == ':')
            {
            if (longToken != "")
                {
                tokens.push_back(longToken);
                longToken = "";
                }
            tokens.push_back( std::string(1,x) );
            }
        else
            {
            longToken += x;
            }
        }
        
#   if 0
    for (std::string s : tokens)
        {
        std::cout << s << std::endl;
        }
#   endif
}
