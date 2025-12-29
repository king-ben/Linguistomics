#include <iomanip>
#include <iostream>
#include <istream>
#include <sstream>
#include <vector>
#include "Msg.hpp"
#include "Node.hpp"
#include "NodeFactory.hpp"
#include "NodeSet.hpp"
#include "Probability.hpp"
#include "RandomVariable.hpp"
#include "Tree.hpp"



Tree::Tree(Tree& t) {

    clone(t);
}

Tree::Tree(std::vector<std::string> tNames, double betaT, RandomVariable* rv, bool ic) {

    buildRandomTree(tNames, betaT, rv, ic);
}

Tree::Tree(std::string treeStr, std::vector<std::string> tNames, double lambda, RandomVariable* rv, bool ic) {

    // read the tree file
    root = NULL;
    numTaxa = 0;
    taxonNames = tNames;
    isClock = ic;

    std::vector<std::string> tokens = tokenizeTreeString(treeStr);
    bool readingBrlen = false;
    Node* p = NULL;
    for (int i=0; i<tokens.size(); i++)
        {
        std::string token = tokens[i];
        //std::cout << "token = " << token << std::endl;
        if (token == "(")
            {
            if (p == NULL)
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
            if (p->getAncestor() != NULL)
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
                for (std::string s : tNames)
                    {
                    if (s == token)
                        break;
                    n++;
                    }
                
                Node* newNode = addNode();
                p->addDescendant(newNode);
                newNode->setAncestor(p);
                newNode->setName(token);
                newNode->setIsLeaf(true);
                newNode->setIndex(n);
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

    // initialize down pass sequence
    initializeDownPassSequence();

    // reindex interior nodes
    int intIdx = numTaxa;
    for (int i=0; i<downPassSequence.size(); i++)
        {
        p = downPassSequence[i];
        if (p->getIsLeaf() == false)
            p->setIndex(intIdx++);
        }

    // check for branch lengths, randomly initializing from prior
    // if they are anot all initialized from the Newick string
    if (isClock == false)
        {
        bool branchLengthsPresent = true;
        for (int i=0; i<downPassSequence.size(); i++)
            {
            p = downPassSequence[i];
            if (p != root && p->getBranchLength() < 0.0001)
                branchLengthsPresent = false;
            }
        if (branchLengthsPresent == false)
            {
            // initialize branch lengths from uniform Dirichlet
            // first, initialize all branch lengths, except the root node
            for (int i=0; i<downPassSequence.size(); i++)
                {
                p = downPassSequence[i];
                if (p != root)
                    p->setBranchLength(Probability::Exponential::rv(rv, lambda));
                else
                    p->setBranchLength(0.0);
                }
            // make certain the two branches incident to the root are the
            // same in length and considered one branch
            std::vector<Node*> rootDes = root->getDescendantsVector();
            if (rootDes.size() != 2)
                Msg::error("Expecting two descendants of the root node");
            double x = rootDes[0]->getBranchLength() + rootDes[1]->getBranchLength();
            rootDes[0]->setBranchLength(x/2.0);
            rootDes[1]->setBranchLength(x/2.0);
            }
        }
    else
        {
        // clock-constrained branch lengths
        initializeClockBranchLengths(rv, lambda);
        }
        
}

Tree::Tree(Tree& t, RbBitSet& taxonMask) {

    makeSubtree(t, taxonMask);
}

Tree::~Tree(void) {

    NodeFactory& nf = NodeFactory::nodeFactory();
    for (int i=0; i<nodes.size(); i++)
        nf.returnToPool(nodes[i]);
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
    newNode->setMyTree(this);
    nodes.push_back(newNode);
    return newNode;
}

void Tree::buildRandomTree(std::vector<std::string> tNames, double lambda, RandomVariable* rv, bool ic) {

    deleteAllNodes();
    
    isClock = ic;
    root = NULL;
    numTaxa = 0;
    taxonNames = tNames;

    // start with a two-species tree
    Node* n1 = addNode();
    Node* n2 = addNode();
    n1->setName(tNames[0]);
    n2->setName(tNames[1]);
    n1->setIsLeaf(true);
    n2->setIsLeaf(true);
    root = addNode();
    root->addDescendant(n1);
    root->addDescendant(n2);
    n1->setAncestor(root);
    n2->setAncestor(root);
    
    for (int i=2; i<tNames.size(); i++)
        {
        Node* newTip = addNode();
        Node* newInt = addNode();
        newTip->setName(tNames[i]);
        newTip->setIsLeaf(true);
        
        Node* p = nodes[rv->uniformRvInt(nodes.size())];
        Node* pAnc = p->getAncestor();
        
        if (pAnc == NULL)
            {
            p->setAncestor(newInt);
            newTip->setAncestor(newInt);
            newInt->addDescendant(p);
            newInt->addDescendant(newTip);
            newInt->setAncestor(NULL);
            root = newInt;
            }
        else
            {
            p->setAncestor(newInt);
            pAnc->removeDescendant(p);
            pAnc->addDescendant(newInt);
            newInt->addDescendant(p);
            newInt->addDescendant(newTip);
            newInt->setAncestor(pAnc);
            newTip->setAncestor(newInt);
            }
        }

    // initialize down pass sequence
    initializeDownPassSequence();

    // reindex interior nodes
    int intIdx = numTaxa;
    for (int i=0; i<downPassSequence.size(); i++)
        {
        Node* p = downPassSequence[i];
        if (p->getIsLeaf() == false)
            p->setIndex(intIdx++);
        }

    // initialize branch lengths
    // first, initialize all branch lengths, except the root node
    if (isClock == false)
        {
        // we do not have a clock constraint, so we choose all of the branch lengths
        // INDEPENDENTLY from a common exponential(lambda) distribution
        for (int i=0; i<downPassSequence.size(); i++)
            {
            Node* p = downPassSequence[i];
            if (p != root)
                p->setBranchLength(Probability::Exponential::rv(rv, lambda));
            else
                p->setBranchLength(0.0);
            }
        }
    else
        {
        // we have a clock constraint
        initializeClockBranchLengths(rv, lambda);
        }
        
    // make certain the two branches incident to the root are the
    // same in length and considered one branch
    std::vector<Node*> rootDes = root->getDescendantsVector();
    if (rootDes.size() != 2)
        Msg::error("Expecting two descendants of the root node");
    double x = rootDes[0]->getBranchLength() + rootDes[1]->getBranchLength();
    rootDes[0]->setBranchLength(x/2.0);
    rootDes[1]->setBranchLength(x/2.0);
    
    print();
}

double Tree::calculateDistance(std::string t1, std::string t2) {

    setAllFlags(false);
    Node* n1 = getLeafNamed(t1);
    Node* p = n1;
    while(p->getAncestor() != NULL)
        {
        p->setFlag(true);
        p = p->getAncestor();
        }
    Node* n2 = getLeafNamed(t2);
    p = n2;
    while(p->getAncestor() != NULL)
        {
        p->setFlag(true);
        p = p->getAncestor();
        }
        
    double d = 0.0;
    bool foundCommonAncestor = false;
    for (int i=0; i<downPassSequence.size(); i++)
        {
        p = downPassSequence[i];
        if (p->getFlag() == true)
            {
            std::set<Node*,CompNode>& des = p->getDescendants()->getNodes();
            int numFlagged = 0;
            for (Node* n : des)
                {
                if (n->getFlag() == true)
                    numFlagged++;
                }
            if (numFlagged > 1)
                foundCommonAncestor = true;
            
            if (foundCommonAncestor == false)
                d += p->getBranchLength();
            }
        }
    
    return d;
}

void Tree::clone(Tree& t) {

    // copy some instance variables
    numTaxa = t.numTaxa;
    taxonNames = t.taxonNames;
    isClock = t.isClock;
    
    // make certain we have the saame number of nodes in each tree
    if (nodes.size() != t.nodes.size())
        {
        deleteAllNodes();
        for (int i=0; i<t.nodes.size(); i++)
            addNode();
        }

    // copy the nodes
    root = nodes[t.root->getOffset()];
    for (int i=0; i<t.nodes.size(); i++)
        {
        Node* pLft = nodes[i];
        Node* pRht = t.nodes[i];
        
        pLft->setIndex( pRht->getIndex() );
        pLft->setIsLeaf( pRht->getIsLeaf() );
        pLft->setName( pRht->getName() );
        pLft->setBranchLength( pRht->getBranchLength() );
        pLft->setNodeTime( pRht->getNodeTime() );
        
        if (pRht->getAncestor() != NULL)
            pLft->setAncestor( nodes[pRht->getAncestor()->getOffset()] );
        else
            pLft->setAncestor(NULL);
            
        pLft->removeDescendants();
        std::set<Node*,CompNode>& rhtNeighbors = pRht->getDescendants()->getNodes();
        for (Node* n : rhtNeighbors)
            pLft->addDescendant( nodes[n->getOffset()] );
        }
            
    // copy the down pass sequence
    downPassSequence.clear();
    for (int i=0; i<t.downPassSequence.size(); i++)
        downPassSequence.push_back( nodes[t.downPassSequence[i]->getOffset()] );
}

void Tree::debugPrint(std::string h) {

    print(h);
    
    std::cout << "nodes: " << std::endl;;
    for (int i=0; i<nodes.size(); i++)
        {
        std::cout << std::setw(4) << nodes[i]->getIndex() << " " << nodes[i] << " < ";
        }

    std::cout << "root: " << root->getIndex() << " " << root << std::endl;
    std::cout << "numTaxa: " << numTaxa << std::endl;
    std::cout << "treeLength: " << getTreeLength() << std::endl;
    for (int i=0; i<taxonNames.size(); i++)
        std::cout << i << " " << taxonNames[i] << std::endl;
}

void Tree::deleteAllNodes(void) {

    NodeFactory& nf = NodeFactory::nodeFactory();
    for (int i=0; i<nodes.size(); i++)
        nf.returnToPool(nodes[i]);
    nodes.clear();
}

Node* Tree::getLeafIndexed(int idx) {

    for (int i=0; i<nodes.size(); i++)
        {
        if (nodes[i]->getIndex() == idx)
            {
            if (nodes[i]->getIsLeaf() == true)
                return nodes[i];
            }
        }
    return NULL;
}

Node* Tree::getLeafNamed(std::string n) {

    for (int i=0; i<nodes.size(); i++)
        {
        if (nodes[i]->getName() == n)
            {
            if (nodes[i]->getIsLeaf() == true)
                return nodes[i];
            }
        }
    return NULL;
}

void Tree::initializeDownPassSequence(void) {

    downPassSequence.clear();
    passDown(root);
    int idx = numTaxa;
    for (int i=0; i<downPassSequence.size(); i++)
        {
        Node* p = downPassSequence[i];
        if (p->getIsLeaf() == false)
            p->setIndex(idx++);
        }
}

std::vector<int> Tree::getAncestorIndices(void) {

    std::vector<int> ancIndices(downPassSequence.size());
    
    for (int i=0; i<downPassSequence.size(); i++)
        {
        Node* p = downPassSequence[i];
        int idx = p->getIndex();
        if (p == root)
            ancIndices[idx] = -1;
        else
            ancIndices[idx] = p->getAncestor()->getIndex();
        }
    
    return ancIndices;
}

std::string Tree::getNewick(int brlenPrecision) {

    std::stringstream ss;
    if (root->getIsLeaf() == true)
        {
        Node* oldRoot = root;
        std::vector<Node*> nbs = root->getDescendantsVector();
        if (nbs.size() > 1)
            Msg::error("Expecting only a single neighbor at the root of the tree");
        Node* newRoot = nbs[0];
        root->setAncestor(newRoot);
        oldRoot->setAncestor(newRoot);
        newRoot->setAncestor(NULL);
        root = newRoot;

        writeTree(root, ss, brlenPrecision);

        newRoot->setAncestor(oldRoot);
        oldRoot->setAncestor(NULL);
        root = oldRoot;
        }
    else
        {
        writeTree(root, ss, brlenPrecision);
        }
    std::string newick = ss.str();
    return newick;
}

int Tree::getTaxonNameIndex(std::string tName) {

    for (int i=0; i<taxonNames.size(); i++)
        {
        if (tName == taxonNames[i])
            return i;
        }
    return -1;
}

double Tree::getTreeLength(void) {

    double tl = 0.0;
    for (int i=0; i<downPassSequence.size(); i++)
        {
        Node* p = downPassSequence[i];
        if (p->getAncestor() != NULL)
            tl += p->getBranchLength();
        }
    return tl;
}

void Tree::initializeClockBranchLengths(RandomVariable* rv, double lam) {

    // check that the tree is rooted and binary
    if (isRooted() == false)
        Msg::error("Clock-constrained trees must be rooted");
    if (isBinary() == false)
        Msg::error("Clock-constrained trees must be binary");
        
    // draw branch times randomly
    for (Node* n : downPassSequence)
        {
        if (n->getIsLeaf() == true)
            {
            n->setNodeTime(0.0);
            }
        else
            {
            double oldestDescendantTime = n->oldestDescendantTime();
            n->setNodeTime( oldestDescendantTime + Probability::Exponential::rv(rv, lam) );
            }
        }
        
#   if 0
    // rescale such that root time is 1.0
    double rootTime = root->getNodeTime();
    for (Node* n : downPassSequence)
        {
        if (n->getIsLeaf() == false)
            {
            double x = n->getNodeTime();
            n->setNodeTime( x/rootTime );
            }
        }
#   endif
        
    // set branch lengths
    for (Node* n : downPassSequence)
        {
        if (n != root)
            {
            n->setBranchLength( n->getAncestor()->getNodeTime() - n->getNodeTime() );
            }
        else
            n->setBranchLength(0.0);
        }
}

bool Tree::isBinary(void) {

    for (int i=0; i<downPassSequence.size(); i++)
        {
        Node* p = downPassSequence[i];
        int ndes = (int)p->numDescendants();
        if (p->getIsLeaf() == true)
            {
            if (ndes != 0)
                return false;
            }
        else
            {
            if (ndes != 2)
                return false;
            }
        }
    return true;
}

bool Tree::isRooted(void) {

    for (Node* n : downPassSequence)
        {
        if (n->getIsLeaf() == false)
            {
            if (n->getDescendants()->size() != 2)
                return false;
            }
        }
    return true;
}

bool Tree::isTaxonPresent(std::string tn) {

    for (int i=0; i<taxonNames.size(); i++)
        {
        if (tn == taxonNames[i])
            return true;
        }
    return false;
}

void Tree::listNodes(Node* p, size_t indent) {

    if (p != NULL)
        {
        std::set<Node*,CompNode> des = p->getDescendants()->getNodes();
        
        for (int i=0; i<indent; i++)
            std::cout << " ";
        std::cout << p->getIndex() << " ( ";
        if (p->getAncestor() != NULL)
            std::cout << "a_" << p->getAncestor()->getIndex() << " ";
        else
        std::cout << "a_NULL ";
        for (Node* n : des)
            {
            std::cout << n->getIndex() << " ";
            }
        std::cout << ") " << std::fixed << std::setprecision(5) << p->getBranchLength() << " ";
    
        if (p->getIsLeaf() == true)
            std::cout << " (" << p->getName() << ")";
            
        if (p == root)
            std::cout << " <- Root";
            
        std::cout << std::endl;

        for (std::set<Node*>::iterator it = des.begin(); it != des.end(); it++)
            {
            listNodes( (*it), indent+3 );
            }
        }
}

void Tree::makeSubtree(Tree& t, RbBitSet& taxonMask) {

    clone(t);
#   if 0
    print("original tree");
#   endif
    
    // prune unsampled taxa
    // 1. mark all nodes that are part of the subtree
    setAllFlags(false);
    std::map<int,int> tipMap;
    for (int i=0, k=0; i<taxonMask.size(); i++)
        {
        if (taxonMask[i] == true)
            {
            Node* p = getLeafIndexed(i);
            if (p == NULL)
                Msg::error("Could not find leaf to mark path to root");
            while (p != NULL)
                {
                p->setFlag(true);
                p = p->getAncestor();
                }
            tipMap.insert( std::make_pair(i,k) );
            k++;
            }
        }
    // 2. remove nodes that are not part of the subtree
    for (int i=0; i<downPassSequence.size(); i++)
        {
        Node* p = downPassSequence[i];
        if (p->getFlag() == false)
            {
            p->removeDescendants();
            Node* a = p->getAncestor();
            if (a != NULL)
                a->removeDescendant(p);
            p->setAncestor(NULL);
            }
        }
    // 3. remove superfluous nodes from the tree
    for (int i=0; i<downPassSequence.size(); i++)
        {
        Node* p = downPassSequence[i];
        if (p->getFlag() == true)
            {
            if (p->numDescendants() == 1)
                {
                std::vector<Node*> des = p->getDescendantsVector();
                Node* d = des[0];
                Node* a = p->getAncestor();
                if (a != NULL)
                    {
                    d->setAncestor(a);
                    a->removeDescendant(p);
                    a->addDescendant(d);
                    d->setBranchLength( d->getBranchLength() + p->getBranchLength() );
                    }
                else
                    {
                    d->setAncestor(NULL);
                    d->setBranchLength(0.0);
                    if (p == root)
                        root = d;
                    }
                p->removeDescendants();
                p->setAncestor(NULL);
                }
            }
        }
    // 4. reinitialize down pass sequence
    initializeDownPassSequence();
    // 5. reindex nodes
    int intIdx = (int)taxonMask.getNumberSetBits();
    for (int i=0; i<downPassSequence.size(); i++)
        {
        Node* p = downPassSequence[i];
        if (p->getIsLeaf() == true)
            {
            std::map<int,int>::iterator it = tipMap.find(p->getIndex());
            if (it == tipMap.end())
                Msg::error("Could not find tip index in tipMat for tip reindexing");
            p->setIndex(it->second);
            }
        else
            {
            p->setIndex(intIdx++);
            }
        }
    // 6. take care of other tree variables
    std::vector<std::string> tempTaxonNames = taxonNames;
    taxonNames.clear();
    for (int i=0; i<tempTaxonNames.size(); i++)
        {
        if (taxonMask[i] == true)
            taxonNames.push_back(tempTaxonNames[i]);
        }
    numTaxa = (int)taxonNames.size();

#   if 0
    print("pruned tree");
#   endif

}

void Tree::newickStream(std::ofstream& strm, int brlenPrecision) {

    if (root->getIsLeaf() == true)
        {
        Node* oldRoot = root;
        std::vector<Node*> nbs = root->getDescendantsVector();
        if (nbs.size() > 1)
            Msg::error("Expecting only a single neighbor at the root of the tree");
        Node* newRoot = nbs[0];
        root->setAncestor(newRoot);
        oldRoot->setAncestor(newRoot);
        newRoot->setAncestor(NULL);
        root = newRoot;

        writeTree(root, strm, brlenPrecision);

        newRoot->setAncestor(oldRoot);
        oldRoot->setAncestor(NULL);
        root = oldRoot;
        }
    else
        {
        writeTree(root, strm, brlenPrecision);
        }
}

std::map<TaxonPair,double,CompTaxonPair> Tree::pairwiseDistances(void) {

    // find all of the taxon pairs for this tree
    std::map<TaxonPair,double,CompTaxonPair> d;
    for (int i=0; i<taxonNames.size(); i++)
        {
        for (int j=i+1; j<taxonNames.size(); j++)
            {
            TaxonPair pair(taxonNames[i], taxonNames[j]);
            d.insert( std::make_pair(pair,0.0) );
            }
        }
        
    // calculate the distances between each pair
    for (std::map<TaxonPair,double,CompTaxonPair>::iterator it = d.begin(); it != d.end(); it++)
        {
        double distance = calculateDistance(it->first.getTaxon1(), it->first.getTaxon2());
        it->second = distance;
        }
    
    return d;
}

void Tree::print(void) {

    listNodes(root, 3);
    std::cout << "DPS: ";
    for (int i=0; i<downPassSequence.size(); i++)
        std::cout << downPassSequence[i]->getIndex() << " ";
    std::cout << std::endl;
}

void Tree::print(std::string header) {

    std::cout << header << std::endl;
    print();
}

void Tree::passDown(Node* p) {

    if (p != NULL)
        {
        std::set<Node*,CompNode>& des = p->getDescendants()->getNodes();
        for (std::set<Node*>::iterator it = des.begin(); it != des.end(); it++)
            {
            passDown( (*it) );
            }
        downPassSequence.push_back(p);
        }
}

Node* Tree::randomNode(RandomVariable* rv) {

    return nodes[rv->uniformRvInt(nodes.size())];
}

void Tree::setAllFlags(bool tf) {

    for (int i=0; i<nodes.size(); i++)
        nodes[i]->setFlag(tf);
}

std::vector<std::string> Tree::tokenizeTreeString(std::string ls) {

    std::vector<std::string> tks;
    
    std::string longToken = "";
    for (int i=0; i<ls.size(); i++)
        {
        char x = ls.at(i);
        if (x ==')' || x == '(' || x == ',' || x == ';' || x == ':')
            {
            if (longToken != "")
                {
                tks.push_back(longToken);
                longToken = "";
                }
            std::string s(1, x);
            tks.push_back( s );
            }
        else
            {
            longToken += x;
            }
        }
#   if 0
    for (std::string s : tks)
        {
        std::cout << s << std::endl;
        }
#   endif
    
    return tks;
}

void Tree::writeTree(Node* p, std::stringstream& ss, int brlenPrecision) {

    if (p != NULL)
        {
        if (p->getIsLeaf() == true)
            {
            ss << p->getIndex()+1;
            //ss << p->getName();
            ss << ":" << std::fixed << std::setprecision(brlenPrecision) << p->getBranchLength();
            }
        else
            {
            ss << "(";
            }
        std::vector<Node*> myDescendants = p->getDescendantsVector();
        for (int i=0; i<(int)myDescendants.size(); i++)
            {
            writeTree(myDescendants[i], ss, brlenPrecision);
            if ( (i + 1) != (int)myDescendants.size() )
                ss << ",";
            }
        if (p->getIsLeaf() == false)
            {
            ss << ")";
            if (p != NULL && p != root)
                ss << ":" << std::fixed << std::setprecision(brlenPrecision) << p->getBranchLength();
            }
        }
}

void Tree::writeTree(Node* p, std::ofstream& strm, int brlenPrecision) {

    if (p != NULL)
        {
        if (p->getIsLeaf() == true)
            {
            strm << p->getIndex()+1;
            //ss << p->getName();
            strm << ":" << std::fixed << std::setprecision(brlenPrecision) << p->getBranchLength();
            }
        else
            {
            strm << "(";
            }
        std::vector<Node*> myDescendants = p->getDescendantsVector();
        for (int i=0; i<(int)myDescendants.size(); i++)
            {
            writeTree(myDescendants[i], strm, brlenPrecision);
            if ( (i + 1) != (int)myDescendants.size() )
                strm << ",";
            }
        if (p->getIsLeaf() == false)
            {
            strm << ")";
            if (p != NULL && p != root)
                strm << ":" << std::fixed << std::setprecision(brlenPrecision) << p->getBranchLength();
            }
        }
}


