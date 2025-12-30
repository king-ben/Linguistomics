#include <algorithm>
#include <iomanip>
#include <iostream>
#include "Msg.hpp"
#include "Node.hpp"
#include "ParameterStatistics.hpp"
#include "Statistics.hpp"
#include "Tree.hpp"



Tree::Tree(void) : root(nullptr), numTaxa(0) {

}

Tree::~Tree(void) {

    deleteAllNodes();
}

Tree::Tree(const Tree& t) : root(nullptr), numTaxa(0) {

    clone(t);
}

Tree& Tree::operator=(const Tree& t) {

    if (this != &t)
        clone(t);
    return *this;
}

Tree::Tree(std::string newickString, std::map<int,std::string> translateMap) : root(nullptr), numTaxa(0) {

    std::vector<std::string> tokens = parseNewickString(newickString);
    
    Node* p = nullptr;
    bool readingBrlen = false;
    for (int i=0; i<tokens.size(); i++)
        {
        std::string token = tokens[i];
        if (token == "(")
            {
            if (p == nullptr)
                {
                root = addNode();
                p = root;
                }
            else
                {
                Node* newNode = addNode();
                p->addNeighbor(newNode);
                newNode->addNeighbor(p);
                newNode->setAncestor(p);
                p = newNode;
                }
            readingBrlen = false;
            }
        else if (token == ")" || token == ",")
            {
            if (p->getAncestor() == nullptr)
                Msg::error("Expecting non null ancestor");
            p = p->getAncestor();
            readingBrlen = false;
            }
        else if (token == ":")
            {
            readingBrlen = true;
            }
        else if (token == ";")
            {
            if (p != root)
                Msg::warning("Expecting to end at the root");
            readingBrlen = false;
            }
        else
            {
            if (readingBrlen == false)
                {
                numTaxa++;
                Node* newNode = addNode();
                newNode->setIsLeaf(true);
                if (translateMap.size() > 0)
                    {
                    int key = atoi(token.c_str());
                    std::map<int,std::string>::iterator it = translateMap.find(key);
                    if (it == translateMap.end())
                        Msg::error("Could not find " + std::to_string(key) + " in the translate table");
                    newNode->setIndex(key-1);
                    newNode->setName(it->second);
                    }
                else
                    {
                    newNode->setIndex(numTaxa-1);
                    newNode->setName(token);
                    }
                p->addNeighbor(newNode);
                newNode->addNeighbor(p);
                newNode->setAncestor(p);
                p = newNode;
                }
            else
                {
                double x = atof(token.c_str());
                p->setBrlen(x);
                }
            readingBrlen = false;
            }
        }
        
    initializeDownPassSequence();
    
    int intIdx = numTaxa;
    for (int i=0; i<downPassSequence.size(); i++)
        {
        Node* pd = downPassSequence[i];
        if (pd->getIsLeaf() == false)
            pd->setIndex(intIdx++);
        }
}

Node* Tree::allocateNode(void) {

    return new Node;
}

Node* Tree::addNode(void) {

    Node* newNode = allocateNode();
    newNode->setOffset(nodes.size());
    nodes.push_back(newNode);
    return newNode;
}

void Tree::clone(const Tree& t) {

    // copy some instance variables
    numTaxa = t.numTaxa;
    
    // make certain we have the same number of nodes in each tree
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
        pLft->setName( pRht->getName() );
        pLft->setBrlen( pRht->getBrlen() );
        
        if (pRht->getAncestor() != nullptr)
            pLft->setAncestor( nodes[pRht->getAncestor()->getOffset()] );
        else
            pLft->setAncestor(nullptr);
            
        pLft->removeAllNeighbors();
        std::set<Node*>& rhtNeighbors = pRht->getNeighbors();
        for (Node* n : rhtNeighbors)
            pLft->addNeighbor( nodes[n->getOffset()] );
        }
            
    // copy the down pass sequence
    downPassSequence.clear();
    for (size_t i=0; i<t.downPassSequence.size(); i++)
        downPassSequence.push_back( nodes[t.downPassSequence[i]->getOffset()] );
}

void Tree::deleteAllNodes(void) {

    for (Node* p : nodes)
        delete p;
    nodes.clear();
    downPassSequence.clear();
    root = nullptr;
}

std::string Tree::getNewick(int brlenPrecision) {

    std::stringstream ss;
    if (root->getIsLeaf() == true)
        {
        Node* oldRoot = root;
        std::vector<Node*> nbs = root->getDescendants();
        if (nbs.size() > 1)
            Msg::error("Expecting only a single neighbor at the root of the tree");
        Node* newRoot = nbs[0];
        root->setAncestor(newRoot);
        oldRoot->setAncestor(newRoot);
        newRoot->setAncestor(nullptr);
        root = newRoot;

        writeTree(root, ss, brlenPrecision);

        newRoot->setAncestor(oldRoot);
        oldRoot->setAncestor(nullptr);
        root = oldRoot;
        }
    else
        {
        writeTree(root, ss, brlenPrecision);
        }
    std::string newick = ss.str();
    return newick;
}

std::map<BitSet,double> Tree::getPartitions(void) {

    std::map<BitSet,double> parts;
    std::vector<BitSet> taxonBipartitions;
    for (int i=0; i<nodes.size(); i++)
        taxonBipartitions.push_back(BitSet(numTaxa));
    for (int i=0; i<numTaxa; i++)
        taxonBipartitions[i].set(i);
    for (int i=0; i<downPassSequence.size(); i++)
        {
        Node* p = downPassSequence[i];
        if (p->getIsLeaf() == false)
            {
            std::vector<Node*> pDes = p->getDescendants();
            for (Node* n : pDes)
                taxonBipartitions[p->getIndex()] |= taxonBipartitions[n->getIndex()];
            }
        parts.insert( std::make_pair(taxonBipartitions[p->getIndex()], p->getBrlen()) );
        }
    
    return parts;
}

void Tree::initializeDownPassSequence(void) {

    downPassSequence.clear();
    passDown(root);
}

void Tree::listNode(Node* p, int indent) {

    if (p != nullptr)
        {
        std::set<Node*>& neighbors = p->getNeighbors();
        
        for (int i=0; i<indent; i++)
            std::cout << " ";
        std::cout << p->getIndex() << " ( ";
        for (Node* n : neighbors)
            {
            if (n == p->getAncestor())
                std::cout << "a_";
            std::cout << n->getIndex() << " ";
            }
        std::cout << ") " << p->getBrlen() << " ";
        std::cout << p->getName() << " ";
        std::cout << std::endl;
        for (Node* n : neighbors)
            {
            if (n != p->getAncestor())
                listNode(n, indent + 3);
            }
        }
}

std::vector<std::string> Tree::parseNewickString(std::string ns) {

    std::vector<std::string> tokens;
    std::string word = "";
    for (int i=0; i<ns.length(); i++)
        {
        char c = ns[i];
        std::string myChar(1, c);
        
        if (c == '(' || c == ')' || c == ',' || c == ':' || c == ';')
            {
            if (word != "")
                tokens.push_back(word);
            word = "";
            tokens.push_back(myChar);
            }
        else
            {
            word += myChar;
            }
        }
    return tokens;
}

void Tree::passDown(Node* p) {

    if (p != nullptr)
        {
        std::set<Node*>& neighbors = p->getNeighbors();
        for (Node* n : neighbors)
            {
            if (n != p->getAncestor())
                passDown(n);
            }
        downPassSequence.push_back(p);
        }
}

void Tree::print(void) {

    std::cout << "Tree (" << numTaxa << " taxa, " << nodes.size() << " nodes):" << std::endl;
    listNode(root, 3);
}

void Tree::writeData(Node* p, std::stringstream& ss, int brlenPrecision) {

    ss << ":" << std::fixed << std::setprecision(brlenPrecision) << p->getBrlen();
}

void Tree::writeTree(Node* p, std::stringstream& ss, int brlenPrecision) {

    if (p != nullptr)
        {
        if (p->getIsLeaf() == true)
            {
            ss << p->getName();
            writeData(p, ss, brlenPrecision);
            }
        else
            {
            ss << "(";
            }

        std::vector<Node*> myDescendants = p->getDescendants();
        for (int i=0; i<(int)myDescendants.size(); i++)
            {
            writeTree(myDescendants[i], ss, brlenPrecision);
            if ( (i + 1) != (int)myDescendants.size() )
                ss << ",";
            }
        if (p->getIsLeaf() == false)
            {
            ss << ")";
            if (p != root)
                writeData(p, ss, brlenPrecision);
            }
        }
}


// ConsensusTree implementation

ConsensusTree::ConsensusTree(void) : Tree() {

}

ConsensusTree::ConsensusTree(std::map<BitSet,ParameterStatistics*>& partitions, int numSamples, std::map<int,std::string>& translateMap) : Tree() {

    // determine the number of taxa from the first partition
    if (partitions.empty())
        Msg::error("No partitions provided for consensus tree");
    numTaxa = static_cast<int>(partitions.begin()->first.size());
    
    // sort partitions by probability (descending) and then by number of bits (descending)
    // this ensures we process larger clades first for proper tree building
    std::vector<std::pair<BitSet, double>> sortedParts;
    for (auto& kv : partitions)
        {
        double prob = static_cast<double>(kv.second->size()) / numSamples;
        if (prob >= 0.5)
            sortedParts.push_back(std::make_pair(kv.first, prob));
        }
    
    // sort by number of set bits descending (larger clades first)
    std::sort(sortedParts.begin(), sortedParts.end(),
        [](const std::pair<BitSet,double>& a, const std::pair<BitSet,double>& b) {
            return a.first.getNumberSetBits() > b.first.getNumberSetBits();
        });

    // add the tip nodes first
    for (int i=0; i<numTaxa; i++)
        {
        NodeConTree* newTip = static_cast<NodeConTree*>(addNode());
        newTip->setIndex(i);
        newTip->setIsLeaf(true);
        newTip->setProb(1.0f);
        if (translateMap.size() > 0)
            {
            int key = i + 1;
            auto it = translateMap.find(key);
            if (it == translateMap.end())
                Msg::error("Could not find " + std::to_string(key) + " in the translate table");
            newTip->setName(it->second);
            }
        }

    // build the tree by adding internal nodes for each partition
    root = nullptr;
    int intIdx = numTaxa;
    
    for (auto& part : sortedParts)
        {
        const BitSet& bs = part.first;
        double prob = part.second;
        
        // Get the statistics for this partition
        auto pit = partitions.find(bs);
        if (pit == partitions.end())
            Msg::error("Could not find partition in map");
        
        ParameterStatistics* stats = pit->second;
        std::pair<float,float> mv = Statistics::getMeanAndVariance(stats->getValues());
        double meanBrlen = mv.first;
        CredibleInterval ci = Statistics::getCredibleInterval(stats->getValues());
        
        size_t numBits = bs.getNumberSetBits();
        
        if (numBits >= 2)
            {
            // This is an internal node (clade)
            NodeConTree* newInt = static_cast<NodeConTree*>(addNode());
            newInt->setIndex(intIdx++);
            newInt->setBrlen(static_cast<float>(meanBrlen));
            newInt->setLowerCi(static_cast<float>(ci.lower));
            newInt->setUpperCi(static_cast<float>(ci.upper));
            newInt->setProb(static_cast<float>(prob));
            
            if (root == nullptr)
                root = newInt;
            
            // Find all taxa in this clade and connect them
            for (int i=0; i<numTaxa; i++)
                {
                if (bs.isSet(i) == true)
                    {
                    Node* taxonNode = nodes[i];
                    Node* currentAncestor = taxonNode->getAncestor();
                    
                    // Disconnect from current ancestor if it exists
                    if (currentAncestor != nullptr)
                        {
                        currentAncestor->removeNeighbor(taxonNode);
                        taxonNode->removeNeighbor(currentAncestor);
                        
                        // Connect current ancestor to new internal node if not already connected
                        if (currentAncestor != newInt)
                            {
                            bool alreadyConnected = false;
                            for (Node* n : newInt->getNeighbors())
                                {
                                if (n == currentAncestor)
                                    {
                                    alreadyConnected = true;
                                    break;
                                    }
                                }
                            if (!alreadyConnected)
                                {
                                currentAncestor->addNeighbor(newInt);
                                newInt->addNeighbor(currentAncestor);
                                newInt->setAncestor(currentAncestor);
                                }
                            }
                        }
                    
                    // Connect taxon to new internal node
                    newInt->addNeighbor(taxonNode);
                    taxonNode->addNeighbor(newInt);
                    taxonNode->setAncestor(newInt);
                    }
                }
            }
        else if (numBits == 1)
            {
            // This is a tip branch - set its branch length and stats
            size_t idx = bs.getFirstSetBit();
            NodeConTree* tipNode = static_cast<NodeConTree*>(nodes[idx]);
            tipNode->setBrlen(static_cast<float>(meanBrlen));
            tipNode->setLowerCi(static_cast<float>(ci.lower));
            tipNode->setUpperCi(static_cast<float>(ci.upper));
            }
        }

    initializeDownPassSequence();
}

Node* ConsensusTree::allocateNode(void) {

    return new NodeConTree;
}

std::string ConsensusTree::getNewick(int brlenPrecision) {

    // Use the parent's getNewick but it will call our writeData
    return Tree::getNewick(brlenPrecision);
}

void ConsensusTree::writeData(Node* p, std::stringstream& ss, int brlenPrecision) {

    NodeConTree* pCon = static_cast<NodeConTree*>(p);
    ss << "[&prob=" << std::fixed << std::setprecision(3) << pCon->getProb();
    ss << ",CI={" << std::setprecision(brlenPrecision) << pCon->getLowerCi() 
       << "," << pCon->getUpperCi() << "}]";
    ss << ":" << std::setprecision(brlenPrecision) << p->getBrlen();
}
