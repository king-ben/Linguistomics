#include <iomanip>
#include <iostream>
#include "Msg.hpp"
#include "Node.hpp"
#include "ParameterStatistics.hpp"
#include "Tree.hpp"



Tree::Tree(std::string newickString, std::map<int,std::string> translateMap) {

    std::vector<std::string> tokens = parseNewickString(newickString);
    
    Node* p = NULL;
    bool readingBrlen = false;
    numTaxa = 0;
    for (int i=0; i<tokens.size(); i++)
        {
        std::string token = tokens[i];
        if (token == "(")
            {
            if (p == NULL)
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
            if (p->getAncestor() == NULL)
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
            {
            pd->setIndex(intIdx++);
            }
        }
}

Tree::Tree(std::map<RbBitSet,ParameterStatistics*>& partitions, std::map<int,std::string> translateMap) {

    // get the number of samples, which should be equal to the maximum number of partitions sampled
    int num = 0;
    numTaxa = -1;
    for (std::map<RbBitSet,ParameterStatistics*>::iterator it = partitions.begin(); it != partitions.end(); it++)
        {
        int n = it->second->size();
        if (n > num)
            num = n;
        if (numTaxa < 0)
            numTaxa = (int)it->first.size();
        }

    std::map<RbBitSet,double,CompBitSetSize> sortedParts;
    for (std::map<RbBitSet,ParameterStatistics*>::iterator it = partitions.begin(); it != partitions.end(); it++)
        {
        double prob = (double)it->second->size()/num;
        if (prob >= 0.5)
            sortedParts.insert( std::make_pair(RbBitSet(it->first),prob) );
        }
        
    for (std::map<RbBitSet,double,CompBitSetSize>::iterator it = sortedParts.begin(); it != sortedParts.end(); it++)
        {
        std::cout << it->first << " " << it->second << std::endl;
        }



    // add the tips
    for (int i=0; i<numTaxa; i++)
        {
        Node* newTip = addNode();
        newTip->setIndex(i);
        newTip->setIsLeaf(true);

        if (translateMap.size() > 0)
            {
            int key = i + 1;
            std::map<int,std::string>::iterator it = translateMap.find(key);
            if (it == translateMap.end())
                Msg::error("Could not find " + std::to_string(key) + " in the translate table");
            newTip->setName(it->second);
            }


        }
    root = NULL;
    int intIdx = numTaxa;
    for (std::map<RbBitSet,double,CompBitSetSize>::iterator it = sortedParts.begin(); it != sortedParts.end(); it++)
        {
        
        std::map<RbBitSet,ParameterStatistics*>::iterator pit = partitions.find(it->first);
        if (pit == partitions.end())
            Msg::error("Could not find partition");
        double v = pit->second->getMean();
        CredibleInterval ci = pit->second->getCredibleInterval();
        
        if (it->first.getNumberSetBits() >= 2)
            {
            Node* newInt = addNode();
            newInt->setIndex(intIdx++);
            newInt->setBrlen(v);
            newInt->setLowerCi(ci.lower);
            newInt->setUpperCi(ci.upper);
            newInt->setProb(it->second);
            if (root == NULL)
                root = newInt;
            for (int i=0; i<numTaxa; i++)
                {
                if (it->first.isSet(i) == true)
                    {
                    newInt->addNeighbor(nodes[i]);
                    nodes[i]->addNeighbor(newInt);
                    Node* a = nodes[i]->getAncestor();
                    nodes[i]->setAncestor(newInt);
                    if (a != NULL)
                        {
                        a->removeNeighbor(nodes[i]);
                        nodes[i]->removeNeighbor(a);
                        a->addNeighbor(newInt);
                        newInt->addNeighbor(a);
                        newInt->setAncestor(a);
                        }
                    }
                }
                
            }
        else
            {
            size_t idx = it->first.getFirstSetBit();
            nodes[idx]->setBrlen(v);
            nodes[idx]->setLowerCi(ci.lower);
            nodes[idx]->setUpperCi(ci.upper);
            nodes[idx]->setProb(it->second);
            }
            
        }

    initializeDownPassSequence();
    
    // assign probabilitie to the branches
    std::vector<RbBitSet> tempPartitions(downPassSequence.size());
    for (Node* p : downPassSequence)
        {
        RbBitSet& bs = tempPartitions[p->getIndex()];
        bs.resize(numTaxa);
        bs.unset();
        if (p->getIsLeaf() == true)
            {
            bs.set(p->getIndex());
            p->setCladeProbability(1.0);
            }
        else
            {
            std::vector<Node*> pDesc = p->getDescendants();
            for (Node* d : pDesc)
                bs |= tempPartitions[d->getIndex()];
            
            std::map<RbBitSet,ParameterStatistics*>::iterator it = partitions.find(bs);
            double prob = 0.0;
            if (it == partitions.end())
                {
                RbBitSet bs2 = bs;
                bs2.flip();
                it = partitions.find(bs2);
                if (it == partitions.end())
                    Msg::error("Could not find branch partition in list of partitions");
                prob = (double)it->second->size()/num;
                p->setCladeProbability(prob);
                }
            prob = (double)it->second->size()/num;
            p->setCladeProbability(prob);
            
            }
        }

}

Node* Tree::addNode(void) {

    Node* newNode = new Node;
    nodes.push_back(newNode);
    return newNode;
}

std::string Tree::getNewick(int brlenPrecision, std::map<RbBitSet,ParameterStatistics*>&) {

    

    return getNewick(brlenPrecision);
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

std::map<RbBitSet,double> Tree::getPartitions(void) {

    std::map<RbBitSet,double> parts;
    std::vector<RbBitSet> taxonBipartitiions;
    for (int i=0; i<nodes.size(); i++)
        taxonBipartitiions.push_back(RbBitSet(numTaxa));
    for (int i=0; i<numTaxa; i++)
        taxonBipartitiions[i].set(i);
    for (int i=0; i<downPassSequence.size(); i++)
        {
        Node* p = downPassSequence[i];
        if (p->getIsLeaf() == false)
            {
            std::vector<Node*> pDes = p->getDescendants();
            for (Node* n : pDes)
                taxonBipartitiions[p->getIndex()] |= taxonBipartitiions[n->getIndex()];
            }
        parts.insert( std::make_pair(taxonBipartitiions[p->getIndex()],p->getBrlen()) );
        }
    
//    for (std::map<RbBitSet,double>::iterator it = parts.begin(); it != parts.end(); it++)
//        std::cout << it->first << " " << it->second << std::endl;
    
    return parts;
}

void Tree::initializeDownPassSequence(void) {

    downPassSequence.clear();
    passDown(root);
}

void Tree::listNode(Node* p, int indent) {

    if (p != NULL)
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
        std::cout << ") " << p->getProb() << " ";
        std::cout << p->getBrlen() << " ";
        
        if (p->getUpperCi() - p->getLowerCi() > 0.00001)
            {
            std::cout << "(" << p->getLowerCi() << ", " << p->getUpperCi() << ") ";
            }
        
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

    if (p != NULL)
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

    listNode(root, 3);
}

void Tree::writeData(Node* p, std::stringstream& ss, int brlenPrecision) {
    ss << "[&Probability=" << p->getCladeProbability() << "]:" << std::fixed << std::setprecision(brlenPrecision) << p->getBrlen();
}

void Tree::writeTree(Node* p, std::stringstream& ss, int brlenPrecision) {

    if (p != NULL)
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
