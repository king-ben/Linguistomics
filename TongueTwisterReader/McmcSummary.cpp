#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include "Alignment.hpp"
#include "AlignmentDistribution.hpp"
#include "Container.hpp"
#include "json.hpp"
#include "McmcSummary.hpp"
#include "Msg.hpp"
#include "Partition.hpp"
#include "RandomVariable.hpp"
#include "Subset.hpp"
#include "Tree.hpp"



McmcSummary::McmcSummary(RandomVariable* r) : conTree(NULL), rv(r) {

    statePartitions = NULL;
}

McmcSummary::~McmcSummary(void) {

    for (int i=0; i<stats.size(); i++)
        delete stats[i];
    if (statePartitions != NULL)
        delete statePartitions;
}

McmcSummary& McmcSummary::operator+=(const McmcSummary& rhs) {

    // combine the ParameterStatistics
    for (int i=0; i<stats.size(); i++)
        {
        std::string thisName = stats[i]->getName();
        ParameterStatistics* rhsStats = rhs.getParameterNamed(thisName);
        if (rhsStats != NULL)
            *(stats[i]) += *rhsStats;
        }
        
    // combine the alignments
    for (int i=0; i<alignments.size(); i++)
        {
        std::string thisName = alignments[i]->getName();
        AlignmentDistribution* rhsAln = rhs.getAlignmentNamed(thisName);
        if (rhsAln != NULL)
            *(alignments[i]) += *rhsAln;
        }
        
    // combine the tree partition information
    for (std::map<RbBitSet,ParameterStatistics*>::const_iterator it = rhs.partitions.begin(); it != rhs.partitions.end(); it++)
        {
        std::map<RbBitSet,ParameterStatistics*>::iterator partitionsIt = this->partitions.find(it->first);
        if (partitionsIt == this->partitions.end())
            {
            ParameterStatistics* newStats = new ParameterStatistics(*(it->second));
            partitions.insert( std::make_pair(it->first,newStats) );
            }
        else
            {
            *(partitionsIt->second) += *(it->second);
            }
        }

    return *this;
}

void McmcSummary::addPartion(std::map<RbBitSet,double>& parts) {

    for (std::map<RbBitSet,double>::iterator it = parts.begin(); it != parts.end(); it++)
        {
        std::map<RbBitSet,ParameterStatistics*>::iterator partitionsIt = partitions.find(it->first);
        if (partitionsIt == partitions.end())
            {
            ParameterStatistics* newStats = new ParameterStatistics;
            newStats->addValue(it->second);
            newStats->setName(it->first.bitString());
            partitions.insert( std::make_pair(it->first,newStats) );
            }
        else
            {
            partitionsIt->second->addValue(it->second);
            }
        }
}

std::vector<std::string> McmcSummary::breakString(std::string str) {

    std::vector<std::string> broken;
    std::string word = "";
    for (int i=0; i<str.length(); i++)
        {
        char c = str[i];
        if (c == ',' || c == ';')
            {
            if (word != "")
                broken.push_back(word);
            broken.push_back( std::string(1,c) );
            word = "";
            }
        else
            {
            word +=  c;
            }
        }
    if (word != "")
        broken.push_back(word);
    
    return broken;
}

void McmcSummary::calculateRates(DoubleMatrix& m) {

    int ns = 0;
    
    // 1. get state frequencies (and the number of states). We look
    // up the apropriate ParameterStatistics objects for this
    std::set<int> observedStates;
    for (int i=0; i<stats.size(); i++)
        {
        ParameterStatistics* s = stats[i];
        if (s->getName()[0] == 'F')
            {
            int st = getFreqElement(s->getName());
            observedStates.insert(st);
            }
        }
    for (int i=0; i<observedStates.size(); i++) // make certain all of the states from 0 to (numStates-1) are found
        {
        std::set<int>::iterator it = observedStates.find(i);
        if (it == observedStates.end())
            Msg::error("Could not find state " + std::to_string(i) + " in frequencies summary");
        }
    ns = (int)observedStates.size();
    if (ns != 0 && ns != numStates)
        Msg::error("Inconsistent number of states (calculateRates)");
    std::vector<double> f(numStates);
    for (int i=0; i<f.size(); i++)
        f[i] = 1.0 / numStates;
    for (int i=0; i<stats.size(); i++)
        {
        ParameterStatistics* s = stats[i];
        if (s->getName()[0] == 'F')
            {
            int st = getFreqElement(s->getName());
            f[st] = s->getMean();
            }
        }
    if (statePartitions != NULL && ns != 0) // check consistencey with a state partition object, if we have one
        {
        int tempNumStates = statePartitions->getNumElements();
        if (tempNumStates != numStates)
            Msg::error("Mismatch in the number of states from the frequencies and the state partitions (calculateRates)");
        }

    int nr = inferNumberOfRates();

    // allocate a vector to hold the exchangeability parameters
    DoubleMatrix r(numStates,numStates);
    for (int i=0; i<numStates; i++)
        {
        for (int j=0; j<numStates; j++)
            r(i,j) = 1.0;
        }
        
    if (nr > 0)
        {
        bool usingCustom = false;
        if (nr > 0 && nr < numStates * (numStates-1)/2)
            usingCustom = true;

        for (int i=0; i<stats.size(); i++)
            {
            ParameterStatistics* s = stats[i];
            if (s->getName()[0] == 'R')
                {
                int r1, r2;
                getRateElements(s->getName(), r1, r2);
                if (usingCustom == true)
                    {
                    Subset* s1 = statePartitions->findSubsetIndexed(r1);
                    Subset* s2 = statePartitions->findSubsetIndexed(r2);
                    for (int x1 : s1->getValues())
                        {
                        for (int x2 : s2->getValues())
                            {
                            double ave = s->getMean();
                            r(x1,x2) = ave;
                            r(x2,x1) = ave;
                            }
                        }
                    
                    }
                else
                    {
                    r(r1-1,r2-1) = s->getMean();
                    r(r2-1,r1-1) = s->getMean();
                    }
                }
            }
        }
        
    for (int i=0; i<numStates; i++)
        {
        for (int j=0; j<numStates; j++)
            {
            if (r(i,j) < 0.0 && i != j)
                Msg::error("Could not find all rate parameters");
            }
        }


    // fill in the rates (Q)
    //DoubleMatrix m(numStates,numStates);
    double averageRate = 0.0;
    for (int i=0; i<numStates; i++)
        {
        double sum = 0.0;
        for (int j=0; j<numStates; j++)
            {
            if (i != j)
                {
                m(i,j) = r(i,j) * f[j];
                sum += m(i,j);
                }
            }
        m(i,i) = -sum;
        averageRate += f[i] * sum;
        }
        
    // rescale so average is one (Q)
    double scaleFactor = 1.0 / averageRate;
    for (int i=0; i<numStates; i++)
        for (int j=0; j<numStates; j++)
            m(i,j) *= scaleFactor;
                        
#   if 1
    double sum = 0.0;
    std::cout << "Frequencies:" << std::endl;
    for (int i=0; i<f.size(); i++)
        {
        std::cout << "F[" << i << "] = " << f[i] << std::endl;
        sum += f[i];
        }
    std::cout << "Frequencies sum = " << sum << std::endl;
    std::cout << "Rates:" << std::endl;
    sum = 0.0;
    for (int i=0; i<numStates; i++)
        {
        for (int j=0; j<numStates; j++)
            {
            std::cout << r(i,j) << " ";
            sum += r(i,j);
            }
        std::cout << std::endl;
        }
    std::cout << "Rates sum = " << sum << std::endl;
    std::cout << "Average Rates:" << std::endl;
    sum = 0.0;
    for (int i=0; i<numStates; i++)
        {
        for (int j=0; j<numStates; j++)
            {
            if (i != j)
                {
                std::cout << m(i,j) << " ";
                sum += m(i,j);
                }
            }
        std::cout << std::endl;
        }
    std::cout << "Average rates  um = " << sum << std::endl;
#   endif

}

void McmcSummary::calculateFrequencies(std::vector<double>& f) {

    std::set<int> observedStates;
    for (int i=0; i<stats.size(); i++)
        {
        ParameterStatistics* s = stats[i];
        if (s->getName()[0] == 'F')
            {
            int st = getFreqElement(s->getName());
            observedStates.insert(st);
            }
        }
    for (int i=0; i<observedStates.size(); i++) // make certain all of the states from 0 to (numStates-1) are found
        {
        std::set<int>::iterator it = observedStates.find(i);
        if (it == observedStates.end())
            Msg::error("Could not find state " + std::to_string(i) + " in frequencies summary");
        }
    int ns = (int)observedStates.size();
    if (ns != 0 && ns != numStates)
        Msg::error("Inconsistency in the number of states (2)");
    f.resize(numStates);
    for (int i=0, n=(int)f.size(); i<n; i++)
        f[i] = 1.0 / numStates;
    for (int i=0; i<stats.size(); i++)
        {
        ParameterStatistics* s = stats[i];
        if (s->getName()[0] == 'F')
            {
            int st = getFreqElement(s->getName());
            f[st] = s->getMean();
            }
        }

}

double McmcSummary::average(std::vector<double>& x) {

    double sum = 0.0;
    int n = (int)x.size();
    for (int i=0; i<n; i++)
        sum += x[i];
    return sum / n;
}

void McmcSummary::calculateAverageRates(DoubleMatrix& m) {

    int ns = 0;
    
    // 1. get state frequencies (and the number of states). We look
    // up the apropriate ParameterStatistics objects for this
    std::set<int> observedStates;
    for (int i=0; i<stats.size(); i++)
        {
        ParameterStatistics* s = stats[i];
        if (s->getName()[0] == 'F')
            {
            int st = getFreqElement(s->getName());
            observedStates.insert(st);
            }
        }
    for (int i=0; i<observedStates.size(); i++) // make certain all of the states from 0 to (numStates-1) are found
        {
        std::set<int>::iterator it = observedStates.find(i);
        if (it == observedStates.end())
            Msg::error("Could not find state " + std::to_string(i) + " in frequencies summary");
        }
    ns = (int)observedStates.size();
    if (ns != 0 && ns != numStates)
        Msg::error("Inconsistency in the number of states (2)");
    std::vector<double> f(numStates);
    for (int i=0, n=(int)f.size(); i<n; i++)
        f[i] = 1.0 / numStates;
    for (int i=0; i<stats.size(); i++)
        {
        ParameterStatistics* s = stats[i];
        if (s->getName()[0] == 'F')
            {
            int st = getFreqElement(s->getName());
            f[st] = s->getMean();
            }
        }
    if (statePartitions != NULL && ns == numStates) // check consistencey with a state partition object, if we have one
        {
        int tempNumStates = statePartitions->getNumElements();
        if (tempNumStates != numStates)
            Msg::error("Mismatch in the number of states from the frequencies and the state partitions (calculateAverageRates)");
        }
        
    int nr = inferNumberOfRates();

    // allocate a vector to hold the exchangeability parameters
    DoubleMatrix r(numStates,numStates);
    for (int i=0; i<numStates; i++)
        {
        for (int j=0; j<numStates; j++)
            r(i,j) = 1.0;
        r(i,i) = 0.0;
        }
        
    if (nr > 0)
        {
        bool usingCustom = false;
        if (nr > 0 && nr < numStates * (numStates-1)/2)
            usingCustom = true;
        for (int i=0; i<stats.size(); i++)
            {
            ParameterStatistics* s = stats[i];
            if (s->getName()[0] == 'R')
                {
                nr++;
                int r1, r2;
                getRateElements(s->getName(), r1, r2);
                if (usingCustom == true)
                    {
                    // custom model
                    Subset* s1 = statePartitions->findSubsetIndexed(r1);
                    Subset* s2 = statePartitions->findSubsetIndexed(r2);
                    for (int x1 : s1->getValues())
                        {
                        for (int x2 : s2->getValues())
                            {
                            double ave = s->getMean();
                            r(x1,x2) = ave;
                            r(x2,x1) = ave;
                            }
                        }
                    
                    }
                else
                    {
                    // GTR
                    r(r1-1,r2-1) = s->getMean();
                    r(r2-1,r1-1) = s->getMean();
                    }
                }
            }
        }
    
        
    for (int i=0; i<numStates; i++)
        {
        for (int j=0; j<numStates; j++)
            {
            if (r(i,j) < 0.0 && i != j && nr != 0)
                Msg::error("Could not find all rate parameters");
            }
        }

    // fill in the rates (Q)
    //DoubleMatrix m(numStates,numStates);
    double averageRate = 0.0;
    for (int i=0; i<numStates; i++)
        {
        double sum = 0.0;
        for (int j=0; j<numStates; j++)
            {
            if (i != j)
                {
                m(i,j) = r(i,j) * f[j];
                sum += m(i,j);
                }
            }
        m(i,i) = -sum;
        averageRate += f[i] * sum;
        }
        
    // rescale so average is one (Q)
    double scaleFactor = 1.0 / averageRate;
    for (int i=0; i<numStates; i++)
        for (int j=0; j<numStates; j++)
            m(i,j) *= scaleFactor;
            
    // now get average rate from i to j (R)
    for (int i=0; i<numStates; i++)
        for (int j=0; j<numStates; j++)
            m(i,j) *= f[i];
            
#   if 1
    double sum = 0.0;
    std::cout << "Frequencies:" << std::endl;
    for (int i=0; i<f.size(); i++)
        {
        std::cout << "F[" << i << "] = " << f[i] << std::endl;
        sum += f[i];
        }
    std::cout << "Frequencies sum = " << sum << std::endl;
    std::cout << "Rates:" << std::endl;
    sum = 0.0;
    for (int i=0; i<numStates; i++)
        {
        for (int j=0; j<numStates; j++)
            {
            std::cout << r(i,j) << " ";
            sum += r(i,j);
            }
        std::cout << std::endl;
        }
    std::cout << "Rates sum = " << sum << std::endl;
    std::cout << "Average Rates:" << std::endl;
    sum = 0.0;
    for (int i=0; i<numStates; i++)
        {
        for (int j=0; j<numStates; j++)
            {
            if (i != j)
                {
                std::cout << m(i,j) << " ";
                sum += m(i,j);
                }
            }
        std::cout << std::endl;
        }
    std::cout << "Average rates sum = " << sum << std::endl;
#   endif
}

int McmcSummary::getFreqElement(std::string str) {

    int n = 0;
    std::string tempStr = "";
    bool readingBetweenSquareBrackets = false;
    for (int i=0; i<str.length(); i++)
        {
        char c = str[i];
        if (c == '[')
            readingBetweenSquareBrackets = true;
        else if (c == ']')
            readingBetweenSquareBrackets = false;
        else
            {
            if (readingBetweenSquareBrackets == true)
                tempStr += std::string(1, c);
            }
        }
    if (tempStr != "")
        n = std::stoi(tempStr);
    return n;
}

void McmcSummary::getRateElements(std::string str, int& r1, int& r2) {

    r1 = 0;
    r2 = 0;

    std::string tempStr = "";
    bool readingFirstNum = false, readingSecondNum = false;
    for (int i=0; i<str.length(); i++)
        {
        char c = str[i];
        if (c == '[')
            {
            readingFirstNum = true;
            readingSecondNum = false;
            }
        else if (c == ']')
            {
            readingFirstNum = false;
            readingSecondNum = false;
            if (tempStr != "")
                r2 = std::stoi(tempStr);
            tempStr = "";
            }
        else if (c == '-')
            {
            readingFirstNum = false;
            readingSecondNum = true;
            if (tempStr != "")
                r1 = std::stoi(tempStr);
            tempStr = "";
            }
        else
            {
            if (readingFirstNum == true || readingSecondNum == true)
                tempStr += std::string(1, c);
            }
        }
}

std::vector<CredibleInterval> McmcSummary::getCredibleIntervals(void) {

    std::vector<CredibleInterval> cis;
    for (int i=0; i<stats.size(); i++)
        cis.push_back( stats[i]->getCredibleInterval() );
    return cis;
}

std::vector<double> McmcSummary::getMeans(void) {

    std::vector<double> means;
    for (int i=0; i<stats.size(); i++)
        means.push_back( stats[i]->getMean() );
    return means;
}

std::string McmcSummary::getCognateName(std::string str) {

    std::vector<std::string> segs;
    std::string word = "";
    for (int i=0; i<str.length(); i++)
        {
        char c = str[i];
        if (c == '/' || c == '\\' || c == '.')
            {
            if (word != "")
                {
                segs.push_back(word);
                word = "";
                }
            }
        else
            {
            word += c;
            }
        }
    if (word != "")
        segs.push_back(word);
        
//    for (int i=0; i<segs.size(); i++)
//        {
//        std::cout << i << " -- " << segs[i] << std::endl;
//        }
    return segs[segs.size()-2];
}

AlignmentDistribution* McmcSummary::getAlignmentNamed(std::string str) const {

    for (AlignmentDistribution* a : alignments)
        {
        if (a->getName() == str)
            return a;
        }
    return NULL;
}

ParameterStatistics* McmcSummary::getParameterNamed(std::string str) const {

    for (ParameterStatistics* s : stats)
        {
        if (s->getName() == str)
            return s;
        }
    return NULL;
}

bool McmcSummary::hasSemicolon(std::string str) {

    std::size_t found = str.find(';');
    if (found != std::string::npos)
        return true;
    return false;
}

int McmcSummary::inferNumberOfRates(void) {

    int numRateClasses = 0;
    for (int i=0; i<stats.size(); i++)
        {
        ParameterStatistics* s = stats[i];
        if (s->getName()[0] == 'R')
            numRateClasses++;
        }
    std::cout << "stats.size() = " << stats.size() << std::endl;
    std::cout << "numRateClasses = " << numRateClasses << std::endl;
    std::cout << "numStates*(numStates-1)/2 = " << numStates*(numStates-1)/2 << std::endl;
    
    if ( statePartitions != NULL && numRateClasses > 0 && numRateClasses != numStates*(numStates-1)/2 ) // check consistencey with a state partition object, if we have one
        {
        int numSubsets = statePartitions->numSubsets();
        int tempNumRates = numSubsets * (numSubsets-1) / 2;
        std::set<Subset*>& subsets = statePartitions->getSubsets();
        for (Subset* s : subsets)
            {
            if (s->getNumElements() > 1)
                tempNumRates++;
            }
        if (tempNumRates != numRateClasses)
            Msg::error("Mismatch in the number of rate classes");
        }
    return numRateClasses;
}

int McmcSummary::inferNumberOfStates(void) {

    // 1. get state frequencies (and the number of states). We look
    // up the apropriate ParameterStatistics objects for this
    std::set<int> observedStates;
    for (int i=0; i<stats.size(); i++)
        {
        ParameterStatistics* s = stats[i];
        if (s->getName()[0] == 'F')
            {
            int st = getFreqElement(s->getName());
            observedStates.insert(st);
            }
        }
    for (int i=0; i<observedStates.size(); i++) // make certain all of the states from 0 to (numStates-1) are found
        {
        std::set<int>::iterator it = observedStates.find(i);
        if (it == observedStates.end())
            Msg::error("Could not find state " + std::to_string(i) + " in frequencies summary");
        }
    int ns = (int)observedStates.size();
    if (statePartitions != NULL && ns != 0) // check consistencey with a state partition object, if we have one
        {
        int tempNumStates = statePartitions->getNumElements();
        if (tempNumStates != ns)
            Msg::error("Mismatch in the number of states from the frequencies and the state partitions (inferNumberOfStates)");
        }
    return ns;
}

std::map<int,std::string> McmcSummary::interpretTranslateString(std::vector<std::string> translateTokens) {

    std::map<int,std::string> translateMap;
    bool readingKey = true;
    std::string key = "";
    std::string val = "";
    for (int i=0; i<translateTokens.size(); i++)
        {
        std::string token = translateTokens[i];
        if (token == "," || token == ";")
            {
            readingKey = true;
            int intKey = atoi(key.c_str());
            translateMap.insert( std::make_pair(intKey,val) );
            }
        else
            {
            if (readingKey == true)
                {
                key = token;
                readingKey = false;
                }
            else
                val = token;
            }
        //std::cout << i << " " << token << std::endl;
        }
//    for (std::map<int,std::string>::iterator it = translateMap.begin(); it != translateMap.end(); it++)
//        std::cout << it->first << " -> " << it->second << std::endl;
        
    return translateMap;
}

std::string McmcSummary::interpretTreeString(std::string str) {

    std::string ns = "";
    bool reading = false;
    for (int i=0; i<str.length(); i++)
        {
        char c = str[i];
        if (c == '(')
            reading = true;
            
        if (reading == true)
            {
            ns += c;
            }
            
        if (c == ';')
            reading = false;
        }
    
    return ns;
}

void McmcSummary::print(void) {

    // print summary of tsv file
    std::vector<CredibleInterval> cis = getCredibleIntervals();
    std::vector<double> means = getMeans();
    std::cout << std::fixed << std::setprecision(4);
    for (int i=0; i<stats.size(); i++)
        {
        std::cout << stats[i]->getName() << ": ";
        std::cout << cis[i].median << " (";
        std::cout << cis[i].lower << ", " << cis[i].upper << ")";
        std::cout << std::endl;
        }
        
    if (statePartitions != NULL)
        printPartitionFreqs();
        
    // print summary of aln files
    for (int i=0; i<alignments.size(); i++)
        {
        alignments[i]->print();
        }
    
    // print summary of tre file
    int n = 0;
    for (auto it = partitions.begin(); it != partitions.end(); it++)
        {
        if (it->second->size() > n)
            n = it->second->size();
        }
    std::cout << std::fixed << std::setprecision(6);
    for (auto it = partitions.begin(); it != partitions.end(); it++)
        {
        auto ci = it->second->getCredibleInterval();
        std::cout << it->second->getName() << " ";
        std::cout << (double)it->second->size() / n << " ";
        std::cout << it->second->getMean() << " (" << ci.lower << ", " << ci.upper << ")";
        std::cout << std::endl;
        }
    conTree->print();
    std::cout << conTree->getNewick(4) << std::endl;
    
    std::vector<int> dist(1000,0);
    int maxSize = 0;
    std::cout << "raw = c(";
    for (int i=0; i<alignments.size(); i++)
        {
        int ciSize = alignments[i]->ciSize();
        if (ciSize > maxSize)
            maxSize = ciSize;
        dist[ciSize]++;

        std::cout << ciSize;
        if (i + 1 != alignments.size())
            std::cout << ",";
        }
    std::cout << ")" << std::endl;

    for (int i=1; i<maxSize; i++)
        std::cout << i << " -- " << dist[i] << std::endl;
    std::cout << "v = c(";
    for (int i=1; i<maxSize; i++)
        {
        std::cout << i;
        if (i + 1 != maxSize)
            std::cout << ",";
        }
    std::cout << ")" << std::endl;
    
    std::cout << "x = c(";
    for (int i=1; i<maxSize; i++)
        {
        std::cout << dist[i];
        if (i + 1 != maxSize)
            std::cout << ",";
        }
    std::cout << ")" << std::endl;
    
    std::cout << "x <- c(";
    for (int i=0; i<alignments.size(); i++)
        {
        Alignment* map = alignments[i]->getMapAlignment();
        std::map<double,double> gaps = map->gapInfo();
        for (std::map<double,double>::iterator it = gaps.begin(); it != gaps.end(); it++)
            std::cout << it->first << ",";
        }
    std::cout << std::endl;
    std::cout << "y <- c(";
    for (int i=0; i<alignments.size(); i++)
        {
        Alignment* map = alignments[i]->getMapAlignment();
        std::map<double,double> gaps = map->gapInfo();
        for (std::map<double,double>::iterator it = gaps.begin(); it != gaps.end(); it++)
            std::cout << it->second << ",";
        }
    std::cout << std::endl;

}

void McmcSummary::printPartitionSet(void) {

    if (statePartitions != NULL)
        statePartitions->print();
}

void McmcSummary::printPartitionFreqs(void) {

    if (statePartitions != NULL)
        {
        std::vector<double> subsetFreqs(statePartitions->numSubsets(), 0.0);
        for (int i=0; i<stats.size(); i++)
            {
            std::string statName = stats[i]->getName();
             if (statName[0] == 'F')
                {
                statName.erase(remove(statName.begin(), statName.end(), 'F'), statName.end());
                statName.erase(remove(statName.begin(), statName.end(), '['), statName.end());
                statName.erase(remove(statName.begin(), statName.end(), ']'), statName.end());
                int idx = atoi(statName.c_str());
                int pIdx = statePartitions->indexOfSubsetWithValue(idx);
                subsetFreqs[pIdx-1] += stats[i]->getMean();
                //std::cout << statName << " " << idx << " " << pIdx << std::endl;
                }
            }
            
        double maxFreq = 0.0;
        for (int i=0; i<subsetFreqs.size(); i++)
            {
            if (subsetFreqs[i] > maxFreq)
                maxFreq = subsetFreqs[i];
            }
        double maxR = sqrt(maxFreq / 3.14);

        for (int i=0; i<subsetFreqs.size(); i++)
            {
            double r = (sqrt(subsetFreqs[i] / 3.14) / maxR) * 100.0;
            std::cout << "Subset " << i+1 << " freqeuncy = " << subsetFreqs[i] << " " << r << std::endl;
            }
            
            
            
            
        
        double maxVal = 0.0;
        for (int i=0; i<stats.size(); i++)
            {
            std::string statName = stats[i]->getName();
             if (statName[0] == 'R')
                {
                if (stats[i]->getMean() > maxVal)
                    maxVal = stats[i]->getMean();
                }
            }
        std::cout << "maxVal = " << maxVal << std::endl;
        maxR = sqrt(maxVal / 3.14);
        for (int i=0; i<stats.size(); i++)
            {
            std::string statName = stats[i]->getName();
             if (statName[0] == 'R')
                {
                double r = (sqrt(stats[i]->getMean() / 3.14) / maxR) * 36.0;
                std::cout << stats[i]->getName() << " " << stats[i]->getMean() << " " << r << std::endl;
                }
            }
           
            

        }
}

void McmcSummary::writeMatrix(std::ofstream& file, DoubleMatrix &m, std::string name) {
    file << name << " = [\n";
    auto rows = m.getNumRows();
    auto cols = m.getNumCols();
    double max = 0;
    for (int i = 0; i < rows; ++i)
        {
        if (i)
            file << ",\n";
        file << "  [";
        for (int j = 0; j < cols; ++j)
            {
            double r = m(i, j);
            if (r > max)
                max = r;
            if (j)
                file << ",";
            file << r;
            }
        file << "]";
        }
    file << "\n];\n";
    file << name << "Max = " << max << ";\n\n";
}

void McmcSummary::writeConsensusTree(std::ofstream& file) {

    file << conTree->getNewick(4) << ";";
}

void McmcSummary::output(std::ofstream& file, int maxAlignment) {
    double cutoff = 0.95;

    file << "// This file is auto-generated by TongueTwisterReader. Do not edit\n\n";
    file << "AlignmentCutoff = " << cutoff << ";\n\n";

    // create a json object
    auto json = nlohmann::json::object();
    
    // output state partition information
    if (statePartitions != NULL)
        json["state_part"] = statePartitions->toJson();

    // output the consensus tree
    json["consensus"]["tree"]   = conTree->getNewick(4);
    json["consensus"]["n_taxa"] = conTree->numTaxa;
    
    // output information on mean and credible interval for all real-valued parameters
    auto jStats = nlohmann::json::array();
    for (int i=0; i<stats.size(); i++)
        {
        std::cout << stats[i] << std::endl;
        nlohmann::json cogStats = stats[i]->toJson();
        jStats.push_back(cogStats);
        }
    json["stats"] = jStats;
    
    for (int i=0; i< jStats.size(); i++)
        {
        std::string name = jStats[i]["cognate"];
        //std::cout << "Name = " << name << std::endl;
        if (name == "Insertion" || name == "Deletion")
            {
            file << name << "RateMean = " << jStats[i]["mean"] << ";\n";
            file << name << "RateLow = " << jStats[i]["lower"] << ";\n";
            file << name << "RateHigh = " << jStats[i]["upper"] << ";\n\n";
            }
        }

    int numRateClasses = inferNumberOfRates();
    std::cout << "numRateClasses = " << numRateClasses << std::endl;
    if (numRateClasses == 0)
        {
        // have the JC69 model
        }
    else if (numStates * (numStates-1) / 2 == numRateClasses)
        {
        // have the GTR model
        }
    else
        {
        // have the custom model
        int numSubsets = statePartitions->numSubsets();

        DoubleMatrix lo(numSubsets,numSubsets);
        DoubleMatrix md(numSubsets,numSubsets);
        DoubleMatrix hi(numSubsets,numSubsets);
        for (int i=0; i<numSubsets; i++)
            {
            for (int j=0; j<numSubsets; j++)
                {
                lo(i,j) = -1.0;
                md(i,j) = -1.0;
                hi(i,j) = -1.0;
                }
            }
        for (int i=0; i<stats.size(); i++)
            {
            ParameterStatistics* s = stats[i];
            if (s->getName()[0] == 'R')
                {
                int r1 = 0, r2 = 0;
                getRateElements(s->getName(), r1, r2);
                CredibleInterval interval = s->getCredibleInterval();
                double mdVal = s->getMean();
                double loVal = interval.lower;
                double hiVal = interval.upper;
                r1--;
                r2--;
                lo(r1,r2) = loVal;
                lo(r2,r1) = loVal;
                md(r1,r2) = mdVal;
                md(r2,r1) = mdVal;
                hi(r1,r2) = hiVal;
                hi(r2,r1) = hiVal;
                std::cout << s->getName() << " " << r1 << " " << r2 << std::endl;
                }
            }
        for (int i=0; i<numSubsets; i++)
            {
            for (int j=0; j<numSubsets; j++)
                {
                if (md(i,j) < 0.0)
                    Msg::error("Failed to initialize all partition rates");
                }
            }
            
        file << "StatClass[][] TransitionStats = [\n";
        for (int i = 0; i < numSubsets; ++i)
            {
            if (i)
                file << ",\n";
            file << "  [";
            for (int j = 0; j < numSubsets; ++j)
                {
                if (j)
                    file << ",";
                file << "new(" << lo(i,j) << "," << md(i,j) << "," << hi(i,j) << ")";
                }
            file << "]";
            }
        file << "\n];\n\n";
        
        
        // print out the information for the subsets
        std::vector<CredibleInterval> ncq = partitionRates();
        file << "StatClass[][] NCQRates = [\n";
        for (int i = 0; i < numSubsets; ++i)
            {
            if (i)
                file << ",\n";
            file << "  [";
            for (int j = 0; j < numSubsets; ++j)
                {
                if (j)
                    file << ",";
                CredibleInterval& ci = ncq[i * numSubsets + j];
                file << "new(" << ci.lower << "," << ci.median << "," << ci.upper << ")";
                }
            file << "]";
            }
        file << "\n];\n\n";
        
        std::vector<CredibleInterval> ncf = partitionFrequencies();
        file << "StatClass[] NCQFreqs = [";
        for (int i = 0; i < numSubsets; ++i)
            {
            if (i)
                file << ",\n";
            CredibleInterval& ci = ncf[i];
            file << "new(" << ci.lower << "," << ci.median << "," << ci.upper << ")";
            }
        file << "];\n\n";
        }
    
    // output average rates of change
    // No credible intervals on this information.
    int ns = inferNumberOfStates();
    if (numStates != ns && ns != 0)
        Msg::error("Inconsistency in the number of states");
    DoubleMatrix aveRates(numStates,numStates);
    calculateAverageRates(aveRates);
    writeMatrix(file, aveRates, "AverageRates");


    DoubleMatrix q(numStates,numStates);
    calculateRates(q);
    writeMatrix(file, q, "QRates");
        
    int numSubsets = statePartitions->numSubsets();
    DoubleMatrix sq(numSubsets,numSubsets);
    double maxValue = 0.0;
    for (int i=0; i<numStates; i++)
        {
        int iss = statePartitions->indexOfSubsetWithValue(i) - 1;
        for (int j=0; j<numStates; j++)
            {
            if (i != j)
                {
                int jss = statePartitions->indexOfSubsetWithValue(j) - 1;
                std::cout << iss << " -> " << jss << std::endl;
                //double rate = q(i,j);
                double rate = aveRates(i,j);
                sq(iss,jss) += rate;
                if (sq(iss,iss) > maxValue)
                    maxValue = sq(iss,iss);
                }
            }
        }
    // temp, for printing circles
    for (int i=0; i<numSubsets; i++)
        {
        for (int j=0; j<numSubsets; j++)
            {
            double r = sq(i,j)/maxValue;
            double d = sqrt(5184.0 * r);
            std::cout << std::fixed << std::setprecision(0) << d << " ";
            }
        std::cout << std::endl;
        }
    std::cout << std::setprecision(6);

    sq.print();
    
    std::vector<double> f;
    calculateFrequencies(f);
    std::vector<double> catFreqs(numSubsets);
    for (int i=0; i<numStates; i++)
        {
        int iss = statePartitions->indexOfSubsetWithValue(i) - 1;
        catFreqs[iss] += f[i];
        }
    for (int i=0; i<numSubsets; i++)
        std::cout << i << " -- " << catFreqs[i] << std::endl;
    
    // output information on the gappiness of the alignment
    int nt1 = alignments[0]->getMapAlignment()->getNumTaxa()+1;
    std::vector<std::vector<double>> gapSpectrumCnt;
    gapSpectrumCnt.resize(3);
    for (int i=0; i<3; i++)
        {
        gapSpectrumCnt[i].resize(nt1);
        for (int j=0; j<nt1; j++)
            gapSpectrumCnt[i][j] = 0.0;
        }
    for (int n=0; n<alignments.size(); n++)
        {
        AlignmentDistribution* cogAlignments = alignments[n];
        Alignment* aln = cogAlignments->getMapAlignment();
        std::vector<double> spectrum = aln->getGapSpectrum();
        if (spectrum.size() == 2)
            {
            gapSpectrumCnt[0][(int)spectrum[0]]++;
            gapSpectrumCnt[1][(int)spectrum[1]]++;
            }
        else if (spectrum.size() == 3)
            {
            gapSpectrumCnt[0][(int)spectrum[0]]++;
            gapSpectrumCnt[1][(int)(round(spectrum[1]))]++;
            gapSpectrumCnt[2][(int)spectrum[2]]++;
            }
        }
    file << "GapSpectrum = [\n";
    for (int i=0; i<3; i++)
        {
        file << "[";
        for (int j=0; j<nt1; j++)
            {
            if (j != 0)
                file << ",";
            file << gapSpectrumCnt[i][j];
            }
        file << "]";
        if (i != 2)
            file << ",";
        file << "\n";
        }
    file << "];\n\n";

    file << "AlignIndexClass[] Alignments = [\n";

        
    // output the credible set for all alignments
    auto sampledAlignments = nlohmann::json::array();
    for (int i=0; i<alignments.size(); i++)
        {
        auto& a = alignments[i];
        
        if (a->size() > 0)
            {
            auto cogAlns = a->toJson(cutoff, maxAlignment, file);
            sampledAlignments.push_back(cogAlns);
            }
        }
    json["cog_alns"] = sampledAlignments;

    file << "];\n\n";

    // output the pooled partition subset frequencies (if a subset is available)
    if (statePartitions != NULL && ns != 0)
        {
        std::map<int,double> subsetFreqs;
        
        for (int i=0; i<stats.size(); i++)
            {
            std::string n = stats[i]->getName();
            if (n[0] == 'F' && n[1] == '[')
                {
                int x = parseNumberFromFreqHeader(n);
                double val = stats[i]->getMean();
                
                //std::cout << n << " " << x << std::endl;
                Subset* s = statePartitions->findSubsetWithValue(x);
                int subsetIndex = s->getIndex();
                if (s == NULL)
                    Msg::error("Could not find subset with value " + std::to_string(x));
                
                std::map<int,double>::iterator it = subsetFreqs.find(subsetIndex);
                if (it == subsetFreqs.end())
                    {
                    subsetFreqs.insert( std::make_pair(subsetIndex, val) );
                    }
                else
                    {
                    it->second += val;
                    }
                
                }
            }
            double sum = 0.0;
            for (std::map<int,double>::iterator it = subsetFreqs.begin(); it != subsetFreqs.end(); it++)
                {
                std::cout << it->first << " " << it->second << std::endl;
                sum += it->second;
                }
            std::cout << "sum = " << sum << std::endl;
        json["part_freqs"] = statePartitions->toFile(subsetFreqs, file);
        }
        
}
    
int McmcSummary::parseNumberFromFreqHeader(std::string str) {

    bool parse = false;
    std::string numStr = "";
    for (int i=0; i<str.length(); i++)
        {
        char c = str[i];
        if (c == '[')
            parse = true;
        else if (c == ']')
            parse = false;
        if (c != '[' && c != ']' && parse == true)
            numStr += std::string(1, c);
        }
    int num = atoi(numStr.c_str());
    return num;
}

std::vector<CredibleInterval> McmcSummary::partitionFrequencies() {

    // make certain that we have a partition
    std::vector<CredibleInterval> ncf;
    if (statePartitions == nullptr)
        return ncf;

    int numSubsets = statePartitions->numSubsets();
    ncf.resize(numSubsets);

    // set up vectors of vectors that will hold the frequencies
    std::vector<std::vector<double>> aves;
    aves.resize(numSubsets);
    std::vector<int> nums;
    nums.resize(numSubsets);
    for (int i=0; i<numSubsets; i++)
        nums[i] = 0;

    // fill in the frequencies for the natural class categories
    for (int i=0; i<numStates; i++)
        {
        int iss = statePartitions->indexOfSubsetWithValue(i) - 1;

        std::string fStr = "F[" + std::to_string(i) + "]";
        ParameterStatistics* f = getParameterNamed(fStr);
        if (f == NULL)
            Msg::error("Could not find parameter for partition frequencies: " + fStr);
                    
        int nF = f->size();
                    
        std::vector<double>& a = aves[iss];
        if (a.size() == 0)
            {
            a.resize(nF);
            for (int n=0; n<nF; n++)
                a[n] = 0.0;
            }
        for (int n=0; n<nF; n++)
            a[n] += (*f)[n];
        nums[iss]++;
        }

    // average
    for (int i=0; i<numSubsets; i++)
        {
        std::string iLabel = partitionOrder[i];
        Subset* ssi = statePartitions->findSubsetLabeled(iLabel);
        int idxI = ssi->getIndex() - 1;
        
        std::vector<double>& a = aves[idxI];
//        int num = nums[idxI];
//        for (int n=0; n<a.size(); n++)
//            a[n] /= num;
        sort(a.begin(), a.end());
        double asize = (int)a.size();
        ncf[i].lower  = a[(size_t)(asize*0.025)];
        ncf[i].upper  = a[(size_t)(0.975*asize)];
        ncf[i].median = average(a); // yeah, it's not the median
        }

    return ncf;
}

std::vector<CredibleInterval> McmcSummary::partitionRates(void) {

    // make certain that we have a partition
    std::vector<CredibleInterval> ncq;
    if (statePartitions == nullptr)
        return ncq;

    int numSubsets = statePartitions->numSubsets();
    ncq.resize(numSubsets*numSubsets);
    
    // figure out the number of samples and check for consistence
    int numSamples = -1;
    for (int i=0; i<stats.size(); i++)
        {
        if (numSamples == -1)
            numSamples = stats[i]->size();
        if (stats[i]->size() != numSamples)
            Msg::error("Inconsistent number of samples taken");
        }
    std::vector<DoubleMatrix*> Q;
    Q.resize(numSamples);
    for (int i=0; i<numSamples; i++)
        Q[i] = new DoubleMatrix(numStates,numStates);
    double** pi = new double*[numSamples];
    for (int i=0; i<numSamples; i++)
        {
        pi[i] = new double[numStates];
        for (int j=0; j<numStates; j++)
            pi[i][j] = 0.0;
        }
    
    // set up the rate matrix and stationary frequencies
    //    get stationary frequencies
    for (int i=0; i<numStates; i++)
        {
//        int iss = statePartitions->indexOfSubsetWithValue(i) - 1;
        std::string fStr = "F[" + std::to_string(i) + "]";
        ParameterStatistics* f = getParameterNamed(fStr);
        for (int n=0; n<numSamples; n++)
            pi[n][i] = (*f)[n];
        }
        
    //    fill in off diagonals
    for (int i=0; i<numStates; i++)
        {
        int iss = statePartitions->indexOfSubsetWithValue(i) - 1;
        
        for (int j=0; j<numStates; j++)
            {
            if (i != j)
                {
                int jss = statePartitions->indexOfSubsetWithValue(j) - 1;
                int ir = iss;
                int jr = jss;
                if (iss > jss)
                    {
                    ir = jss;
                    jr = iss;
                    }

                std::string rStr = "R[" + std::to_string(ir+1) + "-" + std::to_string(jr+1) + "]";
                ParameterStatistics* r = getParameterNamed(rStr);

                if (r == NULL)
                    Msg::error("Could not find parameter for partition rates: " + rStr);
                    
                int nR = r->size();
                if (nR != numSamples)
                    Msg::error("Unequal length r and f vectors");
                    
                for (int n=0; n<nR; n++)
                    (*Q[n])(i,j) = (*r)[n] * pi[n][j];
                }
            }
        }
        
    //    fill in diagonal elements and rescale
    for (int n=0; n<numSamples; n++)
        {
        double averageRate = 0.0;
        for (int i=0; i<numStates; i++)
            {
            double sum = 0.0;
            for (int j=0; j<numStates; j++)
                {
                if (i != j)
                    sum += (*Q[n])(i,j);
                }
            (*Q[n])(i,i) = -sum;
            averageRate += pi[n][i] * sum;
            }
        double factor = 1.0 / averageRate;
        for (int i=0; i<numStates; i++)
            for (int j=0; j<numStates; j++)
                (*Q[n])(i,j) *= factor;
        }
    
    // set up vectors of vectors that will hold the rates
    std::vector<std::vector<std::vector<double>>> aves;
    aves.resize(numSubsets);
    for (int i=0; i<numSubsets; i++)
        aves[i].resize(numSubsets);
    std::vector<std::vector<int>> nums;
    nums.resize(numSubsets);
    for (int i=0; i<numSubsets; i++)
        {
        nums[i].resize(numSubsets);
        for (int j=0; j<numSubsets; j++)
            nums[i][j] = 0;
        }
        
    // fill in the rates for the natural class categories
#   if 1
    for (int i=0; i<numStates; i++)
        {
        int iss = statePartitions->indexOfSubsetWithValue(i) - 1;

        for (int j=0; j<numStates; j++)
            {
            if (i != j)
                {
                int jss = statePartitions->indexOfSubsetWithValue(j) - 1;
                    
                std::vector<double>& a = aves[iss][jss];
                if (a.size() == 0)
                    {
                    a.resize(numSamples);
                    for (int n=0; n<numSamples; n++)
                        a[n] = 0.0;
                    }
                for (int n=0; n<numSamples; n++)
                    a[n] += (*Q[n])(i,j);
                nums[iss][jss]++;
                }
            }
        }
#   else
    for (int i=0; i<numStates; i++)
        {
        int iss = statePartitions->indexOfSubsetWithValue(i) - 1;

        for (int j=0; j<numStates; j++)
            {
            if (i != j)
                {
                int jss = statePartitions->indexOfSubsetWithValue(j) - 1;

                int ir = iss;
                int jr = jss;
                if (iss > jss)
                    {
                    ir = jss;
                    jr = iss;
                    }
                std::string rStr = "R[" + std::to_string(ir+1) + "-" + std::to_string(jr+1) + "]";
                std::string fStr = "F[" + std::to_string(j) + "]";
                ParameterStatistics* r = getParameterNamed(rStr);
                ParameterStatistics* f = getParameterNamed(fStr);
                if (r == NULL || f == NULL)
                    Msg::error("Could not find parameter for partition rates: " + rStr + " " + fStr);
                    
                int nR = r->size();
                int nF = f->size();
                if (nR != nF)
                    Msg::error("Unequal length r and f vectors");
                    
                std::vector<double>& a = aves[iss][jss];
                if (a.size() == 0)
                    {
                    a.resize(nR);
                    for (int n=0; n<nR; n++)
                        a[n] = 0.0;
                    }
                for (int n=0; n<nR; n++)
                    a[n] += (*r)[n] * (*f)[n];
                nums[iss][jss]++;
                }
            }
        }
#   endif
        
    // average
    for (int i=0; i<numSubsets; i++)
        {
        std::string iLabel = partitionOrder[i];
        Subset* ssi = statePartitions->findSubsetLabeled(iLabel);
        int idxI = ssi->getIndex() - 1;
        
        for (int j=0; j<numSubsets; j++)
            {
            std::string jLabel = partitionOrder[j];
            Subset* ssj = statePartitions->findSubsetLabeled(jLabel);
            int idxJ = ssj->getIndex() - 1;
            
            std::vector<double>& a = aves[idxI][idxJ];
            int num   = nums[idxI][idxJ];
            int asize = (int)a.size();
            for (int n=0; n<asize; n++)
                a[n] /= num;
            sort(a.begin(), a.end());
            int k = i * numSubsets + j;
            ncq[k].lower  = a[(size_t)(0.025*asize)];
            ncq[k].upper  = a[(size_t)(0.975*asize)];
            ncq[k].median = average(a); // yeah, it's not the median
            }
        }
        
    for (int i=0; i<Q.size(); i++)
        delete Q[i];
    for (int i=0; i<numSamples; i++)
        delete [] pi[i];
    delete [] pi;

    return ncq;
}

void McmcSummary::readAlnFile(std::string fn, int bi, Partition* prt) {

    std::ifstream ifs(fn);
    nlohmann::json j;
    try
        {
        j = nlohmann::json::parse(ifs);
        }
    catch (nlohmann::json::parse_error& ex)
        {
        Msg::warning("Error parsing JSON file" + fn + " at byte " + std::to_string(ex.byte));
        return;
        }
        
    std::string cognateName = getCognateName(fn);

    auto it = j.find("PartitionSets");
    hasPartitions = true;
    Partition* part = NULL;
    if (it == j.end())
        {
        hasPartitions = false;
        //Msg::error("Could not find partition set in the JSON file");
        }
    if (hasPartitions == true)
        {
        nlohmann::json jsPart = j["PartitionSets"];
        part = new Partition(jsPart);
        }
    if (part == NULL)
        part = prt;

    AlignmentDistribution* dist = new AlignmentDistribution(rv, part);
    alignments.push_back(dist);
    dist->setName(cognateName);
    

    it = j.find("Samples");
    if (it == j.end())
        Msg::error("Could not find sampled alignments in the JSON file");
    nlohmann::json js = j["Samples"];
    for (int i=0; i<js.size(); i++)
        {
        if (i > bi)
            {
            Alignment* aln = new Alignment( js[i]["Data"], taxa );
            dist->addAlignment(aln);
            }
        }
        
    //dist->print();
}

void McmcSummary::readConfigFile(std::string fn) {

    std::ifstream ifs(fn);
    nlohmann::json j;
    try
        {
        j = nlohmann::json::parse(ifs);
        }
    catch (nlohmann::json::parse_error& ex)
        {
        Msg::warning("Error parsing JSON file" + fn + " at byte " + std::to_string(ex.byte));
        return;
        }

    auto jtaxa = j["Taxa"];
    for (int i = 0; i < jtaxa.size(); i++)
        taxa.push_back(jtaxa[i].get<std::string>());

    auto it = j.find("PartitionSets");
    if (it == j.end())
        {
        Msg::warning("Could not find partition set in the JSON configuration file");
        statePartitions = NULL;
        }
    else
        {
        nlohmann::json jsPart = j["PartitionSets"];
        statePartitions = new Partition(jsPart);
        
        partitionOrder.push_back("Long Vowel");
        partitionOrder.push_back("Nasal Vowel");
        partitionOrder.push_back("Diphthong");
        partitionOrder.push_back("Short Vowel");
        partitionOrder.push_back("Nasal Consonant");
        partitionOrder.push_back("Liquid");
        partitionOrder.push_back("Approximant");
        partitionOrder.push_back("Affricate");
        partitionOrder.push_back("Fricatives");
        partitionOrder.push_back("Stops");
        }
    
    if (statePartitions != NULL)
        statePartitions->print();
}

Partition* McmcSummary::readPartition(std::string fn) {

    std::ifstream ifs(fn);
    nlohmann::json j;
    try
        {
        j = nlohmann::json::parse(ifs);
        }
    catch (nlohmann::json::parse_error& ex)
        {
        Msg::warning("Error parsing JSON file" + fn + " at byte " + std::to_string(ex.byte));
        return NULL;
        }


    auto it = j.find("PartitionSets");
    if (it == j.end())
        {
        Msg::warning("Could not find partition set in the JSON configuration file");
        return NULL;
        }
    else
        {
        nlohmann::json jsPart = j["PartitionSets"];
        return new Partition(jsPart);
        }
}

int McmcSummary::readNumStates(std::string fn) {

    std::ifstream ifs(fn);
    nlohmann::json j;
    try
        {
        j = nlohmann::json::parse(ifs);
        }
    catch (nlohmann::json::parse_error& ex)
        {
        Msg::warning("Error parsing JSON file" + fn + " at byte " + std::to_string(ex.byte));
        return 0;
        }

    auto it = j.find("NumberOfStates");
    if (it == j.end())
        {
        Msg::warning("Could not find partition set in the JSON configuration file, \"" + fn + "\"");
        return 0;
        }
    else
        {
        int ns = j["NumberOfStates"];
        return ns;
        }
}

void McmcSummary::readTreFile(std::string fn, int bi) {

	// open the file
	std::ifstream fstrm(fn.c_str());
	if (!fstrm)
        Msg::error("Cannot open file \"" + fn + "\"");

	std::string lineString = "";
	int line = 0;
    bool readingTranslateTable = false, readingTree = false;
    std::vector<std::string> translateTokens;
    std::string treeString = "";
    std::map<int,std::string> translateMap;
    int treeCount = 0;
	while( getline(fstrm, lineString).good() )
		{
        //std::cout << line << " -- " << lineString << std::endl;
		std::istringstream linestream(lineString);
		int ch;
		std::string word = "";
		int wordNum = 0;
		std::string cmdString = "";
		do
			{
			word = "";
			linestream >> word;
            if (word == "")
                continue;
            if (word == "translate")
                {
                readingTranslateTable = true;
                }
            else if (word == "tree")
                {
                readingTree = true;
                }
            else
                {
                if (readingTranslateTable == true)
                    {
                    std::vector<std::string> brokenWord = breakString(word);
                    if (hasSemicolon(word) == true)
                        {
                        for (int i=0; i<brokenWord.size(); i++)
                            translateTokens.push_back(brokenWord[i]);
                        readingTranslateTable = false;
                        translateMap = interpretTranslateString(translateTokens);
                        }
                    else
                        {
                        for (int i=0; i<brokenWord.size(); i++)
                            translateTokens.push_back(brokenWord[i]);
                        }
                    }
                else if (readingTree == true)
                    {
                    if (hasSemicolon(word) == true)
                        {
                        treeString += word;
                        readingTree = false;
                        treeCount++;
                        if (treeCount > bi)
                            {
                            std::string newickStr = interpretTreeString(treeString);
                            Tree t(newickStr, translateMap);
                            std::map<RbBitSet,double> parts = t.getPartitions();
                            addPartion(parts);
                            }
                        treeString = "";
                        }
                    else
                        treeString += word;
                    }
                }

			wordNum++;
            } while ( (ch=linestream.get()) != EOF );
						
		line++;
		}

	// close the file
	fstrm.close();
 
    conTree = new Tree(partitions, translateMap);
}

void McmcSummary::readTsvFile(std::string fn, int bi) {

	// open the file
	std::ifstream fstrm(fn.c_str());
	if (!fstrm)
        Msg::error("Cannot open file \"" + fn + "\"");

	std::string lineString = "";
	int line = 0;
	while( getline(fstrm, lineString).good() )
		{
        //std::cout << line << " -- " << lineString << std::endl;
		std::istringstream linestream(lineString);
		int ch;
		std::string word = "";
		int wordNum = 0;
		std::string cmdString = "";
		do
			{
			word = "";
			linestream >> word;
            if (word == "")
                continue;
			if (line == 0)
				{
                ParameterStatistics* s = new ParameterStatistics;
                s->setName(word);
                stats.push_back(s);
				}
			else
				{
                if (line > bi)
                    {
                    double x = atof(word.c_str());
                    ParameterStatistics* s = stats[wordNum];
                    s->addValue(x);
                    }
				}
			wordNum++;
            } while ( (ch=linestream.get()) != EOF );
						
		line++;
		}

	
	// close the file
	fstrm.close();
}
