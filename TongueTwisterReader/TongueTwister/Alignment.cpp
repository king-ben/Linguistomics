#include <iostream>
#include <iomanip>
#include <istream>
#include <sstream>
#include <fstream>
#include "Alignment.hpp"
#include "Msg.hpp"



Alignment::Alignment(nlohmann::json& js, int ns, std::vector<std::string> canonicalTaxonList) {
    
    //std::cout << j.dump() << std::endl;
    
    matrix = NULL;
    indelMatrix = NULL;
    
    // set the state information
    numStates = ns;
    gapCode = numStates;
    if (gapCode <= 1)
        Msg::error("There must be at least two states in the model");

    // read the word name
    auto it = js.find("Name");
    if (it == js.end())
        Msg::error("Could not find word name in the JSON object");
    name = js["Name"];
    
    // read the number of languages encoded for this word
    it = js.find("Data");
    if (it == js.end())
        Msg::error("Could not find word data in the JSON object");
    numTaxa = (int)js["Data"].size();
    if (numTaxa <= 0)
        Msg::error("Must have at least one taxon in the word");
            
    // set up the taxon mask
    taxonMask.resize(canonicalTaxonList.size());
    for (int i=0; i<taxonMask.size(); i++)
        taxonMask[i] = false;
        
    // read the json data
    numChar = 0;
    for (int n=0; n<numTaxa; n++)
        {
        nlohmann::json jw = js["Data"][n];
        
        // check that there is taxon name information
        it = jw.find("Taxon");
        if (it == jw.end())
            Msg::error("Could not find taxon name in the JSON object");
        std::string tName = jw["Taxon"];
        
        // find the index of the taxon name and set the mask accordingly
        int taxonIdx = 0;
        bool foundTaxon = false;
        for (int j=0; j<canonicalTaxonList.size(); j++)
            {
            if (tName == canonicalTaxonList[j])
                {
                foundTaxon = true;
                break;
                }
            taxonIdx++;
            }
        if (foundTaxon == false)
            Msg::error("Could not find taxon " + tName + " in list of canonical taxa");
        taxonMask[taxonIdx] = true;
        taxonMap.insert( std::make_pair(taxonIdx,n) );
        taxonNames.push_back(tName);
        
        // check that there is segment information
        it = jw.find("Segments");
        if (it == jw.end())
            Msg::error("Could not find segment data in the JSON object");
        std::vector<int> segInfo = jw["Segments"];
        
        // check that the number of word segments is the same for each taxon
        if (n == 0)
            {
            numChar = (int)segInfo.size();
            if (numChar <= 0)
                Msg::error("Must have at least one segment in the word");

            if (matrix == NULL)
                {
                matrix = new int*[numTaxa];
                matrix[0] = new int[numTaxa * numChar];
                for (size_t i=1; i<numTaxa; i++)
                    matrix[i] = matrix[i-1] + numChar;
                for (size_t i=0; i<numTaxa; i++)
                    for (size_t j=0; j<numChar; j++)
                        matrix[i][j] = 0;
                }
            if (indelMatrix == NULL)
                {
                indelMatrix = new bool*[numTaxa];
                indelMatrix[0] = new bool[numTaxa * numChar];
                for (size_t i=1; i<numTaxa; i++)
                    indelMatrix[i] = indelMatrix[i-1] + numChar;
                for (size_t i=0; i<numTaxa; i++)
                    for (size_t j=0; j<numChar; j++)
                        indelMatrix[i][j] = true;
                }
            }
        else
            {
            if (segInfo.size() != numChar)
                Msg::error("Inconsistent segment lengths for word " + name);
            }
            
        //read the segment information
        for (int j=0; j<segInfo.size(); j++)
            {
            if (segInfo[j] == -1)
                matrix[n][j] = gapCode;
            else
                matrix[n][j] = segInfo[j];
            }
        }
        
    // check if the alignment is complete and if the taxa are in the correct order
    bool isCompletelySampled = true;
    for (int i=0; i<taxonMask.size(); i++)
        {
        if (taxonMask[i] == false)
            isCompletelySampled = false;
        }
    if (isCompletelySampled == true)
        {
        for (int i=0; i<numTaxa; i++)
            {
            if (taxonNames[i] != canonicalTaxonList[i])
                Msg::error("Taxa must be in the same order as in the canonical taxon list");
            }
        }
    else
        {
        
        }
        
    // fill in the indel matrix
    for (int i=0; i<numTaxa; i++)
        {
        for (int j=0; j<numChar; j++)
            {
            if (isIndel(i, j) == true)
                indelMatrix[i][j] = false;
            }
        }
            
    std::cout << "   * Word alignment \"" << name << "\" has " << numCompleteTaxa() << " taxa and " << numChar << " syllables" << std::endl;
}

Alignment::~Alignment(void) {

	delete [] matrix[0];
	delete [] matrix;
    delete [] indelMatrix[0];
    delete [] indelMatrix;
}

std::string Alignment::bomLessString(std::string& str) {

    std::stringstream is(str);
    std::string rStr = "";
    char c;
    while (is.get(c))
        {
        if (c != '\xEF' && c != '\xBB' && c != '\xBF')
            rStr += c;
        }
    return rStr;
}

int Alignment::getCharacter(size_t i, size_t j) {

    return matrix[i][j];
}

std::vector<std::vector<int> > Alignment::getIndelMatrix(void) {

    // note that this returns a numSites X numTaxa matrix (i.e.,
    // the rows contain the information for the i-th site while
    // the columns contain the information the j-th taxon)
    std::vector<std::vector<int> > m;
    m.resize(numChar);
    for (int i=0; i<numChar; i++)
        m[i].resize(numTaxa);
    for (int i=0; i<numChar; i++)
        {
        for (int j=0; j<numTaxa; j++)
            {
            if (isIndel(j, i) == true)
                m[i][j] = 0;
            else
                m[i][j] = 1;
            }
        }
    return m;
}

std::vector<int> Alignment::getRawSequence(int txnIdx) {

    std::vector<int> seq;
    for (int j=0; j<numChar; j++)
        {
        if (isIndel(txnIdx, j) == false)
            {
            seq.push_back( matrix[txnIdx][j] );
            }
        }
    return seq;
}

std::vector<std::vector<int> > Alignment::getRawSequenceMatrix(void) {

    std::vector<std::vector<int> > seqMat;
    
    seqMat.resize(numTaxa);
    for (int i=0; i<numTaxa; i++)
        {
        seqMat[i] = getRawSequence(i);
        }
        
    return seqMat;
}

std::string Alignment::getTaxonMaskString(void) {

    std::string str = "";
    for (size_t i=0; i<taxonMask.size(); i++)
        {
        if (taxonMask[i] == false)
            str += "0";
        else
            str += "1";
        }
    return str;
}

bool Alignment::hasAllGapColumn(void) {

    for (int j=0; j<numChar; j++)
        {
        bool hasNonGapChar = false;
        for (int i=0; i<numTaxa; i++)
            {
            if (matrix[i][j] != gapCode)
                hasNonGapChar = true;
            }
        if (hasNonGapChar == false)
            return true;
        }
    return false;
}

bool Alignment::hasBOM(std::string& str) {

    std::stringstream is(str);

    // read the first byte
    char const c0 = (char)is.get();
    if (c0 != '\xEF')
        {
        is.putback(c0);
        return false;
        }

    // read the second byte
    char const c1 = (char)is.get();
    if (c1 != '\xBB')
        {
        is.putback(c1);
        is.putback(c0);
        return false;
        }

    // peek the third byte
    char const c2 = (char)is.peek();
    if (c2 != '\xBF')
        {
        is.putback(c1);
        is.putback(c0);
        return false;
        }

    return true;
}

bool Alignment::isAlignmentComplete(void) {

    for (int i=0; i<taxonMask.size(); i++)
        {
        if (taxonMask[i] == false)
            return false;
        }
    return true;
}

bool Alignment::isIndel(size_t i, size_t j) {

    if (matrix[i][j] == gapCode)
        return true;
    return false;
}

bool Alignment::isInteger(const std::string& str) {

    return str.find_first_not_of("0123456789") == std::string::npos;
}

void Alignment::listTaxa(void) {

	int i = 1;
	for (std::vector<std::string>::iterator p=taxonNames.begin(); p != taxonNames.end(); p++)
		std::cout << std::setw(4) << i++ << " -- " << (*p) << '\n';
}

std::string Alignment::getTaxonName(int i) {

	return taxonNames[i];
}

int Alignment::getTaxonIndex(std::string ns) {

    for (int i=0; i<taxonNames.size(); i++)
        {
        if (ns == taxonNames[i])
            return i;
        }
    return -1;
}

int Alignment::numCompleteTaxa(void) {

    int n = 0;
    for (int i=0; i<taxonMask.size(); i++)
        {
        if (taxonMask[i] == true)
            n++;
        }
    return n;
}

void Alignment::print(void) {

	int** x = matrix;
		
    std::cout << "Name: " << name << std::endl;
    std::cout << "Taxa: ( ";
    for (int i=0; i<taxonNames.size(); i++)
        std::cout << taxonNames[i] << " ";
    std::cout << ")" << std::endl;
    std::cout << "Name: " << name << std::endl;
	std::cout << "        ";
	for (size_t i=0; i<numTaxa; i++)
		std::cout << std::setw(3) << i;
	std::cout << '\n';
	std::cout << "------------------------";
	for (size_t i=0; i<numTaxa; i++)
		std::cout << "---";
	std::cout << '\n';
	for (size_t j=0; j<numChar; j++)
		{
		std::cout << std::setw(4) << j+1 << " -- ";
		for (size_t i=0; i<numTaxa; i++)
            {
            if (x[i][j] == gapCode)
                std::cout << std::setw(3) << "-";
            else
                std::cout << std::setw(3) << x[i][j];
            }
		std::cout << '\n';
		}
    std::cout << std::endl;
}

void Alignment::printIndels(void) {

    bool** x = indelMatrix;
        
    std::cout << "        ";
    for (size_t i=0; i<numTaxa; i++)
        std::cout << std::setw(3) << i;
    std::cout << '\n';
    std::cout << "------------------------";
    for (size_t i=0; i<numTaxa; i++)
        std::cout << "---";
    std::cout << '\n';
    for (size_t j=0; j<numChar; j++)
        {
        std::cout << std::setw(4) << j+1 << " -- ";
        for (size_t i=0; i<numTaxa; i++)
            {
            std::cout << std::setw(3) << x[i][j];
            }
        std::cout << '\n';
        }
    std::cout << std::endl;
}

