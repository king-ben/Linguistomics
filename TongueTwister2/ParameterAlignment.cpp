#include "Alignment.hpp"
#include "BitMask.hpp"
#include "Model.hpp"
#include "Msg.hpp"
#include "ParameterAlignment.hpp"
#include "ParameterIndelRates.hpp"



ParameterAlignment::ParameterAlignment(Model* m, RandomVariable* r, std::string n) : Parameter(m, r, n), indelRates(nullptr) {

    indelRates = modelPtr->findParameter<ParameterIndelRates>();
    if (indelRates == nullptr)
        Msg::error("Could not find insertion/deletion parameter when instantiating alignment");
}

ParameterAlignment::~ParameterAlignment(void) {

    delete alignment[0];
    delete alignment[1];
}

size_t ParameterAlignment::calculateMaximumAlignmentLength(nlohmann::json j) {

    int sum = 0;
    for (size_t i=0; i<j["Data"].size(); i++)
        {
        nlohmann::json jw = j["Data"][i];
                        
        // check that there is segment information
        auto it = jw.find("Segments");
        if (it == jw.end())
            Msg::error("Could not find segment data in the JSON object");
        
        // get the segment information and increment the total number of non-gap states
        // that are observed for this cognate
        std::vector<int> segInfo = jw["Segments"];
        for (size_t j=0; j<segInfo.size(); j++)
            {
            if (segInfo[j] != -1)
                sum++;
            }
        }
    return sum;
}

size_t ParameterAlignment::getNumTaxa(void) {

    return alignment[0]->getNumTaxa();
}

bool ParameterAlignment::initializeFromJson(nlohmann::json jsonObj, size_t ns, std::vector<std::string> ctl) {

    // check that all the necessary keys are in the JSON object
    auto it = jsonObj.find("Name");
    if (it == jsonObj.end())
        Msg::error("Could not find word \"Name\" in the JSON object");
    it =jsonObj.find("Data");
    if (it == jsonObj.end())
        Msg::error("Could not find word \"Data\" in the JSON object");

    // set the state information
    numStates = ns;
    gapCode = numStates;
    if (gapCode <= 1)
        Msg::error("There must be at least two states in the model");
        
    // set the canonical taxon list
    canonicalTaxonList = ctl;

    // read the word name
    std::string name = jsonObj["Name"];
    
    // read the number of languages encoded for this word
    size_t numTaxa = static_cast<int>(jsonObj["Data"].size());
    if (numTaxa <= 0)
        Msg::error("Must have at least one taxon in the word");
                    
    // calculate the maximum alignment length
    size_t maxAlignmentLength = calculateMaximumAlignmentLength(jsonObj);
    size_t allocatedLength = maxAlignmentLength + 1;
    
    // allocate the alignment (allocate the copy after this one is filled out
    alignment[0] = new Alignment(numTaxa, allocatedLength);

    // read the json data
    size_t numChar = 0;
    //taxonMask.resize(canonicalTaxonList.size());
    taxonMask = 0;
    for (size_t n=0; n<numTaxa; n++)
        {
        nlohmann::json jw = jsonObj["Data"][n];
        
        // check that there is taxon name information
        it = jw.find("Taxon");
        if (it == jw.end())
            Msg::error("Could not find taxon name in the JSON object");
        std::string tName = jw["Taxon"];
        
        // find the index of the taxon name and set the mask accordingly
        int taxonIdx = 0;
        bool foundTaxon = false;
        for (size_t j=0; j<canonicalTaxonList.size(); j++)
            {
            if (tName == canonicalTaxonList[j])
                {
                foundTaxon = true;
                break;
                }
            taxonIdx++;
            }
        if (foundTaxon == false)
            Msg::error("Could not find taxon \"" + tName + "\" in list of canonical taxa");
        //taxonMask.set(taxonIdx); // sets the bit at taxonIdx to be true 
        BitMask::setBit(taxonMask, taxonIdx); // sets the bit at taxonIdx to be true 
        taxonNames.push_back(tName);
        
        // check that there is segment information
        it = jw.find("Segments");
        if (it == jw.end())
            Msg::error("Could not find segment data in the JSON object");
        std::vector<int> segInfo = jw["Segments"];
        
        // check that the number of word segments is the same for each taxon
        if (n == 0)
            {
            numChar = static_cast<int>(segInfo.size());
            if (numChar <= 0)
                Msg::error("Must have at least one segment in the word");

            }
        else
            {
            if (segInfo.size() != numChar)
                Msg::error("Inconsistent segment lengths for word " + name);
            }
            
        // read the segment information
        for (size_t j=0; j<segInfo.size(); j++)
            {
            if (segInfo[j] == -1)
                (*alignment[0])(n,j) = static_cast<int>(gapCode);
            else
                (*alignment[0])(n,j) = segInfo[j];
            }
        }
    alignment[0]->setNumSites(numChar);
    
    // check if the alignment is complete and if the taxa are in the correct order
    bool isCompletelySampled = true;
    for (size_t i=0; i<canonicalTaxonList.size(); i++)
        {
        if (BitMask::isSet(taxonMask, i) == false)
            isCompletelySampled = false;
        }
    if (isCompletelySampled == true)
        {
        for (size_t i=0; i<numTaxa; i++)
            {
            if (taxonNames[i] != canonicalTaxonList[i])
                Msg::error("Taxa must be in the same order as in the canonical taxon list");
            }
        }
    else
        {
        
        }
        
    // put in the sentinel values at the end of each row
    for (size_t i=0; i<numTaxa; i++)
        (*alignment[0])(i,numChar) = -1;
        
    // make the scratch alignment as a copy
    alignment[1] = new Alignment(*alignment[0]);
    
    return true;
}

void ParameterAlignment::keep(void) {

    (*alignment[1]) = (*alignment[0]);
}

double ParameterAlignment::lnPriorProbability(void) {

    int len = static_cast<int>(alignment[0]->getNumSites());
    double lambda = indelRates->getInsertionRate();
    double mu = indelRates->getDeletionRate();
    double prob = lambda / mu;
    double lnP = (len-1) * log(prob) + log(1.0-prob);
    return lnP;
}

size_t ParameterAlignment::longestNameLength(void) {

    size_t len = 0;
    for (size_t i=0; i<taxonNames.size(); i++)
        {
        if (taxonNames[i].length() > len)
            len = taxonNames[i].length();
        }
    return len;
}

void ParameterAlignment::print(void) {

    std::cout << "   * " << name << std::endl;
    for (int idx=0; idx<1; idx++)
        {
        std::cout << "     Address:         " << alignment[idx] << std::endl;
        std::cout << "     Num. Taxa:       " << alignment[idx]->getNumTaxa() << std::endl;
        std::cout << "     Num. Sites:      " << alignment[idx]->getNumSites() << std::endl;
        std::cout << "     Max. Num. Sites: " << alignment[idx]->getMaximumNumberOfSites()-1 << std::endl;
        size_t longLen = longestNameLength();
        for (size_t i=0; i<alignment[idx]->getNumTaxa(); i++)
            {
            std::cout << "     " << taxonNames[i] << " ";
            for (size_t j=0; j<longLen-taxonNames[i].length(); j++)
                std::cout << " ";
            for (size_t j=0; j<alignment[idx]->getNumSites(); j++)
                {
                if (static_cast<size_t>((*alignment[idx])(i,j)) == numStates)
                    std::cout << std::setw(3) << "-";
                else
                    std::cout << std::setw(3) << (*alignment[idx])(i,j);
                }
            std::cout << std::endl;
            }
        }
}

void ParameterAlignment::restore(void) {

    (*alignment[0]) = (*alignment[1]);
}
