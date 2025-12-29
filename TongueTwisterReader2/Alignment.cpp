#include <iomanip>
#include <iostream>
#include "Alignment.hpp"
#include "Msg.hpp"
#include "Sequence.hpp"

int Alignment::max_int8 = pow(2.0, sizeof(int8_t) * 8 - 1) - 1;


Alignment::Alignment(nlohmann::json js) {

    for (size_t i=0; i<js.size(); i++)
        {
        const std::string& tName = js[i]["Taxon"];
        const std::vector<int>& rawInfo = js[i]["Segments"];
        for (size_t j=0; j<rawInfo.size(); j++)
            {
            if (rawInfo[j] > max_int8)
                Msg::error("Character state value is too large");
            }
        std::vector<int8_t> tInfo(rawInfo.size());
        for (size_t j=0; j<rawInfo.size(); j++)
            tInfo[j] = static_cast<int8_t>(rawInfo[j]);
        matrix.push_back( Sequence(tName, tInfo) );
        }
    numTaxa = (int)matrix.size();
    numChar = (int)matrix[0].size();
}

bool Alignment::operator==(const Alignment& aln) const {

    if (matrix.size() != aln.matrix.size())
        return false;
    for (size_t i=0; i<matrix.size(); i++)
        {
        if (matrix[i].size() != aln.matrix[i].size())
            return false;
        for (int j=0; j<matrix[i].size(); j++)
            {
            if (matrix[i][j] != aln.matrix[i][j])
                return false;
            }
        }
    return true;
}

void Alignment::from_json(const nlohmann::json& /*j*/) {

    //j.at("name").get_to(p.name);
    //j.at("address").get_to(p.address);
    //j.at("age").get_to(p.age);
}

std::map<double,double> Alignment::gapInfo(void) {

    std::map<double,double> info;
    
    for (int j=0; j<numChar; j++)
        {
        int numGaps = 0;
        for (int i=0; i<numTaxa; i++)
            {
            int charCode = matrix[i][j];
            if (charCode == -1)
                numGaps++;
            }
            
        double fracGaps = (double)numGaps / numTaxa;
        double pos = (double)(j+1) / numChar;
        info.insert( std::make_pair(pos,fracGaps) );
        }
    
    return info;
}

std::vector<double> Alignment::getGapSpectrum(void) {

    std::vector<double> spectrum;
    if (numChar == 1)
        {
        spectrum.resize(1);
        int n = 0;
        for (int i=0; i<numTaxa; i++)
            {
            if (matrix[i][0] == -1)
                n++;
            }
        spectrum[0] = n;
        }
    else if (numChar == 2)
        {
        spectrum.resize(2);
        for (int j=0; j<2; j++)
            {
            int n = 0;
            for (int i=0; i<numTaxa; i++)
                {
                if (matrix[i][j] == -1)
                    n++;
                }
            spectrum[j] = n;
            }
        }
    else
        {
        spectrum.resize(3);
        for (int i=0; i<3; i++)
            spectrum[i] = 0.0;
        for (int j=0; j<numChar; j++)
            {
            int n = 0;
            for (int i=0; i<numTaxa; i++)
                {
                if (matrix[i][j] == -1)
                    n++;
                }
            if (j == 0)
                spectrum[0] = n;
            else if (j == numChar-1)
                spectrum[2] = n;
            else
                spectrum[1] += n;
            }
        spectrum[1] /= (numChar-2);
        }
    
    return spectrum;
}

int Alignment::lengthOfLongestName(void) {

    int maxLen = 0;
    for (int i=0; i<matrix.size(); i++)
        {
        std::string n = matrix[i].getName();
        if (n.length() > maxLen)
            maxLen = (int)n.length();
        }
    return maxLen;
}

void Alignment::print(std::string h) {

    std::cout << h << std::endl;
    print();
}

void Alignment::print(void) {

    int maxLen = lengthOfLongestName();        
    for (int i=0; i<matrix.size(); i++)
        {
        std::string n = matrix[i].getName();
        std::cout << n << " ";
        for (int j=0; j<maxLen-n.length(); j++)
            std::cout << " ";
        for (int j=0; j<matrix[i].size(); j++)
            {
            int charCode = matrix[i][j];
            if (charCode == -1)
                std::cout << " - ";
            else
                std::cout << std::setw(2) << charCode << " ";
            }
        std::cout << std::endl;
        }
}
