#include <cmath>
#include <iomanip>
#include <iostream>
#include "Alignment.hpp"
#include "Msg.hpp"
#include "Sequence.hpp"

int Alignment::max_int8 = static_cast<int>(std::pow(2.0, sizeof(int8_t) * 8 - 1) - 1);



Alignment::Alignment(nlohmann::json js) {

    for (size_t i = 0; i < js.size(); i++)
        {
        const std::string& tName = js[i]["Taxon"];
        const std::vector<int>& rawInfo = js[i]["Segments"];
        
        std::vector<int8_t> tInfo(rawInfo.size());
        for (size_t j = 0; j < rawInfo.size(); j++)
            {
            if (rawInfo[j] > max_int8)
                Msg::error("Character state value is too large");
            tInfo[j] = static_cast<int8_t>(rawInfo[j]);
            }
        matrix.emplace_back(tName, std::move(tInfo));
        }
    numTaxa = static_cast<int>(matrix.size());
    numChar = static_cast<int>(matrix[0].size());
}

bool Alignment::operator==(const Alignment& aln) const {

    if (matrix.size() != aln.matrix.size())
        return false;
    
    for (size_t i = 0; i < matrix.size(); i++)
        {
        if (matrix[i].size() != aln.matrix[i].size())
            return false;
        for (size_t j = 0; j < matrix[i].size(); j++)
            {
            if (matrix[i][j] != aln.matrix[i][j])
                return false;
            }
        }
    return true;
}

std::map<double,double> Alignment::gapInfo(void) const {

    std::map<double,double> info;
    
    for (int j = 0; j < numChar; j++)
        {
        int numGaps = 0;
        for (int i = 0; i < numTaxa; i++)
            {
            int charCode = matrix[i][j];
            if (charCode == -1)
                numGaps++;
            }
            
        double fracGaps = static_cast<double>(numGaps) / numTaxa;
        double pos = static_cast<double>(j + 1) / numChar;
        info.insert(std::make_pair(pos, fracGaps));
        }
    
    return info;
}

std::vector<double> Alignment::getGapSpectrum(void) const {

    std::vector<double> spectrum;
    
    if (numChar == 1)
        {
        spectrum.resize(1);
        int n = 0;
        for (int i = 0; i < numTaxa; i++)
            {
            if (matrix[i][0] == -1)
                n++;
            }
        spectrum[0] = n;
        }
    else if (numChar == 2)
        {
        spectrum.resize(2);
        for (int j = 0; j < 2; j++)
            {
            int n = 0;
            for (int i = 0; i < numTaxa; i++)
                {
                if (matrix[i][j] == -1)
                    n++;
                }
            spectrum[j] = n;
            }
        }
    else
        {
        spectrum.resize(3, 0.0);
        for (int j = 0; j < numChar; j++)
            {
            int n = 0;
            for (int i = 0; i < numTaxa; i++)
                {
                if (matrix[i][j] == -1)
                    n++;
                }
            if (j == 0)
                spectrum[0] = n;
            else if (j == numChar - 1)
                spectrum[2] = n;
            else
                spectrum[1] += n;
            }
        spectrum[1] /= (numChar - 2);
        }
    
    return spectrum;
}

int Alignment::lengthOfLongestName(void) const {

    size_t maxLen = 0;
    for (size_t i = 0; i < matrix.size(); i++)
        {
        const std::string& n = matrix[i].getName();
        if (n.length() > maxLen)
            maxLen = n.length();
        }
    return static_cast<int>(maxLen);
}

void Alignment::print(std::string h) const {

    std::cout << h << std::endl;
    print();
}

void Alignment::print(void) const {

    int maxLen = lengthOfLongestName();
    for (size_t i = 0; i < matrix.size(); i++)
        {
        const std::string& n = matrix[i].getName();
        std::cout << n << " ";
        for (size_t j = 0; j < maxLen - n.length(); j++)
            std::cout << " ";
        for (size_t j = 0; j < matrix[i].size(); j++)
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

nlohmann::json Alignment::toFile(std::ostream& file) const {

    nlohmann::json j = nlohmann::json::array();

    // json
    for (int i = 0; i < numTaxa; i++)
        {
        const Sequence& m = matrix[i];
        nlohmann::json sj;
        sj["language"] = m.getName();
        sj["sequence"] = m.getSequence();
        j.push_back(sj);
        }

    // file output
    bool com = false;
    for (int i = 0; i < numTaxa; i++)
        {
        const Sequence& m = matrix[i];
        if (com)
            file << ",";

        file << "[";
        bool comma = false;
        for (auto si : m.getSequence())
            {
            if (comma)
                file << ",";
            file << static_cast<int>(si);
            comma = true;
            }
        file << "]";

        com = true;
        }

    return j;
}
