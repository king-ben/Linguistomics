#include <algorithm>
#include <cmath>
#include <limits>
#include <regex>
#include <utility>
#include <vector>
#include "Msg.hpp"
#include "Simplex.hpp"



Simplex::Simplex(double bf) : burnFraction(bf) {

}

std::pair<int,int> Simplex::parseBracketedNumbers(const std::string& s) {

    static const std::regex pattern(R"(\[(\d+)(?:-(\d+))?\])");
    std::smatch match;

    if (!std::regex_search(s, match, pattern))
        Msg::error("No bracketed number found");

    int first = std::stoi(match[1]);
    int second = match[2].matched ? std::stoi(match[2]) : -1;

    return { first, second };
}

std::map<float,StatInfo> Simplex::sortByKL(void) {

    std::map<float,StatInfo> res;

    for (size_t n=0; n<values.size(); n++)
        {
        double kl = getKullbackLiebler(n);
        StatInfo val;
        val.index1 = getFirstIndex(n);
        val.index2 = getSecondIndex(n);
        val.mean = getMean(n);
        val.priorMean = getPriorMean(n);
        val.variance = getVariance(n);
        val.lowerCI = getLowerCI(n);
        val.upperCI = getUpperCI(n);
        res.insert( std::make_pair(kl,val) );
//    int             index1;
//    int             index2;
//    std::string     name;
        
        }
        
    return res;
}

void Simplex::valuesAtIndex(size_t idx, std::vector<float>& vec) {

    if (idx >= values[0].size())
        Msg::error("Trying to access an out-of-range value");
    if (dimension != values.size())
        Msg::error("Dimensions issue");
    for (size_t i=0; i<dimension; i++)
        {
        vec[i] = values[i][idx];
        }
}




