#include <iostream>
#include "Container.hpp"
#include "Exchangeabilities.hpp"
#include "McmcSummary.hpp"
#include "Msg.hpp"
#include "ParameterStatistics.hpp"
#include "Statistics.hpp"



Exchangeabilities::Exchangeabilities(McmcSummary& samples, double bf) : Simplex(bf) {

    // collect all of the frequency samples
    std::vector<ParameterStatistics*> rateStats;
    for (size_t i=0; i<samples.size(); i++)
        {
        ParameterStatistics* stat = samples[i];
        std::string parmName = stat->getName();
        if (parmName[0] == 'R' || parmName[0] == 'N')
            rateStats.push_back(stat);
        }
                
    // move the values to this object
    values.resize(rateStats.size());
    fromTo.resize(rateStats.size());
    dimension = values.size();
    for (size_t n=0; n<rateStats.size(); n++)
        {
        std::pair<int,int> bracketValues = parseBracketedNumbers(rateStats[n]->getName());
        fromTo[n] = bracketValues;

        const std::vector<float> parmValues = rateStats[n]->getValues();
        size_t startIdx = static_cast<size_t>(parmValues.size() * burnFraction);
        size_t numSamples = parmValues.size() - startIdx;
        values[n].resize(numSamples);
        
        for (size_t i=startIdx, k=0; i<parmValues.size(); i++)
            values[n][k++] = parmValues[i];
        }
        
    // check that all of the columns have the same number of entries
    size_t numEntries = values[0].size();
    for (size_t i=1; i<values.size(); i++)
        {
        if (values[i].size() != numEntries)
            Msg::error("Unequal number of rate parameters sampled");
        }
        
    // check that each row sums to 1.0
    for (size_t j=0; j<numEntries; j++)
        {
        double sum = 0.0;
        for (size_t i=0; i<dimension; i++)
            sum += values[i][j];
        if (fabs(1.0-sum) > 0.0001)
            Msg::error("Row sum is not 1.0: " + std::to_string(sum));
        }

    // precalculate the mean, variance, anc CI
    mean.resize(dimension);
    priorMean.resize(dimension);
    variance.resize(dimension);
    lowerCI.resize(dimension);
    upperCI.resize(dimension);
    kl.resize(dimension);
    for (size_t n=0; n<dimension; n++)
        {
        double d = Statistics::calcKlToPrior(values[n], 1, dimension);
        std::pair<double,double> mv = Statistics::getMeanAndVariance(values[n], 0.0);
        CredibleInterval ci = Statistics::getCredibleInterval(values[n], 0.0);
        
        mean[n] = mv.first;
        priorMean[n] = 1.0 / dimension;
        variance[n] = mv.second;
        lowerCI[n] = ci.lower;
        upperCI[n] = ci.upper;
        kl[n] = d;
        }
}

double Exchangeabilities::getAverageRate(int r1, int r2) const {

    if (r1 > r2)
        std::swap(r1, r2);
    for (size_t i=0; i<size(); i++)
        {
        if (fromTo[i].first == r1 && fromTo[i].second == r2)
            return mean[i];
        }
    Msg::error("Could not find average rate for groups");
    return 0.0;
}

void Exchangeabilities::print(void) {

    for (size_t i=0; i<dimension; i++)
        {
        std::cout << std::fixed;
        std::cout << fromTo[i].first << "-" << fromTo[i].second << " -- ";
        std::cout << mean[i] << " (" << lowerCI[i] << ", " << upperCI[i] << ") KL=";
        std::cout << kl[i] << " N=" << values[i].size() << std::endl;
        }
}

nlohmann::json Exchangeabilities::toJson(void) {

    nlohmann::json j = nlohmann::json::object();
    
    for (size_t i=0; i<dimension; i++)
        {
        std::string name = "R[" + std::to_string(i) + "]";
        j["statname"] = name;
        j["mean"] = mean[i];
        j["lower"] = lowerCI[i];
        j["upper"] = upperCI[i];
        }

    return j;
}
