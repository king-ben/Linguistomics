#include <iomanip>
#include <iostream>
#include "McmcSummary.hpp"
#include "Msg.hpp"
#include "ParameterStatistics.hpp"
#include "Partition.hpp"
#include "StateFrequencies.hpp"
#include "Statistics.hpp"
#include "Subset.hpp"



StateFrequencies::StateFrequencies(McmcSummary& samples, double bf) : Simplex(bf) {

    // collect all of the frequency samples
    std::vector<ParameterStatistics*> freqStats;
    for (size_t i=0; i<samples.size(); i++)
        {
        ParameterStatistics* stat = samples[i];
        std::string parmName = stat->getName();
        if (parmName[0] == 'F')
            freqStats.push_back(stat);
        }
                
    // move the values to this object
    values.resize(freqStats.size());
    stateValue.resize(freqStats.size());
    dimension = values.size();
    for (size_t n=0; n<freqStats.size(); n++)
        {
        std::pair<int,int> bracketValues = parseBracketedNumbers(freqStats[n]->getName());
        stateValue[n] = bracketValues.first;
        const std::vector<float> parmValues = freqStats[n]->getValues();
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
            Msg::error("Unequal number of frequency parameters sampled");
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

StateFrequencies::StateFrequencies(const StateFrequencies& f, Partition* part) : Simplex(0.0) {

    const SampleVector& fValues = f.getValues();

    dimension = part->numSubsets();
    size_t numSamples = f.numSamples();

    values.resize(dimension);
    for (size_t i=0; i<dimension; i++)
        values[i].resize(numSamples);
    for (size_t i=0; i<dimension; i++)
        for (size_t j=0; j<numSamples; j++)
            values[i][j] = 0.0;
 
    stateValue.resize(dimension);
    for (int i=0; i<dimension; i++)
        stateValue[i] = i+1;
        
    for (int i=0; i<part->numSubsets(); i++)
        {
        Subset* ss = part->findSubsetIndexed(i+1);
        std::set<int>& states = ss->getValues();
        for (int idx : states)
            {
            for (size_t n=0; n<numSamples; n++)
                values[i][n] += fValues[idx][n];
            }
        }

    // check that each row sums to 1.0
    for (size_t j=0; j<numSamples; j++)
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

void StateFrequencies::print(void) {

    std::cout << "   * State Frequencies" << std::endl;
    for (size_t i=0; i<dimension; i++)
        {
        std::cout << "      ";
        std::cout << std::fixed;
        std::cout << stateValue[i] << " -- ";
        std::cout << mean[i] << " (" << lowerCI[i] << ", " << upperCI[i] << ") KL=";
        std::cout << kl[i] << " N=" << values[i].size() << std::endl;
        }
    std::cout << std::endl;
}

nlohmann::json StateFrequencies::toJson(void) {

    nlohmann::json j = nlohmann::json::object();
    
    for (size_t i=0; i<dimension; i++)
        {
        std::string name = "F[" + std::to_string(i) + "]";
        j["statname"] = name;
        j["mean"] = mean[i];
        j["lower"] = lowerCI[i];
        j["upper"] = upperCI[i];
        }

    return j;
}
