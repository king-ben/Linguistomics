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
        std::pair<double,double> mv = Statistics::getMeanAndVariance(values[n], burnFraction);
        CredibleInterval ci = Statistics::getCredibleInterval(values[n], burnFraction);
        
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

void StateFrequencies::print(Partition* part) {

    if (part == nullptr)
        print();

    std::cout << "   * Natural Class Frequencies" << std::endl;
    
    // fill in the new samples
    size_t numSamples = values[0].size();
    std::vector<std::vector<float>> f;
    f.resize(part->numSubsets());
    for (size_t i=0; i<f.size(); i++)
        {
        f[i].resize(numSamples);
        for (size_t j=0; j<numSamples; j++)
            f[i][j] = 0.0;
        }
    for (size_t n=0; n<numSamples; n++)
        {
        for (int i=0; i<dimension; i++)
            {
            int subsetIdx = part->indexOfSubsetWithValue(i);
            f[subsetIdx-1][n] += values[i][n];
            }
        }
    
    size_t len = part->longestLabelLength();
    for (int i=0; i<f.size(); i++)
        {
        Subset* ss = part->findSubsetIndexed(i+1);
        double d = Statistics::calcKlToPrior(f[i], ss->getNumElements(), dimension);
        std::pair<double,double> mv = Statistics::getMeanAndVariance(f[i], burnFraction);
        CredibleInterval ci = Statistics::getCredibleInterval(f[i], burnFraction);

        std::cout << "      " << std::setw(2) << i+1 << " " << ss->getLabel();
        for (size_t j=0; j<len-ss->getLabel().length(); j++)
            std::cout << " ";
        std::cout << " -- ";
        std::cout << std::fixed << std::setprecision(6);
        std::cout << (double)ss->getNumElements() / dimension << " -> " << mv.first << " (" << ci.lower << ", " << ci.upper << ") KL=";
        std::cout << d << " N=" << f[i].size() << std::endl;
        }
    std::cout << std::endl;
}
