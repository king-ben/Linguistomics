#include "IndelRates.hpp"
#include "McmcSummary.hpp"
#include "Msg.hpp"
#include "ParameterStatistics.hpp"



IndelRates::IndelRates(McmcSummary& samples, double bf) : burnFraction(bf) {

    // find lambda and mu
    ParameterStatistics* lambdaStat = nullptr;
    ParameterStatistics* muStat = nullptr;
    std::vector<ParameterStatistics*> freqStats;
    for (size_t i=0; i<samples.size(); i++)
        {
        ParameterStatistics* stat = samples[i];
        std::string parmName = stat->getName();
        if (parmName == "Lambda")
            lambdaStat = stat;
        if (parmName == "Mu")
            muStat = stat;
        }
    if (lambdaStat == nullptr || muStat == nullptr)
        Msg::error("Could not find insertion and/or deletion rates");
        
    // check that both have the same number of values
    if (lambdaStat->size() != muStat->size())
        Msg::error("Unequal number of insertion and deletion samples made");

    // move the values to this object
    size_t startIdx = static_cast<size_t>(lambdaStat->size() * burnFraction);
    const std::vector<float>& lambdaValues = lambdaStat->getValues();
    const std::vector<float>& muValues = muStat->getValues();
    for (size_t i=startIdx; i<lambdaStat->size(); i++)
        {
        std::pair<float,float> x = std::make_pair(lambdaValues[i],muValues[i]);
        values.push_back(x);
        }
        
    // precalculate the mean, variance, anc CI
    std::vector<float> x(values.size());
    std::vector<float> y(values.size());
    for (size_t n=0; n<values.size(); n++)
        {
        x[n] = values[n].first;
        y[n] = values[n].second;
        }
    std::pair<double,double> mvX = Statistics::getMeanAndVariance(x, 0.0);
    std::pair<double,double> mvY = Statistics::getMeanAndVariance(y, 0.0);
    CredibleInterval ciX = Statistics::getCredibleInterval(x, 0.0);
    CredibleInterval ciY = Statistics::getCredibleInterval(y, 0.0);
    mean.first = mvX.first;
    mean.second = mvY.first;
    variance.first = mvX.second;
    variance.second = mvY.second;
    lowerCI.first = ciX.lower;
    lowerCI.second = ciY.lower;
    upperCI.first = ciX.upper;
    upperCI.second = ciY.upper;
}

nlohmann::json IndelRates::toJson(void) {

    nlohmann::json j = nlohmann::json::object();
    
    j["statname"] = "Lambda";
    j["mean"] = mean.first;
    j["lower"] = lowerCI.first;
    j["upper"] = upperCI.first;

    j["statname"] = "Mu";
    j["mean"] = mean.second;
    j["lower"] = lowerCI.second;
    j["upper"] = upperCI.second;

    return j;
}
