#include <iomanip>
#include <iostream>
#include "Update.hpp"
#include "UpdateStatistics.hpp"



void UpdateStatistics::accept(Update* update) {

    uint64_t id = update->getUpdateId();
    UpdateMap::iterator it = info.find(id);
    if (it == info.end())
        {
        UpdateInfo x;
        x.updateIdx = update->getUpdateIdx();
        x.name = update->getUpdateName();
        x.numProposals = 1;
        x.numAcceptances = 1;
        x.updateType = update->getUpdateType();
        x.tuningParm = update->getTuningParameter();
        x.updatePtr = update;
        info.insert( std::make_pair(id,x) );
        }
    else 
        {
        it->second.numProposals++;
        it->second.numAcceptances++;
        }
}

void UpdateStatistics::print(void) {

    // get the length of the longest name
    size_t longestNameLength = 0;
    for (auto [key,val] : info)
        {
        if (val.name.length() > longestNameLength)
            longestNameLength = val.name.length();
        }
        
    // sort by the value.name
    std::vector<std::pair<uint64_t, UpdateInfo>> sortedUpdates;
    sortedUpdates.reserve(info.size());
    for (auto& kv : info)
        sortedUpdates.push_back(kv);
    std::sort(sortedUpdates.begin(), sortedUpdates.end(),[](const std::pair<uint64_t, UpdateInfo>& a, const std::pair<uint64_t, UpdateInfo>& b) {

        return a.second.name < b.second.name;
    });

    std::cout << "   MCMC Update Summary" << std::endl;
    for (auto& kv : sortedUpdates) 
        {
        UpdateInfo& u = kv.second;

        std::cout << "   * ";
        std::cout << u.name << " ";
        for (size_t i=0; i<longestNameLength-u.name.length(); i++)
            std::cout << " ";
        std::cout << std::setw(4) << u.numProposals << " ";
        std::cout << std::setw(4) << u.numAcceptances << " ";
        if (u.numProposals > 0)
            std::cout << std::fixed << std::setprecision(1) << std::setw(5) << (double)(u.numAcceptances * 100.0) / u.numProposals << "% ";
        else 
            std::cout << std::setw(5) << "N/A" << " ";
        std::cout << std::fixed << std::setprecision(3) << "(" << u.tuningParm << ") ";
        std::cout << std::endl;
        }
    std::cout << std::endl;
}

void UpdateStatistics::reject(Update* update) {

    uint64_t id = update->getUpdateId();
    UpdateMap::iterator it = info.find(id);
    if (it == info.end())
        {
        UpdateInfo x;
        x.updateIdx = update->getUpdateIdx();
        x.name = update->getUpdateName();
        x.numProposals = 1;
        x.numAcceptances = 0;
        x.updateType = update->getUpdateType();
        x.tuningParm = update->getTuningParameter();
        x.updatePtr = update;
        info.insert( std::make_pair(id,x) );
        }
    else 
        {
        it->second.numProposals++;
        }
}

void UpdateStatistics::tune(void) {

    for (auto& [key,val] : info)
        {
        if (val.numProposals >= 200)
            {
            std::cout << std::setprecision(6);
            double acceptanceProb = static_cast<double>(val.numAcceptances) / val.numProposals;
            if (val.updateType == probability)
                val.tuningParm = tuneProbability(acceptanceProb, val.tuningParm);
            else if (val.updateType == window)
                val.tuningParm = tuneWindow(acceptanceProb, val.tuningParm);
            else if (val.updateType == simplex)
                val.tuningParm = tuneSimplex(acceptanceProb, val.tuningParm);
            else if (val.updateType == factor)
                val.tuningParm = tuneFactor(acceptanceProb, val.tuningParm);  
                
            val.updatePtr->setUpdateTuningParameter(val.updateIdx, val.tuningParm);
            val.numProposals = 0;
            val.numAcceptances = 0;
            }
        }
}

double UpdateStatistics::tuneSimplex(double acceptanceProb, double currentTuningValue) {

    double targetProb = 0.22;  
      
    double newTuningValue = currentTuningValue;
    if (acceptanceProb > targetProb)
        newTuningValue /= (1.0 + ((acceptanceProb-targetProb)/(1.0 - targetProb)) );
    else
        newTuningValue *= (2.0 - acceptanceProb/targetProb);
    
    if (newTuningValue < 4.0)
        newTuningValue = 4.0;
    else if (newTuningValue > 40000.0)
        newTuningValue = 40000.0;
        
    return newTuningValue;
}

double UpdateStatistics::tuneFactor(double acceptanceProb, double currentTuningValue) {

    double targetProb = 0.44;
    
    double newTuningValue = currentTuningValue;
    if (acceptanceProb > targetProb)
        newTuningValue *= (1.0 + ((acceptanceProb-targetProb)/(1.0 - targetProb)) );
    else
        newTuningValue /= (2.0 - acceptanceProb/targetProb);

    if (newTuningValue < log(1.01))
        newTuningValue = log(1.01);
    else if (newTuningValue > log(1000.0))
        newTuningValue = log(1000.0);
        
    return newTuningValue;
}

double UpdateStatistics::tuneWindow(double acceptanceProb, double currentTuningValue) {

    double targetProb = 0.44;
    
    double newTuningValue = currentTuningValue;
    if (acceptanceProb > targetProb)
        newTuningValue *= (1.0 + ((acceptanceProb-targetProb)/(1.0 - targetProb)) );
    else
        newTuningValue /= (2.0 - acceptanceProb/targetProb);

    if (newTuningValue < 0.0001)
        newTuningValue = 0.0001;
    else if (newTuningValue > 10.0)
        newTuningValue = 10.0;
        
    return newTuningValue;
}

double UpdateStatistics::tuneProbability(double acceptanceProb, double currentTuningValue) {

    double targetProb = 0.22;
    
    double newTuningValue = currentTuningValue;
    if (acceptanceProb > targetProb)
        newTuningValue *= (1.0 + ((acceptanceProb-targetProb)/(1.0 - targetProb)) );
    else
        newTuningValue /= (2.0 - acceptanceProb/targetProb);

    if (newTuningValue < 0.0001)
        newTuningValue = 0.0001;
    else if (newTuningValue > 0.9999)
        newTuningValue = 0.9999;
        
    return newTuningValue;
}


void UpdateStatistics::zeroOut(void) {

    for (auto& [key,val] : info)
        {
        val.numProposals = 0;
        val.numAcceptances = 0;
        }
}
