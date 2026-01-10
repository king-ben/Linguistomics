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
        x.name = update->getUpdateName();
        x.numProposals = 1;
        x.numAcceptances = 1;
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
        x.name = update->getUpdateName();
        x.numProposals = 1;
        x.numAcceptances = 0;
        x.tuningParm = update->getTuningParameter();
        x.updatePtr = update;
        info.insert( std::make_pair(id,x) );
        }
    else 
        {
        it->second.numProposals++;
        }
}

void UpdateStatistics::zeroOut(void) {

    for (auto& [key,val] : info)
        {
        val.numProposals = 0;
        val.numAcceptances = 0;
        }
}
