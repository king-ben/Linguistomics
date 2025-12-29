#include <iomanip>
#include <iostream>
#include "Parameter.hpp"
#include "UpdateInfo.hpp"


UpdateInfo::UpdateInfo(void) {

}

UpdateInfo::~UpdateInfo(void) {

    for (std::map<std::pair<Parameter*,int>,AcceptTries*>::iterator it = info.begin(); it != info.end(); it++)
        delete it->second;
    info.clear();
}

void UpdateInfo::accept(std::pair<Parameter*,int>& key) {

    AcceptTries* at = getUpdateInfo(key);
    at->numTries++;
    at->numAccepts++;
}

AcceptTries* UpdateInfo::getUpdateInfo(std::pair<Parameter*,int>& key) {

    std::map<std::pair<Parameter*,int>,AcceptTries*>::iterator it = info.find(key);
    if (it == info.end())
        {
        AcceptTries* newInfo = new AcceptTries;
        info.insert( std::make_pair(key,newInfo) );
        return newInfo;
        }
    return it->second;
}

void UpdateInfo::print(void) {

    int len = 0;
    for (std::map<std::pair<Parameter*,int>,AcceptTries*>::iterator it = info.begin(); it != info.end(); it++)
        {
        std::string updateName = it->first.first->getUpdateName(it->first.second);
        if (updateName.length() > len)
            len = (int)updateName.length();
        }
        
    for (std::map<std::pair<Parameter*,int>,AcceptTries*>::iterator it = info.begin(); it != info.end(); it++)
        {
        std::string updateName = it->first.first->getUpdateName(it->first.second);
        std::cout << "   * Acceptance rate for update of " << updateName << " ";
        for (int i=0; i<len-updateName.length(); i++)
            std::cout << " ";
        std::cout << "= ";
        if (it->second->numTries > 0)
            std::cout << std::fixed << std::setprecision(1) << ((double)it->second->numAccepts / it->second->numTries) * 100.0 << "%";
        else
            std::cout << "N/A";

        std::cout << std::endl;
        }
}

void UpdateInfo::reject(std::pair<Parameter*,int>& key) {

    AcceptTries* at = getUpdateInfo(key);
    at->numTries++;
}

