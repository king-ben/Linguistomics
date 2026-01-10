#ifndef UpdateStatistics_hpp
#define UpdateStatistics_hpp

#include <string>
#include <unordered_map>
#include "Update.hpp"

struct UpdateInfo {

    size_t      updateIdx;
    std::string name;
    int         numProposals;
    int         numAcceptances;
    double      tuningParm;
    UpdateType  updateType;
    Update*     updatePtr;
};

typedef std::unordered_map<uint64_t, UpdateInfo> UpdateMap;



class UpdateStatistics {

    public:
        void        accept(Update* update);  
        void        print(void);
        void        reject(Update* update);
        void        tune(void);
        void        zeroOut(void); 
    
    private:
    
        double      tuneSimplex(double acceptanceProb, double currentTuningValue);
        double      tuneFactor(double acceptanceProb, double currentTuningValue);
        double      tuneWindow(double acceptanceProb, double currentTuningValue);
        double      tuneProbability(double acceptanceProb, double currentTuningValue);
        UpdateMap   info;
};

#endif
