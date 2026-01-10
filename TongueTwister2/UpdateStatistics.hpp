#ifndef UpdateStatistics_hpp
#define UpdateStatistics_hpp

#include <string>
#include <unordered_map>
class Update;

struct UpdateInfo {

    std::string name;
    int         numProposals;
    int         numAcceptances;
    double      tuningParm;
    Update*     updatePtr;
};

typedef std::unordered_map<uint64_t, UpdateInfo> UpdateMap;



class UpdateStatistics {

    public:
        void        accept(Update* update);  
        void        print(void);
        void        reject(Update* update);
        void        zeroOut(void); 
    
    private:
        UpdateMap   info;
};

#endif
