#ifndef UpdateFrequencies_hpp
#define UpdateFrequencies_hpp

#include <vector>
#include "Update.hpp"
class Model;
class ParameterFrequencies;



class UpdateFrequencies : public Update {

    public:
                                    UpdateFrequencies(void) = delete;
                                    UpdateFrequencies(Model* m, RandomVariable* r, ParameterFrequencies* p);
        double                      getTuningParameter(void) { return updateInfo[lastUpdate].tuningParameter; }
        size_t                      getUpdateIdx(void) { return updateInfo[lastUpdate].updateIdx; }
        uint64_t                    getUpdateId(void) { return updateInfo[lastUpdate].updateHash; }
        UpdateType                  getUpdateType(void) { return updateInfo[lastUpdate].updateType;}
        std::string                 getUpdateName(void) { return updateInfo[lastUpdate].updateName; }
        std::string                 parameterType(void) { return "ParameterFrequencies"; }
        void                        setDependants(void);
        double                      update(void);
        double                      update(double power);
    
    private:
        ParameterFrequencies*       myParm;
        static double               minVal;
        size_t                      numStates;
        std::vector<double>         oldValues;
        std::vector<double>         newValues;
        std::vector<double>         alphaForward;
        std::vector<double>         alphaReverse;
        size_t                      lastUpdate;
        std::vector<std::string>    updateNames;
        std::vector<uint64_t>       updateHashes;
};

#endif
