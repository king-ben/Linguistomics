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
        double                      getTuningParameter(void) { return tuningValues[lastUpdate]; }
        uint64_t                    getUpdateId(void);
        std::string                 getUpdateName(void) { return updateNames[lastUpdate]; }
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
        std::string                 updateNames[4];
        uint64_t                    updateHashes[4];
        double                      tuningValues[4];
};

#endif
