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
        std::string             getUpdateName(void) { return "Stationary Frequencies Update"; }
        std::string             parameterType(void) { return "ParameterFrequencies"; }
        void                    setDependants(void);
        double                  update(void);
        double                  update(double power);
        double                  updateFromPrior(void);
    
    private:
        double                  update(int k);
        ParameterFrequencies*   myParm;
        static double           minVal;
        size_t                  numStates;
        std::vector<double>     oldValues;
        std::vector<double>     newValues;
        std::vector<double>     alphaForward;
        std::vector<double>     alphaReverse;
};

#endif
