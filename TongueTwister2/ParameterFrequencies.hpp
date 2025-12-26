#ifndef ParameterFrequencies_hpp
#define ParameterFrequencies_hpp

#include <vector>
#include "Parameter.hpp"



class ParameterFrequencies : public Parameter {

    public:
                                    ParameterFrequencies(void) = delete;
                                    ParameterFrequencies(Model* m, RandomVariable* r, std::string n, size_t ns);
                                   ~ParameterFrequencies(void);
        std::vector<double>&        getAlpha(void) { return alpha; }
        std::vector<double>&        getFrequencies(void) { return freqs[0]; }
        const std::vector<double>&  getFrequencies(void) const { return freqs[0]; }
        size_t                      getNumStates(void) { return numStates; }
        void                        keep(void);
        double                      lnPriorProbability(void);
        void                        print(void);
        void                        restore(void);
    
    private:
        size_t                      numStates;
        std::vector<double>         freqs[2];
        std::vector<double>         alpha;
        bool                        isPriorFlat;
};

#endif
