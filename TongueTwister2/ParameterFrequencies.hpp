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
        std::vector<double>&        getFrequencies(size_t idx) { return freqs[idx]; }
        const std::vector<double>&  getFrequencies(size_t idx) const { return freqs[idx]; }
        size_t                      getNumStates(void) { return numStates; }
        void                        keep(void);
        double                      lnPriorProbability(void);
        void                        print(void);
        void                        restore(void);
    
    private:
                                    // accessed during likelihood calculations
        size_t                      numStates;         // 8 bytes
        
                                    // vectors (24 bytes each)
        std::vector<double>         freqs[2];          // 48 bytes
        std::vector<double>         alpha;             // 24 bytes
};

#endif
