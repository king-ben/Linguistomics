#ifndef ParameterEquilibirumFrequencies_hpp
#define ParameterEquilibirumFrequencies_hpp

#include <vector>
#include "Parameter.hpp"
class Model;



class ParameterEquilibirumFrequencies : public Parameter {
    public:
                                        ParameterEquilibirumFrequencies(void) = delete;
                                        ParameterEquilibirumFrequencies(const ParameterEquilibirumFrequencies& pr) = delete;
                                        ParameterEquilibirumFrequencies(RandomVariable* r, Model* m, std::string n, int ns);
                                       ~ParameterEquilibirumFrequencies(void);
        void                            accept(void);
        void                            fillParameterValues(double* x, int& start, int maxNumValues);
        std::string                     getHeader(void);
        std::string                     getJsonString(void);
        int                             getNumValues(void) { return numStates; }
        std::vector<double>&            getValue(void) { return freqs[0]; }
        std::string                     getString(void);
        char*                           getCString(void);
        std::string                     getUpdateName(int idx);
        double                          lnPriorProbability(void);
        void                            print(void);
        void                            reject(void);
        double                          update(int iter);
        
    private:
        std::vector<int>                randomlyChooseIndices(int k, int n);
        double                          updateFromPrior(void);
        int                             numStates;
        std::vector<double>             freqs[2];
        std::vector<double>             alpha;
        static double                   minVal;
        std::vector<double>             oldValues;
        std::vector<double>             newValues;
        std::vector<double>             alphaForward;
        std::vector<double>             alphaReverse;
};

#endif
