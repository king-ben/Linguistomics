#ifndef ParameterExchangabilityRates_hpp
#define ParameterExchangabilityRates_hpp

#include <vector>
#include "Parameter.hpp"
class Model;



class ParameterExchangabilityRates : public Parameter {

    public:
                                        ParameterExchangabilityRates(void) = delete;
                                        ParameterExchangabilityRates(const ParameterExchangabilityRates& pr) = delete;
                                        ParameterExchangabilityRates(RandomVariable* r, Model* m, std::string n, int ns);
                                        ParameterExchangabilityRates(RandomVariable* r, Model* m, std::string n, int ns, std::vector<std::string> labs);
                                       ~ParameterExchangabilityRates(void);
        void                            accept(void);
        void                            fillParameterValues(double* x, int& start, int maxNumValues);
        std::string                     getHeader(void);
        std::string                     getJsonString(void);
        int                             getNumValues(void) { return numRates; }
        std::string                     getString(void);
        char*                           getCString(void);
        std::string                     getUpdateName(int idx);
        std::vector<double>&            getValue(void) { return rates[0]; }
        double                          lnPriorProbability(void);
        void                            print(void);
        void                            reject(void);
        void                            setEqual(void);
        double                          update(int iter);
        
    private:
        std::vector<int>                randomlyChooseIndices(int k, int n);
        double                          updateFromPrior(void);
        int                             numStates;
        int                             numRates;
        std::vector<double>             rates[2];
        std::vector<double>             alpha;
        std::vector<std::string>        rateLabels;
        static double                   minVal;
        std::vector<double>             oldValues;
        std::vector<double>             newValues;
        std::vector<double>             alphaForward;
        std::vector<double>             alphaReverse;
};

#endif
