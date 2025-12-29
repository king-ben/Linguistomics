#ifndef ParameterIndelRates_hpp
#define ParameterIndelRates_hpp

#include <vector>
#include "Parameter.hpp"
class Model;



class ParameterIndelRates : public Parameter {

    public:
                                        ParameterIndelRates(void) = delete;
                                        ParameterIndelRates(const ParameterIndelRates& pa) = delete;
                                        ParameterIndelRates(RandomVariable* r, Model* m, std::string n, double slen, double insLam, double delLam);
                                       ~ParameterIndelRates(void);
        void                            accept(void);
        void                            fillParameterValues(double* x, int& start, int maxNumValues);
        double                          getDeletionRate(void);
        double                          getExpectedSequenceLength(void);
        std::string                     getJsonString(void);
        std::string                     getHeader(void);
        double                          getInsertionRate(void);
        int                             getNumValues(void) { return 2; }
        std::string                     getString(void);
        char*                           getCString(void);
        std::string                     getUpdateName(int idx);
        double                          lnPriorProbability(void);
        void                            print(void);
        void                            reject(void);
        double                          update(int iter);
        
    protected:
        double                          expectedEpsilon(double slen);
        double                          updateFromPrior(void);
        double                          insertionLambda;
        double                          deletionLambda;
        double                          insertionRate[2];
        double                          deletionRate[2];
};

#endif
