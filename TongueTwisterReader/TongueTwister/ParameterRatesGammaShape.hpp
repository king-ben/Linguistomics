#ifndef ParameterRatesGammaShape_hpp
#define ParameterRatesGammaShape_hpp

#include <string>
#include <vector>
#include "Parameter.hpp"
class Model;
class RandomVariable;



class ParameterRatesGammaShape : public Parameter {

    public:
                                ParameterRatesGammaShape(void) = delete;
                                ParameterRatesGammaShape(const ParameterRatesGammaShape& p) = delete;
                                ParameterRatesGammaShape(RandomVariable* r, Model* m, std::string n, double ep, int nc);
                               ~ParameterRatesGammaShape(void);
        void                    accept(void);
        void                    fillParameterValues(double* x, int& start, int maxNumValues);
        std::string             getHeader(void);
        std::string             getJsonString(void);
        int                     getNumValues(void) { return 1; }
        std::vector<double>&    getRates(void) { return rates[0]; }
        std::string             getString(void);
        char*                   getCString(void);
        std::string             getUpdateName(int idx);
        double                  lnPriorProbability(void);
        void                    print(void);
        void                    reject(void);
        double                  update(int iter);
    
    protected:
        double                  alpha[2];
        double                  expPriorVal;
        int                     numCategories;
        std::vector<double>     rates[2];
};

#endif
