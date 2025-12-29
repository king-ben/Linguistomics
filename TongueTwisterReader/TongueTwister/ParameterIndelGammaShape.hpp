#ifndef ParameterIndelGammaShape_hpp
#define ParameterIndelGammaShape_hpp


#include <string>
#include <vector>
#include "Parameter.hpp"
class Model;
class RandomVariable;



class ParameterIndelGammaShape : public Parameter {

    public:
                                ParameterIndelGammaShape(void) = delete;
                                ParameterIndelGammaShape(const ParameterIndelGammaShape& p) = delete;
                                ParameterIndelGammaShape(RandomVariable* r, Model* m, std::string n, double ep, int nc);
                               ~ParameterIndelGammaShape(void);
        void                    accept(void);
        void                    fillParameterValues(double* x, int& start, int maxNumValues);
        std::string             getHeader(void);
        std::string             getJsonString(void);
        std::string             getString(void);
        char*                   getCString(void);
        int                     getNumValues(void) { return 1; }
        std::vector<double>&    getRates(void) { return rates[0]; }
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
