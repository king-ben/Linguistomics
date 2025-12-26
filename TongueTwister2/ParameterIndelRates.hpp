#ifndef ParameterIndelRates_hpp
#define ParameterIndelRates_hpp

#include "Parameter.hpp"
class Model;



class ParameterIndelRates : public Parameter {

    public:
                    ParameterIndelRates(void) = delete;
                    ParameterIndelRates(Model* m, RandomVariable* r, std::string n, double lamLam, double muLam);
        double      getDeletionLambda(void) { return deletionLambda; }
        double      getDeletionRate(void) { return deletionRate[0]; }
        double      getInsertionLambda(void) { return insertionLambda; }
        double      getInsertionRate(void) { return insertionRate[0]; }
        void        keep(void);
        double      lnPriorProbability(void);
        double      lnPriorProbability(double lambda, double mu);
        void        print(void);
        void        restore(void);
        void        setDeletionRate(double x) { deletionRate[0] = x; }
        void        setInsertionRate(double x) { insertionRate[0] = x; }
        
    protected:
        double      insertionLambda;
        double      deletionLambda;
        double      insertionRate[2];
        double      deletionRate[2];
    
};

#endif
