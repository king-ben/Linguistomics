#ifndef SubstitutionModel_hpp
#define SubstitutionModel_hpp

#include <cstddef>
#include "Container.hpp"
class MathCache;
class ParameterFrequencies;
class RateMatrix;



class Ctmc {

    public:
                                    Ctmc(void) = delete;
                                    Ctmc(size_t ns);
        virtual                    ~Ctmc(void) = default;
        virtual MathCache*          allocateMathCache(void) = 0;
        bool                        checkTransitionProbabilities(DoubleMatrix* m, double thresshold) const;
        virtual void                computeTransitionProbs(double branchLength, DoubleMatrix* output, MathCache* cache=nullptr) const = 0;
        size_t                      getNumStates(void) const { return numStates; }
        
    protected:
        size_t                      numStates;
};

class CtmcJC69 : public Ctmc {

    public:
                                    CtmcJC69(void) = delete;
                                    CtmcJC69(size_t ns);
                                   ~CtmcJC69(void);
        MathCache*                  allocateMathCache(void) override { return nullptr; }
        void                        computeTransitionProbs(double branchLength, DoubleMatrix* output, MathCache* cache=nullptr) const override;
};

class CtmcF81 : public Ctmc {

    public:
                                    CtmcF81(void) = delete;
                                    CtmcF81(size_t ns, ParameterFrequencies* f);
                                   ~CtmcF81(void);
        MathCache*                  allocateMathCache(void) override { return nullptr; }
        void                        computeTransitionProbs(double branchLength, DoubleMatrix* output, MathCache* cache=nullptr) const override;
        
    private:
        ParameterFrequencies*       freqsParm;
};

class CtmcGeneral : public Ctmc {

    public:
                                    CtmcGeneral(void) = delete;
                                    CtmcGeneral(size_t ns, RateMatrix* r);
                                   ~CtmcGeneral(void);
        MathCache*                  allocateMathCache(void) override;
        void                        computeTransitionProbs(double branchLength, DoubleMatrix* output, MathCache* cache) const override;
        
    private:
        double                      factorial(size_t x) const;
        int                         logBase2Plus1(double x) const;
        int                         setQvalue(double tol) const;
        RateMatrix*                 rateMatrix;
        size_t                      largestFactorial;
        double*                     factorialTable;
};


#endif
