#ifndef TransitionProbabilityCalculator_hpp
#define TransitionProbabilityCalculator_hpp

#include "Threads.hpp"
class Ctmc;
class MathCache;



class TransitionProbabilityCalculator : public TransitionProbabilityTask {

    public:
                                TransitionProbabilityCalculator(void) = delete;
                                TransitionProbabilityCalculator(Ctmc* m);
                               ~TransitionProbabilityCalculator(void);
        void                    set(double bl, DoubleMatrix* out);
    
    protected:
        void                    computeTransitionProbabilities(void);
        
    private:
        Ctmc*                   model;
        MathCache*              cache;
        DoubleMatrix*           output;
        double                  branchLength;
};

#endif
