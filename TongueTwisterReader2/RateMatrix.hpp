#ifndef RateMatrix_hpp
#define RateMatrix_hpp

#include <vector>
#include "Container.hpp"



class RateMatrix {

    public:
                                            RateMatrix(void) = delete;
                                            RateMatrix(size_t ns);
        virtual                            ~RateMatrix(void);
        size_t                              getNumSamples(void) { return matrices.size(); }
        size_t                              getNumSamples(void) const { return matrices.size(); }
        size_t                              getNumStates(void) { return numStates; }
        size_t                              getNumStates(void) const { return numStates; }
        std::vector<DoubleMatrix*>&         getRateMatrices(void) { return matrices; }
        const std::vector<DoubleMatrix*>&   getRateMatrices(void) const { return matrices; }
        DoubleMatrix&                       getMean(void) { return *mean; }
        DoubleMatrix&                       getLowerCI(void) { return *lowerCI; }
        DoubleMatrix&                       getUpperCI(void) { return * upperCI; }
    
    protected:
        void                                calculateMeanAndVariance(void);
        size_t                              numStates;
        std::vector<DoubleMatrix*>          matrices;
        DoubleMatrix*                       mean;
        DoubleMatrix*                       variance;
        DoubleMatrix*                       lowerCI;
        DoubleMatrix*                       upperCI;
};

#endif
