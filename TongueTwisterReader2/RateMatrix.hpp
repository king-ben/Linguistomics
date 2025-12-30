#ifndef RateMatrix_hpp
#define RateMatrix_hpp

#include <vector>
#include "Container.hpp"



class RateMatrix {

    public:
                                    RateMatrix(void) = delete;
                                    RateMatrix(size_t ns);
        virtual                    ~RateMatrix(void);
        std::vector<DoubleMatrix*>& getRateMatrices(void) { return matrices; }
        DoubleMatrix&               getMean(void) { return *mean; }
    
    protected:
        void                        calculateMeanAndVariance(void);
        size_t                      numStates;
        std::vector<DoubleMatrix*>  matrices;
        DoubleMatrix*               mean;
        DoubleMatrix*               variance;
        DoubleMatrix*               lowerCI;
        DoubleMatrix*               upperCI;
};

#endif
