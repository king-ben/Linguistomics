#ifndef RateMatrix_hpp
#define RateMatrix_hpp

#include <vector>
#include "Container.hpp"

class RateMatrixHelper;


class RateMatrix {

    public:
                                RateMatrix(void);
                               ~RateMatrix(void);
        std::vector<double>&    getEquilibriumFrequencies(void) { return equilibriumFrequencies[activeMatrix]; }
        DoubleMatrix&           getRateMatrix(void) { return *Q[activeMatrix]; }
        void                    flipActiveValues(void);
        void                    initialize(int d, RateMatrixHelper* h);
        void                    updateRateMatrix(std::vector<double>& rates, std::vector<double>& f);
        
    private:
                                RateMatrix(const RateMatrix&) = delete;
        RateMatrixHelper*       rateMatrixHelper;
        int                     numStates;
        int                     activeMatrix;
        DoubleMatrix*           Q[2];
        std::vector<double>     equilibriumFrequencies[2];
        bool                    isInitialized;
};

#endif
