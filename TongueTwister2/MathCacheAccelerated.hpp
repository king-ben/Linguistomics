#ifndef MathCacheAccelerated_hpp
#define MathCacheAccelerated_hpp

#include "Container.hpp"
#include "MathCache.hpp"

#ifndef ACCELERATE_NEW_LAPACK
#define ACCELERATE_NEW_LAPACK
#endif
#include <Accelerate/Accelerate.h>



class MathCacheAccelerated : public MathCache {

    public:
                            MathCacheAccelerated(void);
                           ~MathCacheAccelerated(void);
        
        void                computeMatrixExponential(const DoubleMatrix& Q, double t, DoubleMatrix& P);
        
    private:
        void                blasMultiply(DoubleMatrix& A, DoubleMatrix& B, DoubleMatrix& C);
        double              infinityNorm(DoubleMatrix& A);
};

#endif
