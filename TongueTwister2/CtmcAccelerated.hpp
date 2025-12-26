#ifndef CtmcAccelerated_hpp
#define CtmcAccelerated_hpp

#include "Ctmc.hpp"

class RateMatrix;
class MathCacheAccelerated;



class CtmcAccelerated : public Ctmc {

    public:
                                    CtmcAccelerated(void) = delete;
                                    CtmcAccelerated(size_t ns, RateMatrix* r);
                                   ~CtmcAccelerated(void);
        MathCache*                  allocateMathCache(void) override;
        void                        computeTransitionProbs(double branchLength, DoubleMatrix* output, MathCache* cache) const override;
        
        const double*               getQData(void) const;
        
    private:
        RateMatrix*                 rateMatrix;
};

#endif
