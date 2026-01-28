#ifndef RateMatrixNaturalClass_hpp
#define RateMatrixNaturalClass_hpp

#include "RateMatrix.hpp"
class Exchangeabilities;
class Partition;
class StateFrequencies;



class RateMatrixNaturalClass : public RateMatrix{

    public:
                        RateMatrixNaturalClass(void) = delete;
                        RateMatrixNaturalClass(const StateFrequencies& f, const Exchangeabilities& r, Partition* part);
    
    private:
        DoubleMatrix    ncRates;
};

#endif
