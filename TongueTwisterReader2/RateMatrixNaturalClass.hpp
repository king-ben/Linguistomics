#ifndef RateMatrixNaturalClass_hpp
#define RateMatrixNaturalClass_hpp

#include "RateMatrix.hpp"
class Exchangeabilities;
class StateFrequencies;



class RateMatrixNaturalClass : public RateMatrix{

    public:
                    RateMatrixNaturalClass(void) = delete;
                    RateMatrixNaturalClass(const StateFrequencies& f, const Exchangeabilities& r);
    
    private:
};

#endif
