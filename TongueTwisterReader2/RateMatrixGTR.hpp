#ifndef RateMatrixGTR_hpp
#define RateMatrixGTR_hpp

#include "RateMatrix.hpp"
class Exchangeabilities;
class StateFrequencies;



class RateMatrixGTR : public RateMatrix {

    public:
                    RateMatrixGTR(void) = delete;
                    RateMatrixGTR(const StateFrequencies& f, const Exchangeabilities& r);
    
    private:
};

#endif
