#ifndef RateMatrixAverage_hpp
#define RateMatrixAverage_hpp

#include "RateMatrix.hpp"
class StateFrequencies;



class RateMatrixAverage : public RateMatrix{

    public:
                    RateMatrixAverage(void) = delete;
                    RateMatrixAverage(size_t ns);
                    RateMatrixAverage(const RateMatrix& Q, const StateFrequencies& f);
    
    private:
};

#endif
