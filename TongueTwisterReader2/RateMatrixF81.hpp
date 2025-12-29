#ifndef RateMatrixF81_hpp
#define RateMatrixF81_hpp

#include "RateMatrix.hpp"
class StateFrequencies;



class RateMatrixF81 : public RateMatrix{

    public:
                    RateMatrixF81(void) = delete;
                    RateMatrixF81(const StateFrequencies& f);
    
    private:
};

#endif
