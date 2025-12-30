#ifndef RateMatrix_hpp
#define RateMatrix_hpp

#include <vector>
#include "Container.hpp"
class ParameterExchangeabilities;
class ParameterFrequencies;
class Partition;
class RateMatrixHelper;



class RateMatrix {

    public:
                                    RateMatrix(void) = delete;
                                    RateMatrix(size_t ns);
        virtual                    ~RateMatrix(void);
        virtual RateMatrixHelper*   getHelper(void) = 0;
        DoubleMatrix&               getQ(void) { return *Q[activeMatrix]; }
        size_t                      getNumStates(void) const { return numStates; }
        void                        flipActiveValues(void);
        void                        keep(void);
        void                        print(void);
        void                        restore(void);
        virtual void                updateRateMatrix(void) = 0;
        
    protected:
                                    RateMatrix(const RateMatrix&) = delete;
        DoubleMatrix*               Q[2];
        size_t                      numStates;
        size_t                      activeMatrix;
};



class RateMatrixGTR : public RateMatrix {

    public:
                                    RateMatrixGTR(void) = delete;
                                    RateMatrixGTR(size_t ns, ParameterExchangeabilities* e, ParameterFrequencies* f);
        RateMatrixHelper*           getHelper(void) { return nullptr; }
        void                        updateRateMatrix(void);
    
    private:
        ParameterExchangeabilities* exchParm;
        ParameterFrequencies*       freqParm;
};



class RateMatrixCustom : public RateMatrix {

    public:
                                    RateMatrixCustom(void) = delete;
                                    RateMatrixCustom(size_t ns, ParameterExchangeabilities* e, ParameterFrequencies* f, Partition* p);
                                   ~RateMatrixCustom(void);
        RateMatrixHelper*           getHelper(void) { return rateMatrixHelper; }
        void                        updateRateMatrix(void);
    
    private:
        ParameterExchangeabilities* exchParm;
        ParameterFrequencies*       freqParm;
        RateMatrixHelper*           rateMatrixHelper;
        
                                    // cached pointers for inner loop (avoid repeated vector access)
        const int*                  rateIndexMap;       // points into helper
};

#endif
