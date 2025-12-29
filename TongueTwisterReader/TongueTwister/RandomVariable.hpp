#ifndef RandomVariable_hpp
#define RandomVariable_hpp

#include <cstdint>
#include <vector>

class RandomVariable {
    static const int N = 624;

    public:
                    RandomVariable(void);
                    RandomVariable(uint32_t seed);
        double      uniformRv(void);
        int         uniformRvInt(size_t scale);


    protected:
        uint32_t    extractU32(void);
        void        initialize(uint32_t seed);
        void        twist(void);
        uint16_t    index;
        uint32_t    mt[N];
};

#endif
