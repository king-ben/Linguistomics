#ifndef RandomVariable_hpp
#define RandomVariable_hpp

#include <random>



class RandomVariable {

    public:
                                RandomVariable(void);
                                RandomVariable(uint32_t seed);
                                RandomVariable(const RandomVariable& r) = delete;
        double                  uniformRv(void);

    private:
        std::mt19937*           rng;
};

#endif
