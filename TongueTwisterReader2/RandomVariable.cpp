#include "RandomVariable.hpp"



RandomVariable::RandomVariable(void) {

    std::random_device rd;
    std::seed_seq seq{ rd(), rd(), rd(), rd(), rd(), rd(), rd(), rd() };
    rng = new std::mt19937(seq);
}

RandomVariable::RandomVariable(uint32_t seed) {

    if (seed == 0)
        {
        std::random_device rd;
        std::seed_seq seq{ rd(), rd(), rd(), rd(), rd(), rd(), rd(), rd() };
        rng = new std::mt19937(seq);
        }
    else 
        {
        rng = new std::mt19937(seed);
        }
}

double RandomVariable::uniformRv(void) {

    std::uniform_real_distribution<double> u01(0.0, 1.0);
    return u01(*rng);
}
