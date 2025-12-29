#include <cmath>
#include <iostream>
#include <limits>
#include "Msg.hpp"
#include "RandomVariable.hpp"

enum {
    M = 397,
    R = 31,
    A = 0x9908B0DF,
    F = 1812433253,
    U = 11,
    S = 7,
    B = 0x9D2C5680,
    T = 15,
    C = 0xEFC60000,
    L = 18,
    MASK_LOWER = (1ull << R) - 1,
    MASK_UPPER = (1ull << R)
};

RandomVariable::RandomVariable(void) {

    uint32_t seed = (uint32_t)time(NULL);
    initialize(seed);
}

RandomVariable::RandomVariable(uint32_t seed) {

    initialize(seed);
}

uint32_t RandomVariable::extractU32(void) {

    int i = index;
    if (index >= N)
        {
        twist();
        i = 0;
        }

    uint32_t y = mt[i];
    index = (uint16_t)(i + 1);

    y ^= (mt[i] >> U);
    y ^= (y << S) & B;
    y ^= (y << T) & C;
    y ^= (y >> L);

    return y;
}

void RandomVariable::initialize(uint32_t seed) {
    memset(mt, 0, sizeof(mt));
    mt[0] = seed;
    for (uint32_t i=1; i<N; i++)
        {
        mt[i] = (F * (mt[i - 1] ^ (mt[i - 1] >> 30)) + i);
        }
    index = N;
    
    for (size_t i=0; i<10000; i++)
        extractU32();
}

void RandomVariable::twist(void) {

    for (uint32_t i=0; i<N; i++)
        {
        uint32_t x = (mt[i] & MASK_UPPER) + (mt[(i + 1) % N] & MASK_LOWER);
        uint32_t xA = x >> 1;

        if ( x & 0x1 )
            xA ^= A;

        mt[i] = mt[(i + M) % N] ^ xA;
        }
    index = 0;
}

double RandomVariable::uniformRv(void) {

    return (double)extractU32() / UINT32_MAX;
}

int RandomVariable::uniformRvInt(size_t scale) {

    return (int) (((double)extractU32() / UINT32_MAX) * (double)scale);
}
