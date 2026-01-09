#include <algorithm>
#include <cstdint>
#include "RandomVariable.hpp"
#include "Update.hpp"



Update::Update(Model* m, RandomVariable* r) : 
    model(m), rng(r), updatedParameter(nullptr), 
    rateMatrixNeedsUpdate(false), allTiprobsNeedUpdate(false), 
    singleBranchChanged(false), changedBranchLength(0.0) {
    
}

void Update::clearDependencyFlags(void) {

    rateMatrixNeedsUpdate = false;
    allTiprobsNeedUpdate = false;
    singleBranchChanged = false;
    changedBranchLength = 0.0;
}

double Update::priorSampleProb(double power) {

    /* Hill, A. V. 1910. The possible effects of the aggregation of the molecules of 
          hæmoglobin on its dissociation curves. Journal of Physiology, 40, iv–vii. */

    double h = 4.0;       // cliffiness (large h -> steeper cliff)
    double beta50 = 0.02; // cliff position (at what power is the probability 50/50)
    power = std::clamp(power, 0.0, 1.0); // keeps power in bounds

    // beta50 must be > 0, h must be > 0 for sane behavior
    if (beta50 <= 0.0) 
        return 0.0;
    if (h <= 0.0)      
        return 0.0;

    if (power == 0.0) 
        return 1.0;

    const double x = power / beta50;
    const double xh = std::pow(x, h);
    return 1.0 / (1.0 + xh);
}


uint64_t Update::rotl64(uint64_t x, int b) {

    return (x << b) | (x >> (64 - b));
}

uint64_t Update::load64_le(const unsigned char* p) {

    uint64_t v = 0;
    for (int i = 0; i < 8; i++)
        v |= (uint64_t)p[i] << (8 * i);
    return v;
}

uint64_t Update::load_tail_le(const unsigned char* p, size_t n) {

    uint64_t v = 0;
    for (size_t i = 0; i < n; i++)
        v |= (uint64_t)p[i] << (8 * i);
    return v;
}

void Update::sip_round(uint64_t& v0, uint64_t& v1, uint64_t& v2, uint64_t& v3) {

    v0 += v1; v1 = rotl64(v1, 13); v1 ^= v0; v0 = rotl64(v0, 32);
    v2 += v3; v3 = rotl64(v3, 16); v3 ^= v2;
    v0 += v3; v3 = rotl64(v3, 21); v3 ^= v0;
    v2 += v1; v1 = rotl64(v1, 17); v1 ^= v2; v2 = rotl64(v2, 32);
}

uint64_t Update::siphash24(const unsigned char* data, size_t len, uint64_t k0, uint64_t k1) {

    uint64_t v0 = 0x736f6d6570736575ULL ^ k0;
    uint64_t v1 = 0x646f72616e646f6dULL ^ k1;
    uint64_t v2 = 0x6c7967656e657261ULL ^ k0;
    uint64_t v3 = 0x7465646279746573ULL ^ k1;

    const unsigned char* end = data + (len & ~static_cast<size_t>(7));

    for (const unsigned char* p = data; p != end; p += 8) {
        uint64_t m = load64_le(p);
        v3 ^= m;
        sip_round(v0, v1, v2, v3);
        sip_round(v0, v1, v2, v3);
        v0 ^= m;
    }

    uint64_t b = (uint64_t)len << 56;
    size_t tail = len & 7;
    if (tail != 0)
        b |= load_tail_le(end, tail);

    v3 ^= b;
    sip_round(v0, v1, v2, v3);
    sip_round(v0, v1, v2, v3);
    v0 ^= b;

    v2 ^= 0xff;
    sip_round(v0, v1, v2, v3);
    sip_round(v0, v1, v2, v3);
    sip_round(v0, v1, v2, v3);
    sip_round(v0, v1, v2, v3);

    return v0 ^ v1 ^ v2 ^ v3;
}

uint64_t Update::hashUpdateName(const std::string& name) {

    // fixed 128-bit key -> deterministic across runs/machines.
    const uint64_t k0 = 0x0706050403020100ULL;
    const uint64_t k1 = 0x0f0e0d0c0b0a0908ULL;
    return siphash24(reinterpret_cast<const unsigned char*>(name.data()), name.size(), k0, k1);
}
