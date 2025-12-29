#include "McmcPhase.hpp"



McmcPhase::McmcPhase(int fbl, int pbl, int tl, int bl, int sl, int sf) {

    firstBurnLength = fbl;
    preBurnLength   = pbl;
    tuneLength      = tl;
    burnLength      = bl;
    sampleLength    = sl;
    sampleFrequency = sf;
    phases          = { "preburn", "tune", "burn", "sample" };
}

int McmcPhase::getLengthForPhase(std::string ph) {

    if (ph == "preburn")
        return preBurnLength;
    else if (ph == "tune")
        return tuneLength;
    else if (ph == "burn")
        return burnLength;
    else if (ph == "sample")
        return sampleLength;
    return 0;
}
