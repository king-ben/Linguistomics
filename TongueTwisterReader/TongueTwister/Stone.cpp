#include "Stone.h"



Stone::Stone(double x) {

    power = x;
    samples.clear();
}

double Stone::maxValue(void) {

    double mx = 0.0;
    for (int i=0; i<samples.size(); i++)
        {
        if (samples[i] > mx)
            mx = samples[i];
        }
    return mx;
}
