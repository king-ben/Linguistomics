#include <cmath>
#include <iostream>
#include "Msg.hpp"
#include "SteppingStones.h"
#include "Stone.h"


SteppingStones::SteppingStones(std::vector<double> pows) {

    for (int i=0; i<pows.size(); i++)
        stones.push_back( new Stone(pows[i]) );
}

SteppingStones::~SteppingStones(void) {

    for (int i=0; i<stones.size(); i++)
        delete stones[i];
}

void SteppingStones::addSample(int powIdx, double lnL) {

    stones[powIdx]->addSample(lnL);
}

void SteppingStones::clearSamples(int powIdx) {

    stones[powIdx]->clearSamples();
}

int SteppingStones::numSamplesPerStone(void) {

    int ns = stones[0]->numSamples();
    for (int i=1; i<stones.size(); i++)
        {
        if (stones[i]->numSamples() != ns)
            {
            for (int j=0; j<stones.size(); j++)
                std::cout << "Stone[" << j << "]: " << stones[j]->numSamples() << std::endl;
            Msg::error("Uneven number of samples");
            }
        }
    return ns;
}

double SteppingStones::marginalLikelihood(void) {

    std::vector<double> powers;
    for (int i=0; i<numStones(); i++)
        powers.push_back( stones[i]->getPower() );
        
    double marginal = 0.0;
    for (size_t i=1; i<numStones(); i++)
        {
        int samplesPerPath = numSamplesPerStone();
        double max = stones[i]->maxValue();
        
        // mean( exp(samples-max)^(beta[k-1]-beta[k]) )     or
        // mean( exp( (samples-max)*(beta[k-1]-beta[k]) ) )
        double mean = 0.0;
        std::vector<double>& likelihoodSamples = stones[i]->getSamples();
        for (size_t j=0; j<samplesPerPath; j++)
            {
            mean += exp( (likelihoodSamples[j] - max) * (powers[i-1]-powers[i]) ) / samplesPerPath;
            }
        
        marginal += log(mean) + (powers[i-1] - powers[i]) * max;
        }

    return marginal;
}
