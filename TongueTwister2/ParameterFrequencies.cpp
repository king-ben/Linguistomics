#include "ParameterFrequencies.hpp"
#include "Probability.hpp"
#include "RandomVariable.hpp"



ParameterFrequencies::ParameterFrequencies(Model* m, RandomVariable* r, std::string n, size_t ns) : 
    Parameter(m, r, n), numStates(ns) {

    alpha.resize(numStates);
    for (size_t i=0; i<numStates; i++)
        alpha[i] = 1.0;
    isPriorFlat = true;
    
    freqs[0].resize(numStates);
    freqs[1].resize(numStates);
    Probability::Dirichlet::rv(rng, alpha, freqs[0]);
    freqs[1] = freqs[0];
}

ParameterFrequencies::~ParameterFrequencies(void) {

}

void ParameterFrequencies::keep(void) {

    freqs[1] = freqs[0];
}

double ParameterFrequencies::lnPriorProbability(void) {

    if (isPriorFlat == true)
        return 0.0;
    return Probability::Dirichlet::lnPdf(alpha, freqs[0]);
}

void ParameterFrequencies::print(void) {

}

void ParameterFrequencies::restore(void) {

    freqs[0] = freqs[1];
}
