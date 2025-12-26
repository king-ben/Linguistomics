#include "ParameterExchangeabilities.hpp"
#include "Probability.hpp"
#include "RandomVariable.hpp"



ParameterExchangeabilities::ParameterExchangeabilities(Model* m, RandomVariable* r, std::string n, size_t ns) : 
    Parameter(m, r, n), numStates(ns) {

    numRates = numStates * (numStates-1) / 2;
    
    alpha.resize(numRates);
    for (size_t i=0; i<numRates; i++)
        alpha[i] = 1.0;
    isPriorFlat = true;
    
    rates[0].resize(numRates);
    rates[1].resize(numRates);
    Probability::Dirichlet::rv(rng, alpha, rates[0]);
    rates[1] = rates[0];
}

ParameterExchangeabilities::ParameterExchangeabilities(Model* m, RandomVariable* r, std::string n, size_t ns, size_t ng) : 
    Parameter(m, r, n), numStates(ns) {

    numRates = (ng * (ng-1) / 2) + ng;
    
    alpha.resize(numRates);
    for (size_t i=0; i<numRates; i++)
        alpha[i] = 1.0;
    
    rates[0].resize(numRates);
    rates[1].resize(numRates);
    Probability::Dirichlet::rv(rng, alpha, rates[0]);
    rates[1] = rates[0];
}

ParameterExchangeabilities::~ParameterExchangeabilities(void) {

}

void ParameterExchangeabilities::keep(void) {

    rates[1] = rates[0];
}

double ParameterExchangeabilities::lnPriorProbability(void) {

    if (isPriorFlat == true)
        return 0.0;
    return Probability::Dirichlet::lnPdf(alpha, rates[0]);
}

void ParameterExchangeabilities::print(void) {

}

void ParameterExchangeabilities::restore(void) {

    rates[0] = rates[1];
}
