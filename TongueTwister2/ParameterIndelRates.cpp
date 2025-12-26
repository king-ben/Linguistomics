#include "ParameterIndelRates.hpp"
#include "Probability.hpp"
#include "RandomVariable.hpp"



ParameterIndelRates::ParameterIndelRates(Model* m, RandomVariable* r, std::string n, double lamLam, double muLam) : 
    Parameter(m, r, n), insertionLambda(lamLam), deletionLambda(muLam) {

    do {
        insertionRate[0] = Probability::Exponential::rv(rng, insertionLambda);
        deletionRate[0]  = Probability::Exponential::rv(rng, deletionLambda);
        } while(insertionRate[0] >= deletionRate[0]);
    insertionRate[1] = insertionRate[0];
    deletionRate[1]  = deletionRate[0];
}

void ParameterIndelRates::keep(void) {

    insertionRate[1] = insertionRate[0];
    deletionRate[1]  = deletionRate[0];
}

double ParameterIndelRates::lnPriorProbability(void) {

    return lnPriorProbability(getInsertionRate(), getDeletionRate());
}

double ParameterIndelRates::lnPriorProbability(double lambda, double mu) {

    double lnP = log(insertionLambda) - insertionLambda*lambda + log(deletionLambda) - deletionLambda*mu - 
                 log(1.0 - deletionLambda / (insertionLambda + deletionLambda));
    return lnP;
}

void ParameterIndelRates::print(void) {

}

void ParameterIndelRates::restore(void) {

    insertionRate[0] = insertionRate[1];
    deletionRate[0]  = deletionRate[1];
}
