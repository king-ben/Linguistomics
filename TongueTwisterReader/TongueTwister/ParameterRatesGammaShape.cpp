#include <cmath>
#include <iomanip>
#include <iostream>
#include "Model.hpp"
#include "Msg.hpp"
#include "ParameterRatesGammaShape.hpp"
#include "Probability.hpp"
#include "RandomVariable.hpp"
#include "TransitionProbabilities.hpp"



ParameterRatesGammaShape::ParameterRatesGammaShape(RandomVariable* r, Model* m, std::string n, double ep, int nc) : Parameter(r, m, n) {

    //parmId = RatesGammaShapeParm;

    std::cout << "   * Setting up gamma shape parameter for site rates " << std::endl;

    updateChangesRateMatrix = false;
    
    expPriorVal = ep;
    numCategories = nc;
    
    alpha[0] = Probability::Exponential::rv(rv, expPriorVal);
    alpha[1] = alpha[0];
    
    rates[0].resize(numCategories);
    rates[1].resize(numCategories);
    Probability::Gamma::discretization(rates[0], alpha[0], alpha[0], numCategories, false);
    rates[1] = rates[0];

    parmStrLen = 20;
    parmStr = new char[parmStrLen];
}

ParameterRatesGammaShape::~ParameterRatesGammaShape(void) {

    delete [] parmStr;
}

void ParameterRatesGammaShape::accept(void) {

    alpha[1] = alpha[0];
    rates[1] = rates[0];
}

void ParameterRatesGammaShape::fillParameterValues(double* x, int& start, int maxNumValues) {

    if (start + 1 > maxNumValues)
        Msg::error("Exceeding bounds when filling parameter value array");
    int j = start;
    x[j++] = alpha[0];
    start = j;
}

std::string ParameterRatesGammaShape::getHeader(void) {

    std::string str = "RatesAlpha";
    str += '\t';
    return str;
}

std::string ParameterRatesGammaShape::getJsonString(void) {

    std::string str = "";
    
    return str;
}

std::string ParameterRatesGammaShape::getString(void) {

    std::string str = std::to_string(alpha[0]);
    str += '\t';
    return str;
}

char* ParameterRatesGammaShape::getCString(void) {

    snprintf(parmStr, parmStrLen, "%1.6lf\t", alpha[0]);
    return parmStr;
}

std::string ParameterRatesGammaShape::getUpdateName(int) {

    return "gamma shape parameter";
}

double ParameterRatesGammaShape::lnPriorProbability(void) {

    return log(expPriorVal) - expPriorVal * alpha[0];
}

void ParameterRatesGammaShape::print(void) {

    std::cout << std::fixed << std::setprecision(6);
    std::cout << alpha[0] << " " << std::endl;
}

void ParameterRatesGammaShape::reject(void) {

    alpha[0] = alpha[1];
    rates[0] = rates[1];
    modelPtr->flipActiveLikelihood();
}

double ParameterRatesGammaShape::update(int) {

    lastUpdateType.first = this;
    lastUpdateType.second = 0;
    
    double tuning = log(4.0);

    double newAlpha = alpha[0] * exp(tuning*(rv->uniformRv()-0.5));
    double lnP = log(newAlpha) - log(alpha[0]);
    
    alpha[0] = newAlpha;
    Probability::Gamma::discretization(rates[0], alpha[0], alpha[0], numCategories, false);

    updateChangesTransitionProbabilities = true;
    TransitionProbabilities* tip = modelPtr->getTransitionProbabilities();
    tip->flipActive();
    tip->setNeedsUpdate(true);
    tip->setTransitionProbabilities();

    modelPtr->setUpdateLikelihood();
    modelPtr->flipActiveLikelihood();

    return lnP;
}
