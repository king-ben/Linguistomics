#include <cmath>
#include <iomanip>
#include <iostream>
#include "Model.hpp"
#include "Msg.hpp"
#include "ParameterIndelGammaShape.hpp"
#include "Probability.hpp"
#include "RandomVariable.hpp"
#include "TransitionProbabilities.hpp"



ParameterIndelGammaShape::ParameterIndelGammaShape(RandomVariable* r, Model* m, std::string n, double ep, int nc) : Parameter(r, m, n) {

    //parmId = IndelGammaShapeParm;

    std::cout << "   * Setting up gamma shape parameter for insertion/deletion rates " << std::endl;

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

ParameterIndelGammaShape::~ParameterIndelGammaShape(void) {

    delete [] parmStr;
}
void ParameterIndelGammaShape::accept(void) {

    alpha[1] = alpha[0];
    for (int i=0; i<numCategories; i++)
        rates[1][i] = rates[0][i];
}

void ParameterIndelGammaShape::fillParameterValues(double* x, int& start, int maxNumValues) {

    if (start + 1 > maxNumValues)
        Msg::error("Exceeding bounds when filling parameter value array");
    int j = start;
    x[j++] = alpha[0];
    start = j;
}

std::string ParameterIndelGammaShape::getHeader(void) {

    std::string str = "IndelAlpha";
    str += '\t';
    return str;
}

std::string ParameterIndelGammaShape::getJsonString(void) {

    std::string str = "";
    
    return str;
}

std::string ParameterIndelGammaShape::getString(void) {

    std::string str = std::to_string(alpha[0]);
    str += '\t';
    return str;
}

char* ParameterIndelGammaShape::getCString(void) {

    snprintf(parmStr, parmStrLen, "%1.6lf\t", alpha[0]);
    return parmStr;
}

std::string ParameterIndelGammaShape::getUpdateName(int) {

    return "indel shape parameter";
}

double ParameterIndelGammaShape::lnPriorProbability(void) {

    return log(expPriorVal) - expPriorVal * alpha[0];
}

void ParameterIndelGammaShape::print(void) {

    std::cout << std::fixed << std::setprecision(6);
    std::cout << alpha[0] << " " << std::endl;
}

void ParameterIndelGammaShape::reject(void) {

    alpha[0] = alpha[1];
    for (int i=0; i<numCategories; i++)
        rates[0][i] = rates[1][i];
    modelPtr->flipActiveLikelihood();
}

double ParameterIndelGammaShape::update(int) {

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
