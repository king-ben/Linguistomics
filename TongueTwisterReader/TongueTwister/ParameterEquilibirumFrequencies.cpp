#include <iomanip>
#include <iostream>
#include "Model.hpp"
#include "Msg.hpp"
#include "ParameterEquilibirumFrequencies.hpp"
#include "Probability.hpp"
#include "RandomVariable.hpp"
#include "RateMatrix.hpp"
#include "TransitionProbabilities.hpp"

double ParameterEquilibirumFrequencies::minVal = 0.001;



ParameterEquilibirumFrequencies::ParameterEquilibirumFrequencies(RandomVariable* r, Model* m, std::string n, int ns) : Parameter(r, m, n) {

    //parmId = EquilibriumFrequenciesParm;

    std::cout << "   * Setting up equilibrium frequencies parameter " << std::endl;

    updateChangesRateMatrix = true;

    numStates = ns;
    freqs[0].resize(numStates);
    freqs[1].resize(numStates);
    alpha.resize(numStates);
    for (int i=0; i<numStates; i++)
        alpha[i] = 1.0;
    Probability::Dirichlet::rv(rv, alpha, freqs[0]);
    for (int i=0; i<numStates; i++)
        freqs[0][i] = 1.0 / numStates;
    Probability::Helper::normalize(freqs[0], minVal);
    freqs[1] = freqs[0];
    
    parmStrLen = numStates*Model::PrintNumberSize;
    parmStr = new char[parmStrLen];
}

ParameterEquilibirumFrequencies::~ParameterEquilibirumFrequencies(void) {

    delete [] parmStr;
}

void ParameterEquilibirumFrequencies::accept(void) {

    for (int i=0; i<numStates; i++)
        freqs[1][i] = freqs[0][i];
}

void ParameterEquilibirumFrequencies::fillParameterValues(double* x, int& start, int maxNumValues) {

    if (start + numStates > maxNumValues)
        Msg::error("Exceeding bounds when filling parameter value array");
    int j = start;
    for (int i=0; i<numStates; i++)
        x[j++] = freqs[0][i];
    start = j;
}

std::string ParameterEquilibirumFrequencies::getJsonString(void) {

    std::string str = "";
    for (int i=0; i<numStates; i++)
        {
        str += std::to_string(freqs[0][i]);
        str += '\t';
        }
    return str;
}

std::string ParameterEquilibirumFrequencies::getHeader(void) {

    std::string str = "";
    for (int i=0; i<numStates; i++)
        {
        str += "F[";
        str += std::to_string(i);
        str += "]";
        str += '\t';
        }
    return str;
}

std::string ParameterEquilibirumFrequencies::getString(void) {

    std::string str = "";
    for (int i=0; i<numStates; i++)
        {
        str += std::to_string(freqs[0][i]);
        str += '\t';
        }
    return str;
}

char* ParameterEquilibirumFrequencies::getCString(void) {
    char* p = parmStr;
    for (int i=0; i<numStates; i++)
        p += snprintf(p, parmStrLen - (p - parmStr), "%1.6lf\t", freqs[0][i]);
    *p = '\0';
    return parmStr;
}

std::string ParameterEquilibirumFrequencies::getUpdateName(int) {

    return "equilibrium frequencies";
}

double ParameterEquilibirumFrequencies::lnPriorProbability(void) {

    return Probability::Helper::lnGamma(numStates-1);
}

void ParameterEquilibirumFrequencies::print(void) {

    std::cout << "[ ";
    std::cout << std::fixed << std::setprecision(6);
    for (int i=0; i<numStates; i++)
        {
        std::cout << freqs[0][i] << " ";
        }
    std::cout << "]" << std::endl;
}

std::vector<int> ParameterEquilibirumFrequencies::randomlyChooseIndices(int k, int n) {

    std::vector<int> possibleIndices(n);
    for (int i=0; i<n; i++)
        possibleIndices[i] = i;
        
    std::vector<int> indices(k);
    for (int i=0; i<k; i++)
        {
        int pos = rv->uniformRvInt(n-i);
        indices[i] = possibleIndices[pos];
        possibleIndices[pos] = possibleIndices[n-i-1];
        }
    
    return indices;
}

void ParameterEquilibirumFrequencies::reject(void) {

    for (int i=0; i<numStates; i++)
        freqs[0][i] = freqs[1][i];
    modelPtr->flipActiveLikelihood();

}

double ParameterEquilibirumFrequencies::update(int) {

    lastUpdateType.first  = this;
    lastUpdateType.second = 0;

    int k = 1;
    double bounce = 1.0;
    
    double lnP = 0.0;
    if (k == 1)
        {
        // resize vectors for move
        oldValues.resize(2, 0.0);
        newValues.resize(2, 0.0);
        alphaForward.resize(2, 0.0);
        alphaReverse.resize(2, 0.0);
        
        // parameterize dirichlet for forward move
        double alpha0 = 100.0;
        int indexToUpdate = rv->uniformRvInt(numStates);
        oldValues[0] = freqs[0][indexToUpdate];
        oldValues[1] = 1.0 - oldValues[0];
        alphaForward[0] = oldValues[0] * alpha0 + bounce;
        alphaForward[1] = oldValues[1] * alpha0 + bounce;
        if (alphaForward[0] < 0.0)
            Msg::error("Negative alpha[0] in k=1 proposal for exchangability rates");
        if (alphaForward[1] < 0.0)
            Msg::error("Negative alpha[1] in k=1 proposal for exchangability rates");
        
        // propose new values
        bool err = Probability::Dirichlet::rv(rv, alphaForward, newValues);
        if (err == true)
            Msg::error("Problem proposing dirichlet values in k=1 for exchangability rates");
        double factor = newValues[1] / oldValues[1];
        for (int i=0; i<numStates; i++)
            freqs[0][i] *= factor;
        freqs[0][indexToUpdate] = newValues[0];
        Probability::Helper::normalize(newValues, minVal);

        // parameterize dirichlet for reverse move
        alphaReverse[0] = newValues[0] * alpha0 + bounce;
        alphaReverse[1] = newValues[1] * alpha0 + bounce;

        // calculate proposal probability
        lnP = Probability::Dirichlet::lnPdf(alphaReverse, oldValues) - Probability::Dirichlet::lnPdf(alphaForward, newValues);
        lnP += (numStates - 2) * log(factor); // Jacobian
        }
    else if (k > 1 && k < numStates)
        {
        double alpha0 = 1000.0;
        std::vector<int> indicesToUpdate = randomlyChooseIndices(k, numStates);
        std::map<size_t,size_t> mapper;
        for (size_t i=0; i<indicesToUpdate.size(); i++)
            mapper.insert( std::make_pair(indicesToUpdate[i], i) );
            
        oldValues.resize(k+1, 0.0);
        newValues.resize(k+1, 0.0);
        alphaForward.resize(k+1, 0.0);
        alphaReverse.resize(k+1, 0.0);
        
        for (size_t i=0; i<numStates; i++)
            {
            std::map<size_t,size_t>::iterator it = mapper.find(i);
            if (it != mapper.end())
                oldValues[it->second] += freqs[0][it->first];
            else
                oldValues[k] += freqs[0][i];
            }
        
        for (size_t i=0; i<k+1; i++)
            alphaForward[i] = oldValues[i] * alpha0 + bounce;
        
        // draw a new value for the reduced vector
        bool err = false;
        do
            {
            err = Probability::Dirichlet::rv(rv, alphaForward, newValues);
            if (err == true)
                Msg::warning("Trying to recover from Dirichlet random variable error in state frequencies update (k=intermediate)");
            } while (err == true);
        Probability::Helper::normalize(newValues, minVal);
        
        // fill in the Dirichlet parameters for the reverse probability calculations
        for (int i=0; i<k+1; i++)
            alphaReverse[i] = newValues[i] * alpha0 + bounce;
        
        // fill in the full vector
        double factor = newValues[k] / oldValues[k];
        for (size_t i=0; i<numStates; i++)
            {
            std::map<size_t,size_t>::iterator it = mapper.find(i);
            if (it != mapper.end())
                freqs[0][i] = newValues[it->second];
            else
                freqs[0][i] *= factor;
            }
        
        // Hastings ratio
        lnP = Probability::Dirichlet::lnPdf(alphaReverse, oldValues) - Probability::Dirichlet::lnPdf(alphaForward, newValues);
        lnP += (numStates - k - 1) * log(factor); // Jacobian
        }
    else
        {
        // update all of the rates
        double alpha0 = 100000.0;
        oldValues.resize(numStates, 0.0);
        newValues.resize(numStates, 0.0);
        alphaForward.resize(numStates, 0.0);
        alphaReverse.resize(numStates, 0.0);

        for (int i=0; i<numStates; i++)
            {
            oldValues[i] = freqs[0][i];
            alphaForward[i] = oldValues[i] * alpha0 + bounce;
            }
        
        bool err = false;
        do
            {
            err = Probability::Dirichlet::rv(rv, alphaForward, newValues);
            if (err == true)
                Msg::warning("Trying to recover from Dirichlet random variable error in state frequencies update (k=N)");
            } while (err == true);
        Probability::Helper::normalize(newValues, minVal);
        
        for (int i=0; i<numStates; i++)
            alphaReverse[i] = newValues[i] * alpha0 + bounce;

        lnP = Probability::Dirichlet::lnPdf(alphaReverse, oldValues) - Probability::Dirichlet::lnPdf(alphaForward, newValues);
        
        freqs[0] = newValues;
        }
    
    // update the rate matrix and transition probabilities
    RateMatrix* rmat = modelPtr->getRateMatrix();
    rmat->flipActiveValues();
    rmat->updateRateMatrix(modelPtr->getExchangabilityRates(), freqs[0]);

    updateChangesTransitionProbabilities = true;
    TransitionProbabilities* tip = modelPtr->getTransitionProbabilities();
    tip->flipActive();
    tip->setNeedsUpdate(true);
    tip->setTransitionProbabilities();

    modelPtr->setUpdateLikelihood();
    modelPtr->flipActiveLikelihood();

    // quick check that the exchangability rates are still valid
    double sum = 0.0;
    bool isBad = false;
    for (int i=0; i<freqs[0].size(); i++)
        {
        sum += freqs[0][i];
        if (freqs[0][i] < 0.0)
            isBad = true;
        }
    if (isBad == true || fabs(1.0-sum) > 0.000001)
        lnP = -1e100;

    return lnP;
}

double ParameterEquilibirumFrequencies::updateFromPrior(void) {

    //lastUpdateType = "random equilibrium frequencies";
    lastUpdateType.first = this;
    lastUpdateType.second = 1;

    // draw from the prior distribution, which is a flat Dirichlet distribution
    Probability::Dirichlet::rv(rv, alpha, freqs[0]);

    // update the rate matrix and transition probabilities
    RateMatrix* rmat = modelPtr->getRateMatrix();
    rmat->flipActiveValues();
    rmat->updateRateMatrix(modelPtr->getExchangabilityRates(), freqs[0]);

    updateChangesTransitionProbabilities = true;
    TransitionProbabilities* tip = modelPtr->getTransitionProbabilities();
    tip->flipActive();
    tip->setNeedsUpdate(true);
    tip->setTransitionProbabilities();

    return 0.0;
}
