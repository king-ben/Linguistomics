#include <iomanip>
#include <iostream>
#include "Model.hpp"
#include "Msg.hpp"
#include "ParameterExchangabilityRates.hpp"
#include "Probability.hpp"
#include "RandomVariable.hpp"
#include "RateMatrix.hpp"
#include "TransitionProbabilities.hpp"

double ParameterExchangabilityRates::minVal = 0.000001;



ParameterExchangabilityRates::ParameterExchangabilityRates(RandomVariable* r, Model* m, std::string n, int ns) : Parameter(r, m, n) {

    //parmId = ExchangabilityRatesParm;

    std::cout << "   * Setting up GTR exchangeability rates parameter " << std::endl;

    updateChangesRateMatrix = true;

    numStates = ns;
    numRates = numStates * (numStates-1) / 2;
    for (int i=0; i<numStates; i++)
        {
        for (int j=i+1; j<numStates; j++)
            {
            std::string str = std::to_string(i) + "-" + std::to_string(j);
            rateLabels.push_back(str);
            }
        }
            
    rates[0].resize(numRates);
    rates[1].resize(numRates);
    alpha.resize(numRates);
    for (int i=0; i<numRates; i++)
        alpha[i] = 1.0;
    Probability::Dirichlet::rv(rv, alpha, rates[0]);
    for (int i=0; i<numRates; i++)
        rates[0][i] = 1.0 / numRates;
    Probability::Helper::normalize(rates[0], minVal);
    rates[1] = rates[0];
    
    parmStrLen = numRates* Model::PrintNumberSize;
    parmStr = new char[parmStrLen];
}

ParameterExchangabilityRates::ParameterExchangabilityRates(RandomVariable* r, Model* m, std::string n, int ns, std::vector<std::string> labs) : Parameter(r, m, n) {

    std::cout << "   * Setting up custom exchangeability rates parameter " << std::endl;

    updateChangesRateMatrix = true;

    numStates = ns;
    rateLabels = labs;
    numRates = (int)rateLabels.size();
    
    rates[0].resize(numRates);
    rates[1].resize(numRates);
    alpha.resize(numRates);
    for (int i=0; i<numRates; i++)
        alpha[i] = 1.0;
    Probability::Dirichlet::rv(rv, alpha, rates[0]);
    rates[1] = rates[0];

    parmStrLen = numRates* Model::PrintNumberSize;
    parmStr = new char[parmStrLen];
}

ParameterExchangabilityRates::~ParameterExchangabilityRates(void) {

    delete [] parmStr;
}

void ParameterExchangabilityRates::accept(void) {

    for (int i=0; i<numRates; i++)
        rates[1][i] = rates[0][i];
}

void ParameterExchangabilityRates::fillParameterValues(double* x, int& start, int maxNumValues) {

    if (start + numRates > maxNumValues)
        Msg::error("Exceeding bounds when filling parameter value array");
    int j = start;
    for (int i=0; i<numRates; i++)
        x[j++] = rates[0][i];
    start = j;
}

std::string ParameterExchangabilityRates::getHeader(void) {

    std::string str = "";
    for (int i=0; i<rateLabels.size(); i++)
        {
        std::string s = "R[";
        s += rateLabels[i];
        s += "]";
        s += '\t';
        str += s;
        }
    return str;
}

std::string ParameterExchangabilityRates::getJsonString(void) {

    std::string str = "";
    
    return str;
}

std::string ParameterExchangabilityRates::getString(void) {

    std::string str = "";
    for (int i=0; i<numRates; i++)
        {
        str += std::to_string(rates[0][i]);
        str += '\t';
        }
    return str;
}

char* ParameterExchangabilityRates::getCString(void) {

    char* p = parmStr;
    for (int i=0; i<numRates; i++)
        p += snprintf(p, parmStrLen - (p - parmStr), "%1.6lf\t", rates[0][i]);
    *p = '\0';
    return parmStr;
}

std::string ParameterExchangabilityRates::getUpdateName(int idx) {

    if (idx == 0)
        return "exchangeability rates (k=1)";
    else if (idx == 1)
        return "exchangeability rates (k=10)";
    else if (idx == 2)
        return "exchangeability rates (All)";
    else
        return "random exchangeability rates";
}

double ParameterExchangabilityRates::lnPriorProbability(void) {

    return Probability::Helper::lnGamma(numRates-1);
}

void ParameterExchangabilityRates::print(void) {

    std::cout << "[ ";
    std::cout << std::fixed << std::setprecision(6);
    for (int i=0; i<rates[0].size(); i++)
        {
        std::cout << rates[0][i] << " ";
        }
    std::cout << "]" << std::endl;
}

std::vector<int> ParameterExchangabilityRates::randomlyChooseIndices(int k, int n) {

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

void ParameterExchangabilityRates::reject(void) {

    for (int i=0; i<numRates; i++)
        rates[0][i] = rates[1][i];
    modelPtr->flipActiveLikelihood();
}

void ParameterExchangabilityRates::setEqual(void) {

    double val = (double)1.0 / numRates;
    double sum = 0.0;
    for (int i=0; i<numRates; i++)
        {
        rates[0][i] = val;
        sum += rates[0][i];
        }
    for (int i=0; i<numRates; i++)
        rates[0][i] /= sum;
    rates[1] = rates[0];
}

double ParameterExchangabilityRates::update(int) {

    double bounce = 1.0;
    
    // choose the number of rates to update
    int k = 1;
    double u = rv->uniformRv();
    if (u <= 0.5)
        k = 1;
    else if (u > 0.5 && u <= 0.9)
        k = 10;
    else
        k = numRates;
        
    k = 1;
            
    // update the rates
    double lnP = 0.0;
    if (k == 1)
        {
        lastUpdateType.first = this;
        lastUpdateType.second = 0;

        // resize vectors for move
        oldValues.resize(2, 0.0);
        newValues.resize(2, 0.0);
        alphaForward.resize(2, 0.0);
        alphaReverse.resize(2, 0.0);

        // parameterize dirichlet for forward move
        double alpha0 = 200.0;
        int indexToUpdate = rv->uniformRvInt(numRates);
        oldValues[0] = rates[0][indexToUpdate];
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
        for (int i=0; i<numRates; i++)
            rates[0][i] *= factor;
        rates[0][indexToUpdate] = newValues[0];
        Probability::Helper::normalize(rates[0], minVal);

        // parameterize dirichlet for reverse move
        alphaReverse[0] = newValues[0] * alpha0 + bounce;
        alphaReverse[1] = newValues[1] * alpha0 + bounce;

        // calculate proposal probability
        lnP = Probability::Dirichlet::lnPdf(alphaReverse, oldValues) - Probability::Dirichlet::lnPdf(alphaForward, newValues);
        lnP += (numRates - 2) * log(factor); // Jacobian
        }
    else if (k > 1 && k < numRates)
        {
        lastUpdateType.first = this;
        lastUpdateType.second = 1;

        // resize vectors for move
        oldValues.resize(k+1, 0.0);
        newValues.resize(k+1, 0.0);
        alphaForward.resize(k+1, 0.0);
        alphaReverse.resize(k+1, 0.0);

        // parameterize dirichlet for forward move
        double alpha0 = 1000.0;
        std::vector<int> indicesToUpdate = randomlyChooseIndices(k, numRates);
        std::map<size_t,size_t> mapper;
        for (size_t i=0; i<indicesToUpdate.size(); i++)
            mapper.insert( std::make_pair(indicesToUpdate[i], i) );
        for (size_t i=0; i<numRates; i++)
            {
            std::map<size_t,size_t>::iterator it = mapper.find(i);
            if (it != mapper.end())
                oldValues[it->second] += rates[0][it->first];
            else
                oldValues[k] += rates[0][i];
            }
        for (size_t i=0; i<k+1; i++)
            {
            alphaForward[i] = oldValues[i] * alpha0 + bounce;
            if (alphaForward[i] < 0.0)
                Msg::error("Negative alpha[" + std::to_string(i)  + ")] in k=" + std::to_string(k) + " proposal for exchangability rates");
            }
        
        // propose new values
        bool err = Probability::Dirichlet::rv(rv, alphaForward, newValues);
        if (err == true)
            Msg::error("Problem proposing dirichlet values in k=" + std::to_string(k) + " for exchangability rates");
        double factor = newValues[k] / oldValues[k];
        for (size_t i=0; i<numRates; i++)
            {
            std::map<size_t,size_t>::iterator it = mapper.find(i);
            if (it != mapper.end())
                rates[0][i] = newValues[it->second];
            else
                rates[0][i] *= factor;
            }
        Probability::Helper::normalize(rates[0], minVal);

        // parameterize dirichlet for reverse move
        for (size_t i=0; i<k+1; i++)
            alphaReverse[i] = newValues[i] * alpha0 + bounce;

        // calculate proposal probability
        lnP = Probability::Dirichlet::lnPdf(alphaReverse, oldValues) - Probability::Dirichlet::lnPdf(alphaForward, newValues);
        lnP += (numRates - k - 1) * log(factor); // Jacobian
        }
    else
        {
        lastUpdateType.first = this;
        lastUpdateType.second = 2;

        // resize vectors for move
        oldValues.resize(numRates);
        newValues.resize(numRates);
        alphaForward.resize(numRates);
        alphaReverse.resize(numRates);

        // parameterize dirichlet for forward move
        double alpha0 = 4000.0;
        for (int i=0; i<numRates; i++)
            {
            oldValues[i] = rates[0][i];
            alphaForward[i] = oldValues[i] * alpha0 + bounce;
            if (alphaForward[i] < 0.0)
                Msg::error("Negative alpha[" + std::to_string(i)  + ")] in k=" + std::to_string(numRates) + " proposal for exchangability rates");
            }
        
        // propose new values
        bool err = Probability::Dirichlet::rv(rv, alphaForward, newValues);
        if (err == true)
            Msg::error("Problem proposing dirichlet values in k=" + std::to_string(numRates) + " for exchangability rates");
        Probability::Helper::normalize(newValues, minVal);
        
        // parameterize dirichlet for reverse move
        for (int i=0; i<numRates; i++)
            {
            rates[0][i] = newValues[i];
            alphaReverse[i] = newValues[i] * alpha0 + bounce;
            }

        // calculate proposal probability
        lnP = Probability::Dirichlet::lnPdf(alphaReverse, oldValues) - Probability::Dirichlet::lnPdf(alphaForward, newValues);
        }
    
    // update the rate matrix and transition probabilities
    RateMatrix* rmat = modelPtr->getRateMatrix();
    rmat->flipActiveValues();
    rmat->updateRateMatrix(rates[0], modelPtr->getEquilibriumFrequencies());

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
    for (int i=0; i<rates[0].size(); i++)
        {
        sum += rates[0][i];
        if (rates[0][i] < 0.0)
            isBad = true;
        }
    if (isBad == true || fabs(1.0-sum) > 0.000001)
        lnP = -1e100;

    return lnP;
}

double ParameterExchangabilityRates::updateFromPrior(void) {

    lastUpdateType.first = this;
    lastUpdateType.second = 3;

    // draw from the prior distribution, which is a flat Dirichlet distribution
    Probability::Dirichlet::rv(rv, alpha, rates[0]);

    // update the rate matrix and transition probabilities
    RateMatrix* rmat = modelPtr->getRateMatrix();
    rmat->flipActiveValues();
    rmat->updateRateMatrix(rates[0], modelPtr->getEquilibriumFrequencies());

    updateChangesTransitionProbabilities = true;
    TransitionProbabilities* tip = modelPtr->getTransitionProbabilities();
    tip->flipActive();
    tip->setNeedsUpdate(true);
    tip->setTransitionProbabilities();

    return 0.0;
}
