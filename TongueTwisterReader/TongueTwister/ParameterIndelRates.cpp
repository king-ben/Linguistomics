#include <cmath>
#include <iomanip>
#include <iostream>
#include "Model.hpp"
#include "Msg.hpp"
#include "ParameterIndelRates.hpp"
#include "Probability.hpp"
#include "RandomVariable.hpp"
#include "TransitionProbabilities.hpp"



ParameterIndelRates::ParameterIndelRates(RandomVariable* r, Model* m, std::string n, double /*slen*/, double insLam, double delLam) : Parameter(r, m, n) {
    
    //parmId = IndelRatesParm;

    std::cout << "   * Setting up insertion/deletion rates parameter " << std::endl;

    updateChangesRateMatrix = false;
    insertionLambda = insLam;
    deletionLambda  = delLam;
    do {
        insertionRate[0] = Probability::Exponential::rv(rv, insertionLambda);
        deletionRate[0]  = Probability::Exponential::rv(rv, deletionLambda);
        } while(insertionRate[0] >= deletionRate[0]);
    insertionRate[1] = insertionRate[0];
    deletionRate[1]  = deletionRate[0];

    parmStrLen = 40;
    parmStr = new char[parmStrLen];
}

ParameterIndelRates::~ParameterIndelRates(void) {

}

void ParameterIndelRates::accept(void) {

    insertionRate[1] = insertionRate[0];
    deletionRate[1]  = deletionRate[0];
}

double ParameterIndelRates::expectedEpsilon(double slen) {

    double eps = 0.01;
    double increment = 0.001;
    bool increase = true;
    int nSwitches = 0;
    bool stopLoop = false;
    do {
        double expVal = eps / (1.0 - eps);

        if (increase == true && expVal < slen)
            {
            //std::cout << " 1 " << eps << " " << expVal << " " << slen << std::endl;
            increment *= 0.5;
            increase = false;
            nSwitches++;
            eps += increment;
            }
        else if (increase == true && expVal > slen)
            {
            //std::cout << " 2 " << eps << " " << expVal << " " << slen << std::endl;
            eps -= increment;
            }
        else if (increase == false && expVal < slen)
            {
            //std::cout << " 3 " << eps << " " << expVal << " " << slen << std::endl;
            eps += increment;
            }
        else
            {
            //std::cout << " 4 " << eps << " " << expVal << " " << slen << std::endl;
            increment *= 0.5;
            increase = true;
            nSwitches++;
            eps -= increment;
            }
    if (fabs(expVal - slen) < 0.00001)
        stopLoop = true;
    if (nSwitches > 100)
        stopLoop = true;
        
    } while (stopLoop == false);
    return eps;
}


void ParameterIndelRates::fillParameterValues(double* x, int& start, int maxNumValues) {

    if (start + 2 > maxNumValues)
        Msg::error("Exceeding bounds when filling parameter value array");
    int j = start;
    x[j++] = getInsertionRate();
    x[j++] = getDeletionRate();
    start = j;
}

double ParameterIndelRates::getDeletionRate(void) {

    return deletionRate[0];
}

double ParameterIndelRates::getExpectedSequenceLength(void) {

    double lambda = getInsertionRate();
    double mu = getDeletionRate();
    double x = lambda / mu;
    return x / (1.0 - x);
}

std::string ParameterIndelRates::getHeader(void) {

    std::string str = "Insertion";
    str += '\t';
    str += "Deletion";
    str += '\t';
    return str;
}

double ParameterIndelRates::getInsertionRate(void) {

    return insertionRate[0];
}

std::string ParameterIndelRates::getJsonString(void) {

    std::string str = "";
    
    return str;
}

std::string ParameterIndelRates::getString(void) {

    std::string str = std::to_string(getInsertionRate());
    str += '\t';
    str += std::to_string(getDeletionRate());
    str += '\t';
    return str;
}

char* ParameterIndelRates::getCString(void) {

    snprintf(parmStr, parmStrLen, "%1.6lf\t%1.6lf", getInsertionRate(), getDeletionRate());
    return parmStr;
}

std::string ParameterIndelRates::getUpdateName(int idx) {

    if (idx == 0)
        return "indel rates";
    else
        return "random indel rates";
}

double ParameterIndelRates::lnPriorProbability(void) {

    double lambda = getInsertionRate();
    double mu = getDeletionRate();
    double lnP = log(insertionLambda) - insertionLambda*lambda + log(deletionLambda) - deletionLambda*mu - log(1.0 - deletionLambda / (insertionLambda + deletionLambda));
    return lnP;
}

void ParameterIndelRates::print(void) {

    std::cout << std::fixed << std::setprecision(6);
    std::cout << "[ " << getInsertionRate() << " " << getDeletionRate() << " ]" << std::endl;
}

void ParameterIndelRates::reject(void) {

    insertionRate[0] = insertionRate[1];
    deletionRate[0]  = deletionRate[1];
    modelPtr->flipActiveLikelihood();
}

double ParameterIndelRates::update(int) {

    lastUpdateType.first = this;
    lastUpdateType.second = 0;

    // initialize some variables
    double window = 0.03;
    double theta1 = insertionRate[0];
    double theta2 = deletionRate[0];
    double lowerLimit = 0.0;
    double upperLimit = 10.0;

    // propose new values using rejection sampling
    double theta1Prime, theta2Prime;
    do {
        theta1Prime = theta1 + (rv->uniformRv() - 0.5) * window; // lambda (x)
        theta2Prime = theta2 + (rv->uniformRv() - 0.5) * window; // mu (y)
        } while(theta1Prime <= lowerLimit || theta1Prime >= theta2Prime || theta2Prime >= upperLimit);
    insertionRate[0] = theta1Prime;
    deletionRate[0]  = theta2Prime;
        
    // forward move proposal probability
    // get vertices
    std::pair<double, double> vertex1 = std::make_pair(theta1 - 0.5*window, theta2 - 0.5*window);
    std::pair<double, double> vertex2 = std::make_pair(theta1 - 0.5*window, theta2 + 0.5*window);
    std::pair<double, double> vertex3 = std::make_pair(theta1 + 0.5*window, theta2 + 0.5*window);
    std::pair<double, double> vertex4 = std::make_pair(theta1 + 0.5*window, theta2 - 0.5*window);
    
    // number of intersections
    int outOfBounds = 0;
    if (vertex1.first < lowerLimit)     // mu-axis intersects window
        outOfBounds += 2;               // vertex 2 necessarily is out of bounds if vertex 1 is
    if (vertex4.first > vertex4.second) // lambda=mu intersects window
        outOfBounds++;
    
    // probability of area
    double pForward = 0.0;
    double r;
    if (outOfBounds == 0)
        pForward = 1.0 / (pow(window,2.0));
    if (outOfBounds == 1)
        {
        // subtract triangular region defined by vertex 4
        r = pow(abs(vertex4.second - vertex4.first),2) / 2;
        pForward = 1.0 / (pow(window, 2.0) - r);
        }
    if (outOfBounds == 2)
        {
        // subtract rectangular region defined by vertex 1 & 2
        r = abs(vertex1.first - lowerLimit) * window;
        pForward = 1.0/(pow(window,2.0) - r);
        }
    if (outOfBounds == 3)
        {
        // subtract both triangular and rectuangular regions
        double r1 = pow(abs(vertex4.second - vertex4.first),2.0) * 0.5;
        double r2 = abs(vertex1.first) * window;
        if (vertex1.second >= lowerLimit)
            pForward = 1.0 / (pow(window,2.0) - r1 - r2);
        else
            {
            // both regions overlap
            double rOverlap = pow(abs(vertex1.second - lowerLimit),2) * 0.5;
            pForward = 1.0 / (pow(window,2.0) - r1 - r2 + rOverlap);
            }
        }
    
    // reverse move proposal probability
    // get vertices
    vertex1 = std::make_pair(theta1Prime - 0.5*window, theta2Prime - 0.5*window);
    vertex2 = std::make_pair(theta1Prime - 0.5*window, theta2Prime + 0.5*window);
    vertex3 = std::make_pair(theta1Prime + 0.5*window, theta2Prime + 0.5*window);
    vertex4 = std::make_pair(theta1Prime + 0.5*window, theta2Prime - 0.5*window);
    
    // number of intersections
    outOfBounds = 0;
    if (vertex1.first < lowerLimit)     // mu-axis intersects window
        outOfBounds += 2;               // vertex 2 necessarily is out of bounds if vertex 1 is
    if (vertex4.first > vertex4.second) // lambda=mu intersects window
        outOfBounds++;
    
    // probability of area
    double pReverse = 0.0;
    if (outOfBounds == 0)
        pReverse = 1.0/(pow(window,2.0));
    if (outOfBounds == 1)
        {
        // subtract triangular region defined by vertex 4
        r = pow(abs(vertex4.second - vertex4.first),2.0) * 0.5;
        pReverse = 1.0 / (pow(window,2.0) - r);
        }
    if (outOfBounds == 2)
        {
        // subtract rectangular region defined by vertex 1 & 2
        r = abs(vertex1.first - lowerLimit) * window;
        pReverse = 1.0 / (pow(window,2.0) - r);
        }
    if (outOfBounds == 3)
        {
        // subtract both triangular and rectuangular regions
        double r1 = pow(abs(vertex4.second - vertex4.first),2.0) * 0.5;
        double r2 = abs(vertex1.first) * window;
        if (vertex1.second >= lowerLimit)
            pReverse = 1.0 / (pow(window,2.0) - r1 - r2);
        else
            { // both regions overlap
            double rOverlap = pow(abs(vertex1.second - lowerLimit),2.0) * 0.5;
            pReverse = 1.0/(pow(window,2.0) - r1 - r2 + rOverlap);
            }
        }
        
    updateChangesTransitionProbabilities = false;
    modelPtr->setUpdateLikelihood();
    modelPtr->flipActiveLikelihood();
    
    return log(pReverse) - log(pForward);
}

double ParameterIndelRates::updateFromPrior(void) {

    lastUpdateType.first = this;
    lastUpdateType.second = 1;

    return 0.0;
}
