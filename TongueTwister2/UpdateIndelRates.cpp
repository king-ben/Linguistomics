#include <cmath>
#include "ParameterIndelRates.hpp"
#include "Probability.hpp"
#include "RandomVariable.hpp"
#include "UpdateIndelRates.hpp"



UpdateIndelRates::UpdateIndelRates(Model* m, RandomVariable* r, ParameterIndelRates* p) : Update(m, r), myParm(p) {

    updateId = hashUpdateName(getUpdateName());
    tuningParameter = 0.03;
}

void UpdateIndelRates::setDependants(void) {

    clearDependencyFlags();
    
    updatedParameter      = myParm;
    rateMatrixNeedsUpdate = false;
    allTiprobsNeedUpdate  = false;
    singleBranchChanged   = false;
}

double UpdateIndelRates::update(void) {

    // initialize some variables
    double window = tuningParameter;
    double theta1 = myParm->getInsertionRate();
    double theta2 = myParm->getDeletionRate();
    double lowerLimit = 0.0;
    double upperLimit = 10.0;

    // propose new values using rejection sampling
    double theta1Prime, theta2Prime;
    do {
        theta1Prime = theta1 + (rng->uniformRv() - 0.5) * window; // lambda (x)
        theta2Prime = theta2 + (rng->uniformRv() - 0.5) * window; // mu (y)
    } while(theta1Prime <= lowerLimit || theta1Prime >= theta2Prime || theta2Prime >= upperLimit);
    myParm->setInsertionRate(theta1Prime);
    myParm->setDeletionRate(theta2Prime);
        
    // forward move proposal probability
    // get vertices
    std::pair<double, double> vertex1 = std::make_pair(theta1 - 0.5*window, theta2 - 0.5*window);
    std::pair<double, double> vertex2 = std::make_pair(theta1 - 0.5*window, theta2 + 0.5*window);
    std::pair<double, double> vertex3 = std::make_pair(theta1 + 0.5*window, theta2 + 0.5*window);
    std::pair<double, double> vertex4 = std::make_pair(theta1 + 0.5*window, theta2 - 0.5*window);
    (void)vertex2; (void)vertex3; // suppress unused warnings
    
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
        r = pow(fabs(vertex4.second - vertex4.first),2) / 2;
        pForward = 1.0 / (pow(window, 2.0) - r);
        }
    if (outOfBounds == 2)
        {
        // subtract rectangular region defined by vertex 1 & 2
        r = fabs(vertex1.first - lowerLimit) * window;
        pForward = 1.0/(pow(window,2.0) - r);
        }
    if (outOfBounds == 3)
        {
        // subtract both triangular and rectangular regions
        double r1 = pow(fabs(vertex4.second - vertex4.first),2.0) * 0.5;
        double r2 = fabs(vertex1.first) * window;
        if (vertex1.second >= lowerLimit)
            pForward = 1.0 / (pow(window,2.0) - r1 - r2);
        else
            {
            // both regions overlap
            double rOverlap = pow(fabs(vertex1.second - lowerLimit),2) * 0.5;
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
        r = pow(fabs(vertex4.second - vertex4.first),2.0) * 0.5;
        pReverse = 1.0 / (pow(window,2.0) - r);
        }
    if (outOfBounds == 2)
        {
        // subtract rectangular region defined by vertex 1 & 2
        r = fabs(vertex1.first - lowerLimit) * window;
        pReverse = 1.0 / (pow(window,2.0) - r);
        }
    if (outOfBounds == 3)
        {
        // subtract both triangular and rectangular regions
        double r1 = pow(fabs(vertex4.second - vertex4.first),2.0) * 0.5;
        double r2 = fabs(vertex1.first) * window;
        if (vertex1.second >= lowerLimit)
            pReverse = 1.0 / (pow(window,2.0) - r1 - r2);
        else
            { // both regions overlap
            double rOverlap = pow(fabs(vertex1.second - lowerLimit),2.0) * 0.5;
            pReverse = 1.0/(pow(window,2.0) - r1 - r2 + rOverlap);
            }
        }
        
    setDependants();
           
    return log(pReverse) - log(pForward);
}

double UpdateIndelRates::update(double power) {

    if (rng->uniformRv() < priorSampleProb(power))
        return updateFromPrior();
    return update();
}

double UpdateIndelRates::updateFromPrior(void) {

    // get parameters of prior
    double insertionLambda = myParm->getDeletionLambda();
    double deletionLambda = myParm->getDeletionLambda();

    // propose new insertion (lambda) and deletion (mu) rates
    double lambda, mu;
    do {
        lambda = Probability::Exponential::rv(rng, insertionLambda);
        mu     = Probability::Exponential::rv(rng, deletionLambda);
        } while(lambda >= mu);
        
    // calculate the proposal probability
    double lnForwardProposal = myParm->lnPriorProbability(lambda, mu);
    double lnReverseProposal = myParm->lnPriorProbability(myParm->getInsertionRate(), myParm->getDeletionRate());
        
    // set the parameters
    myParm->setInsertionRate(lambda);
    myParm->setDeletionRate(mu);

    // and update dependents
    setDependants();

    return lnReverseProposal - lnForwardProposal;
}
