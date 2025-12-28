#include <cmath>
#include "MathCache.hpp"
#include "Msg.hpp"
#include "ParameterFrequencies.hpp"
#include "RateMatrix.hpp"
#include "Ctmc.hpp"

#undef CHECK_TIPROBS



#pragma mark Base Class

Ctmc::Ctmc(size_t ns) : numStates(ns) {

}

bool Ctmc::checkTransitionProbabilities(DoubleMatrix* m, double thresshold) const {

    for (size_t i=0; i<numStates; i++)
        {
        double sum = 0.0;
        for (size_t j=0; j<numStates; j++)
            {
            sum += (*m)(i,j);
            if ((*m)(i,j) < 0.0)
                return false;
            }
        if (fabs(1.0 - sum) > thresshold)
            return false;
        }
    return true;
}

#pragma mark JC69

CtmcJC69::CtmcJC69(size_t ns) : Ctmc(ns) {

}

CtmcJC69::~CtmcJC69(void) {

}

void CtmcJC69::computeTransitionProbs(double branchLength, DoubleMatrix* output, MathCache*) const {

    const double n = static_cast<double>(numStates);
    const double oneOverN = 1.0 / n;
    const double expTerm = std::exp(-n * branchLength / (n - 1.0));
    
    double* p = output->begin();
    double* end = p + (numStates * numStates);

    // fill entire matrix (compiler should do something clever here)
    const double pDiff = oneOverN * (1.0 - expTerm);
    while (p < end)
        *p++ = pDiff;

    // walk the diagonal
    const double pSame = oneOverN + (1.0 - oneOverN) * expTerm;
    p = output->begin();
    size_t stride = numStates + 1;
    for (size_t i=0; i<numStates; i++)
        {
        *p = pSame;
        p += stride;
        }
}

#pragma mark F81

CtmcF81::CtmcF81(size_t ns, ParameterFrequencies* f) : Ctmc(ns), freqsParm(f) {

}

CtmcF81::~CtmcF81(void) {

}

void CtmcF81::computeTransitionProbs(double branchLength, DoubleMatrix* output, MathCache*) const {

    std::vector<double>& freqs = freqsParm->getFrequencies();
    double* f = freqs.data();
    double sumSq = 0.0;
    for (size_t i=0; i<numStates; i++)
        {
        sumSq += (*f) * (*f);
        f++;
        } 
    const double mu = 1.0 / (1.0 - sumSq);
    const double expTerm = std::exp(-mu * branchLength);
    const double oneMinusExpTerm = 1.0 - expTerm;

    double* p = output->begin();
    f = freqs.data();
    for (size_t i=0; i<numStates; i++)
        {
        for (size_t j=0; j<numStates; j++)
            {
            if (i != j)
                (*p) = oneMinusExpTerm * f[j];
            else
                (*p) = expTerm + oneMinusExpTerm * f[j];
            p++;
            }
        }
}

#pragma mark General Model

CtmcGeneral::CtmcGeneral(size_t ns, RateMatrix* r) : Ctmc(ns), rateMatrix(r) {

    // initialize the factorial table for the setQvalue function
    largestFactorial = 50;
    factorialTable = new double[largestFactorial+1];
    factorialTable[0] = 1.0;
    for (size_t i=1; i<=largestFactorial; i++)
        factorialTable[i] = factorialTable[i-1] * i;
}

CtmcGeneral::~CtmcGeneral(void) {

    delete [] factorialTable;
}

MathCache* CtmcGeneral::allocateMathCache(void) {

    return new MathCache;
}

/* The method approximates the matrix exponential, P = e^A, using
   the algorithm 11.3.1, described in:

   Golub, G. H., and C. F. Van Loan. 1996. Matrix Computations, Third Edition.
      The Johns Hopkins University Press, Baltimore, Maryland.

   The method has the advantage of error control. The error is controlled by
   setting qValue appropriately (using the function SetQValue). 
   
   The scaling-and-squaring approach:
   1. Choose j such that ||A||/2^j < 0.5 (or some threshold)
   2. Compute R = pade(A / 2^j) ≈ exp(A / 2^j)
   3. Square R exactly j times: exp(A) = R^(2^j) */
void CtmcGeneral::computeTransitionProbs(double branchLength, DoubleMatrix* output, MathCache* cache) const {

    // check for invalid branch length
    if (branchLength < 0.0)
        Msg::error("Cannot compute transition probabilities for negative branch length: " + std::to_string(branchLength));

    DoubleMatrix& Q = rateMatrix->getQ();
    
    size_t size = output->getNumRows();
    DoubleMatrix* d  = cache->pushMatrix(size);
    DoubleMatrix* n  = cache->pushMatrix(size);
    DoubleMatrix* a  = cache->pushMatrix(size);
    DoubleMatrix* x  = cache->pushMatrix(size);
    DoubleMatrix* cx = cache->pushMatrix(size);

    // A is the matrix Q * v and p = exp(a)
    Q.multiply(branchLength, *a);

    // set up identity matrices
    d->setIdentity();
    n->setIdentity();
    x->setIdentity();

    // Compute the infinity norm (maximum absolute row sum) for proper scaling.
    // This replaces the incorrect maxDiagonal() which returns 0 for rate matrices
    // because all diagonal elements are negative.
    double infinityNorm = 0.0;
    for (size_t i = 0; i < size; i++)
        {
        double rowSum = 0.0;
        for (size_t jj = 0; jj < size; jj++)
            rowSum += fabs((*a)(i, jj));
        if (rowSum > infinityNorm)
            infinityNorm = rowSum;
        }
    
    int y = logBase2Plus1(infinityNorm);
    int j = y < 0 ? 0 : y;

    a->divideByPowerOfTwo(j);
    
    double c = 1.0;
    int qValue = setQvalue(1e-7);  // Note: using 1e-7 for better accuracy
    for (int k=1; k<=qValue; k++)
        {
        c = c * (qValue - k + 1.0) / ((2.0 * qValue - k + 1.0) * k);

        // X = AX
        cache->multiply(*a, *x);

        // N = N + cX
        x->multiply(c, *cx);
        n->add(*cx);

        // D = D + (-1)^k*cX 
        int negativeFactor = (k % 2 == 0 ? 1 : -1);
        if (negativeFactor == -1)
            cx->multiply(negativeFactor);
        d->add(*cx);
        }

    cache->gaussianElimination(*d, *n, *output);
    
    // Square the result j times to get exp(A) = R^(2^j)
    // Note: The original code incorrectly used (j+1) instead of 2^j
    if (j > 0)
        cache->power(*output, 1 << j);

    for (auto p = output->begin(), end = output->end(); p < end; ++p)
        *p = fabs(*p);

    cache->popMatrix(5);
   
#   if defined(CHECK_TIPROBS)
    if (checkTransitionProbabilities(output, 0.00001) == false)
        Msg::error("Problem calculating transition probabilities");
#   endif
}

int CtmcGeneral::logBase2Plus1(double x) const {

	int j = 0;
	while(x > 1.0 - 1.0e-07)
		{
		x /= 2.0;
		j++;
		}
	return (j);
}

int CtmcGeneral::setQvalue(double tol) const {
	
	double x = pow(2.0, 3.0 - (0 + 0)) * factorial(0) * factorial(0) / (factorial(0+0) * factorial(0+0+1));
	int qV = 0;
	while (x > tol)
		{
		qV++;
		x = pow(2.0, 3.0 - (qV + qV)) * factorial(qV) * factorial(qV) / (factorial(qV+qV) * factorial(qV+qV+1));
		}
	return (qV);
}

double CtmcGeneral::factorial(size_t x) const {

	if (x > largestFactorial)
        Msg::error("Factorial exceeds limits");
    return factorialTable[x];
}
