#include <cmath>
#include <iomanip>
#include <iostream>
#include "Alignment.hpp"
#include "MatrixMath.hpp"
#include "Model.hpp"
#include "Msg.hpp"
#include "Node.hpp"
#include "ParameterTree.hpp"
#include "RateMatrix.hpp"
#include "Threads.hpp"
#include "TransitionProbabilities.hpp"
#include "Tree.hpp"
#include "UserSettings.hpp"

double factorial(int x);
int logBase2Plus1(double x);
int setQvalue(double tol);



class TransitionProbabilitiesTask : public ThreadTask {

    public:
    
        TransitionProbabilitiesTask() {
            numStates  = 0;
            Tree       = NULL;
            Q          = NULL;
            P          = NULL;
        }
        
        void init(Tree* tree, DoubleMatrix* q, TransitionProbabilitiesInfo& info, int activeProbs) {
        
            numStates  = info.numStates;
            Tree       = tree;
            Q          = q;
            P          = info.probs[activeProbs];
        }

        /* The method approximates the matrix exponential, P = e^A, using
           the algorithm 11.3.1, described in:

           Golub, G. H., and C. F. Van Loan. 1996. Matrix Computations, Third Edition.
              The Johns Hopkins University Press, Baltimore, Maryland.

           The method has the advantage of error control. The error is controlled by
           setting qValue appropriately (using the function SetQValue). */
        void computeMatrixExponential(MathCache& cache, int qValue, double v, DoubleMatrix* probs) {
        
            assert(probs->getNumRows() == probs->getNumCols());
            auto size = probs->getNumRows();

            auto d  = cache.pushMatrix(size);
            auto n  = cache.pushMatrix(size);
            auto a  = cache.pushMatrix(size);
            auto x  = cache.pushMatrix(size);
            auto cx = cache.pushMatrix(size);

            // A is the matrix Q * v and p = exp(a)
            Q->multiply(v, *a);

            // set up identity matrices
            d->setIdentity();
            n->setIdentity();
            x->setIdentity();

            auto maxAValue = a->maxDiagonal();
            int y = logBase2Plus1(maxAValue);
            int j = y < 0 ? 0 : y;

            a->divideByPowerOfTwo(j);
            
            double c = 1.0;
            for (int k = 1; k <= qValue; k++)
                {
                c = c * (qValue - k + 1.0) / ((2.0 * qValue - k + 1.0) * k);

                /* X = AX */
                cache.multiply(*a, *x);

                /* N = N + cX */
                x->multiply(c, *cx);
                n->add(*cx);

                /* D = D + (-1)^k*cX */
                int negativeFactor = (k % 2 == 0 ? 1 : -1);
                if (negativeFactor == -1)
                    cx->multiply(negativeFactor);
                d->add(*cx);
                }

            cache.gaussianElimination(*d, *n, *probs);
            if (j > 0)
              cache.power(*probs, j+1);

            for (auto p = probs->begin(), end = probs->end(); p < end; ++p)
                *p = fabs(*p);

            cache.popMatrix(5);
        }

        virtual void Run(MathCache& cache) {
        
            int qValue = setQvalue(10e-7);
            std::vector<Node*>& traversalSeq = Tree->getDownPassSequence();
            for (int n = 0; n < traversalSeq.size(); n++)
                {
                Node* p = traversalSeq[n];
                if (p->getAncestor() != NULL)
                    {
                    int pIdx = p->getIndex();
                    double v = p->getBranchLength();
                    computeMatrixExponential(cache, qValue, v, P[pIdx]);
                    }
                }
            }

    private:
    
        int             numStates;
        Tree*           Tree;
        DoubleMatrix*   Q;
        DoubleMatrix**  P;
};



TransitionProbabilities::TransitionProbabilities(void) {

    modelPtr          = NULL;
    numNodes          = 0;
    numRateCategories = 1;
    numStates         = 0;
    substitutionModel = jc69;
    threadPool        = NULL;
    isInitialized     = false;
    needsUpdate       = true;
    activeProbs       = 0;
}

TransitionProbabilities::~TransitionProbabilities(void) {

    for (std::map<RbBitSet,TransitionProbabilitiesInfo>::iterator it = transProbs.begin(); it != transProbs.end(); it++)
        {
        for (int s=0; s<2; s++)
            delete [] it->second.probs[s];
        }
}

bool TransitionProbabilities::areTransitionProbabilitiesValid(double tolerance) {

    bool allGood = true;
    for (std::map<RbBitSet,TransitionProbabilitiesInfo>::iterator it = transProbs.begin(); it != transProbs.end(); it++)
        {
        for (int n=0; n<it->second.numMatrices; n++)
            {
            for (int i=0; i<numStates; i++)
                {
                double sum = 0.0;
                for (int j=0; j<numStates; j++)
                    sum += (*it->second.probs[activeProbs][n])(i,j);
                if (fabs(1.0 - sum) > tolerance)
                    allGood = false;
                }
            }
        }
    return allGood;
}

void TransitionProbabilities::flipActive(void) {

    activeProbs ^= 1;
}

DoubleMatrix** TransitionProbabilities::getTransitionProbabilities(RbBitSet& bs) {

    std::map<RbBitSet,TransitionProbabilitiesInfo>::iterator it = transProbs.find(bs);
    if (it == transProbs.end())
        Msg::error("Could not find transition probability vector for mask " + bs.bitString());
    return it->second.probs[activeProbs];
}

DoubleMatrix& TransitionProbabilities::getTransitionProbabilities(RbBitSet& bs, int nodeIdx) {

    std::map<RbBitSet,TransitionProbabilitiesInfo>::iterator it = transProbs.find(bs);
    if (it == transProbs.end())
        Msg::error("Could not find transition probability vector for mask " + bs.bitString());
    return *it->second.probs[activeProbs][nodeIdx];
}

void TransitionProbabilities::initialize(Model* m, ThreadPool* pool, std::vector<Alignment*>& alns, int nn, int ns, int sm) {

    if (isInitialized == true)
        {
        Msg::warning("Transition probabilities can only be initialized once");
        return;
        }
                
    UserSettings& settings = UserSettings::userSettings();
    modelPtr = m;
    threadPool = pool;
    numNodes = nn;
    numStates = ns;
    substitutionModel = sm;
    numRateCategories = settings.getNumRateCategories();
    std::cout << "   * Number of states = " << numStates << std::endl;
    std::cout << "   * Number of gamma rate categories = " << numRateCategories << std::endl;

    for (int i=0; i<alns.size(); i++)
        {
        // get the taxon mask
        std::vector<bool>& bm = alns[i]->taxonMask;
        RbBitSet mask(bm);

        // if the mask is not found in the map, insert it
        std::map<RbBitSet,TransitionProbabilitiesInfo>::iterator it = transProbs.find(mask);
        if (it == transProbs.end())
            {
            Tree* t = modelPtr->getTree(mask);
            if (t == NULL)
                Msg::error("Could not find tree for mask " + mask.bitString() + " when initializing transition probabilities");
            int nNodes = t->getNumNodes();
            TransitionProbabilitiesInfo info;
            info.numMatrices = nNodes;
            info.numStates = numStates;

            for (int s=0; s<2; s++)
                {
                info.probs[s] = new DoubleMatrix*[nNodes];
                for (int mi=0; mi<nNodes; mi++)
                    {
                    info.probs[s][mi] = new DoubleMatrix(numStates, numStates);
                    info.probs[s][mi]->setIdentity();
                    }
                }

            transProbs.insert( std::make_pair(mask,info) );
            }
        }
        
    stationaryFreqs[0].resize(numStates);
    stationaryFreqs[1].resize(numStates);
        
    isInitialized = true;
    
    //print();
}

void TransitionProbabilities::print(void) {

    std::cout << std::fixed << std::setprecision(5);
    for (std::map<RbBitSet,TransitionProbabilitiesInfo>::iterator it = transProbs.begin(); it != transProbs.end(); it++)
        {
        for (int n=0; n<it->second.numMatrices; n++)
            {
            std::cout << "Transition probabilities for node " << n << " (" << it->first.bitString() << ")" << std::endl;
            it->second.probs[activeProbs][n]->print();
            }
        }
}

void TransitionProbabilities::getStationaryFrequencies(std::vector<double>& f) {

    for (int i=0; i<numStates; i++)
        f[i] = stationaryFreqs[activeProbs][i];
}

void TransitionProbabilities::setTransitionProbabilities(void) {

    if (needsUpdate == false)
        return;
        
    if (substitutionModel == jc69)
        {
        // the transition probabilities can be calculated analytically (and quickly)
        setTransitionProbabilitiesJc69();
        }
    else
        {
        // the transition probabilities are calculated using either the Eigen
        // system of the rate matrix or using the Pade approximation
        setTransitionProbabilitiesUsingPadeMethod();
        }
        
    needsUpdate = false;
}

void TransitionProbabilities::setTransitionProbabilitiesJc69(void) {

    // calculate transition probabilities under the Jukes-Cantor (1969) model
    for (std::map<RbBitSet,TransitionProbabilitiesInfo>::iterator it = transProbs.begin(); it != transProbs.end(); it++)
        {
        Tree* t = modelPtr->getTree(it->first);
        if (t == NULL)
            Msg::error("Could not find tree for mask " + it->first.bitString());
                
        DoubleMatrix** probs = it->second.probs[activeProbs];
        
        std::vector<Node*>& traversalSeq = t->getDownPassSequence();
        for (int n=0; n<traversalSeq.size(); n++)
            {
            Node* p = traversalSeq[n];
            DoubleMatrix& tp = *probs[p->getIndex()];
            double v = p->getBranchLength();
            
            double x = -((double)numStates/(numStates-1));
            double pChange = (1.0/numStates) - (1.0/numStates) * exp(x * v);
            double pNoChange = (1.0/numStates) + ((double)(numStates-1)/numStates) * exp(x * v);
            tp.fill(pChange);
            for (int i=0; i<numStates; i++)
                tp(i,i) = pNoChange;
                
            }
        
        }
            
    double sf = 1.0 / numStates;
    for (int i=0; i<numStates; i++)
        stationaryFreqs[activeProbs][i] = sf;
}

void TransitionProbabilities::setTransitionProbabilitiesUsingPadeMethod(void) {

    RateMatrix* rmat = modelPtr->getRateMatrix();
    DoubleMatrix& Q = rmat->getRateMatrix();

    auto tasks = new TransitionProbabilitiesTask[transProbs.size()];
    auto task = tasks;
    
    // update the main tree
    for (auto it = transProbs.begin(); it != transProbs.end(); it++)
        {
        Tree* t = modelPtr->getTree(it->first);
        if (t == NULL)
            Msg::error("Could not find tree for mask " + it->first.bitString());
        task->init(t, &Q, it->second, activeProbs);
        threadPool->PushTask(task);
        ++task;
        }
        
    threadPool->Wait();
    delete[] tasks;
        
    auto& rmatFreqs = rmat->getEquilibriumFrequencies();
    for (int i=0; i<numStates; i++)
        stationaryFreqs[activeProbs][i] = rmatFreqs[i];
}

int logBase2Plus1(double x) {

	int j = 0;
	while(x > 1.0 - 1.0e-07)
		{
		x /= 2.0;
		j++;
		}
	return (j);
}

int setQvalue(double tol) {
	
	double x = pow(2.0, 3.0 - (0 + 0)) * factorial(0) * factorial(0) / (factorial(0+0) * factorial(0+0+1));
	int qV = 0;
	while (x > tol)
		{
		qV++;
		x = pow(2.0, 3.0 - (qV + qV)) * factorial(qV) * factorial(qV) / (factorial(qV+qV) * factorial(qV+qV+1));
		}
	return (qV);
}

double factorial(int x) {
	
	double fac = 1.0;
	for (int i=0; i<x; i++)
		fac *= (i+1);
	return (fac);
}
