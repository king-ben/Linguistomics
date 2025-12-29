#include <iomanip>
#include "SiteLikelihood.hpp"



SiteLikelihood::SiteLikelihood(int nn, int ns) {

    // set dimensions
    numNodes = nn;
    numStates = ns;
    numStates1 = numStates+1;
        
    probsH = new double*[numNodes];
    probsH[0] = new double[numNodes * numStates];
    probsI = new double* [numNodes];
    probsI[0] = new double[numNodes * numStates1];
    for (int i=1; i<numNodes; i++)
        {
        probsH[i] = probsH[i-1] + numStates;
        probsI[i] = probsI[i-1] + numStates1;
        zeroOutH(i);
        zeroOutI(i);
        }
}

SiteLikelihood::~SiteLikelihood(void) {

    delete [] probsH[0];
    delete [] probsH;
    delete [] probsI[0];
    delete [] probsI;
}

void SiteLikelihood::print(void) {

    std::cout << std::fixed << std::setprecision(0);
    
    std::cout << "probsH vector:" << std::endl;
    for (int i=0; i<numNodes; i++)
        {
        std::cout << std::setw(3) << i << " -- ";
        for (int j=0; j<numStates; j++)
            std::cout << probsH[i][j] << " ";
        std::cout << std::endl;
        }

    std::cout << "probsI vector:" << std::endl;
    for (int i=0; i<numNodes; i++)
        {
        std::cout << std::setw(3) << i << " -- ";
        for (int j=0; j<numStates1; j++)
            std::cout << probsI[i][j] << " ";
        std::cout << std::endl;
        }
}

void SiteLikelihood::zeroOutH(void) {

    for (int i=0; i<numNodes; i++)
        memset( probsH[i], 0, numStates*sizeof(double) );
}

void SiteLikelihood::zeroOutI(void) {

    for (int i=0; i<numNodes; i++)
        memset( probsI[i], 0, numStates1*sizeof(double) );
}

void SiteLikelihood::zeroOutH(int idx) {

    memset( probsH[idx], 0, numStates*sizeof(double) );
}

void SiteLikelihood::zeroOutI(int idx) {

    memset( probsI[idx], 0, numStates1*sizeof(double) );
}
