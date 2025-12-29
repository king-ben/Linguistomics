#include <iostream>
#include "IndelMatrix.hpp"
#include "Msg.hpp"


IndelMatrix::IndelMatrix(int nt, int maxNs) : IntMatrix(maxNs, nt) {
    numSites = 0;
    numTaxa = nt;
    maxNumSites = maxNs;
}

void IndelMatrix::print(void) {

    std::cout << this << std::endl;
    for (int i=0; i<numSites; i++)
        {
        for (int j=0; j<numTaxa; j++)
            std::cout << (*this)(i,j);
        std::cout << std::endl;
        }
}

void IndelMatrix::setNumSites(int ns) {

    if (ns > maxNumSites)
        {
        maxNumSites *= 2;
        IntMatrix::create(maxNumSites, numTaxa);
        Msg::warning("Increasing the number of sites for an IndelMatrix");
        }
    numSites = ns;
}
