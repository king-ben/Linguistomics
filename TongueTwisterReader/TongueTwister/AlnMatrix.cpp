#include <iomanip>
#include "AlnMatrix.hpp"
#include "Msg.hpp"



AlnMatrix::AlnMatrix(int nr, int maxNs) : IntMatrix(nr, maxNs) {
    numSites = 0;
    numRows = nr;
    maxNumSites = maxNs;
}

AlnMatrix& AlnMatrix::operator=(const AlnMatrix& rhs) {

    if (this != &rhs)
        {
        numSites = rhs.numSites;
        maxNumSites = rhs.maxNumSites;
        IntMatrix::copy(rhs);
        }
    return *this;
}

bool AlnMatrix::operator==(const AlnMatrix& m) const {

    if (numRows != m.numRows || numCols != m.numCols || numSites != m.numSites)
        return false;
    return IntMatrix::operator==(m);
}

void AlnMatrix::print(void) {

    std::cout << "Address:    " << this << std::endl;
    std::cout << "Num. Taxa:  " << numRows << std::endl;
    std::cout << "Num. Sites: " << numSites << std::endl;
    for (int i=0; i<numRows; i++)
        {
        for (int j=0; j<numSites; j++)
            std::cout << std::setw(3) << (*this)(i,j);
        std::cout << std::endl;
        }
}

void AlnMatrix::setNumSites(int ns) {

    if (ns > maxNumSites)
        {
        maxNumSites *= 2;
        IntMatrix::create(numRows, maxNumSites);
        Msg::warning("Increasing the number of sites for an AlnMatrix");
        }
    numSites = ns;
}
