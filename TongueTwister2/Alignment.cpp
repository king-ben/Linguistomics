#include <iomanip>
#include "Alignment.hpp"
#include "Msg.hpp"



Alignment::Alignment(const Alignment& a) : IntMatrix(a.numTaxa, a.maxNumSites) {

    this->numTaxa = a.numTaxa;
    this->numSites = a.numSites;
    this->maxNumSites = a.maxNumSites;
    IntMatrix::copy(a);
}

Alignment::Alignment(size_t nr, size_t maxNs) : IntMatrix(nr, maxNs) {

    numTaxa = nr;
    numSites = 0;
    maxNumSites = maxNs;
}

Alignment& Alignment::operator=(const Alignment& rhs) {

    if (this != &rhs)
        {
        numTaxa = rhs.numTaxa;
        numSites = rhs.numSites;
        maxNumSites = rhs.maxNumSites;
        IntMatrix::copy(rhs);
        }
    return *this;
}

bool Alignment::operator==(const Alignment& m) const {

    if (numRows != m.numRows || numCols != m.numCols || numSites != m.numSites)
        return false;
    return IntMatrix::operator==(m);
}

void Alignment::print(void) {

    std::cout << "Address:    " << this << std::endl;
    std::cout << "Num. Taxa:  " << numRows << std::endl;
    std::cout << "Num. Sites: " << numSites << std::endl;
    for (size_t i=0; i<numRows; i++)
        {
        for (size_t j=0; j<numSites; j++)
            std::cout << std::setw(3) << (*this)(i,j);
        std::cout << std::endl;
        }
}

void Alignment::setNumSites(size_t ns) {

    if (ns >= maxNumSites)
        Msg::error("Exceeding maximum alignment length");
    numSites = ns;
}
