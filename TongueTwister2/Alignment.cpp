#include <iomanip>
#include "Alignment.hpp"
#include "Msg.hpp"



Alignment::Alignment(const Alignment& a) : IntMatrix(a.numTaxa, a.maxNumSegments) {

    this->numTaxa = a.numTaxa;
    this->numSegments = a.numSegments;
    this->maxNumSegments = a.maxNumSegments;
    IntMatrix::copy(a);
}

Alignment::Alignment(size_t nr, size_t maxNs) : IntMatrix(nr, maxNs) {

    numTaxa = nr;
    numSegments = 0;
    maxNumSegments = maxNs;
}

Alignment& Alignment::operator=(const Alignment& rhs) {

    if (this != &rhs)
        {
        numTaxa = rhs.numTaxa;
        numSegments = rhs.numSegments;
        maxNumSegments = rhs.maxNumSegments;
        IntMatrix::copy(rhs);
        }
    return *this;
}

bool Alignment::operator==(const Alignment& m) const {

    if (numRows != m.numRows || numCols != m.numCols || numSegments != m.numSegments)
        return false;
    return IntMatrix::operator==(m);
}

void Alignment::print(void) {

    std::cout << "Address:       " << this << std::endl;
    std::cout << "Num. Taxa:     " << numRows << std::endl;
    std::cout << "Num. Segments: " << numSegments << std::endl;
    for (size_t i=0; i<numRows; i++)
        {
        for (size_t j=0; j<numSegments; j++)
            std::cout << std::setw(3) << (*this)(i,j);
        std::cout << std::endl;
        }
}

void Alignment::setNumSegments(size_t ns) {

    if (ns >= maxNumSegments)
        Msg::error("Exceeding maximum alignment length");
    numSegments = ns;
}
