#include <iomanip>
#include <iostream>
#include "JphMatrix.hpp"



JphMatrix::JphMatrix(void) {
    
    nRows = 0;
    nCols = 0;
    m = NULL;
}

JphMatrix::JphMatrix(const JphMatrix& copy) {

    std::cout << "copying matrix through copy constructor" << std::endl;
    clone(copy);
}

JphMatrix::JphMatrix(size_t nr, size_t nc) {

    initialize(nr, nc);
}

JphMatrix::~JphMatrix(void) {

    if (m != NULL)
        {
        delete [] m[0];
        delete [] m;
        }
}

JphMatrix& JphMatrix::operator=(JphMatrix& rhs) {

    std::cout << "copying matrix through assignment operator" << std::endl;
    if (this != &rhs)
        clone(rhs);
    return *this;
}

void JphMatrix::clone(const JphMatrix& copy) {

    if (nRows != copy.nRows || nCols != copy.nCols)
        {
        if (m == NULL)
            {
            delete [] m[0];
            delete [] m;
            }
        initialize(copy.nRows, copy.nCols);
        }
    double* pLft = m[0];
    double* pRht = copy.m[0];
    for (int i=0; i<nRows*nCols; ++i)
        *(pLft++) = *(pRht++);
}

void JphMatrix::initialize(size_t nr, size_t nc) {

    nRows = nr;
    nCols = nc;
    m = new double*[nRows];
    m[0] = new double[nRows * nCols];
    for (int i=1; i<nRows; i++)
        m[i] = m[i-1] + nCols;
    for (size_t i=0; i<nRows; i++)
        for (size_t j=0; j<nCols; j++)
            m[i][j] = 0.0;
}

void JphMatrix::print(int prec) {

    std::cout << std::fixed << std::setprecision(prec);
    std::cout << this << std::endl;
    for (int i=0; i<nRows; i++)
        {
        double sum = 0.0;
        for (int j=0; j<nCols; j++)
            {
            sum += m[i][j];
            std::cout << m[i][j] << " ";
            }
        std::cout << " [" << sum << "]" << std::endl;
        }
}

void JphMatrix::print(std::string s, int prec) {

    std::cout << s << std::endl;
    print(prec);
}

void JphMatrix::printSummary(void) {

    std::cout << "Matrix address:      " << this << std::endl;
    std::cout << "  values address:    " << m << " (begin)" << std::endl;
    std::cout << "  values address:    " << &m[nRows*nCols-1] << " (end)" << std::endl;
    std::cout << "  Number of rows:    " << nRows << std::endl;
    std::cout << "  Number of columns: " << nCols << std::endl;
}
