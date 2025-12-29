#ifndef JphMatrix_hpp
#define JphMatrix_hpp

#include <iostream>
#include <string>

class JphMatrix {

    public:
                        explicit JphMatrix(void);
                        explicit JphMatrix(size_t nr, size_t nc);
                        explicit JphMatrix(const JphMatrix& copy);
                       ~JphMatrix(void);
        JphMatrix&      operator=(JphMatrix& rhs);
        double*         operator[](int i) { return m[i]; }
        const double*   operator[](int i) const { return m[i]; }
        size_t          getNumRows(void) { return nRows; }
        size_t          getNumRows(void) const { return nRows; }
        size_t          getNumCols(void) { return nCols; }
        void            initialize(size_t nr, size_t nc);
        void            print(int prec);
        void            print(std::string s, int prec);
        void            printSummary(void);
    
    private:
        void            clone(const JphMatrix& copy);
        double**        m;
        size_t          nRows;
        size_t          nCols;
};

typedef JphMatrix DoubleMatrix;

#endif
