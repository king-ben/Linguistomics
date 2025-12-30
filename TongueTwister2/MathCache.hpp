#ifndef MathCache_hpp
#define MathCache_hpp

#include "Container.hpp"


class MathCache {

    public:
                        MathCache(void);
        virtual        ~MathCache(void);
        void            backSubstitutionRow(DoubleMatrix& U, double* b);
        void            computeLandU(DoubleMatrix& A, DoubleMatrix& L, DoubleMatrix& U);
        void            forwardSubstitutionRow(DoubleMatrix& L, double* b);
        void            gaussianElimination(DoubleMatrix& A, DoubleMatrix& B, DoubleMatrix& X);
        void            multiply(DoubleMatrix& A, DoubleMatrix& B);
        void            power(DoubleMatrix& m, int power);

        DoubleMatrix*   pushMatrix(size_t rows, size_t columns);
        DoubleMatrix*   pushMatrix(size_t size);
        void            popMatrix(size_t n);
        DoubleArray*    pushArray(size_t size);
        void            popArray();

        static const    size_t bufferSize = 16;

                        // place first so they're at known offset and likely in same cache line
        size_t          matrixCount;
        size_t          arrayCount;
        
                        // indexed into after count check
        DoubleMatrix    matrixBuffer[bufferSize];
        DoubleArray     arrayBuffer[bufferSize];
};

#endif
