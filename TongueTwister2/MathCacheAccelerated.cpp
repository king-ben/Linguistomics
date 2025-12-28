#include "MathCacheAccelerated.hpp"
#include <cmath>

static const double PADE_B[14] = {
    64764752532480000.0, 32382376266240000.0, 7771770303897600.0,
    1187353796428800.0, 129060195264000.0, 10559470521600.0,
    670442572800.0, 33522128640.0, 1323241920.0, 40840800.0,
    960960.0, 16380.0, 182.0, 1.0
};

static const double THETA_13 = 5.371920351148152;



MathCacheAccelerated::MathCacheAccelerated(void) : MathCache() {

}

MathCacheAccelerated::~MathCacheAccelerated(void) {

}

void MathCacheAccelerated::blasMultiply(DoubleMatrix& A, DoubleMatrix& B, DoubleMatrix& C) {

    int n = static_cast<int>(A.getNumRows());
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                n, n, n, 1.0, A.begin(), n, B.begin(), n, 0.0, C.begin(), n);
}

double MathCacheAccelerated::infinityNorm(DoubleMatrix& A) {

    size_t n = A.getNumRows();
    double maxSum = 0.0;
    
    for (size_t i = 0; i < n; i++)
        {
        double rowSum = 0.0;
        for (size_t j = 0; j < n; j++)
            rowSum += std::fabs(A(i, j));
        if (rowSum > maxSum)
            maxSum = rowSum;
        }
    return maxSum;
}

void MathCacheAccelerated::computeMatrixExponential(const DoubleMatrix& Q, double t, DoubleMatrix& P) {

    size_t n = Q.getNumRows();
    
    // allocate working matrices from base class buffer
    DoubleMatrix* A  = pushMatrix(n);
    DoubleMatrix* A2 = pushMatrix(n);
    DoubleMatrix* A4 = pushMatrix(n);
    DoubleMatrix* A6 = pushMatrix(n);
    DoubleMatrix* U  = pushMatrix(n);
    DoubleMatrix* V  = pushMatrix(n);
    DoubleMatrix* temp1 = pushMatrix(n);
    DoubleMatrix* temp2 = pushMatrix(n);
    
    // A = Q * t 
    Q.multiply(t, *A);
    
    // compute scaling factor using infinity norm
    double normA = infinityNorm(*A);
    int s = 0;
    if (normA > THETA_13)
        {
        s = static_cast<int>(std::ceil(std::log2(normA / THETA_13)));
        double scale = 1.0 / (1L << s);
        A->multiply(scale);
        }
    
    // compute A^2, A^4, A^6 using BLAS
    blasMultiply(*A, *A, *A2);
    blasMultiply(*A2, *A2, *A4);
    blasMultiply(*A2, *A4, *A6);
    
    // temp1 = b13*A6 + b11*A4 + b9*A2 + b7*I
    temp1->setIdentity();
    temp1->multiply(PADE_B[7]);
    for (size_t i = 0; i < n; i++)
        {
        for (size_t j = 0; j < n; j++)
            {
            (*temp1)(i,j) += PADE_B[9] * (*A2)(i,j) + PADE_B[11] * (*A4)(i,j) + PADE_B[13] * (*A6)(i,j);
            }
        }
    
    // temp2 = A6 * temp1
    blasMultiply(*A6, *temp1, *temp2);
    
    // temp2 += b5*A4 + b3*A2 + b1*I
    for (size_t i = 0; i < n; i++)
        {
        for (size_t j = 0; j < n; j++)
            {
            (*temp2)(i,j) += PADE_B[5] * (*A4)(i,j) + PADE_B[3] * (*A2)(i,j);
            }
        (*temp2)(i,i) += PADE_B[1];
        }
    
    // U = A * temp2
    blasMultiply(*A, *temp2, *U);
    
    // temp1 = b12*A6 + b10*A4 + b8*A2 + b6*I
    temp1->setIdentity();
    temp1->multiply(PADE_B[6]);
    for (size_t i = 0; i < n; i++)
        {
        for (size_t j = 0; j < n; j++)
            {
            (*temp1)(i,j) += PADE_B[8] * (*A2)(i,j) + PADE_B[10] * (*A4)(i,j) + PADE_B[12] * (*A6)(i,j);
            }
        }
    
    // V = A6 * temp1
    blasMultiply(*A6, *temp1, *V);
    
    // V += b4*A4 + b2*A2 + b0*I
    for (size_t i = 0; i < n; i++)
        {
        for (size_t j = 0; j < n; j++)
            {
            (*V)(i,j) += PADE_B[4] * (*A4)(i,j) + PADE_B[2] * (*A2)(i,j);
            }
        (*V)(i,i) += PADE_B[0];
        }
    
    // compute V-U and V+U for solve
    for (size_t i = 0; i < n; i++)
        {
        for (size_t j = 0; j < n; j++)
            {
            double u = (*U)(i,j);
            double v = (*V)(i,j);
            (*temp1)(i,j) = v - u;
            (*temp2)(i,j) = v + u;
            }
        }
    
    // solve (V - U) * P = (V + U) using base class method 
    gaussianElimination(*temp1, *temp2, P);
    
    // squaring phase: P = P^(2^s) 
    for (int sq = 0; sq < s; sq++)
        {
        blasMultiply(P, P, *temp1);
        P.copy(*temp1);
        }
    
    // absolute value cleanup
    for (size_t i = 0; i < n; i++)
        {
        for (size_t j = 0; j < n; j++)
            {
            P(i,j) = std::fabs(P(i,j));
            }
        }
    
    // release working matrices
    popMatrix(8);
}
