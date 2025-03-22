#pragma once
#include "MatrixView.hpp"

void recursiveMatMulImpl(MatrixView A, MatrixView B, MatrixView C)
{
    if (A.row_count() == 0 || A.col_count() == 0 || B.col_count() == 0) // Empty matrices
        return;

    int n = A.row_count(), m = A.col_count(), p = B.col_count();

    if (n == m && m == p && p == 1) // base case
    {
        C(0, 0) += A(0, 0) * B(0, 0);
        return;
    }

    MatrixView A11 = A.getSubMatrix(0, n / 2, 0, m / 2);
    MatrixView A12 = A.getSubMatrix(0, n / 2, m / 2, m);
    MatrixView A21 = A.getSubMatrix(n / 2, n, 0, m / 2);
    MatrixView A22 = A.getSubMatrix(n / 2, n, m / 2, m);

    MatrixView B11 = B.getSubMatrix(0, m / 2, 0, p / 2);
    MatrixView B12 = B.getSubMatrix(0, m / 2, p / 2, p);
    MatrixView B21 = B.getSubMatrix(m / 2, m, 0, p / 2);
    MatrixView B22 = B.getSubMatrix(m / 2, m, p / 2, p);

    MatrixView C11 = C.getSubMatrix(0, n / 2, 0, p / 2);
    MatrixView C12 = C.getSubMatrix(0, n / 2, p / 2, p);
    MatrixView C21 = C.getSubMatrix(n / 2, n, 0, p / 2);
    MatrixView C22 = C.getSubMatrix(n / 2, n, p / 2, p);

    recursiveMatMulImpl(A11, B11, C11);
    recursiveMatMulImpl(A12, B21, C11);
    recursiveMatMulImpl(A11, B12, C12);
    recursiveMatMulImpl(A12, B22, C12);
    recursiveMatMulImpl(A21, B11, C21);
    recursiveMatMulImpl(A22, B21, C21);
    recursiveMatMulImpl(A21, B12, C22);
    recursiveMatMulImpl(A22, B22, C22);
}

void recursiveMatMul(MatrixView A, MatrixView B, MatrixView C)
{
    C.clear();
    
    if (A.col_count() != B.row_count() || B.col_count() != C.col_count() || A.row_count() != C.row_count()) // Invalid matrix dimensions
        return;
        
    recursiveMatMulImpl(A, B, C);
}
