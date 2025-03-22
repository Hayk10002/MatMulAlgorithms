#include "MatrixView.hpp"

void naiveMatMul(MatrixView A, MatrixView B, MatrixView C)
{
    for (int i = 0; i < C.row_count(); i++)
        for (int j = 0; j < C.col_count(); j++)
        {
            C(i, j) = 0;
            for (int k = 0; k < A.col_count(); k++)
                C(i, j) += A(i, k) * B(k, j);
        }
}

void betterIterativeMatMul(MatrixView A, MatrixView B, MatrixView C)
{
    C.clear();
    for (int i = 0; i < C.row_count(); i++)
        for (int k = 0; k < A.col_count(); k++)
            for (int j = 0; j < C.col_count(); j++)
                C(i, j) += A(i, k) * B(k, j);
}