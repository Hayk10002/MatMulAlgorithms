#include "MatrixView.hpp"

void naiveMatMul(MatrixView A, MatrixView B, MatrixView C)
{   
    int n = A.row_count(), m = A.col_count(), p = B.col_count();
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < p; j++)
        {
            C(i, j) = 0;
            for (int k = 0; k < m; k++)
                C(i, j) += A(i, k) * B(k, j);
        }
    }
}

void naiveCacheFriendlyMatMul(MatrixView A, MatrixView B, MatrixView C)
{
    C.clear();
    int n = A.row_count(), m = A.col_count(), p = B.col_count();
    for (int i = 0; i < n; i++)
        for (int k = 0; k < m; k++)
            for (int j = 0; j < p; j++)
                C(i, j) += A(i, k) * B(k, j);
}