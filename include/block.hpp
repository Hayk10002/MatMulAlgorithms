#pragma once

#include "getL1CacheSize.hpp"
#include "MatrixView.hpp"

void blockMatMul(MatrixView A, MatrixView B, MatrixView C, MatMulMode mode, int block_size)
{
    if (A.col_count() != B.row_count() || B.col_count() != C.col_count() || A.row_count() != C.row_count()) // Invalid matrix dimensions
        return;

    if (mode == MatMulMode::Overwrite)
        C.clear();
    
    int n = A.row_count(), m = A.col_count(), p = B.col_count();
    
    for (int i = 0; i < n; i += block_size)
        for (int k = 0; k < m; k += block_size)
            for (int j = 0; j < p; j += block_size)
                for (int i1 = i; i1 < i + block_size && i1 < n; i1++)
                    for (int k1 = k; k1 < k + block_size && k1 < m; k1++)
                        for (int j1 = j; j1 < j + block_size && j1 < p; j1++)
                            C(i1, j1) += A(i1, k1) * B(k1, j1);
}

void cacheFriendlyBlockMatMul(MatrixView A, MatrixView B, MatrixView C, MatMulMode mode)
{
    blockMatMul(A, B, C, mode, sqrt(getL1CacheSize() / sizeof(int) / 4));
}