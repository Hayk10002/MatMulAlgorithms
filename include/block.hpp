#pragma once

#include "getL1CacheSize.hpp"
#include "MatrixView.hpp"

void blockMatMul(MatrixView A, MatrixView B, MatrixView C, int block_size = sqrt(getL1CacheSize() / sizeof(int) / 4))
{
    if (A.col_count() != B.row_count() || B.col_count() != C.col_count() || A.row_count() != C.row_count()) // Invalid matrix dimensions
        return;

    C.clear();

    for (int i = 0; i < C.row_count(); i += block_size)
    {
        for (int j = 0; j < C.col_count(); j += block_size)
        {
            for (int k = 0; k < A.col_count(); k += block_size)
            {
                for (int i1 = i; i1 < i + block_size && i1 < C.row_count(); i1++)
                {
                    for (int j1 = j; j1 < j + block_size && j1 < C.col_count(); j1++)
                    {
                        for (int k1 = k; k1 < k + block_size && k1 < A.col_count(); k1++)
                        {
                            C(i1, j1) += A(i1, k1) * B(k1, j1);
                        }
                    }
                }
            }
        }
    }
}