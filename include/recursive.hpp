#pragma once
#include "MatrixView.hpp"

#include <vector>

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

    MatrixView A11 = A.getSubMatrix(0    , n / 2, 0    , m / 2);
    MatrixView A12 = A.getSubMatrix(0    , n / 2, m / 2, m    );
    MatrixView A21 = A.getSubMatrix(n / 2, n    , 0    , m / 2);
    MatrixView A22 = A.getSubMatrix(n / 2, n    , m / 2, m    );

    MatrixView B11 = B.getSubMatrix(0    , m / 2, 0    , p / 2);
    MatrixView B12 = B.getSubMatrix(0    , m / 2, p / 2, p    );
    MatrixView B21 = B.getSubMatrix(m / 2, m    , 0    , p / 2);
    MatrixView B22 = B.getSubMatrix(m / 2, m    , p / 2, p    );

    MatrixView C11 = C.getSubMatrix(0    , n / 2, 0    , p / 2);
    MatrixView C12 = C.getSubMatrix(0    , n / 2, p / 2, p    );
    MatrixView C21 = C.getSubMatrix(n / 2, n    , 0    , p / 2);
    MatrixView C22 = C.getSubMatrix(n / 2, n    , p / 2, p    );

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

void StrassenMatMul(MatrixView A, MatrixView B, MatrixView C)
{
    int s = A.row_count();
    if ((s & (s - 1)) != 0 ||
        A.col_count() != s || B.row_count() != s || B.col_count() != s || C.row_count() != s || C.col_count() != s) // Invalid matrix dimensions
        return;

    if (s == 1)
    {
        C(0, 0) = A(0, 0) * B(0, 0);
        return;
    }

    static std::vector<int> buffer{};
    static int current_buffer_start = 0;

    int potential_new_size = 5 * (s * s - 1) / 3; 
    buffer.reserve(potential_new_size);
    if (buffer.size() < potential_new_size) 
        buffer.resize(potential_new_size);

    MatrixView A11 = A.getSubMatrix(0    , s / 2, 0    , s / 2);
    MatrixView A12 = A.getSubMatrix(0    , s / 2, s / 2, s    );
    MatrixView A21 = A.getSubMatrix(s / 2, s    , 0    , s / 2);
    MatrixView A22 = A.getSubMatrix(s / 2, s    , s / 2, s    );

    MatrixView B11 = B.getSubMatrix(0    , s / 2, 0    , s / 2);
    MatrixView B12 = B.getSubMatrix(0    , s / 2, s / 2, s    );
    MatrixView B21 = B.getSubMatrix(s / 2, s    , 0    , s / 2);
    MatrixView B22 = B.getSubMatrix(s / 2, s    , s / 2, s    );

    MatrixView C11 = C.getSubMatrix(0    , s / 2, 0    , s / 2);
    MatrixView C12 = C.getSubMatrix(0    , s / 2, s / 2, s    );
    MatrixView C21 = C.getSubMatrix(s / 2, s    , 0    , s / 2);
    MatrixView C22 = C.getSubMatrix(s / 2, s    , s / 2, s    );

    MatrixView x({buffer.begin() + current_buffer_start + 0 * s * s / 4, size_t(s * s / 4)}, s / 2);
    MatrixView y({buffer.begin() + current_buffer_start + 1 * s * s / 4, size_t(s * s / 4)}, s / 2);
    MatrixView u({buffer.begin() + current_buffer_start + 2 * s * s / 4, size_t(s * s / 4)}, s / 2);
    MatrixView v({buffer.begin() + current_buffer_start + 3 * s * s / 4, size_t(s * s / 4)}, s / 2);
    MatrixView w({buffer.begin() + current_buffer_start + 4 * s * s / 4, size_t(s * s / 4)}, s / 2);

    current_buffer_start += 5 * s * s / 4;

    add(A11, A22, x);
    add(B11, B22, y);
    StrassenMatMul(x, y, u);
    add(A21, A22, x);
    StrassenMatMul(x, B11, C21);
    sub(B12, B22, x);
    StrassenMatMul(A11, x, C12);
    sub(B21, B11, x);
    StrassenMatMul(A22, x, v);
    add(A11, A12, x);
    StrassenMatMul(x, B22, w);
    sub(A21, A11, x);
    add(B11, B12, y);
    StrassenMatMul(x, y, C22);
    sub(A12, A22, x);
    add(B21, B22, y);
    StrassenMatMul(x, y, C11);
    C11 += u;
    C11 += v;
    C11 -= w;
    C22 += u;
    C22 += C12;
    C22 -= C21;
    C12 += w;
    C21 += v;

    current_buffer_start -= 5 * s * s / 4;
}
