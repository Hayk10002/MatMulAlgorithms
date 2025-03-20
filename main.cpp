#include <iostream>
#include <format>
#include <vector>
#include <span>

struct MatrixView
{
private:
    int row_size;
    std::span<int> data;

    int row_start;
    int row_end;
    int col_start;
    int col_end;
    
public:
    MatrixView(std::span<int> data, int row_size, int row_start, int row_end, int col_start, int col_end)
        : data(data), row_size(row_size), row_start(row_start), row_end(row_end), col_start(col_start), col_end(col_end)
    {
    }

    MatrixView(std::span<int> data, int row_size): 
        MatrixView(data, row_size, 0, data.size() / row_size, 0, row_size)
    {
    }

    int& operator()(int row, int col)
    {
        return data[row_size * (row_start + row) + col_start + col];
    }

    int row_count() const
    {
        return row_end - row_start;
    }

    int col_count() const
    {
        return col_end - col_start;
    }

    MatrixView getSubMatrix(int row_start, int row_end, int col_start, int col_end)
    {
        return MatrixView(data, row_size, 
            this->row_start + row_start, 
            this->row_start + row_end, 
            this->col_start + col_start, 
            this->col_start + col_end);
    }
};

void naiveMatMul(MatrixView A, MatrixView B, MatrixView C)
{
    for (int i = 0; i < C.row_count(); i++)
    {
        for (int j = 0; j < C.col_count(); j++)
        {
            C(i, j) = 0;
            for (int k = 0; k < A.col_count(); k++)
            {
                C(i, j) += A(i, k) * B(k, j);
            }
        }
    }
}

void recursiveMatMulImpl(MatrixView A, MatrixView B, MatrixView C)
{
    if (A.col_count() != B.row_count() || B.col_count() != C.col_count() || A.row_count() != C.row_count()) // Invalid matrix dimensions
        return;    

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
    for (int i = 0; i < C.row_count(); i++)
        for (int j = 0; j < C.col_count(); j++)
            C(i, j) = 0;
    recursiveMatMulImpl(A, B, C);
}

int main()
{
    std::vector<int> A = {1, 2, 3, 4, 5, 6};
    std::vector<int> B = {7, 8, 9, 10, 11, 12};
    std::vector<int> C(9, 0);

    MatrixView A_view(A, 2);
    MatrixView B_view(B, 3);
    MatrixView C_view(C, 3);

    naiveMatMul(A_view, B_view, C_view);

    for (int i = 0; i < C.size(); i++)
    {
        std::cout << C[i] << " ";
    }
    std::cout << std::endl;

    recursiveMatMul(A_view, B_view, C_view);

    for (int i = 0; i < C.size(); i++)
    {
        std::cout << C[i] << " ";
    }
    std::cout << std::endl;

    return 0;
}