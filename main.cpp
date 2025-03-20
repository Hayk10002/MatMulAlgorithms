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

void recursiveMatMul(MatrixView A, MatrixView B, MatrixView C)
{
    if (A.row_count() == 1 && A.col_count() == 1 && B.row_count() == 1 && B.col_count() == 1)
    {
        C(0, 0) += A(0, 0) * B(0, 0);
    }
    else
    {
        int row_mid = A.row_count() / 2;
        int col_mid = A.col_count() / 2;

        MatrixView A11 = A.getSubMatrix(0, row_mid, 0, col_mid);
        MatrixView A12 = A.getSubMatrix(0, row_mid, col_mid, A.col_count());
        MatrixView A21 = A.getSubMatrix(row_mid, A.row_count(), 0, col_mid);
        MatrixView A22 = A.getSubMatrix(row_mid, A.row_count(), col_mid, A.col_count());

        MatrixView B11 = B.getSubMatrix(0, row_mid, 0, col_mid);
        MatrixView B12 = B.getSubMatrix(0, row_mid, col_mid, B.col_count());
        MatrixView B21 = B.getSubMatrix(row_mid, B.row_count(), 0, col_mid);
        MatrixView B22 = B.getSubMatrix(row_mid, B.row_count(), col_mid, B.col_count());

        MatrixView C11 = C.getSubMatrix(0, row_mid, 0, col_mid);
        MatrixView C12 = C.getSubMatrix(0, row_mid, col_mid, C.col_count());
        MatrixView C21 = C.getSubMatrix(row_mid, C.row_count(), 0, col_mid);
        MatrixView C22 = C.getSubMatrix(row_mid, C.row_count(), col_mid, C.col_count());

        recursiveMatMul(A11, B11, C11);
        recursiveMatMul(A12, B21, C11);
        recursiveMatMul(A11, B12, C12);
        recursiveMatMul(A12, B22, C12);
        recursiveMatMul(A21, B11, C21);
        recursiveMatMul(A22, B21, C21);
        recursiveMatMul(A21, B12, C22);
        recursiveMatMul(A22, B22, C22);
    }
}


int main()
{
    std::vector<int> A = {1, 2, 3, 4, 5, 6};
    std::vector<int> B = {7, 8, 9, 10, 11, 12};
    std::vector<int> C(4, 0);

    MatrixView A_view(A, 2);
    MatrixView B_view(B, 2);
    MatrixView C_view(C, 2);

    naiveMatMul(A_view, B_view, C_view);

    for (int i = 0; i < C.size(); i++)
    {
        std::cout << C[i] << " ";
    }
    std::cout << std::endl;

    return 0;
}