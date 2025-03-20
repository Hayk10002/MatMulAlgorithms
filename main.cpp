#include <iostream>
#include <format>
#include <vector>
#include <span>

struct MatrixView
{
    std::span<int> data;
    int row_size;
    int rows_start;
    int rows_end;
    int cols_start;
    int cols_end;

    MatrixView(std::span<int> data, int row_size, int rows_start, int rows_end, int cols_start, int cols_end)
        : data(data), row_size(row_size), rows_start(rows_start), rows_end(rows_end), cols_start(cols_start), cols_end(cols_end)
    {
    }

    MatrixView(std::span<int> data, int row_size): 
        MatrixView(data, row_size, 0, data.size() / row_size, 0, row_size)
    {
    }

    int& operator()(int row, int col)
    {
        return data[row_size * (rows_start + row) + cols_start + col];
    }

    int row_count() const
    {
        return rows_end - rows_start;
    }

    int col_count() const
    {
        return cols_end - cols_start;
    }

    MatrixView getSubMatrix(int row_start, int row_end, int col_start, int col_end)
    {
        return MatrixView(data, row_size, row_start, row_end, col_start, col_end);
    }
};

void naiveMatMul(MatrixView A, MatrixView B, MatrixView C)
{
    for (int i = 0; i < C.rows_end - C.rows_start; i++)
    {
        for (int j = 0; j < C.cols_end - C.cols_start; j++)
        {
            C(i, j) = 0;
            for (int k = 0; k < A.cols_end - A.cols_start; k++)
            {
                C(i, j) += A(i, k) * B(k, j);
            }
        }
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