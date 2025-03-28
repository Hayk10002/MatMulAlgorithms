#pragma once
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

    int operator()(int row, int col) const 
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

    void clear()
    {
        for (int i = 0; i < row_count(); i++)
            for (int j = 0; j < col_count(); j++)
                (*this)(i, j) = 0;
    }

    MatrixView& operator=(MatrixView O)
    {
        for (int i = 0; i < row_count() && i < O.row_count(); i++)
            for (int j = 0; j < col_count() && i < O.col_count(); j++)
                (*this)(i, j) = O(i, j);

        return *this;
    }

    MatrixView& operator+=(MatrixView O)
    {
        for (int i = 0; i < row_count() && i < O.row_count(); i++)
            for (int j = 0; j < col_count() && i < O.col_count(); j++)
                (*this)(i, j) += O(i, j);

        return *this;
    }

    MatrixView& operator-=(MatrixView O)
    {
        for (int i = 0; i < row_count() && i < O.row_count(); i++)
            for (int j = 0; j < col_count() && i < O.col_count(); j++)
                (*this)(i, j) -= O(i, j);

        return *this;
    }
};

void add(MatrixView A, MatrixView B, MatrixView C) // C = A + B
{
    int row_count = std::min({A.row_count(), B.row_count(), C.row_count()});
    int col_count = std::min({A.col_count(), B.col_count(), C.col_count()});
    for (int i = 0; i < row_count; i++)
            for (int j = 0; j < col_count; j++)
                C(i, j) = A(i, j) + B(i, j);
}

void sub(MatrixView A, MatrixView B, MatrixView C) // C = A - B
{
    int row_count = std::min({A.row_count(), B.row_count(), C.row_count()});
    int col_count = std::min({A.col_count(), B.col_count(), C.col_count()});
    for (int i = 0; i < row_count; i++)
            for (int j = 0; j < col_count; j++)
                C(i, j) = A(i, j) - B(i, j);
}