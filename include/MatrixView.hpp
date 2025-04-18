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
    MatrixView() = default;
    MatrixView(std::span<int> data, int row_size, int row_start, int row_end, int col_start, int col_end)
        : data(data), row_size(row_size), row_start(row_start), row_end(row_end), col_start(col_start), col_end(col_end)
    {
    }

    MatrixView(std::span<int> data, int row_size): 
        MatrixView(data, row_size, 0, row_size ? data.size() / row_size : 0, 0, row_size)
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

    MatrixView& clone_from(MatrixView O)
    {
        for (int i = 0; i < row_count() && i < O.row_count(); i++)
            for (int j = 0; j < col_count() && i < O.col_count(); j++)
                (*this)(i, j) = O(i, j);

        return *this;
    }

    MatrixView& add_eq(MatrixView O)
    {
        for (int i = 0; i < row_count() && i < O.row_count(); i++)
            for (int j = 0; j < col_count() && i < O.col_count(); j++)
                (*this)(i, j) += O(i, j);

        return *this;
    }

    MatrixView& rem_eq(MatrixView O)
    {
        for (int i = 0; i < row_count() && i < O.row_count(); i++)
            for (int j = 0; j < col_count() && i < O.col_count(); j++)
                (*this)(i, j) -= O(i, j);

        return *this;
    }

    bool is_equal(const MatrixView& other) const
    {
        if (row_count() != other.row_count() || col_count() != other.col_count())
            return false;
        for (int i = 0; i < row_count(); i++)
            for (int j = 0; j < col_count(); j++)
                if ((*this)(i, j) != other(i, j))
                    return false;
        return true;
    }

    bool is_same_view(const MatrixView& other) const
    {
        return row_size == other.row_size &&
               data.data() == other.data.data() &&
               row_start == other.row_start &&
               row_end == other.row_end &&
               col_start == other.col_start &&
               col_end == other.col_end;
    }

    friend class std::hash<MatrixView>;
};

namespace std
{
    template<>
    struct hash<MatrixView>
    {
        std::size_t operator()(const MatrixView& mv) const
        {
            std::size_t h = 0;

            auto hash_combine = [&h](auto const& value)
            {
                h ^= std::hash<std::remove_cvref_t<decltype(value)>>{}(value) + 0x9e3779b9 + (h << 6) + (h >> 2);
            };

            hash_combine(mv.row_size);
            hash_combine(reinterpret_cast<std::uintptr_t>(mv.data.data()));
            hash_combine(mv.row_start);
            hash_combine(mv.row_end);
            hash_combine(mv.col_start);
            hash_combine(mv.col_end);

            return h;
        }
    };

    template<>
    struct hash<std::tuple<MatrixView, MatrixView, MatrixView>>
    {
        std::size_t operator()(const std::tuple<MatrixView, MatrixView, MatrixView>& tup) const
        {
            std::size_t seed = 0;

            auto hash_combine = [&seed](const auto& value)
            {
                using T = std::remove_cvref_t<decltype(value)>;
                seed ^= std::hash<T>{}(value) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
            };

            hash_combine(std::get<0>(tup));
            hash_combine(std::get<1>(tup));
            hash_combine(std::get<2>(tup));

            return seed;
        }
    };
}

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

enum class MatMulMode
{
    Overwrite, // C  = A * B
    Add        // C += A * B
};