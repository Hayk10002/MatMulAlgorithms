#pragma once
#include <functional>
#include <vector>
#include <cmath>
#include <set>
#include <string>
#include "MatrixView.hpp"
#include "getL1CacheSize.hpp"

#undef min
#undef max

std::set<std::string> logs{};

class MatrixMultiplier
{
    struct Multiplier
    {
        using PreconditionType = std::function<bool(int, int, int)>;
        using MultiplierType = std::function<void(const MatrixMultiplier&, MatrixView, MatrixView, MatrixView)>;
        PreconditionType precondition;
        MultiplierType multiplier;

        Multiplier(MultiplierType multiplier) : Multiplier([](int, int, int) { return true; }, multiplier) {}
        Multiplier(PreconditionType precondition, MultiplierType multiplier) : precondition(precondition), multiplier(multiplier) {}
    };

    std::vector<Multiplier> multipliers{};

    static bool valid_for_multiplying(MatrixView A, MatrixView B, MatrixView C)
    {
        return A.col_count() == B.row_count() && B.col_count() == C.col_count() && A.row_count() == C.row_count();
    }

public:
    void naive_iterative(MatrixView A, MatrixView B, MatrixView C) const
    { 
        for (int i = 0; i < C.row_count(); i++)
            for (int j = 0; j < C.col_count(); j++)
            {
                C(i, j) = 0;
                for (int k = 0; k < A.col_count(); k++)
                    C(i, j) += A(i, k) * B(k, j);
            }
    }

    void better_iterative(MatrixView A, MatrixView B, MatrixView C) const
    { 
        C.clear();
        for (int i = 0; i < C.row_count(); i++)
            for (int k = 0; k < A.col_count(); k++)
                for (int j = 0; j < C.col_count(); j++)
                    C(i, j) += A(i, k) * B(k, j);
    }

    struct BlockedMultiplier
    {
        int block_size;
        void operator()(const MatrixMultiplier& mult, MatrixView A, MatrixView B, MatrixView C)
        {
            C.clear();

            for (int i = 0; i < C.row_count(); i += block_size)            
                for (int j = 0; j < C.col_count(); j += block_size)                
                    for (int k = 0; k < A.col_count(); k += block_size)                    
                        mult(A.getSubMatrix(i, std::min(i + block_size, C.row_count()), k, std::min(k + block_size, A.col_count())),
                            B.getSubMatrix(k, std::min(k + block_size, A.col_count()), j, std::min(j + block_size, C.col_count())),
                            C.getSubMatrix(i, std::min(i + block_size, C.row_count()), j, std::min(j + block_size, C.col_count())));
        }
    };

    void recursive(MatrixView A, MatrixView B, MatrixView C) const
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
    
        (*this)(A11, B11, C11);
        (*this)(A12, B21, C11);
        (*this)(A11, B12, C12);
        (*this)(A12, B22, C12);
        (*this)(A21, B11, C21);
        (*this)(A22, B21, C21);
        (*this)(A21, B12, C22);
        (*this)(A22, B22, C22);
    }

    void strassen(MatrixView A, MatrixView B, MatrixView C) const
    {   
        int s = A.row_count();

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
        (*this)(x, y, u);
        add(A21, A22, x);
        (*this)(x, B11, C21);
        sub(B12, B22, x);
        (*this)(A11, x, C12);
        sub(B21, B11, x);
        (*this)(A22, x, v);
        add(A11, A12, x);
        (*this)(x, B22, w);
        sub(A21, A11, x);
        add(B11, B12, y);
        (*this)(x, y, C22);
        sub(A12, A22, x);
        add(B21, B22, y);
        (*this)(x, y, C11);
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

    static MatrixMultiplier one_strategy(Multiplier::MultiplierType multiplier)
    {
        MatrixMultiplier result{};
        result.multipliers.push_back(Multiplier{multiplier});
        return result;
    }

    static MatrixMultiplier add_strategy(Multiplier::PreconditionType precondition, Multiplier::MultiplierType strategy, const MatrixMultiplier& multiplier)
    {
        MatrixMultiplier result{multiplier};
        result.multipliers.insert(result.multipliers.begin(), Multiplier{precondition, strategy});
        return result;
    }

    static MatrixMultiplier into_blocks_then(int block_size, const MatrixMultiplier& multiplier)
    {
        return add_strategy([block_size](int n, int m, int p){ return n > block_size && m > block_size && p > block_size; },
                            BlockedMultiplier{block_size},
                            multiplier);
    }

    static MatrixMultiplier strassen_then(Multiplier::PreconditionType until, const MatrixMultiplier& multiplier)
    {
        return add_strategy([until](int n, int m, int p){ return n == m && m == p && (n & (n - 1)) == 0 && until(n, m, p); },
                            &MatrixMultiplier::strassen,
                            multiplier);
    }

    static MatrixMultiplier recursive_then(Multiplier::PreconditionType until, const MatrixMultiplier& multiplier)
    {
        return add_strategy(until, &MatrixMultiplier::recursive, multiplier);
    }

    static MatrixMultiplier naive_iterative_mutliplier;
    static MatrixMultiplier better_iterative_mutliplier;
    static MatrixMultiplier full_recursive_mutliplier;
    static MatrixMultiplier cache_aware_blocked_multiplier;
    static MatrixMultiplier hybrid_multiplier(int N, int M, int P)
    {
        int max_power_of_2_less_than_NMP = 1 << (int)log2(std::min({N, M, P}));
        return  into_blocks_then(max_power_of_2_less_than_NMP,
                strassen_then ([](int n, int m, int p){ return n * m + m * p + n * p > getL1CacheSize() / sizeof(int); },
                recursive_then([](int n, int m, int p){ return n * m + m * p + n * p > getL1CacheSize() / sizeof(int); },
                better_iterative_mutliplier
        )));
    }

    void operator()(MatrixView A, MatrixView B, MatrixView C) const
    {
        for (auto& multiplier : multipliers)
            if (valid_for_multiplying(A, B, C) &&
                multiplier.precondition(A.row_count(), A.col_count(), B.col_count()))
            {
                multiplier.multiplier(*this, A, B, C);
                return;
            }
    }
};

MatrixMultiplier MatrixMultiplier::naive_iterative_mutliplier{one_strategy(&MatrixMultiplier::naive_iterative)};
MatrixMultiplier MatrixMultiplier::better_iterative_mutliplier{one_strategy(&MatrixMultiplier::better_iterative)};
MatrixMultiplier MatrixMultiplier::full_recursive_mutliplier{one_strategy(&MatrixMultiplier::recursive)};
MatrixMultiplier MatrixMultiplier::cache_aware_blocked_multiplier(
    MatrixMultiplier::into_blocks_then(sqrt(getL1CacheSize() / sizeof(int) / 4),
    MatrixMultiplier::better_iterative_mutliplier));
