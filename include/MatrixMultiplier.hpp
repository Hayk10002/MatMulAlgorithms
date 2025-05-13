#pragma once
#include <functional>
#include <vector>
#include <cmath>
#include <set>
#include <atomic>
#include <thread>
#include <variant>
#include <unordered_set>
#include "MatrixView.hpp"
#include "getL1CacheSize.hpp"

#undef min
#undef max

class MatrixMultiplier
{
    struct Multiplier
    {
        using PreconditionTypeWithSizes = std::function<bool(int, int, int)>;
        using PreconditionTypeWithViews = std::function<bool(MatrixView, MatrixView, MatrixView)>;
        using PreconditionType = std::variant<PreconditionTypeWithSizes, PreconditionTypeWithViews>;
        using MultiplierType = std::function<void(const MatrixMultiplier&, MatrixView, MatrixView, MatrixView, MatMulMode)>;
        PreconditionType precondition;
        MultiplierType multiplier;

        Multiplier(MultiplierType multiplier) : Multiplier([](int, int, int) { return true; }, multiplier) {}
        Multiplier(PreconditionType precondition, MultiplierType multiplier) : precondition(precondition), multiplier(multiplier) {}

        bool can_call_precondition_with_sizes() const
        {
            return std::holds_alternative<PreconditionTypeWithSizes>(precondition);
        }
        bool can_call_precondition_with_views() const
        {
            return std::holds_alternative<PreconditionTypeWithViews>(precondition);
        }
        bool precondition_with_sizes(int n, int m, int p) const
        {
            return std::get<PreconditionTypeWithSizes>(precondition)(n, m, p);
        }
        bool precondition_with_views(MatrixView A, MatrixView B, MatrixView C) const
        {
            return std::get<PreconditionTypeWithViews>(precondition)(A, B, C);
        }

    };

    std::vector<Multiplier> multipliers{};

    static bool valid_for_multiplying(MatrixView A, MatrixView B, MatrixView C)
    {
        return A.col_count() == B.row_count() && B.col_count() == C.col_count() && A.row_count() == C.row_count();
    }

public:
    void naive_iterative(MatrixView A, MatrixView B, MatrixView C, MatMulMode mode) const
    { 
        if (mode == MatMulMode::Overwrite)
            C.clear();

        int n = A.row_count(), m = A.col_count(), p = B.col_count();
        for (int i = 0; i < n; i++)
            for (int j = 0; j < p; j++)
                for (int k = 0; k < m; k++)
                    C(i, j) += A(i, k) * B(k, j);
    }

    void naive_cache_friendly_iterative(MatrixView A, MatrixView B, MatrixView C, MatMulMode mode) const
    {
        if (mode == MatMulMode::Overwrite)
            C.clear();

        int n = A.row_count(), m = A.col_count(), p = B.col_count();
        for (int i = 0; i < n; i++)
            for (int k = 0; k < m; k++)
                for (int j = 0; j < p; j++)
                    C(i, j) += A(i, k) * B(k, j);
    }

    struct BlockedMultiplier
    {
        int block_size;
        void operator()(const MatrixMultiplier& mult, MatrixView A, MatrixView B, MatrixView C, MatMulMode mode)
        {
            if (mode == MatMulMode::Overwrite)
                C.clear();

            int n = A.row_count(), m = A.col_count(), p = B.col_count();

            for (int i = 0; i < n; i += block_size)            
                for (int k = 0; k < m; k += block_size)                    
                    for (int j = 0; j < p; j += block_size)                
                        mult(A.getSubMatrix(i, std::min(i + block_size, n), k, std::min(k + block_size, m)),
                            B.getSubMatrix(k, std::min(k + block_size, m), j, std::min(j + block_size, p)),
                            C.getSubMatrix(i, std::min(i + block_size, n), j, std::min(j + block_size, p)), MatMulMode::Add);
        }
    };

    void recursive(MatrixView A, MatrixView B, MatrixView C, MatMulMode mode) const
    {
        if (A.row_count() == 0 || A.col_count() == 0 || B.col_count() == 0) // Empty matrices
            return;

        if (mode == MatMulMode::Overwrite)
            C.clear();
    
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
    
        (*this)(A11, B11, C11, MatMulMode::Add);
        (*this)(A12, B21, C11, MatMulMode::Add);
        (*this)(A11, B12, C12, MatMulMode::Add);
        (*this)(A12, B22, C12, MatMulMode::Add);
        (*this)(A21, B11, C21, MatMulMode::Add);
        (*this)(A22, B21, C21, MatMulMode::Add);
        (*this)(A21, B12, C22, MatMulMode::Add);
        (*this)(A22, B22, C22, MatMulMode::Add);
    }

    struct MultithreadedRecursiveMultiplier
    {
        void operator()(const MatrixMultiplier& mult, MatrixView A, MatrixView B, MatrixView C, MatMulMode mode)
        {
            if (A.row_count() == 0 || A.col_count() == 0 || B.col_count() == 0) // Empty matrices
            return;

            if (mode == MatMulMode::Overwrite)
                C.clear();
    
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
        
            std::vector<std::thread> threads;
            std::vector<std::array<MatrixView, 3>> parts_to_do_in_this_thread;
            threads.reserve(3);
            threads.push_back(std::thread([&]() 
            { 
                mult(A11, B11, C11, MatMulMode::Add); 
                mult(A12, B21, C11, MatMulMode::Add); 
            }));
            threads.push_back(std::thread([&]() 
            { 
                mult(A11, B12, C12, MatMulMode::Add);
                mult(A12, B22, C12, MatMulMode::Add);
            }));
            threads.push_back(std::thread([&]() 
            { 
                mult(A21, B11, C21, MatMulMode::Add); 
                mult(A22, B21, C21, MatMulMode::Add); 
            }));

            mult(A21, B12, C22, MatMulMode::Add); 
            mult(A22, B22, C22, MatMulMode::Add); 

            for (auto& thread : threads) thread.join();
        }
    };

    void strassen(MatrixView A, MatrixView B, MatrixView C, MatMulMode mode) const
    {   
        int s = A.row_count();

        if (s == 1)
        {
            if (mode == MatMulMode::Overwrite) C(0, 0)  = A(0, 0) * B(0, 0);
            else                               C(0, 0) += A(0, 0) * B(0, 0);
            return;
        }

        std::vector<int> D_vec;
        MatrixView D;
        if (mode == MatMulMode::Add)
        {
            D_vec.resize(s * s);
            D = MatrixView(D_vec, s);
            std::swap(C, D);
        }

        std::vector<int> buffer(5 * s * s / 4);

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

        MatrixView x({buffer.begin() + 0 * s * s / 4, size_t(s * s / 4)}, s / 2);
        MatrixView y({buffer.begin() + 1 * s * s / 4, size_t(s * s / 4)}, s / 2);
        MatrixView u({buffer.begin() + 2 * s * s / 4, size_t(s * s / 4)}, s / 2);
        MatrixView v({buffer.begin() + 3 * s * s / 4, size_t(s * s / 4)}, s / 2);
        MatrixView w({buffer.begin() + 4 * s * s / 4, size_t(s * s / 4)}, s / 2);

        add(A11, A22, x);
        add(B11, B22, y);
        (*this)(x, y, u, MatMulMode::Overwrite);
        add(A21, A22, x);
        (*this)(x, B11, C21, MatMulMode::Overwrite);
        sub(B12, B22, x);
        (*this)(A11, x, C12, MatMulMode::Overwrite);
        sub(B21, B11, x);
        (*this)(A22, x, v, MatMulMode::Overwrite);
        add(A11, A12, x);
        (*this)(x, B22, w, MatMulMode::Overwrite);
        sub(A21, A11, x);
        add(B11, B12, y);
        (*this)(x, y, C22, MatMulMode::Overwrite);
        sub(A12, A22, x);
        add(B21, B22, y);
        (*this)(x, y, C11, MatMulMode::Overwrite);
        C11.add_eq(u);
        C11.add_eq(v);
        C11.rem_eq(w);
        C22.add_eq(u);
        C22.add_eq(C12);
        C22.rem_eq(C21);
        C12.add_eq(w);
        C21.add_eq(v);

        if (mode == MatMulMode::Add) D.add_eq(C);
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

    static MatrixMultiplier possibly_multithreaded(Multiplier::PreconditionTypeWithSizes until, const MatrixMultiplier& multiplier)
    {
        return add_strategy(until, 
                            MultithreadedRecursiveMultiplier{},
                            multiplier);
    }

    static MatrixMultiplier into_blocks_then(int block_size, const MatrixMultiplier& multiplier)
    {
        return add_strategy([block_size](int n, int m, int p){ return n > block_size && m > block_size && p > block_size; },
                            BlockedMultiplier{block_size},
                            multiplier);
    }

    static MatrixMultiplier strassen_then(Multiplier::PreconditionTypeWithSizes until, const MatrixMultiplier& multiplier)
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
    static MatrixMultiplier naive_cache_friendly_mutliplier;
    static MatrixMultiplier full_recursive_mutliplier;
    static MatrixMultiplier cache_aware_blocked_multiplier;
    static MatrixMultiplier hybrid_multiplier(int N, int M, int P)
    {
        int max_power_of_2_less_than_NMP = 1 << (int)log2(std::min({N, M, P}));
        return  into_blocks_then(max_power_of_2_less_than_NMP,
                strassen_then ([](int n, int m, int p){ return n * m + m * p + n * p > getL1CacheSize() / sizeof(int); },
                recursive_then([](int n, int m, int p){ return n * m + m * p + n * p > getL1CacheSize() / sizeof(int); },
                naive_cache_friendly_mutliplier
        )));
    }

    static MatrixMultiplier multithreaded_hybrid_multiplier(int N, int M, int P)
    {
        int max_power_of_2_less_than_NMP = 1 << (int)log2(std::min({N, M, P}));
        return  
                possibly_multithreaded([N](int n, int, int){ return n > (N / 4 + 1); },
                into_blocks_then(max_power_of_2_less_than_NMP / 4,
                strassen_then ([](int n, int m, int p){ return n * m + m * p + n * p > getL1CacheSize() / sizeof(int); },
                recursive_then([](int n, int m, int p){ return n * m + m * p + n * p > getL1CacheSize() / sizeof(int); },
                naive_cache_friendly_mutliplier
        ))));
    }

    void operator()(MatrixView A, MatrixView B, MatrixView C, MatMulMode mode) const
    {
        for (auto& multiplier : multipliers)
            if (valid_for_multiplying(A, B, C) &&
                ((multiplier.can_call_precondition_with_sizes() && multiplier.precondition_with_sizes(A.row_count(), A.col_count(), B.col_count())) ||
                 (multiplier.can_call_precondition_with_views() && multiplier.precondition_with_views(A, B, C))))
            {
                multiplier.multiplier(*this, A, B, C, mode);
                return;
            }
    }
};

MatrixMultiplier MatrixMultiplier::naive_iterative_mutliplier{one_strategy(&MatrixMultiplier::naive_iterative)};
MatrixMultiplier MatrixMultiplier::naive_cache_friendly_mutliplier{one_strategy(&MatrixMultiplier::naive_cache_friendly_iterative)};
MatrixMultiplier MatrixMultiplier::full_recursive_mutliplier{one_strategy(&MatrixMultiplier::recursive)};
MatrixMultiplier MatrixMultiplier::cache_aware_blocked_multiplier(
    MatrixMultiplier::into_blocks_then(sqrt(getL1CacheSize() / sizeof(int) / 4),
    MatrixMultiplier::naive_cache_friendly_mutliplier));
