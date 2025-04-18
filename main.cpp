#include <iostream>
#include <format>
#include <vector>
#include <span>
#include <cmath>
#include <chrono>
#include <algorithm>
#include <functional>
#include <random>

#include "include/MatrixView.hpp"
#include "include/iterative.hpp"
#include "include/block.hpp"
#include "include/recursive.hpp"
#include "include/MatrixMultiplier.hpp"

int generateRandomNumber(int min = 0, int max = 100) {
    static std::random_device rd;
    static std::mt19937 gen(rd());
    std::uniform_int_distribution<int> dis(min, max);
    return dis(gen);
}

template<class F>
auto time(F f, int N, int M, int P)
{
    std::vector<int> A(N * M);
    std::vector<int> B(M * P);
    std::vector<int> C(N * P);

    std::vector<int> E(N * P, 0);

    std::generate(A.begin(), A.end(), [](){ return generateRandomNumber(); });
    std::generate(B.begin(), B.end(), [](){ return generateRandomNumber(); });
    std::generate(C.begin(), C.end(), [](){ return generateRandomNumber(); });

    MatrixView A_view(A, M);
    MatrixView B_view(B, P);
    MatrixView C_view(C, P);

    MatrixView E_view(E, P);
    naiveCacheFriendlyMatMul(A_view, B_view, E_view, MatMulMode::Add);

    auto start = std::chrono::high_resolution_clock::now();

    f(A_view, B_view, C_view, MatMulMode::Overwrite);

    auto t_overwrite = std::chrono::high_resolution_clock::now() - start;
    bool overwrite_check = C == E;

    std::ranges::for_each(E, [](int& x){ x *= 2; });

    start = std::chrono::high_resolution_clock::now();

    f(A_view, B_view, C_view, MatMulMode::Add);

    auto t_add = std::chrono::high_resolution_clock::now() - start;
    bool add_check = C == E;

    return std::pair{ 
        std::pair{ t_overwrite, overwrite_check }, 
        std::pair{ t_add, add_check} 
    };
}

MatrixMultiplier recursive_until_size(int size)
{
    return MatrixMultiplier::recursive_then([size](int n, int m, int p){ return n < size || m < size || p < size; }, MatrixMultiplier::naive_cache_friendly_mutliplier);
}

MatrixMultiplier strassen_until_size(int size)
{
    return MatrixMultiplier::strassen_then([size](int n, int m, int p){ return n < size || m < size || p < size; }, MatrixMultiplier::naive_cache_friendly_mutliplier);
}

template<class F>
void print_test(std::string_view name, F f, int N, int M, int P)
{
    static auto to_ms = [](auto t){ return std::chrono::duration_cast<std::chrono::milliseconds>(t); };
    auto res = time(f, N, M, P);
    std::cout << std::format("{}: OWT: {:>7} {:>5}, ADD: {:>7} {:>5}\n", name, to_ms(res.first.first), res.first.second, to_ms(res.second.first), res.second.second);
}

int main(int argc, char* argv[])
{
    if (argc != 4)
    {
        std::cerr << "Usage: " << argv[0] << " <A_row_count> <A_col_count> <B_col_count>.";
        return 1;
    }

    std::cout << "CTEST_FULL_OUTPUT\n";

    int N = std::stoi(argv[1]);
    int M = std::stoi(argv[2]);
    int P = std::stoi(argv[3]);

    
    if (N <= 1024 && M <= 1024 && P <= 1024)
    {

    print_test("Naive                  (MatrixMultiplier)", MatrixMultiplier::naive_iterative_mutliplier              , N, M, P);
    print_test("Naive                                    ", naiveMatMul                                               , N, M, P);
    
    }
    print_test("Cache-friendly naive   (MatrixMultiplier)", MatrixMultiplier::naive_cache_friendly_mutliplier         , N, M, P);
    print_test("Cache-friendly naive                     ", naiveCacheFriendlyMatMul                                  , N, M, P);
    print_test("Cache-aware-blocked    (MatrixMultiplier)", MatrixMultiplier::cache_aware_blocked_multiplier          , N, M, P);
    print_test("Cache-aware-blocked                      ", cacheFriendlyBlockMatMul                                  , N, M, P);
    print_test("Recursive until size 4 (MatrixMultiplier)", recursive_until_size(4)                                   , N, M, P);
    print_test("Recursive until size 8 (MatrixMultiplier)", recursive_until_size(8)                                   , N, M, P);
    print_test("Recursive until size 16(MatrixMultiplier)", recursive_until_size(16)                                  , N, M, P);
    print_test("Recursive until size 32(MatrixMultiplier)", recursive_until_size(32)                                  , N, M, P);
    print_test("Recursive until size 64(MatrixMultiplier)", recursive_until_size(64)                                  , N, M, P);
    if (N == M && M == P && (N & N - 1) == 0)
    {
    
    print_test("Strassen until size 4  (MatrixMultiplier)", strassen_until_size(4)                                    , N, M, P);
    print_test("Strassen until size 8  (MatrixMultiplier)", strassen_until_size(8)                                    , N, M, P);
    print_test("Strassen until size 16 (MatrixMultiplier)", strassen_until_size(16)                                   , N, M, P);
    print_test("Strassen until size 32 (MatrixMultiplier)", strassen_until_size(32)                                   , N, M, P);
    print_test("Strassen until size 64 (MatrixMultiplier)", strassen_until_size(64)                                   , N, M, P);

    }
    print_test("Hybrid                 (MatrixMultiplier)", MatrixMultiplier::hybrid_multiplier(N, M, P)              , N, M, P);
    print_test("Multithreaded hybrid   (MatrixMultiplier)", MatrixMultiplier::multithreaded_hybrid_multiplier(N, M, P), N, M, P);

    return 0;
}