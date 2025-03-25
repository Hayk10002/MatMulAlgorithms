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

    std::generate(A.begin(), A.end(), [](){ return generateRandomNumber(); });
    std::generate(B.begin(), B.end(), [](){ return generateRandomNumber(); });
    std::generate(C.begin(), C.end(), [](){ return generateRandomNumber(); });

    MatrixView A_view(A, M);
    MatrixView B_view(B, P);
    MatrixView C_view(C, P);

    auto start = std::chrono::high_resolution_clock::now();

    f(A_view, B_view, C_view);

    return std::chrono::high_resolution_clock::now() - start;
}

MatrixMultiplier recursive_until_size(int size)
{
    return MatrixMultiplier::recursive_then([size](int n, int m, int p){ return n < size || m < size || p < size; }, MatrixMultiplier::naive_cache_friendly_mutliplier);
}

MatrixMultiplier strassen_until_size(int size)
{
    return MatrixMultiplier::strassen_then([size](int n, int m, int p){ return n < size || m < size || p < size; }, MatrixMultiplier::naive_cache_friendly_mutliplier);
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

    std::cout << std::format("Naive                  (MatrixMultiplier): {:<8}\n", time(MatrixMultiplier::naive_iterative_mutliplier,      N, M, P));
    std::cout << std::format("Naive                                    : {:<8}\n", time(naiveMatMul,                                       N, M, P));
    
    }
    std::cout << std::format("Cache-friendly naive   (MatrixMultiplier): {:<8}\n", time(MatrixMultiplier::naive_cache_friendly_mutliplier, N, M, P));
    std::cout << std::format("Cache-friendly naive                     : {:<8}\n", time(naiveCacheFriendlyMatMul,                          N, M, P));
    std::cout << std::format("Cache-aware-blocked    (MatrixMultiplier): {:<8}\n", time(MatrixMultiplier::cache_aware_blocked_multiplier,  N, M, P));
    std::cout << std::format("Cache-aware-blocked                      : {:<8}\n", time(cacheFriendlyBlockMatMul,                          N, M, P));
    std::cout << std::format("Recursive until size 4 (MatrixMultiplier): {:<8}\n", time(recursive_until_size(4),                           N, M, P));
    std::cout << std::format("Recursive until size 8 (MatrixMultiplier): {:<8}\n", time(recursive_until_size(8),                           N, M, P));
    std::cout << std::format("Recursive until size 16(MatrixMultiplier): {:<8}\n", time(recursive_until_size(16),                          N, M, P));
    std::cout << std::format("Recursive until size 32(MatrixMultiplier): {:<8}\n", time(recursive_until_size(32),                          N, M, P));
    std::cout << std::format("Recursive until size 64(MatrixMultiplier): {:<8}\n", time(recursive_until_size(64),                          N, M, P));
    if (N == M && M == P && (N & N - 1) == 0)
    {
    
    std::cout << std::format("Strassen until size 4  (MatrixMultiplier): {:<8}\n", time(strassen_until_size(4),                            N, M, P));
    std::cout << std::format("Strassen until size 8  (MatrixMultiplier): {:<8}\n", time(strassen_until_size(8),                            N, M, P));
    std::cout << std::format("Strassen until size 16 (MatrixMultiplier): {:<8}\n", time(strassen_until_size(16),                           N, M, P));
    std::cout << std::format("Strassen until size 32 (MatrixMultiplier): {:<8}\n", time(strassen_until_size(32),                           N, M, P));
    std::cout << std::format("Strassen until size 64 (MatrixMultiplier): {:<8}\n", time(strassen_until_size(64),                           N, M, P));     

    }
    std::cout << std::format("Hybrid                 (MatrixMultiplier): {:<8}\n", time(MatrixMultiplier::hybrid_multiplier(N, M, P),      N, M, P));

    return 0;
}