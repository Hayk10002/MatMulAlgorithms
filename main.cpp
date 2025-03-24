#include <iostream>
#include <format>
#include <vector>
#include <span>
#include <cmath>

#include "include/MatrixView.hpp"
#include "include/iterative.hpp"
#include "include/block.hpp"
#include "include/recursive.hpp"
#include "include/MatrixMultiplier.hpp"
#include <functional>

int main()
{
    std::vector<int> A(10000, 4);
    std::vector<int> B(10000, 6);
    std::vector<int> C(15625);

    MatrixView A_view(A, 80);
    MatrixView B_view(B, 125);
    MatrixView C_view(C, 125);

    naiveMatMul(A_view, B_view, C_view);

    std::cout << C[0] << std::endl;

    blockMatMul(A_view, B_view, C_view);

    std::cout << C[0] << std::endl;

    recursiveMatMul(A_view, B_view, C_view);

    std::cout << C[0] << std::endl;

    StrassenMatMul(A_view, B_view, C_view);

    std::cout << C[0] << std::endl;

    MatrixMultiplier::hybrid_multiplier(125, 80, 125)(A_view, B_view, C_view);

    for (auto log : logs) std::cout << log;

    return 0;
}