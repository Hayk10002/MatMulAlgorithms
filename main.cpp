#include <iostream>
#include <format>
#include <vector>
#include <span>
#include <cmath>

#include "include/MatrixView.hpp"
#include "include/iterative.hpp"
#include "include/block.hpp"
#include "include/recursive.hpp"
#include <functional>

int main()
{
    std::vector<int> A = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
    std::vector<int> B = {7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
    std::vector<int> C(64, 0);

    MatrixView A_view(A, 8);
    MatrixView B_view(B, 8);
    MatrixView C_view(C, 8);

    naiveMatMul(A_view, B_view, C_view);

    for (int i = 0; i < C.size(); i++)
    {
        std::cout << C[i] << " ";
    }
    std::cout << std::endl;

    blockMatMul(A_view, B_view, C_view);

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

    StrassenMatMul(A_view, B_view, C_view);

    for (int i = 0; i < C.size(); i++)
    {
        std::cout << C[i] << " ";
    }
    std::cout << std::endl;

    return 0;
}