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
    std::vector<int> A = {1, 2, 3, 4, 5, 6};
    std::vector<int> B = {7, 8, 9, 10, 11, 12};
    std::vector<int> C(9, 0);

    MatrixView A_view(A, 2);
    MatrixView B_view(B, 3);
    MatrixView C_view(C, 3);

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

    return 0;
}