cmake_minimum_required(VERSION 3.10)

project(MatMulAlgorithms)

set(CMAKE_CXX_STANDARD 23)
add_executable(main main.cpp)

# Add the main executable as a test
enable_testing()
add_test(NAME MainTest COMMAND main --mult all without naive naive_MatMul --sizes default with 
    10_10_10 
    5_20_10 
    20_10_20 
    20_30_40 
    50_50_50 
    100_100_100 
    200_500_500 
    1000_1000_1000 
    2000_2000_2000 
    3000_3000_3000)