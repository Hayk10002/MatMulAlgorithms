cmake_minimum_required(VERSION 3.10)

project(MatMulAlgorithms)

set(CMAKE_CXX_STANDARD 23)
add_executable(main main.cpp)

# Add the main executable as a test
enable_testing()
add_test(NAME MainTest1 COMMAND main 10 10 10)
add_test(NAME MainTest2 COMMAND main 5 20 10)
add_test(NAME MainTest3 COMMAND main 20 10 20)
add_test(NAME MainTest4 COMMAND main 16 16 16)
add_test(NAME MainTest5 COMMAND main 20 30 40)
add_test(NAME MainTest6 COMMAND main 50 50 50)
add_test(NAME MainTest7 COMMAND main 64 64 64)
add_test(NAME MainTest8 COMMAND main 100 100 100)
add_test(NAME MainTest9 COMMAND main 128 128 128)
add_test(NAME MainTest10 COMMAND main 200 500 500)
add_test(NAME MainTest11 COMMAND main 512 512 512)
add_test(NAME MainTest12 COMMAND main 1000 1000 1000)
add_test(NAME MainTest13 COMMAND main 1024 1024 1024)
add_test(NAME MainTest14 COMMAND main 2000 2000 2000)
add_test(NAME MainTest15 COMMAND main 2048 2048 2048)
add_test(NAME MainTest16 COMMAND main 3000 3000 3000)