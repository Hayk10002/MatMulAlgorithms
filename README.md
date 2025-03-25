# MatMulAlgorithms

## Table of Contents
- [Introduction](#introduction)
- [Build and Run](#build-and-run)
- [Possible Output](#possible-output)
- [Benchmarked algorithms](#benchmarked-algorithms)
- [Benchmark results](#benchmark-results)

## Introduction
Here are implemented several sequential algorithms for multiplying two matrices. Algorithms like naive looping, looping in right order to be cache-friendly, dividing into blocks in cache-aware manner, divide and conquer recursive and Strassen's algorithms which are cache-oblivious, and a hybrid algorithm, that uses different algorithms based on the sizes of matrices. 

## Build and Run
To clone and run this project, you'll need [Git](https://git-scm.com) and [CMake](https://cmake.org/) installed on your computer. From your command line:

```bash
# Clone this repository
$ git clone https://github.com/Hayk10002/MatMulAlgorithms

# Go into the repository
$ cd MatMulAlgorithms

# Generate the build files
$ cmake -DCMAKE_BUILD_TYPE=Release -S . -B build

# Build the project
$ cmake --build build --config Release

# Then, run the executable generated in the `build` directory with the appropriate sizes to run the test.
$ your/path/to/exe/main.exe {A_row_count} {A_col_count} {B_col_count} # other dimensions are calculated based on the requirements for the matrix multiplication.
# example - .../main.exe 1000 500 800
# this multiplies two matrices of size 1000x500, 500x800, to get a matrix of size 1000x800.
```

## Possible Output
(for inputs 1024 1024 1024)

```
Naive                  (MatrixMultiplier): 2680739600ns
Naive                                    : 2742736300ns
Cache-friendly naive   (MatrixMultiplier): 214720300ns
Cache-friendly naive                     : 223352700ns
Cache-aware-blocked    (MatrixMultiplier): 284380500ns
Cache-aware-blocked                      : 392763600ns
Recursive until size 4 (MatrixMultiplier): 214823800ns
Recursive until size 8 (MatrixMultiplier): 221695200ns
Recursive until size 16(MatrixMultiplier): 219338200ns
Recursive until size 32(MatrixMultiplier): 229016900ns
Recursive until size 64(MatrixMultiplier): 232560400ns
Strassen until size 4  (MatrixMultiplier): 285143900ns
Strassen until size 8  (MatrixMultiplier): 309053100ns
Strassen until size 16 (MatrixMultiplier): 234409600ns
Strassen until size 32 (MatrixMultiplier): 220563500ns
Strassen until size 64 (MatrixMultiplier): 213096900ns
Hybrid                 (MatrixMultiplier): 172983600ns
```

## Benchmarked Algorithms

Almost all algorithms have time complexity of \(O(n^3)\), but in real live this can be very misleading.

### Naive

Naive matrix multiplication involves three nested loops to compute the product of two matrices. This algorithm is straightforward but not optimized for modern CPU caches.

### Cache friendly

This algorithm is the naive one but the two inner loops are swapped, greatly reducing the amount of cache misses.
The effect is very big. Computing time for large matrices reduces by more or less 10 times.  

### Cache aware

This algorithm divides the starting matrices into submatrices and insures that all elements needed for multiplication of submatrices can simultaniously be fetched and stored into the L1 cache, trying to utilize the cache maximally. 
The benchmarks show though that the speedup is not as much as could be expected.

### Recursive

This algorithm divides the matrices into 4 submatrices recursively, until some size, and then runs the cache friendy naive algorithm on them. It does not divide the matrices until size 1x1x1, because the benchmarking process showed that the time in that case skyrockets too high to benchmark the algorithm for sufficiently large matrices.

### Strassen's algorithm

It is a famous divide and conquer algorithm for matrix multiplication that uses a mathematical trick to reduce the number of simple mutliplications from n^3 to n^2.8, bringing down the time complexity of the algorithm to \(O(n^2.8)\), but it has several overheads for small sizes, and becomes useful only for large sizes. And it has the requirements for all matrices to be powers of 2.
This algorithm like the recursive one, also becomes very inefficient if used to divide the matrices until size 1x1x1, so it also stops dividing the matrices until some specified size.

### Hybrid approach

This algorithm splits the matrices into 4 smaller ones, the top left matrix size is picked to be the highest power of two smaller than the dimenstions of the matrices. The top left matrices are multiplied using Strassen's algorithm, and other multiplications are using the recursive method. Both Strassen's algorithm and recursive one stop subdividing the matrices at the moment when the whole multiplication of the submatrices can be done in the L1 cache, and then multiplies them with the cache-friendly naive method. 

## Benchmark results

| Algorithm | 1000x1000 | 1024x1024 | 2000x2000 | 2048x2048 | 3000x3000 |
| :-------- | :-------: | :-------: | :-------: | :-------: | :-------: |
| Naive | 0.938s | 2.680s | (too long) | (too long) | (too long) |
| Naive cache-friendly | 0.205s | 0.214s | 2.80s | 3.00s | 8.46s |
| Blocked (cache-aware) | 0.254s | 0.284s | 1.99s | 2.43s | 6.84s |
| Recursive until size 4 | 0.202s | 0.214s | 2.64s | 3.14s | 8.37s |
| Recursive until size 8 | 0.195s | 0.221s | 2.74s | 2.89s | 8.27s |
| Recursive until size 16 | 0.215s | 0.219s | 2.67s | 2.88s | 8.40s |
| Recursive until size 32 | 0.215s | 0.229s | 2.68s | 2.88s | 9.10s |
| Recursive until size 64 | 0.217s | 0.233s | 2.64s | 2.88s | 8.62s |
| Strassen until size 4 | - | 0.285s | - | 2.88s | - |
| Strassen until size 8 | - | 0.309s | - | 2.92s | - |
| Strassen until size 16 | - | 0.234s | - | 2.90s | - |
| Strassen until size 32 | - | 0.221s | - | 2.92s | - |
| Strassen until size 64 | - | 0.213s | - | 2.89s | - |
| Hybrid | 0.236s | 0.173s | 1.87s | 1.26s | 6.46s |

