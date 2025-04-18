# MatMulAlgorithms

## Table of Contents
- [Introduction](#introduction)
- [Build and Run](#build-and-run)
- [Possible Output](#possible-output)
- [What is done](#what-is-done)
- [Benchmarked algorithms](#benchmarked-algorithms)
- [Benchmark results](#benchmark-results)

## Introduction
Here are implemented several sequential and one multithreaded algorithms for multiplying two matrices. Algorithms like naive looping, looping in right order to be cache-friendly, dividing into blocks in cache-aware manner, divide and conquer recursive and Strassen's algorithms which are cache-oblivious, and a hybrid algorithm, and a multithreaded algorithm each thread using the hybrid algorithm, that uses different algorithms based on the sizes of matrices. 

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
Naive                  (MatrixMultiplier): OWT:  2933ms  true, ADD:  2676ms  true
Naive                                    : OWT:  2615ms  true, ADD:  2604ms  true
Cache-friendly naive   (MatrixMultiplier): OWT:   222ms  true, ADD:   227ms  true
Cache-friendly naive                     : OWT:   227ms  true, ADD:   223ms  true
Cache-aware-blocked    (MatrixMultiplier): OWT:   306ms  true, ADD:   303ms  true
Cache-aware-blocked                      : OWT:   344ms  true, ADD:   355ms  true
Recursive until size 4 (MatrixMultiplier): OWT:   248ms  true, ADD:   250ms  true
Recursive until size 8 (MatrixMultiplier): OWT:   255ms  true, ADD:   235ms  true
Recursive until size 16(MatrixMultiplier): OWT:   221ms  true, ADD:   241ms  true
Recursive until size 32(MatrixMultiplier): OWT:   234ms  true, ADD:   254ms  true
Recursive until size 64(MatrixMultiplier): OWT:   323ms  true, ADD:   294ms  true
Strassen until size 4  (MatrixMultiplier): OWT:   246ms  true, ADD:   228ms  true
Strassen until size 8  (MatrixMultiplier): OWT:   278ms  true, ADD:   291ms  true
Strassen until size 16 (MatrixMultiplier): OWT:   231ms  true, ADD:   238ms  true
Strassen until size 32 (MatrixMultiplier): OWT:   242ms  true, ADD:   224ms  true
Strassen until size 64 (MatrixMultiplier): OWT:   240ms  true, ADD:   229ms  true
Hybrid                 (MatrixMultiplier): OWT:   199ms  true, ADD:   203ms  true
Multithreaded hybrid   (MatrixMultiplier): OWT:   166ms  true, ADD:   189ms  true
```

## What is done

Every algorithm has two modes, one that overrides the output matrix (OWT), and one that adds the result of the multiplication to the output matrix (ADD). Each mode is benchmarked and then the result is checked with the naive cache friendly result for correctness (the output "true" indicates that the multiplication is correct).

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
| Naive                     | 0.988s | 2.615s | (too long) | (too long) | (too long) |
| Naive cache-friendly      | 0.207s | 0.222s | 2.60s | 2.79s | 8.42s |
| Blocked (cache-aware)     | 0.253s | 0.306s | 2.12s | 2.60s | 7.20s |
| Recursive until size 4    | 0.208s | 0.248s | 2.65s | 2.85s | 8.61s |
| Recursive until size 8    | 0.220s | 0.255s | 2.66s | 2.86s | 8.59s |
| Recursive until size 16   | 0.214s | 0.221s | 2.60s | 2.81s | 9.15s |
| Recursive until size 32   | 0.223s | 0.234s | 2.58s | 2.81s |21.91s |
| Recursive until size 64   | 0.210s | 0.323s | 2.58s | 2.99s | 8.63s |
| Strassen until size 4     | -      | 0.246s | -     | 2.99s | -     |
| Strassen until size 8     | -      | 0.278s | -     | 2.84s | -     |
| Strassen until size 16    | -      | 0.231s | -     | 3.19s | -     |
| Strassen until size 32    | -      | 0.242s | -     | 2.92s | -     |
| Strassen until size 64    | -      | 0.240s | -     | 2.81s | -     |
| Hybrid                    | 0.268s | 0.199s | 1.97s | 1.32s | 6.67s |
| Multithreaded hybrid      | 0.126s | 0.166s | 0.97s | 1.97s | 3.30s |
