# CacheAwareMergeSort

## Table of Contents
- [Introduction](#introduction)
- [Build and Run](#build-and-run)
- [Possible Output](#possible-output)
- [How does this work](#how-does-this-work)

## Introduction
The purpose of this project is to implement and test a cache-aware merge sort algorithm and compare it's performance against a normal merge sort algorithm.

## Build and Run
To clone and run this project, you'll need [Git](https://git-scm.com) and [CMake](https://cmake.org/) installed on your computer. From your command line:

```bash
# Clone this repository
$ git clone https://github.com/Hayk10002/CacheAwareMergeSort

# Go into the repository
$ cd CacheAwareMergeSort

# Generate the build files
$ cmake -DCMAKE_BUILD_TYPE=Release -S . -B build

# Build the project
$ cmake --build build --config Release

# Then, run the executable generated in the `build` directory with the vector size and iteration count to run the test.
$ your/path/to/exe/main.exe {vector_size} {iteration_count}
# example - .../main.exe 1000000 100
```

## Possible Output
(for vector size 1000000 and iteration count 10)

```
L1 Cache Size: 49152 bytes
Chunk Size: auto  Time: 549ms
Chunk Size: 6144  Time: 540ms
```

## How does this work
The project measures the performance of a cache-aware merge sort algorithm by considering the L1 cache size. The algorithm divides the data into chunks that are multiples of half L1 cache trying to minimize cache misses and improve performance.

The code benchmarks the time taken to sort a vector of integers using specified chunk sizes and a normal merge sort.
The outcome almost certainly is that it doesn't matter if the divisions of the data are done according to L1 cache size, because the merge sort already is cache-oblivious, meaning that it maximizes cache utility regardless of the cache size. 
