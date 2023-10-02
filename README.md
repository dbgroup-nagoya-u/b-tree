# B<sup>+</sup>-trees

[![Ubuntu-20.04](https://github.com/dbgroup-nagoya-u/b-tree/actions/workflows/unit_tests.yaml/badge.svg)](https://github.com/dbgroup-nagoya-u/b-tree/actions/workflows/unit_tests.yaml)

This repository contains open source implementations of B<sup>+</sup>-trees for research use. The purpose of this implementation is to reproduce B<sup>+</sup>-trees and measure its performance.

- [Build](#build)
    - [Prerequisites](#prerequisites)
    - [Build Options](#build-options)
    - [Build and Run Unit Tests](#build-and-run-unit-tests)
- [Usage](#usage)
    - [Linking by CMake](#linking-by-cmake)
    - [Implemented B+-Tree Variants](#implemented-b-tree-variants)
    - [Read/Write APIs](#readwrite-apis)
- [Acknowledgments](#acknowledgments)

## Build

**Note**: this is a header only library. You can use this without pre-build.

### Prerequisites

```bash
sudo apt update && sudo apt install -y build-essential cmake
```

### Build Options

#### Tuning Parameters

- `B_TREE_PAGE_SIZE`: Page size in bytes (default `1024`).
- `B_TREE_MAX_DELETED_SPACE_SIZE`: Invoking a clean-up-operation if the deleted space size of a node exceeds this threshold (default `${B_TREE_PAGE_SIZE} / 4`).
- `B_TREE_MIN_FREE_SPACE_SIZE`: Invoking a split operation if the free space size of a node after cleaning up becomes lower than this threshold (default `${B_TREE_PAGE_SIZE} / 8`).
- `B_TREE_MIN_USED_SPACE_SIZE`: Invoking a merge operation if the used space size of a node becomes lower than this threshold (default `${B_TREE_PAGE_SIZE} / 16`).
- `B_TREE_MAX_VARLEN_DATA_SIZE`: The expected maximum size of a variable-length data (default `128`).

#### Build Options for Unit Testing

- `B_TREE_BUILD_TESTS`: Build unit tests for this library if `ON` (default `OFF`).
- `B_TREE_TEST_BUILD_PML`: Build tests for PML-based B+trees if `ON` (default `OFF`).
- `B_TREE_TEST_BUILD_PSL`: Build tests for PSL-based B+trees if `ON` (default `OFF`).
- `B_TREE_TEST_BUILD_OML`: Build tests for OML-based B+trees if `ON` (default `OFF`).
- `B_TREE_TEST_BUILD_OSL`: Build tests for OSL-based B+trees if `ON` (default `OFF`).
- `DBGROUP_TEST_THREAD_NUM`: The maximum number of threads to perform unit tests (default `8`).
- `DBGROUP_TEST_RANDOM_SEED`: A fixed seed value to reproduce unit tests (default `0`).
- `DBGROUP_TEST_EXEC_NUM`: The number of executions per a thread (default `1E5`).
- `DBGROUP_TEST_OVERRIDE_MIMALLOC`: Override entire memory allocation with mimalloc (default `OFF`).
    - NOTE: we use `find_package(mimalloc 1.7 REQUIRED)` to link mimalloc.

### Build and Run Unit Tests

```bash
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release -DB_TREE_BUILD_TESTS=ON \
  -DB_TREE_TEST_BUILD_PML=ON -DB_TREE_TEST_BUILD_PSL=ON \
  -DB_TREE_TEST_BUILD_OML=ON -DB_TREE_TEST_BUILD_OSL=ON ..
make -j
ctest -C Release
```

## Usage

### Linking by CMake

1. Download the files in any way you prefer (e.g., `git submodule`).

    ```bash
    cd <your_project_workspace>
    mkdir external
    git submodule add https://github.com/dbgroup-nagoya-u/b-tree.git external/b-tree
    ```

1. Add this library to your build in `CMakeLists.txt`.

    ```cmake
    add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/external/b-tree")

    add_executable(
      <target_bin_name>
      [<source> ...]
    )
    target_link_libraries(
      <target_bin_name> PRIVATE
      dbgroup::b_tree
    )
    ```

### Implemented B<sup>+</sup>-Tree Variants

We implement four variants of B<sup>+</sup>-trees.

1. Pessimistic Multi-level Locking (PML): using pessimistic locking with multi-level SMOs (i.e., lock coupling[^1]).
2. Pessimistic Single-level Locking (PSL): using pessimistic locking with single-level SMOs (i.e., B<sup>link</sup>-tree[^2]).
3. Optimistic Multi-level Locking (OML): using optimistic locking with multi-level SMOs (i.e., optimistic lock coupling[^3]).
4. Optimistic Single-level Locking (OSL): using optimistic locking with single-level SMOs (i.e., OLFIT[^4]).

Our `::dbgroup::index::b_tree::BTree` uses OSL by default, but you can select a preferred one as follows:

1. `::dbgroup::index::b_tree::BTreePML`
2. `::dbgroup::index::b_tree::BTreePSL`
3. `::dbgroup::index::b_tree::BTreeOML`
4. `::dbgroup::index::b_tree::BTreeOSL`

Note that our B<sup>+</sup>-trees optimize their node layout with fixed-length data (i.e., non-string keys). If you want to control the node layout, please use the following aliases:

- `::dbgroup::index::b_tree::BTreePMLVarLen`
- `::dbgroup::index::b_tree::BTreePMLFixLen`
- `::dbgroup::index::b_tree::BTreePSLVarLen`
- `::dbgroup::index::b_tree::BTreePSLFixLen`
- `::dbgroup::index::b_tree::BTreeOMLVarLen`
- `::dbgroup::index::b_tree::BTreeOMLFixLen`
- `::dbgroup::index::b_tree::BTreeOSLVarLen`
- `::dbgroup::index::b_tree::BTreeOSLFixLen`

### Read/Write APIs

We provide the same read/write APIs for the implemented indexes. See [here](https://github.com/dbgroup-nagoya-u/index-benchmark/wiki/Common-APIs-for-Index-Implementations) for common APIs and usage examples.

## Acknowledgments

This work is based on results from project JPNP16007 commissioned by the New Energy and Industrial Technology Development Organization (NEDO), and it was supported partially by KAKENHI (JP20K19804, JP21H03555, and JP22H03594).

[^1]: [Theodore Johnson and Dennis Sasha, "The Performance of Current B-Tree Algorithms," ACM TODS Vol. 18, No. 1, pp. 51-101, 1993.](https://doi.org/10.1145/151284.151286)
[^2]: [Philip L. Lehman and S. Bing Yao, "Efficient Locking for Concurrent Operations on B-Trees," ACM TODS Vol. 6, No. 4, pp. 650-670, 1981.](https://doi.org/10.1145/319628.319663)
[^3]: [Viktor Leis, Florian Scheibner, Alfons Kemper, and Thomas Neumann. "The ART of Practical Synchronization," In Proc. DaMoN, Article No. 3, 2016.](https://doi.org/10.1145/2933349.2933352)
[^4]: [Sang K. Cha, Sangyong Hwang, Kihong Kim, and Keunjoo Kwon, "Cache-Conscious Concurrency Control of Main-Memory Indexes on Shared-Memory Multiprocessor Systems," In Proc. VLDB, 181-190, 2001.](https://dl.acm.org/doi/10.5555/645927.672375)
