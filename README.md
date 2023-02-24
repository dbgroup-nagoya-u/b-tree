# B<sup>+</sup>-trees

![Ubuntu-20.04](https://github.com/dbgroup-nagoya-u/b-tree/workflows/Ubuntu-20.04/badge.svg?branch=main)

This repository contains open source implementations of B<sup>+</sup>-trees for research use. The purpose of this implementation is to reproduce B<sup>+</sup>-trees and measure its performance.

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

- `B_TREE_BUILD_TESTS`: Building unit tests for this library if `ON` (default `OFF`).
- `DBGROUP_TEST_BUILD_MULTI_THREAD_TESTS`: Building multi-threading tests if `ON` (default `ON`).
- `DBGROUP_TEST_THREAD_NUM`: The maximum number of threads to perform unit tests (default `8`).
- `DBGROUP_TEST_RANDOM_SEED`: A fixed seed value to reproduce unit tests (default `0`).
- `DBGROUP_TEST_EXEC_NUM`: The number of executions per a thread (default `1E5`).
- `DBGROUP_TEST_OVERRIDE_MIMALLOC`: Override entire memory allocation with mimalloc (default `OFF`).
    - NOTE: we use `find_package(mimalloc 1.7 REQUIRED)` to link mimalloc.

### Build and Run Unit Tests

```bash
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release -DB_TREE_BUILD_TESTS=ON ..
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

### Read/Write APIs

We provide the same read/write APIs for the reproduced indexes. See [here](https://github.com/dbgroup-nagoya-u/index-benchmark/wiki/Common-APIs-for-Index-Implementations) for common APIs and usage examples.
