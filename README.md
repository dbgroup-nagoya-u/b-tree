# B<sup>+</sup>-trees

[![Ubuntu 24.04](https://github.com/dbgroup-nagoya-u/b-tree/actions/workflows/ubuntu_24.yaml/badge.svg)](https://github.com/dbgroup-nagoya-u/b-tree/actions/workflows/ubuntu_24.yaml) [![Ubuntu 22.04](https://github.com/dbgroup-nagoya-u/b-tree/actions/workflows/ubuntu_22.yaml/badge.svg)](https://github.com/dbgroup-nagoya-u/b-tree/actions/workflows/ubuntu_22.yaml) [![macOS](https://github.com/dbgroup-nagoya-u/b-tree/actions/workflows/mac.yaml/badge.svg)](https://github.com/dbgroup-nagoya-u/b-tree/actions/workflows/mac.yaml)

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

### Prerequisites

```bash
sudo apt update && sudo apt install -y build-essential cmake
```

### Build Options

#### Tuning Parameters

- `B_TREE_PAGE_SIZE`: Page size in bytes (default `1024`).
- `B_TREE_MAX_VARLEN_DATA_SIZE`: The expected maximum size of a variable-length data (default `32`).

#### Build Options for Unit Testing

Please refer to [index-fixtures](https://github.com/dbgroup-nagoya-u/index-fixtures?tab=readme-ov-file#build-options) for more build options.

- `B_TREE_BUILD_TESTS`: Build unit tests for this library if `ON` (default `OFF`).

### Build and Run Unit Tests

```bash
mkdir build && cd build
cmake .. \
  -DCMAKE_BUILD_TYPE=Release \
  -DB_TREE_BUILD_TESTS=ON
cmake --build . --parallel --config Release
ctest -C Release
```

## Usage

### Linking by CMake

Add this library to your build in `CMakeLists.txt`.

```cmake
FetchContent_Declare(
    b-tree
    GIT_REPOSITORY "https://github.com/dbgroup-nagoya-u/b-tree.git"
    GIT_TAG "<commit_tag_you_want_to_use>"
)
FetchContent_MakeAvailable(b-tree)

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

1. Pessimistic Single-level Locking (PSL): using pessimistic locking with single-level SMOs (i.e., B<sup>link</sup>-tree[^2]).
2. Optimistic Single-level Locking (OSL): using optimistic locking with single-level SMOs (i.e., OLFIT[^4]).

Our `::dbgroup::index::b_tree::BTree` uses OSL by default, but you can select a preferred one as follows:

1. `::dbgroup::index::b_tree::BTreePSL`
2. `::dbgroup::index::b_tree::BTreeOSL`

### Read/Write APIs

We provide the same read/write APIs for the implemented indexes. See [here](https://github.com/dbgroup-nagoya-u/index-benchmark/wiki/Common-APIs-for-Index-Implementations) for common APIs and usage examples.

## Acknowledgments

This work is based on results from project JPNP16007 commissioned by the New Energy and Industrial Technology Development Organization (NEDO), and it was supported partially by KAKENHI (JP20K19804, JP21H03555, and JP22H03594).

<!-- cSpell:disable -->
[^1]: [Theodore Johnson and Dennis Sasha, "The Performance of Current B-Tree Algorithms," ACM TODS Vol. 18, No. 1, pp. 51-101, 1993.](https://doi.org/10.1145/151284.151286)
[^2]: [Philip L. Lehman and S. Bing Yao, "Efficient Locking for Concurrent Operations on B-Trees," ACM TODS Vol. 6, No. 4, pp. 650-670, 1981.](https://doi.org/10.1145/319628.319663)
[^3]: [Viktor Leis, Florian Scheibner, Alfons Kemper, and Thomas Neumann. "The ART of Practical Synchronization," In Proc. DaMoN, Article No. 3, 2016.](https://doi.org/10.1145/2933349.2933352)
[^4]: [Sang K. Cha, Sangyong Hwang, Kihong Kim, and Keunjoo Kwon, "Cache-Conscious Concurrency Control of Main-Memory Indexes on Shared-Memory Multiprocessor Systems," In Proc. VLDB, 181-190, 2001.](https://dl.acm.org/doi/10.5555/645927.672375)
<!-- cSpell:enable -->
