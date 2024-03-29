
/*
 * Copyright 2023 Database Group, Nagoya University
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#ifndef B_TREE_UTILITY_HPP
#define B_TREE_UTILITY_HPP

// C++ standard libraries
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <cstring>

// external sources
#include "memory/utility.hpp"

namespace dbgroup::index::b_tree
{
/*######################################################################################
 * Global enum and constants
 *####################################################################################*/

/**
 * @brief Return codes for B+trees.
 *
 */
enum ReturnCode {
  kSuccess = 0,
  kKeyNotExist = -2,
  kKeyExist,
};

/// the default time interval for garbage collection [us].
constexpr size_t kDefaultGCTime = 10000;  // 10 ms

/// the default number of worker threads for garbage collection.
constexpr size_t kDefaultGCThreadNum = 1;

/// a flag for indicating pessimistic locking.
constexpr bool kPessimisticLock = false;

/// a flag for indicating optimistic locking.
constexpr bool kOptimisticLock = true;

/// a flag for indicating multi-layer locking.
constexpr bool kMultiLayerLock = false;

/// a flag for indicating multi-layer locking.
constexpr bool kSingleLayerLock = true;

/// a flag for indicating optimized page layouts for fixed-length data.
constexpr bool kOptimizeForFixLenData = false;

/*######################################################################################
 * Utility functions
 *####################################################################################*/

/**
 * @tparam T a target class.
 * @retval true if a target class is variable-length data.
 * @retval false if a target class is fixed-length data.
 */
template <class T>
constexpr auto
IsVarLenData()  //
    -> bool
{
  return false;
}

/**
 * @brief Use binary data as variable-length data.
 *
 */
template <>
constexpr auto
IsVarLenData<char *>()  //
    -> bool
{
  return true;
}

/**
 * @brief Use binary data as variable-length data.
 *
 */
template <>
constexpr auto
IsVarLenData<std::byte *>()  //
    -> bool
{
  return true;
}

/*######################################################################################
 * Tuning parameters for B+trees
 *####################################################################################*/

/// The default page size of each node.
constexpr size_t kPageSize = B_TREE_PAGE_SIZE;

/// The maximum size of deleted space size for invoking split.
constexpr size_t kMaxDeletedSpaceSize = B_TREE_MAX_DELETED_SPACE_SIZE;

/// The minimum size of free space size for invoking split.
constexpr size_t kMinFreeSpaceSize = B_TREE_MIN_FREE_SPACE_SIZE;

/// The minimum size of used space size for invoking merge.
constexpr size_t kMinUsedSpaceSize = B_TREE_MIN_USED_SPACE_SIZE;

/// The maximum size of variable-length data.
constexpr size_t kMaxVarLenDataSize = B_TREE_MAX_VARLEN_DATA_SIZE;

/// The page size of virtual memory address.
constexpr size_t kVMPageSize = 4096;

/// Assumes that one cache line is represented by 64 bytes.
constexpr size_t kCacheLineSize = 64;

}  // namespace dbgroup::index::b_tree

#endif
