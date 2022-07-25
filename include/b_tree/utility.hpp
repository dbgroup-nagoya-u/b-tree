
/*
 * Copyright 2021 Database Group, Nagoya University
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

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <cstring>

namespace dbgroup::index::b_tree
{
/*######################################################################################
 * Global constants
 *####################################################################################*/

/// Assumes that one word is represented by 8 bytes.
constexpr size_t kWordSize = sizeof(uintptr_t);

/// Header length in bytes.
constexpr size_t kHeaderLength = 32;

/*######################################################################################
 * Utility enum and classes
 *####################################################################################*/

/**
 * @brief Return codes for B_TREE.
 *
 */
enum ReturnCode {
  kSuccess = 0,
  kKeyNotExist = -2,
  kKeyExist,
};

/**
 * @brief Compare binary keys as CString. The end of every key must be '\\0'.
 *
 */
struct CompareAsCString {
  constexpr auto
  operator()(const void *a, const void *b) const noexcept  //
      -> bool
  {
    if (a == nullptr) return false;
    if (b == nullptr) return true;
    return strcmp(static_cast<const char *>(a), static_cast<const char *>(b)) < 0;
  }
};

/**
 * @tparam T a target class.
 * @retval true if a target class is variable-length data.
 * @retval false if a target class is static-length data.
 */
template <class T>
constexpr auto
IsVariableLengthData()  //
    -> bool
{
  static_assert(std::is_trivially_copyable_v<T>);
  return false;
}

/*######################################################################################
 * Tuning parameters for B-tree
 *####################################################################################*/

/// The default page size of each node
constexpr size_t kPageSize = B_TREE_PAGE_SIZE;

/// The maximum size of deleted space size for invoking split
constexpr size_t kMaxDeletedSpaceSize = B_TREE_MAX_DELETED_SPACE_SIZE;

/// The minimum size of free space size for invoking split
constexpr size_t kMinFreeSpaceSize = B_TREE_MIN_FREE_SPACE_SIZE;

/// The minimum size of used space size for invoking merge
constexpr size_t kMinUsedSpaceSize = B_TREE_MIN_USED_SPACE_SIZE;

/// The maximum size of variable-length data
constexpr size_t kMaxVarDataSize = B_TREE_MAX_VARLEN_DATA_SIZE;

// Check whether the specified page size is valid
static_assert(kPageSize % kWordSize == 0);
static_assert(kMaxVarDataSize * 2 < kPageSize);

}  // namespace dbgroup::index::b_tree

#endif
