
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

#ifndef B_TREE_COMPONENT_COMMON_HPP
#define B_TREE_COMPONENT_COMMON_HPP

// C++ standard libraries
#include <chrono>
#include <cstring>
#include <functional>
#include <memory>

// local sources
#include "b_tree/utility.hpp"

namespace dbgroup::index::b_tree::component
{

/*######################################################################################
 * Internal enum and constants
 *####################################################################################*/

/**
 * @brief Internal return codes for representing results of node modification.
 *
 */
enum NodeRC {
  kCompleted = 0,
  kKeyNotInserted = -7,
  kKeyAlreadyDeleted,
  kKeyAlreadyInserted,
  kNeedSplit,
  kNeedMerge,
  kAbortMerge,
  kNeedRetry,
};

/// an expected maximum height of a tree.
constexpr size_t kExpectedTreeHeight = 8;

/// leave free space for later modifications.
constexpr size_t kNodeCapacityForBulkLoading = kPageSize * 0.9;

/// a sleep time for retrying.
constexpr auto kRetryWait = std::chrono::microseconds{10};

/// a flag for indicating leaf nodes.
constexpr uint32_t kLeafFlag = 0;

/// a flag for indicating internal nodes.
constexpr uint32_t kInnerFlag = 1;

/// a flag for indicating closed-interval.
constexpr bool kClosed = true;

/// the chacheline size for memory alignment.
constexpr std::align_val_t kCacheAlignVal = static_cast<std::align_val_t>(kCacheLineSize);

/// The alignment size for internal pages.
constexpr size_t kPageAlign = kPageSize < kVMPageSize ? kPageSize : kVMPageSize;

/*######################################################################################
 * Internal utility classes
 *####################################################################################*/

/**
 * @brief A dummy struct for representing internal pages.
 *
 */
struct alignas(kPageAlign) Page : public ::dbgroup::memory::DefaultTarget {
  // reuse pages
  static constexpr bool kReusePages = true;

  /// @brief A dummy member variable to ensure the page size.
  uint8_t dummy[kPageSize];
};

/*######################################################################################
 * Internal utility functions
 *####################################################################################*/

/**
 * @brief Shift a memory address by byte offsets.
 *
 * @param addr an original address.
 * @param offset an offset to shift.
 * @return a shifted address.
 */
constexpr auto
ShiftAddr(  //
    const void *addr,
    const int64_t offset)  //
    -> void *
{
  return static_cast<std::byte *>(const_cast<void *>(addr)) + offset;
}

}  // namespace dbgroup::index::b_tree::component

#endif
