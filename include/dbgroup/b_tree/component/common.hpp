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

#ifndef B_TREE_DBGROUP_B_TREE_COMPONENT_COMMON_HPP_
#define B_TREE_DBGROUP_B_TREE_COMPONENT_COMMON_HPP_

// C++ standard libraries
#include <cstddef>
#include <cstdint>
#include <cstring>

// external libraries
#include <dbgroup/constants.hpp>
#include <dbgroup/index/utility.hpp>
#include <dbgroup/memory/utility.hpp>

// local sources
#include "dbgroup/b_tree/utility.hpp"

namespace dbgroup::index::b_tree::component
{
/*############################################################################*
 * Internal constants
 *############################################################################*/

/// @brief The size of record metadata.
constexpr uint32_t kMetaSize = kWordSize;

/// @brief The size of pointers.
constexpr uint32_t kPtrSize = kWordSize;

/// @brief The length of node header.
constexpr uint32_t kHeaderSize = 32;

/*############################################################################*
 * Internal types
 *############################################################################*/

/**
 * @brief Internal return codes for representing results of node modification.
 *
 */
enum NodeRC {
  kCompleted = 0,
  kKeyNotInserted = -100,
  kKeyAlreadyDeleted,
  kKeyAlreadyInserted,
  kNeedSplit,
  kNeedMerge,
  kAbortMerge,
  kNeedRetry,
};

/*############################################################################*
 * Global utility functions
 *############################################################################*/

/**
 * @tparam T A target class.
 * @return The maximum size for storing a given class.
 */
template <class T>
constexpr auto
MaxSize() noexcept  //
    -> uint32_t
{
  return IsVarLenData<T>() ? kMaxVarDataSize : sizeof(T);
}

}  // namespace dbgroup::index::b_tree::component

#endif  // B_TREE_DBGROUP_B_TREE_COMPONENT_COMMON_HPP_
