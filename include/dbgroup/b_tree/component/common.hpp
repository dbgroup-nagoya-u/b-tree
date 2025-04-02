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
#include <cstring>

// external libraries
#include "dbgroup/constants.hpp"
#include "dbgroup/index/utility.hpp"
#include "dbgroup/memory/utility.hpp"

// local sources
#include "dbgroup/b_tree/utility.hpp"

namespace dbgroup::index::b_tree::component
{
/*############################################################################*
 * Internal constants
 *############################################################################*/

/// @brief The size of record metadata.
constexpr size_t kMetaSize = kWordSize;

/// @brief The alignment size for internal pages.
constexpr size_t kPageAlign = kPageSize < kVMPageSize ? kPageSize : kVMPageSize;

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

/**
 * @brief A dummy struct for representing internal pages.
 *
 */
struct alignas(kPageAlign) Page : public ::dbgroup::memory::DefaultTarget {
  // filling zeros in reclaimed pages
  using T = ZeroFilling<kPageSize>;

  // reuse pages
  static constexpr bool kReusePages = true;

  /// @brief A dummy member variable to ensure the page size.
  std::byte dummy[kPageSize]{};
};

}  // namespace dbgroup::index::b_tree::component

#endif  // B_TREE_DBGROUP_B_TREE_COMPONENT_COMMON_HPP_
