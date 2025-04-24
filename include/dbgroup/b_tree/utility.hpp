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

#ifndef B_TREE_DBGROUP_B_TREE_UTILITY_HPP_
#define B_TREE_DBGROUP_B_TREE_UTILITY_HPP_

// C++ standard libraries
#include <cstddef>

namespace dbgroup::index::b_tree
{
/*############################################################################*
 * Global constants
 *############################################################################*/

/// @brief The default time interval for garbage collection [us].
constexpr size_t kDefaultGCTime = 10000;  // 10 ms

/// @brief The default number of worker threads for garbage collection.
constexpr size_t kDefaultGCThreadNum = 1;

/// @brief A flag for using single-level SMO.
constexpr bool kSingleLevelSMOFlag = true;

/// @brief A flag for using multi-level SMO.
constexpr bool kMultiLevelSMOFlag = false;

/// @brief A flag for using optimistic concurrency controls.
constexpr bool kOCCFlag = true;

/// @brief A flag for using pessimistic concurrency controls.
constexpr bool kPCCFlag = false;

/*############################################################################*
 * Tuning parameters for B+trees
 *############################################################################*/

/// @brief The default page size of each node.
constexpr size_t kPageSize = (B_TREE_PAGE_SIZE);

/// @brief The maximum size of variable-length data.
constexpr size_t kMaxVarDataSize = (B_TREE_MAX_VARLEN_DATA_SIZE);

/*############################################################################*
 * Global types
 *############################################################################*/

/**
 * @brief A struct for representing internal pages.
 *
 */
using Page = ::dbgroup::memory::PageTarget<kPageSize>;

}  // namespace dbgroup::index::b_tree

#endif  // B_TREE_DBGROUP_B_TREE_UTILITY_HPP_
