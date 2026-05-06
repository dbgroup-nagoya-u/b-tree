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
constexpr size_t kDefaultGCTime = 100;  // 100 ms

/// @brief The default number of worker threads for garbage collection.
constexpr size_t kDefaultGCThreadNum = 1;

/// @brief A bit mask for using all bits.
constexpr uint32_t kNoMask = ~0U;

/// @brief A bit mask for extracting insert/delete versions.
constexpr uint32_t kInsDelMask = 0xFFFF'F000U;

/// @brief A bit mask for extracting SMO versions.
constexpr uint32_t kSMOMask = 0xFF00'0000U;

/*############################################################################*
 * Tuning parameters for B+trees
 *############################################################################*/

/// @brief The maximum size of variable-length data.
constexpr size_t kMaxVarDataSize = (B_TREE_MAX_VARLEN_DATA_SIZE);

}  // namespace dbgroup::index::b_tree

#endif  // B_TREE_DBGROUP_B_TREE_UTILITY_HPP_
