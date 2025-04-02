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
#include <chrono>
#include <cstddef>

namespace dbgroup::index::b_tree
{
/*############################################################################*
 * Global enum and constants
 *############################################################################*/

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

/*############################################################################*
 * Tuning parameters for B+trees
 *############################################################################*/

/// The default page size of each node.
constexpr size_t kPageSize = (B_TREE_PAGE_SIZE);

/// The maximum size of deleted space size for invoking split.
constexpr size_t kMaxDeletedSpaceSize = (B_TREE_MAX_DELETED_SPACE_SIZE);

/// The minimum size of free space size for invoking split.
constexpr size_t kMinFreeSpaceSize = (B_TREE_MIN_FREE_SPACE_SIZE);

/// The minimum size of used space size for invoking merge.
constexpr size_t kMinUsedSpaceSize = (B_TREE_MIN_USED_SPACE_SIZE);

/// The maximum size of variable-length data.
constexpr size_t kMaxVarLenDataSize = (B_TREE_MAX_VARLEN_DATA_SIZE);

/// @brief A sleep time for backoff [us].
constexpr std::chrono::microseconds kBackOffTime{B_TREE_BACKOFF_TIME};

}  // namespace dbgroup::index::b_tree

#endif  // B_TREE_DBGROUP_B_TREE_UTILITY_HPP_
