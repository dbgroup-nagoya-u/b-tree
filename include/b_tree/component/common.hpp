
/*
 * Copyright 2022 Database Group, Nagoya University
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

#include <cstring>
#include <memory>

#include "b_tree/utility.hpp"

namespace dbgroup::index::b_tree::component
{

/*######################################################################################
 * Internal enum and constants
 *####################################################################################*/

/**
 * @brief Internal return codes to represent results of node modification.
 *
 */
enum NodeRC {
  kCompleted = 0,
  kHasSpace,
  kKeyNotInserted = -6,
  kKeyAlreadyDeleted,
  kKeyAlreadyInserted,
  kNeedCleanUp,
  kNeedSplit,
  kNeedMerge,
};

/*######################################################################################
 * Internal utility functions
 *####################################################################################*/

/**
 * @brief Shift a memory address by byte offsets.
 *
 * @param addr an original address.
 * @param offset an offset to shift.
 * @return void* a shifted address.
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
