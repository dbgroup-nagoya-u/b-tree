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
#include <type_traits>

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
 * Global type aliases
 *############################################################################*/

/// @brief A data class without a pointer for variable-length data.
template <class T>
using WOPtr = std::conditional_t<IsVarLenData<T>(), std::remove_pointer_t<T>, T>;

/**
 * @brief A function for merging a payload with a given payload or delta.
 *
 * The arguments are the current payload and its length, a given payload or
 * delta and its length, and a destination buffer of `MaxPaySize<Payload>()`
 * bytes, and the return value is the length of a merged payload. Since this
 * function may be called multiple times due to retries, it must not have
 * side effects.
 *
 * @tparam Payload A payload class.
 * @tparam T A class of given payloads or deltas.
 */
template <class Payload, class T>
using Merger = size_t (*)(const WOPtr<Payload>* old_pay,
                          size_t old_len,
                          const WOPtr<T>* delta,
                          size_t delta_len,
                          WOPtr<Payload>* out_pay);

/*############################################################################*
 * Global utility functions
 *############################################################################*/

/**
 * @tparam T A target key class.
 * @return The maximum size for storing a given key class.
 */
template <class T>
constexpr auto
MaxKeySize() noexcept  //
    -> uint32_t
{
  return IsVarLenData<T>() ? kMaxVarKeySize : sizeof(T);
}

/**
 * @tparam T A target payload class.
 * @return The maximum size for storing a given payload class.
 */
template <class T>
constexpr auto
MaxPaySize() noexcept  //
    -> uint32_t
{
  return IsVarLenData<T>() ? kMaxVarPaySize : sizeof(T);
}

}  // namespace dbgroup::index::b_tree::component

#endif  // B_TREE_DBGROUP_B_TREE_COMPONENT_COMMON_HPP_
