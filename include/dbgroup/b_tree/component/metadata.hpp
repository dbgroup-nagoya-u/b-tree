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

#ifndef B_TREE_DBGROUP_B_TREE_COMPONENT_METADATA_HPP_
#define B_TREE_DBGROUP_B_TREE_COMPONENT_METADATA_HPP_

// C++ standard libraries
#include <cstddef>
#include <cstdint>

// external libraries
#include <dbgroup/constants.hpp>

// local sources
#include "dbgroup/b_tree/component/common.hpp"

namespace dbgroup::index::b_tree::component
{
/**
 * @brief A struct for representing record metadata.
 *
 */
struct Metadata {
  /*##########################################################################*
   * Public constructors and assignment operators
   *##########################################################################*/

  /**
   * @brief Construct an empty object.
   *
   */
  constexpr Metadata() noexcept = default;

  /**
   * @brief Construct a new insert/update-delta metadata object.
   *
   * @param offset The offset of a record position.
   * @param key_len The length of a key.
   * @param rec_len The length of a record.
   */
  constexpr Metadata(  //
      const size_t offset,
      const size_t key_len,
      const size_t rec_len) noexcept
      : offset{offset}, key_len{key_len}, rec_len{rec_len}
  {
  }

  constexpr Metadata(const Metadata &) = default;
  constexpr Metadata(Metadata &&) = default;

  constexpr auto operator=(const Metadata &) -> Metadata & = default;
  constexpr auto operator=(Metadata &&) -> Metadata & = default;

  /*##########################################################################*
   * Public destructors
   *##########################################################################*/

  /**
   * @brief Destroy the Metadata object.
   *
   */
  ~Metadata() = default;

  /*##########################################################################*
   * Public APIs
   *##########################################################################*/

  /**
   * @return The offset to a payload in bytes.
   */
  [[nodiscard]]
  constexpr auto
  GetPayOff() const noexcept  //
      -> size_t
  {
    return offset + key_len;
  }

  /**
   * @return The size of a payload in bytes.
   */
  [[nodiscard]]
  constexpr auto
  GetPayLen() const noexcept  //
      -> size_t
  {
    return rec_len - key_len;
  }

  /*##########################################################################*
   * Public member variables
   *##########################################################################*/

  /// @brief A flag for indicating whether a record is deleted.
  uint64_t deleted : 1 {};

  /// @brief An offset to a corresponding record.
  uint64_t offset : 31 {};

  /// @brief Length of a key in a corresponding record.
  uint64_t key_len : 16 {};

  /// @brief The total length of a corresponding record.
  uint64_t rec_len : 16 {};
};

/*############################################################################*
 * Static assertions
 *############################################################################*/

// The size of metadata is 8 byte.
static_assert(sizeof(Metadata) == kMetaSize);

}  // namespace dbgroup::index::b_tree::component

#endif  // B_TREE_DBGROUP_B_TREE_COMPONENT_METADATA_HPP_
