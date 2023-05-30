
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

#ifndef B_TREE_COMPONENT_METADATA_HPP
#define B_TREE_COMPONENT_METADATA_HPP

// local sources
#include "b_tree/component/common.hpp"

namespace dbgroup::index::b_tree::component
{
/**
 * @brief A struct for representing record metadata.
 *
 */
struct Metadata {
  /*####################################################################################
   * Public constructors and assignment operators
   *##################################################################################*/

  /**
   * @brief Construct an empty object.
   *
   */
  constexpr Metadata() : is_deleted{}, offset{} {};

  /**
   * @brief Construct a new metadata object.
   *
   */
  constexpr Metadata(  //
      const size_t offset,
      const size_t key_length,
      const size_t record_length)
      : is_deleted{0},
        offset{static_cast<uint32_t>(offset)},
        key_len{static_cast<uint16_t>(key_length)},
        rec_len{static_cast<uint16_t>(record_length)}
  {
  }

  constexpr Metadata(const Metadata &) = default;
  constexpr Metadata(Metadata &&) = default;

  constexpr auto operator=(const Metadata &) -> Metadata & = default;
  constexpr auto operator=(Metadata &&) -> Metadata & = default;

  /*####################################################################################
   * Public destructors
   *##################################################################################*/

  /**
   * @brief Destroy the Metadata object.
   *
   */
  ~Metadata() = default;

  /*####################################################################################
   * Public operators
   *##################################################################################*/

  /**
   * @param comp a metadata to be compared.
   * @retval true if a given metadata is equivalent with this one.
   * @retval false otherwise.
   */
  constexpr auto
  operator==(const Metadata &comp) const  //
      -> bool
  {
    return offset == comp.offset       //
           && key_len == comp.key_len  //
           && rec_len == comp.rec_len;
  }

  /**
   * @param comp a metadata to be compared.
   * @retval true if a given metadata is different from this one.
   * @retval false otherwise.
   */
  constexpr auto
  operator!=(const Metadata &comp) const  //
      -> bool
  {
    return !(*this == comp);
  }

  /*####################################################################################
   * Internal member variables
   *##################################################################################*/

  /// a flag for indicating a corresponding record is deleted.
  uint32_t is_deleted : 1;

  /// an offset to a corresponding record.
  uint32_t offset : 31;

  /// the length of a key in a corresponding record.
  uint16_t key_len{};

  /// the total length of a corresponding record.
  uint16_t rec_len{};
};

}  // namespace dbgroup::index::b_tree::component

#endif
