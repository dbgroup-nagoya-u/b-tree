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

#ifndef B_TREE_DBGROUP_B_TREE_COMPONENT_INNER_NODE_HPP_
#define B_TREE_DBGROUP_B_TREE_COMPONENT_INNER_NODE_HPP_

// C++ standard libraries
#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <utility>

// external libraries
#include <dbgroup/constants.hpp>
#include <dbgroup/index/utility.hpp>
#include <dbgroup/lock/utility.hpp>

// local sources
#include "dbgroup/b_tree/component/common.hpp"
#include "dbgroup/b_tree/component/metadata.hpp"
#include "dbgroup/b_tree/component/node_header.hpp"

namespace dbgroup::index::b_tree::component
{
/**
 * @brief A class for representing inner nodes in B+trees.
 *
 * @tparam Key A target key class.
 * @tparam Comp A comparator class for keys.
 * @tparam Lock A lock class for concurrency controls.
 * @note Each record in inner nodes has a pointer to a child node as its
 * payload. The first record in each node is a dummy record for the lowest
 * child, so searching for a key starts from the second record.
 */
template <class Key,  //
          class Comp,
          ::dbgroup::lock::OptimisticallyLockable Lock>
class InnerNode : public NodeHeader<Key, Comp, Lock>
{
  /*##########################################################################*
   * Type aliases
   *##########################################################################*/

  using Header = NodeHeader<Key, Comp, Lock>;

  // header members used in this class
  using Header::kMaxKeySize;
  using Header::leftmost_;
  using Header::meta_arr_;
  using Header::offset_;
  using Header::page_size_;
  using Header::rec_num_;
  using Header::sib_node_;
  using Header::usage_;
  using typename Header::SIXGuard;
  using typename Header::XGuard;

 public:
  /*##########################################################################*
   * Public constructors and assignment operators
   *##########################################################################*/

  /**
   * @brief Construct a new root node.
   *
   * @param page_size The page size for each node.
   * @param align_val The alignment size for each node.
   * @param level The level where this node is.
   * @param key A separator key.
   * @param key_len The length of a separator key.
   * @param l_child A left child node.
   * @param r_child A right child node.
   */
  InnerNode(  //
      const uint32_t page_size,
      const uint32_t align_val,
      const uint32_t level,
      const Key& key,
      const uint32_t key_len,
      const Header* l_child,
      const Header* r_child) noexcept  //
      : Header{page_size, align_val, level}
  {
    assert(level > 0);

    const uint32_t l_offset = page_size_ - kPtrSize;
    const uint32_t r_offset = page_size_ - 2 * kPtrSize;
    offset_ = r_offset - key_len;
    rec_num_ = 2;
    usage_ = 2 * (kMetaSize + kPtrSize) + key_len;

    std::memcpy(ShiftAddr(this, l_offset), &l_child, kPtrSize);
    meta_arr_[0] = Metadata{l_offset, 0, kPtrSize};

    std::memcpy(ShiftAddr(this, r_offset), &r_child, kPtrSize);
    std::memcpy(ShiftAddr(this, offset_), GetSrcAddr(key), key_len);
    meta_arr_[1] = Metadata{offset_, key_len, key_len + kPtrSize};
  }

  InnerNode(const InnerNode&) = delete;
  InnerNode(InnerNode&&) = delete;

  auto operator=(const InnerNode&) -> InnerNode& = delete;
  auto operator=(InnerNode&&) -> InnerNode& = delete;

  /*##########################################################################*
   * Public destructors
   *##########################################################################*/

  /**
   * @brief Destroy the node object.
   *
   */
  ~InnerNode() = default;

  /*##########################################################################*
   * Public getters
   *##########################################################################*/

  /**
   * @param pos The position of a target child.
   * @return A child node in a given position.
   */
  [[nodiscard]]
  auto
  GetChild(                              //
      const int32_t pos) const noexcept  //
      -> Header*
  {
    Header* child{};
    this->CopyPayloadTo(pos, &child);
    return child;
  }

  /**
   * @brief Get the position of a search key using binary search.
   *
   * @param key A search key.
   * @retval 1st: true if a target key is found.
   * @retval 1st: false otherwise.
   * @retval 2nd: The position of a found record.
   * @note If no specified key is in this node, this returns the minimum
   * position greater than the specified key.
   */
  [[nodiscard]]
  auto
  SearchRecord(                       //
      const Key& key) const noexcept  //
      -> std::pair<bool, int32_t>
  {
    return this->SearchRecordImpl(key, 1);
  }

  /**
   * @brief Get the position of a search key using linear search.
   *
   * @param[in] key A search key.
   * @param[in,out] pos The position of a target record.
   * @param[out] found true if a target record is found.
   * @retval true if the found record has been deleted.
   * @retval false otherwise.
   */
  auto
  LinearSearch(  //
      const Key& key,
      int32_t& pos,
      bool& found) const noexcept  //
      -> bool
  {
    const int32_t rec_num = rec_num_;
    if (pos >= rec_num) {
      if (rec_num == 1) {
        pos = 1;
        found = false;
        return false;
      }
      pos = rec_num - 1;
    }
    return this->LinearSearchImpl(key, pos, found, 1);
  }

  /**
   * @brief Search a child node with a given key using binary search.
   *
   * @param key A search key.
   * @return A child node.
   */
  [[nodiscard]]
  auto
  SearchChild(                        //
      const Key& key) const noexcept  //
      -> Header*
  {
    auto [found, pos] = SearchRecord(key);
    pos -= static_cast<int32_t>(!found);
    return GetChild(pos);
  }

  /*##########################################################################*
   * Write APIs
   *##########################################################################*/

  /**
   * @brief Insert a given down link.
   *
   * @param pos A position where a new record is inserted.
   * @param key A separator key to be written.
   * @param key_len The length of a separator key.
   * @param child A child node to be written.
   * @retval true if a record is written.
   * @retval false if this node should be split before inserting a record.
   */
  void
  InsertChild(  //
      const int32_t pos,
      const Key& key,
      const size_t key_len,
      const Header* child,
      SIXGuard& six_grd) noexcept
  {
    const auto total_len = key_len + kPtrSize + kMetaSize;
    if (!this->HasSpace(total_len)) return;

    // inner nodes do not have deleted records, so the position is not changed
    this->CleanUpIfNeeded(total_len, six_grd);
    this->InsertRecord(pos, key, key_len, &child, kPtrSize, six_grd);
  }

  /**
   * @brief Delete a down link physically using a given position.
   *
   * @param pos A position where a target record is.
   * @retval 1st: The deleted child node.
   * @retval 2nd: true if this node should be merged.
   * @retval 2nd: false otherwise.
   */
  auto
  DeleteChild(  //
      const int32_t pos,
      SIXGuard& six_grd) noexcept  //
      -> Header*
  {
    assert(pos > 0);

    auto* child = GetChild(pos);
    auto& meta = meta_arr_[pos];
    const auto rec_len = meta.rec_len;
    auto&& x_grd = six_grd.UpgradeToX();
    std::memmove(&meta, ShiftAddr(&meta, kMetaSize), kMetaSize * (rec_num_ - pos - 1UL));
    --rec_num_;
    usage_ -= rec_len + kMetaSize;
    VerIncrement<kInsDelMask>(x_grd);
    if (this->NeedMerge()) {
      six_grd = x_grd.DowngradeToSIX();
    }
    return child;
  }

  /*##########################################################################*
   * Public bulkload APIs
   *##########################################################################*/

  /**
   * @brief Bulkload down links into this node.
   *
   * @tparam BulkIter An iterator class of entries.
   * @param iter The begin position of target records.
   * @param iter_end The end position of target records.
   */
  template <class BulkIter>
  void
  Bulkload(  //
      BulkIter& iter,
      const BulkIter& iter_end) noexcept
  {
    const uint32_t page_size = page_size_;
    const auto min_free = std::max(kMaxKeySize, static_cast<uint32_t>(page_size / k8));
    this->BulkloadImpl(iter, iter_end, page_size - min_free);
  }
};

}  // namespace dbgroup::index::b_tree::component

#endif  // B_TREE_DBGROUP_B_TREE_COMPONENT_INNER_NODE_HPP_
