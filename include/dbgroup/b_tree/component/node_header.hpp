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

#ifndef B_TREE_DBGROUP_B_TREE_COMPONENT_NODE_HEADER_HPP_
#define B_TREE_DBGROUP_B_TREE_COMPONENT_NODE_HEADER_HPP_

// C++ standard libraries
#include <bit>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <new>
#include <sstream>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>

// our libraries
#include <dbgroup/constants.hpp>
#include <dbgroup/index/utility.hpp>
#include <dbgroup/lock/utility.hpp>
#include <dbgroup/memory/utility.hpp>

// local sources
#include "dbgroup/b_tree/component/common.hpp"
#include "dbgroup/b_tree/component/metadata.hpp"
#include "dbgroup/b_tree/utility.hpp"

// we intentionally use a zero-length array for record metadata
#pragma GCC diagnostic ignored "-Warray-bounds"
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic ignored "-Wzero-length-bounds"
#endif

namespace dbgroup::index::b_tree::component
{
/**
 * @brief A class for representing the common header of B+tree nodes.
 *
 * @tparam Key A target key class.
 * @tparam Comp A comparator class for keys.
 * @tparam Lock A lock class for concurrency controls.
 * @note Both inner and leaf nodes inherit this class, and traversal code
 * accesses either kind of node through this class. Since both kinds of nodes
 * use the same slotted-page layout, this class also provides the common
 * utilities for manipulating records.
 */
template <class Key,  //
          class Comp,
          ::dbgroup::lock::OptimisticallyLockable Lock>
class NodeHeader
{
 protected:
  /*##########################################################################*
   * Type aliases
   *##########################################################################*/

  /// @brief A type for allocating variable-length data.
  using KeyWOPtr = std::remove_pointer_t<Key>;

  using SIXGuard = typename Lock::SIXGuard;
  using XGuard = typename Lock::XGuard;

 public:
  /*##########################################################################*
   * Public constructors and assignment operators
   *##########################################################################*/

  /**
   * @brief Construct an empty node.
   *
   * @param page_size The page size for each node.
   * @param align_val The alignment size for each node.
   * @param level The level where this node is.
   */
  constexpr NodeHeader(  //
      const uint32_t page_size,
      const uint32_t align_val,
      const uint32_t level) noexcept
      : page_size_{page_size}
      , align_val_{align_val}
      , offset_{page_size}
      , level_{static_cast<uint8_t>(level)}
  {
    static_assert(                          //
        sizeof(NodeHeader) <= kHeaderSize,  //
        "Header regions must be in 32 bytes.");
  }

  /**
   * @brief Construct split nodes.
   *
   * @param[in] l_node A left (source) node.
   * @param[in,out] l_grd The SIX-lock guard of the left node, which is released
   * after splitting.
   * @param[out] key A separator key between the split nodes.
   * @param[out] key_len The length of the separator key.
   * @note This is a split right node.
   * @note The left node is X-locked only while the split-left records are
   * copied back to it.
   */
  explicit NodeHeader(  //
      NodeHeader* const l_node,
      SIXGuard& l_grd,
      Key& key,
      uint32_t& key_len)
      : page_size_{l_node->page_size_}
      , align_val_{l_node->align_val_}
      , offset_{static_cast<uint32_t>(page_size_ - l_node->hk_len_)}
      , level_{l_node->level_}
      , leftmost_{false}
      , sib_node_{l_node->sib_node_}
  {
    auto* const tmp = AllocNode();
    tmp->offset_ = page_size_ - kMaxKeySize;

    // copy split-left records to a temporary region
    const double sep_ratio = l_node->leftmost_ && sib_node_     ? 0.25  // NOLINTBEGIN
                             : !l_node->leftmost_ && !sib_node_ ? 0.75
                                                                : 0.5;  // NOLINTEND
    const uint32_t sep_size = (l_node->usage_ - l_node->hk_len_) * sep_ratio;
    int32_t pos = 0;
    tmp->CopyLowRecord(l_node, pos);
    tmp->template CopyRecords<kSplit>(l_node, pos, l_node->rec_num_, sep_size);
    const auto sep_meta = l_node->meta_arr_[pos];
    tmp->CopyHighKey(ShiftAddr(l_node, sep_meta.offset), sep_meta.key_len);
    std::tie(key, key_len) = tmp->GetHighKey();  // keep the separator key

    // copy split-right records to this node
    CopyHighKey(ShiftAddr(l_node, offset_), l_node->hk_len_);
    CopyRecords(l_node, pos, l_node->rec_num_);

    // copy the split-left records to the left node
    auto&& x_grd = l_grd.UpgradeToX();
    l_node->CopyFrom(tmp);
    l_node->sib_node_ = this;
    VerIncrement<kSMOMask>(x_grd);
    x_grd = {};

    FreeNode(tmp);
  }

  NodeHeader(const NodeHeader&) = delete;
  NodeHeader(NodeHeader&&) = delete;

  auto operator=(const NodeHeader&) -> NodeHeader& = delete;
  auto operator=(NodeHeader&&) -> NodeHeader& = delete;

  /*##########################################################################*
   * Public getters
   *##########################################################################*/

  /**
   * @return The level where this node is.
   * @note This parameter is immutable (readable without locks).
   */
  [[nodiscard]]
  constexpr auto
  Level() const noexcept  //
      -> uint32_t
  {
    return level_;
  }

  /**
   * @return The number of records in this node.
   */
  [[nodiscard]]
  constexpr auto
  RecNum() const noexcept  //
      -> uint32_t
  {
    return rec_num_;
  }

  /**
   * @return The current memory usage.
   */
  [[nodiscard]]
  constexpr auto
  Usage() const noexcept  //
      -> uint32_t
  {
    return kHeaderSize + usage_;
  }

  /**
   * @retval true if this node has been removed.
   * @retval false otherwise.
   */
  [[nodiscard]]
  constexpr auto
  Removed() const noexcept  //
      -> bool
  {
    return removed_;
  }

  /**
   * @return A sibling node.
   */
  [[nodiscard]]
  constexpr auto
  GetSibNode() const noexcept  //
      -> NodeHeader*
  {
    return sib_node_;
  }

  /**
   * @retval 1st: A separator key.
   * @retval 2nd: The length of a separator key.
   */
  [[nodiscard]]
  auto
  GetHighKey() const  //
      -> std::pair<Key, uint32_t>
  {
    assert(hk_len_ > 0);

    const auto* const src_addr = ShiftAddr(this, page_size_ - hk_len_);
    if constexpr (IsVarLenData<Key>()) {
      auto* const key = ::dbgroup::memory::Allocate<KeyWOPtr>(hk_len_);
      std::memcpy(key, src_addr, hk_len_);
      return {key, hk_len_};
    } else {
      Key key{};
      std::memcpy(&key, src_addr, sizeof(Key));
      return {key, sizeof(Key)};
    }
  }

  /**
   * @param key A search key.
   * @param reverse A flag for searching in reverse order.
   * @retval true if this node can include a given key.
   * @retval false otherwise.
   */
  [[nodiscard]]
  auto
  Include(  //
      const Key& key,
      const bool reverse = false) const noexcept  //
      -> bool
  {
    if (hk_len_ == 0) return true;

    Key h_key;
    const auto* const src_addr = ShiftAddr(this, page_size_ - hk_len_);
    if constexpr (IsVarLenData<Key>()) {
      alignas(alignof(KeyWOPtr)) thread_local std::byte tls_key[kMaxKeySize];
      h_key = std::bit_cast<Key>(&tls_key);
      std::memcpy(h_key, src_addr, hk_len_);
    } else {
      std::memcpy(&h_key, src_addr, sizeof(Key));
    }
    return kLess(key, h_key) || (reverse && !kLess(h_key, key));
  }

  /**
   * @return A guard instance for this node.
   */
  [[nodiscard]]
  auto
  GetGuard()  //
      -> Lock::OptGuard
  {
    return lock_.GetVersion();
  }

  /*##########################################################################*
   * Public SMO APIs
   *##########################################################################*/

  /**
   * @retval true if this node is a redundant root node.
   * @retval false otherwise.
   */
  [[nodiscard]]
  constexpr auto
  IsRemovable() const noexcept  //
      -> bool
  {
    return level_ > 0 && leftmost_ && !sib_node_ && rec_num_ == 1;
  }

  /**
   * @return The SIX guard of a right sibling node.
   */
  auto
  GetMergeableSib()  //
      -> SIXGuard
  {
    SIXGuard six_grd{};
    if (sib_node_) {
      auto&& grd = sib_node_->GetGuard();
      do {
        const auto merged_usage = usage_ - hk_len_ + sib_node_->usage_;
        assert(merged_usage < 2 * page_size_);
        if (kHeaderSize + merged_usage >= page_size_ * 0.75) break;  // NOLINT
        six_grd = grd.TryLockSIX();
      } while (!six_grd);
    }
    return six_grd;
  }

  /**
   * @brief Merge a given node into this node.
   *
   * @param[in,out] six_grd The SIX-lock guard of this node, which is released
   * after merging.
   * @param[in] r_node A right-merged (to be removed) node.
   * @param[in,out] r_grd The SIX-lock guard of a right node, which is released
   * after merging.
   * @note A given node must be in the same level as this node.
   */
  void
  Merge(  //
      SIXGuard& six_grd,
      NodeHeader* const r_node,
      SIXGuard& r_grd)
  {
    assert(level_ == r_node->level_);

    auto* const tmp = AllocNode();
    tmp->offset_ = page_size_ - r_node->hk_len_;

    int32_t pos = 0;
    tmp->CopyHighKey(ShiftAddr(r_node, tmp->offset_), r_node->hk_len_);
    tmp->CopyLowRecord(this, pos);
    tmp->CopyRecords(this, pos, rec_num_);
    pos = 0;
    tmp->CopyRecords(r_node, pos, r_node->rec_num_);
    {
      auto&& x_grd = six_grd.UpgradeToX();
      CopyFrom(tmp);
      sib_node_ = r_node->sib_node_;
      r_node->Remove(r_grd, this);
      VerIncrement<kSMOMask>(x_grd);
    }

    FreeNode(tmp);
  }

  /**
   * @brief Set a remove flag to this node.
   *
   * @param[in,out] six_grd The SIX-lock guard of this node, which is released
   * after setting the flag.
   * @param[in] next A side link for traversing to the next target node.
   */
  void
  Remove(  //
      SIXGuard& six_grd,
      NodeHeader* const next = nullptr)
  {
    auto&& x_grd = six_grd.UpgradeToX();
    removed_ = true;
    sib_node_ = next;
    VerIncrement<kSMOMask>(x_grd);
  }

  /*##########################################################################*
   * Public bulkload APIs
   *##########################################################################*/

  /**
   * @brief Link this node with a right sibling node.
   *
   * @param sib_node A right sibling node.
   * @param key A highest key.
   * @param key_len The length of a highest key.
   */
  void
  LinkSiblingNode(  //
      NodeHeader* const sib_node,
      const Key& key,
      const size_t key_len) noexcept
  {
    sib_node_ = sib_node;
    hk_len_ = key_len;
    usage_ += key_len;
    std::memcpy(ShiftAddr(this, page_size_ - key_len), GetSrcAddr(key), key_len);
  }

  /**
   * @brief Link border nodes between partial trees.
   *
   * @param l_node The top border node in a left tree.
   * @param r_node The top border node in a right tree.
   */
  static void
  LinkVerticalBorderNodes(  //
      NodeHeader* l_node,
      NodeHeader* r_node) noexcept
  {
    if (!l_node) return;

    while (true) {
      // link the border nodes
      const auto meta = r_node->meta_arr_[0];
      const auto key_len = meta.key_len;
      l_node->sib_node_ = r_node;
      l_node->hk_len_ = key_len;
      l_node->usage_ += key_len;
      auto* const hk_addr = ShiftAddr(l_node, l_node->page_size_ - key_len);
      std::memcpy(hk_addr, ShiftAddr(r_node, meta.offset), key_len);
      if (l_node->level_ == 0) return;  // all the border nodes are linked

      // go down to the lower level
      l_node->CopyPayloadTo(l_node->rec_num_ - 1, &l_node);
      r_node->CopyPayloadTo(0, &r_node);
    }
  }

  /**
   * @brief Remove the leftmost keys from the leftmost nodes.
   *
   * @param node A root node.
   */
  static void
  RemoveLeftmostKeys(  //
      NodeHeader* node) noexcept
  {
    while (true) {
      node->leftmost_ = true;
      if (node->level_ == 0) break;

      auto& meta = node->meta_arr_[0];
      meta = Metadata{meta.GetPayOff(), 0, kWordSize};
      node->CopyPayloadTo(0, &node);
    }
  }

  /*##########################################################################*
   * Public utilities
   *##########################################################################*/

  /**
   * @return The explanation of this instance.
   */
  [[nodiscard]]
  auto
  ToStr() const  //
      -> std::string
  {
    std::stringstream ss{};
    ss << "address    : " << this << "\n"
       << "  level    : " << std::to_string(level_) << "\n"
       << "  sib_node : " << sib_node_ << "\n"
       << "  rec_num  : " << std::to_string(rec_num_) << "\n"
       << "  usage    : " << std::to_string(usage_) << "\n"
       << "  offset   : " << std::to_string(offset_) << "\n"
       << "  removed  : " << removed_ << "\n"
       << "  leftmost : " << leftmost_ << "\n";
    if (hk_len_ > 0) {
      Key h_key;
      const auto* const src_addr = ShiftAddr(this, page_size_ - hk_len_);
      if constexpr (IsVarLenData<Key>()) {
        alignas(alignof(KeyWOPtr)) thread_local std::byte tls_key[kMaxKeySize];
        h_key = std::bit_cast<Key>(&tls_key);
        std::memcpy(h_key, src_addr, hk_len_);
      } else {
        std::memcpy(&h_key, src_addr, sizeof(Key));
      }

      if constexpr (std::is_arithmetic_v<Key>) {
        ss << "  h_key    : " << std::to_string(h_key) << "\n";
      } else {
        ss << "  h_key    : " << h_key << "\n";
      }
    }
    return ss.str();
  }

 protected:
  /*##########################################################################*
   * Protected destructors
   *##########################################################################*/

  /**
   * @brief Destroy the node object.
   *
   */
  ~NodeHeader() = default;

  /*##########################################################################*
   * Internal constants
   *##########################################################################*/

  /// @brief A comparator instance.
  static constexpr Comp kLess{};

  /// @brief A flag for indicating whether the current SMO is splitting.
  static constexpr bool kSplit = true;

  /// @brief The maximum size of keys.
  static constexpr uint32_t kMaxKeySize = MaxKeySize<Key>();

  /*##########################################################################*
   * Internal utilities
   *##########################################################################*/

  /**
   * @brief Get the position of a search key using binary search.
   *
   * @param key A search key.
   * @param begin_pos The begin position of the search space.
   * @retval 1st: true if a target key is found.
   * @retval 1st: false otherwise.
   * @retval 2nd: The position of a found record.
   * @note If no specified key is in this node, this returns the minimum
   * position greater than the specified key.
   */
  [[nodiscard]]
  auto
  SearchRecordImpl(  //
      const Key& key,
      const int32_t begin_pos) const noexcept  //
      -> std::pair<bool, int32_t>
  {
    auto found = false;
    auto pos = begin_pos;
    int32_t end_pos = rec_num_ - 1;
    while (pos <= end_pos) {
      const int32_t cur_pos = (pos + end_pos) / 2;
      const auto meta = meta_arr_[cur_pos];
      const auto key_len = meta.key_len;
      if (meta.offset + key_len > page_size_ || key_len > kMaxKeySize) [[unlikely]] {
        break;  // read corrupted data due to OCC
      }

      const auto& index_key = GetKey(meta);
      if (kLess(key, index_key)) {
        end_pos = cur_pos - 1;
      } else if (kLess(index_key, key)) {
        pos = cur_pos + 1;
      } else {  // the target key has been found
        pos = cur_pos;
        found = true;
        break;
      }
    }
    return {found, pos};
  }

  /**
   * @brief Get the position of a search key using linear search.
   *
   * @param[in] key A search key.
   * @param[in,out] pos The position of a target record.
   * @param[out] found true if a target record is found.
   * @param[in] begin_pos The begin position of the search space.
   * @retval true if the found record has been deleted.
   * @retval false otherwise.
   * @note A given position must be in [`begin_pos`, `rec_num_`).
   */
  auto
  LinearSearchImpl(  //
      const Key& key,
      int32_t& pos,
      bool& found,
      const int32_t begin_pos) const noexcept  //
      -> bool
  {
    auto meta = meta_arr_[pos];
    const auto& idx_key = GetKey(meta);
    if (kLess(key, idx_key)) {
      found = false;
      while (pos > begin_pos) {
        meta = meta_arr_[--pos];
        const auto& idx_key = GetKey(meta);
        if (!kLess(key, idx_key)) {
          found = !kLess(idx_key, key);
          pos += static_cast<int32_t>(!found);
          break;
        }
      }
    } else if (kLess(idx_key, key)) {
      found = false;
      const int32_t rec_num = rec_num_;
      while (++pos < rec_num) {
        meta = meta_arr_[pos];
        const auto& idx_key = GetKey(meta);
        if (!kLess(idx_key, key)) {
          found = !kLess(key, idx_key);
          break;
        }
      }
    } else {
      found = true;
    }
    return meta.deleted;
  }

  /**
   * @param meta Record metadata.
   * @return The key in the target record.
   */
  [[nodiscard]]
  auto
  GetKey(                                  //
      const Metadata meta) const noexcept  //
      -> Key
  {
    Key key;
    const auto* const src_addr = ShiftAddr(this, meta.offset);
    if constexpr (IsVarLenData<Key>()) {
      alignas(alignof(KeyWOPtr)) thread_local std::byte tls_key[kMaxKeySize];
      key = std::bit_cast<Key>(&tls_key);
      std::memcpy(key, src_addr, meta.key_len);
    } else {
      std::memcpy(&key, src_addr, sizeof(Key));
    }
    return key;
  }

  /**
   * @param[in] pos A target record position.
   * @param[out] out_pay An instance for storing an output payload.
   * @param[in] pay_len The length of a target payload.
   * @note We need `pay_len` because OCC may read corrupted record metadata.
   */
  void
  CopyPayloadTo(  //
      const int32_t pos,
      void* const out_pay,
      const size_t pay_len = kPtrSize) const noexcept
  {
    const auto offset = meta_arr_[pos].GetPayOff();
    if (offset + pay_len > page_size_) [[unlikely]] {
      return;  // read corrupted data due to OCC
    }
    std::memcpy(out_pay, ShiftAddr(this, offset), pay_len);
  }

  /**
   * @param total_len The length of a new record including its metadata.
   * @retval true if this node has enough space for a new record.
   * @retval false if this node should be split before inserting a record.
   */
  [[nodiscard]]
  constexpr auto
  HasSpace(                                   //
      const size_t total_len) const noexcept  //
      -> bool
  {
    return Usage() + total_len <= page_size_;
  }

  /**
   * @brief Clean up this node if its free space is fragmented.
   *
   * @param[in] total_len The length of a new record including its metadata.
   * @param[in,out] six_grd The SIX-lock guard of this node, which is
   * temporarily upgraded to an X lock if this node is cleaned up.
   * @retval true if this node has been cleaned up.
   * @retval false otherwise.
   * @note The positions of records are changed if this node includes deleted
   * records.
   */
  auto
  CleanUpIfNeeded(  //
      const size_t total_len,
      SIXGuard& six_grd)  //
      -> bool
  {
    if (kHeaderSize + kMetaSize * rec_num_ + total_len <= offset_) return false;
    CleanUp(six_grd);
    return true;
  }

  /**
   * @brief Insert a given key/payload pair into a given position.
   *
   * @param[in] pos A position where a new record is inserted.
   * @param[in] key A target key to be written.
   * @param[in] key_len The length of a target key.
   * @param[in] payload A target payload to be written.
   * @param[in] pay_len The length of a target payload.
   * @param[in,out] six_grd The SIX-lock guard of this node, which is released
   * after inserting.
   * @note This node must have enough contiguous free space.
   */
  void
  InsertRecord(  //
      const int32_t pos,
      const Key& key,
      const size_t key_len,
      const void* const payload,
      const size_t pay_len,
      SIXGuard& six_grd) noexcept
  {
    const auto rec_len = key_len + pay_len;
    const auto mov_size = kMetaSize * (rec_num_ - pos);
    auto& meta = meta_arr_[pos];

    auto&& x_grd = six_grd.UpgradeToX();
    offset_ -= rec_len;
    std::memcpy(ShiftAddr(this, offset_), GetSrcAddr(key), key_len);
    std::memcpy(ShiftAddr(this, offset_ + key_len), payload, pay_len);
    std::memmove(ShiftAddr(&meta, kMetaSize), &meta, mov_size);
    meta = Metadata{offset_, key_len, rec_len};
    ++rec_num_;
    usage_ += rec_len + kMetaSize;
    VerIncrement<kInsDelMask>(x_grd);
  }

  /**
   * @retval true if this node should be merged.
   * @retval false otherwise.
   */
  [[nodiscard]]
  constexpr auto
  NeedMerge() const noexcept  //
      -> bool
  {
    return Usage() < page_size_ / k8;
  }

  /**
   * @brief Bulkload records into this node.
   *
   * @tparam BulkIter An iterator class of entries.
   * @param[in,out] iter The begin position of target records.
   * @param[in] iter_end The end position of target records.
   * @param[in] max_usage The maximum usage of this node.
   */
  template <class BulkIter>
  void
  BulkloadImpl(  //
      BulkIter& iter,
      const BulkIter& iter_end,
      const uint32_t max_usage) noexcept
  {
    offset_ = page_size_ - kMaxKeySize;
    for (; iter < iter_end; ++iter) {
      const auto& [key, payload, key_len, pay_len] = ParseEntry(*iter);
      const uint32_t rec_len = key_len + pay_len;
      const uint32_t total_len = rec_len + kMetaSize;
      if (usage_ + total_len > max_usage) break;

      offset_ -= rec_len;
      std::memcpy(ShiftAddr(this, offset_), GetSrcAddr(key), key_len);
      std::memcpy(ShiftAddr(this, offset_ + key_len), GetSrcAddr(payload), pay_len);
      meta_arr_[rec_num_++] = Metadata{offset_, key_len, rec_len};
      usage_ += total_len;
    }
    leftmost_ = false;
  }

  /**
   * @return A temporary node.
   */
  [[nodiscard]]
  auto
  AllocNode()  //
      -> NodeHeader*
  {
    auto* const page = ::operator new(page_size_, static_cast<std::align_val_t>(align_val_));
    return new (page) NodeHeader{page_size_, align_val_, level_};
  }

  /**
   * @param page A page for a temporary node.
   */
  void
  FreeNode(  //
      void* const page)
  {
    ::operator delete(page, static_cast<std::align_val_t>(align_val_));
  }

  /**
   * @brief Clean up this node.
   *
   * @param[in,out] six_grd The SIX-lock guard of this node, which is
   * temporarily upgraded to an X lock while copying back records.
   */
  void
  CleanUp(  //
      SIXGuard& six_grd)
  {
    auto* const tmp = AllocNode();
    tmp->offset_ = page_size_ - hk_len_;

    int32_t pos = 0;
    tmp->CopyHighKey(ShiftAddr(this, tmp->offset_), hk_len_);
    tmp->CopyLowRecord(this, pos);
    tmp->CopyRecords(this, pos, rec_num_);

    auto&& x_grd = six_grd.UpgradeToX();
    CopyFrom(tmp);
    VerIncrement<kInsDelMask>(x_grd);
    six_grd = x_grd.DowngradeToSIX();

    FreeNode(tmp);
  }

  /**
   * @brief Copy records and their relevant data from a temporary node.
   *
   * @param[in] src A temporary node.
   */
  void
  CopyFrom(  //
      const NodeHeader* const src) noexcept
  {
    rec_num_ = src->rec_num_;
    offset_ = src->offset_;
    usage_ = src->usage_;
    hk_len_ = src->hk_len_;
    std::memcpy(meta_arr_, src->meta_arr_, kMetaSize * rec_num_);
    std::memcpy(ShiftAddr(this, offset_), ShiftAddr(src, offset_), page_size_ - offset_);
  }

  /**
   * @brief Copy a highest key from a given address.
   *
   * @param[in] src The address of a source key.
   * @param[in] key_len The length of a key.
   */
  void
  CopyHighKey(  //
      const void* const src,
      const size_t key_len) noexcept
  {
    std::memcpy(ShiftAddr(this, page_size_ - key_len), src, key_len);
    usage_ += key_len;
    hk_len_ = key_len;
  }

  /**
   * @brief Copy a record that includes the lowest key in a node.
   *
   * @param[in] src A source node page.
   * @param[in,out] pos The position of the lowest record.
   */
  void
  CopyLowRecord(  //
      const NodeHeader* const src,
      int32_t& pos) noexcept
  {
    if (src->rec_num_ == 0) return;

    const auto meta = src->meta_arr_[pos++];
    const auto rec_len = meta.rec_len;
    offset_ -= rec_len;
    std::memcpy(ShiftAddr(this, offset_), ShiftAddr(src, meta.offset), rec_len);
    meta_arr_[rec_num_++] = Metadata{meta.deleted, offset_, meta.key_len, rec_len};
    usage_ += rec_len + kMetaSize;
  }

  /**
   * @brief Copy records from a given node.
   *
   * @param[in] src A source node page.
   * @param[in,out] pos The begin position of target records.
   * @param[in] end_pos The end position of target records.
   * @param[in] sep_size The desired size of a left node.
   */
  template <bool kSearchSplitPos = false>
  void
  CopyRecords(  //
      const NodeHeader* const src,
      int32_t& pos,
      const int32_t end_pos,
      [[maybe_unused]] const uint32_t sep_size = 0) noexcept
  {
    for (; pos < end_pos; ++pos) {
      const auto meta = src->meta_arr_[pos];
      if (meta.deleted) continue;

      const auto rec_len = meta.rec_len;
      const auto total_len = rec_len + kMetaSize;
      if constexpr (kSearchSplitPos) {
        if (usage_ + (total_len / 2) > sep_size) break;
      }
      offset_ -= rec_len;
      std::memcpy(ShiftAddr(this, offset_), ShiftAddr(src, meta.offset), rec_len);
      meta_arr_[rec_num_++] = Metadata{offset_, meta.key_len, rec_len};
      usage_ += total_len;
      assert(kHeaderSize + kMetaSize * rec_num_ <= offset_);
    }
  }

  /*##########################################################################*
   * Static assertions
   *##########################################################################*/

  static_assert(  //
      sizeof(Lock) <= kWordSize,
      "The size of a lock class must be smaller than 8 bytes.");

  /*##########################################################################*
   * Internal member variables
   *##########################################################################*/

  /// @brief The page size for each node.
  uint32_t page_size_ : 17 {};

  /// @brief The alignment size for each node.
  uint32_t align_val_ : 15 {};

  /// @brief The total byte length of a record block.
  uint32_t offset_ : 17 {};

  /// @brief The number of records in this node.
  uint32_t rec_num_ : 15 {};

  /// @brief The total usage of this node (without a header).
  uint32_t usage_ : 17 {};

  /// @brief The length of a highest key.
  uint32_t hk_len_ : 15 {};

  /// @brief A level where this node is.
  uint8_t level_{};

  /// @brief A flag for indicating this node is removed from a tree.
  bool removed_{};

  /// @brief A flag for indicating the leftmost node.
  bool leftmost_{true};

  /// @brief A lock for concurrency controls.
  Lock lock_{};

  /// @brief A right sibling node.
  NodeHeader* sib_node_{};

  /// @brief An actual data block (it starts with record metadata).
  Metadata meta_arr_[0];
};

}  // namespace dbgroup::index::b_tree::component

#endif  // B_TREE_DBGROUP_B_TREE_COMPONENT_NODE_HEADER_HPP_
