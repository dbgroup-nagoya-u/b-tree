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

#ifndef B_TREE_DBGROUP_B_TREE_COMPONENT_NODE_HPP_
#define B_TREE_DBGROUP_B_TREE_COMPONENT_NODE_HPP_

// C++ standard libraries
#include <algorithm>
#include <bit>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <optional>
#include <sstream>
#include <tuple>
#include <type_traits>
#include <utility>

// external sources
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
 * @brief A class for representing nodes in B+trees.
 *
 * @tparam Key A target key class.
 * @tparam Comp A comparator class for keys.
 * @tparam Lock A lock class for concurrency controls.
 */
template <class Key,  //
          class Comp,
          ::dbgroup::lock::OptimisticallyLockable Lock>
class Node
{
  /*##########################################################################*
   * Type aliases
   *##########################################################################*/

  using ScanKey = std::optional<std::tuple<Key, size_t, bool>>;
  using SIXGuard = typename Lock::SIXGuard;

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
  constexpr explicit Node(  //
      const uint32_t page_size,
      const uint32_t align_val,
      const uint32_t level = 0) noexcept
      : page_size_{page_size}
      , align_val_{align_val}
      , level_{static_cast<uint8_t>(level)}
  {
    static_assert(                    //
        sizeof(Node) <= kHeaderSize,  //
        "Header regions must be in 32 bytes.");
  }

  /**
   * @brief Construct a new root node.
   *
   * @param level The level where this node is.
   * @param key A separator key.
   * @param key_len The length of a separator key.
   * @param l_child A left child node.
   * @param r_child A right child node.
   */
  Node(  //
      const uint32_t level,
      const Key& key,
      const uint32_t key_len,
      const Node* l_child,
      const Node* r_child) noexcept  //
      : page_size_{l_child->page_size_}
      , align_val_{l_child->align_val_}
      , offset_{static_cast<uint32_t>(page_size_ - (2 * kPtrSize + key_len))}
      , rec_num_{2}
      , usage_{static_cast<uint32_t>(2 * (kMetaSize + kPtrSize) + key_len)}
      , level_{static_cast<uint8_t>(level)}
  {
    const uint32_t l_offset = page_size_ - kPtrSize;
    const uint32_t r_offset = page_size_ - 2 * kPtrSize;

    std::memcpy(ShiftAddr(this, l_offset), &l_child, kPtrSize);
    meta_arr_[0] = Metadata{l_offset, 0, kPtrSize};

    const auto offset = r_offset - key_len;
    std::memcpy(ShiftAddr(this, r_offset), &r_child, kPtrSize);
    std::memcpy(ShiftAddr(this, offset), GetSrcAddr(key), key_len);
    meta_arr_[1] = Metadata{offset, key_len, key_len + kPtrSize};
  }

  /**
   * @brief Construct split nodes.
   *
   * @param[in,out] l_node A left (source) node.
   * @note This is a split right node.
   */
  explicit Node(  //
      Node* l_node)
      : page_size_{l_node->page_size_}
      , align_val_{l_node->align_val_}
      , offset_{static_cast<uint32_t>(page_size_ - l_node->hk_len_)}
      , level_{l_node->level_}
      , leftmost_{false}
      , sib_node_{l_node->sib_node_}
  {
    auto* tmp = AllocNode();
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

    // copy split-right records to this node
    CopyHighKey(ShiftAddr(l_node, offset_), l_node->hk_len_);
    CopyRecords(l_node, pos, l_node->rec_num_);

    // copy the split-left records to the left node
    l_node->CopyFrom(tmp);
    l_node->sib_node_ = this;

    FreeNode(tmp);
  }

  Node(const Node&) = delete;
  Node(Node&&) = delete;

  auto operator=(const Node&) -> Node& = delete;
  auto operator=(Node&&) -> Node& = delete;

  /*##########################################################################*
   * Public destructors
   *##########################################################################*/

  /**
   * @brief Destroy the node object.
   *
   */
  ~Node() = default;

  /*##########################################################################*
   * Public getters
   *##########################################################################*/

  /**
   * @return The level where this node is.
   * @note This parameter is immutable (readable without locks).
   */
  [[nodiscard]]
  constexpr auto
  GetLevel() const noexcept  //
      -> uint32_t
  {
    return level_;
  }

  /**
   * @return The number of records in this node.
   */
  [[nodiscard]]
  constexpr auto
  GetRecNum() const noexcept  //
      -> uint32_t
  {
    return rec_num_;
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
      -> Node*
  {
    return sib_node_;
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
      while (true) {
        const auto merged_usage = usage_ - hk_len_ + sib_node_->usage_;
        if (kHeaderSize + merged_usage >= page_size_ * 0.75) break;  // NOLINT
        six_grd = grd.TryLockSIX();
        if (six_grd) break;
      }
    }
    return six_grd;
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

    const auto* src_addr = ShiftAddr(this, page_size_ - hk_len_);
    if constexpr (IsVarLenData<Key>()) {
      auto* key = ::dbgroup::memory::Allocate<KeyWOPtr>(hk_len_);
      std::memcpy(key, src_addr, hk_len_);
      return {key, hk_len_};
    } else {
      Key key{};
      std::memcpy(&key, src_addr, sizeof(Key));
      return {key, sizeof(Key)};
    }
  }

  /**
   * @return The lowest key.
   */
  [[nodiscard]]
  auto
  GetLowKey() const  //
      -> Key
  {
    assert(rec_num_ > 0);

    Key key;
    const auto meta = meta_arr_[0];
    const auto* src_addr = ShiftAddr(this, meta.offset);
    if constexpr (IsVarLenData<Key>()) {
      const auto key_len = meta.key_len;
      key = ::dbgroup::memory::Allocate<KeyWOPtr>(key_len);
      std::memcpy(key, src_addr, key_len);
    } else {
      std::memcpy(&key, src_addr, sizeof(Key));
    }
    return key;
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
    const auto* src_addr = ShiftAddr(this, page_size_ - hk_len_);
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
    auto found = false;
    auto pos = static_cast<int32_t>(level_ > 0);
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
   * @brief Get the position of a search key using binary search.
   *
   * @param key A search key.
   * @retval 1st: true if a target key is found.
   * @retval 1st: false otherwise.
   * @retval 2nd: true if a target record is deleted.
   * @retval 2nd: false otherwise.
   * @retval 3rd: The position of a found record.
   * @note If no specified key is in this node, this returns the minimum
   * position greater than the specified key.
   */
  [[nodiscard]]
  auto
  CheckUniqueness(                    //
      const Key& key) const noexcept  //
      -> std::tuple<bool, bool, int32_t>
  {
    const auto [found, pos] = SearchRecord(key);
    return {found, meta_arr_[pos].deleted, pos};
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
      -> Node*
  {
    auto [found, pos] = SearchRecord(key);
    pos -= static_cast<int32_t>(!found);

    Node* child{};
    CopyPayloadTo(pos, &child);
    return child;
  }

  /**
   * @tparam Guard A desired guard type.
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
   * Read APIs
   *##########################################################################*/

  /**
   * @param pos A target record position.
   * @return The payload in the target position.
   */
  template <class Payload>
  [[nodiscard]]
  auto
  GetPayload(                            //
      const int32_t pos) const noexcept  //
      -> Payload
  {
    Payload out_pay{};
    const auto offset = meta_arr_[pos].GetPayOff();
    if (offset + sizeof(Payload) <= page_size_) [[likely]] {
      std::memcpy(&out_pay, ShiftAddr(this, offset), sizeof(Payload));
    }
    return out_pay;
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
      void* out_pay,
      const size_t pay_len = kPtrSize) const noexcept
  {
    const auto offset = meta_arr_[pos].GetPayOff();
    if (offset + pay_len > page_size_) [[unlikely]] {
      return;  // read corrupted data due to OCC
    }
    std::memcpy(out_pay, ShiftAddr(this, offset), pay_len);
  }

  /**
   * @brief Copy a record if exist.
   *
   * @param[in] pos A position to be copied.
   * @param[in,out] key A destination address for a key.
   * @param[in,out] payload A destination address for a payload.
   * @retval true if a target record is active (not deleted).
   * @retval false otherwise.
   */
  [[nodiscard]]
  auto
  CopyRecordTo(  //
      const int32_t pos,
      void* key,
      void* payload) const noexcept  //
      -> bool
  {
    const auto meta = meta_arr_[pos];
    if (meta.deleted) return false;

    std::memcpy(key, ShiftAddr(this, meta.offset), meta.key_len);
    std::memcpy(payload, ShiftAddr(this, meta.GetPayOff()), meta.GetPayLen());
    return true;
  }

  /**
   * @brief Get the end position of records for scanning and check it has been finished.
   *
   * @tparam kForward A flag for indicating forward or backward order.
   * @param[in] end_key A scan-end key.
   * @param[out] is_end A flag for indicating the last scan node.
   * @param[out] end_pos The end position of this scan operation.
   */
  template <bool kForward>
  void
  SearchEndPosition(  //
      const ScanKey& end_key,
      bool& is_end,
      int32_t& end_pos) const noexcept
  {
    if (!end_key) {
      if constexpr (kForward) {
        is_end = !sib_node_;
        end_pos = rec_num_;
      } else {
        is_end = leftmost_;
        end_pos = static_cast<int32_t>(rec_num_ == 0);
      }
      return;
    }

    const auto& [key, _, closed] = *end_key;
    if constexpr (kForward) {
      is_end = !sib_node_ || Include(key);
      end_pos = rec_num_;
    } else {
      is_end = leftmost_ || (rec_num_ > 0 && !kLess(key, GetKey(meta_arr_[0])));
      end_pos = 0;
    }
    if (is_end) {
      const auto [found, deleted, pos] = CheckUniqueness(key);
      if constexpr (kForward) {
        end_pos = pos + static_cast<int32_t>(found && !deleted && closed);
      } else {
        end_pos = pos + static_cast<int32_t>(!found || deleted || !closed);
      }
    }
  }

  /*##########################################################################*
   * Write APIs
   *##########################################################################*/

  /**
   * @brief Insert a given kay/payload pair.
   *
   * @param pos A position where a new record is inserted.
   * @param found A flag for indicating the same (deleted) key is found.
   * @param key A target key to be written.
   * @param key_len The length of a target key.
   * @param payload A target payload to be written.
   * @param pay_len The length of a target payload.
   * @retval true if a record is written.
   * @retval false if this node should be split before inserting a record.
   */
  auto
  Insert(  //
      int32_t pos,
      bool found,
      const Key& key,
      const size_t key_len,
      const void* payload,
      const size_t pay_len) noexcept  //
      -> bool
  {
    const auto rec_len = key_len + pay_len;
    const auto total_len = rec_len + kMetaSize;
    if (found) {  // reuse the deleted record
      auto& meta = meta_arr_[pos];
      std::memcpy(ShiftAddr(this, meta.GetPayOff()), payload, pay_len);
      meta.deleted = false;
      usage_ += pos > 0 ? total_len : 0;
    } else {  // insert a new record
      if (kHeaderSize + usage_ + total_len > page_size_) return false;
      if (kHeaderSize + kMetaSize * rec_num_ + total_len > offset_) {
        CleanUp();
        pos = SearchRecord(key).second;
      }
      const auto mov_size = kMetaSize * (rec_num_ - pos);
      auto& meta = meta_arr_[pos];
      offset_ -= rec_len;
      std::memcpy(ShiftAddr(this, offset_), GetSrcAddr(key), key_len);
      std::memcpy(ShiftAddr(this, offset_ + key_len), payload, pay_len);
      std::memmove(ShiftAddr(&meta, kMetaSize), &meta, mov_size);
      meta = Metadata{offset_, key_len, rec_len};
      ++rec_num_;
      usage_ += total_len;
    }
    return true;
  }

  /**
   * @brief Update a record using a given kay/payload pair.
   *
   * @param pos A position where a target record is.
   * @param payload A target payload to be written.
   * @param merger A function for merging payloads.
   * @return The payload before updating.
   */
  template <class Payload, class Delta>
  auto
  Update(  //
      int32_t pos,
      const Delta& payload,
      Payload (*merger)(const Payload&, const Delta&)) noexcept  //
      -> Payload
  {
    Payload out_pay{};
    auto* pay_addr = ShiftAddr(this, meta_arr_[pos].GetPayOff());
    std::memcpy(&out_pay, pay_addr, sizeof(Payload));
    if (merger) {
      const auto& merged = merger(out_pay, payload);
      std::memcpy(pay_addr, &merged, sizeof(Payload));
    } else {
      std::memcpy(pay_addr, &payload, sizeof(Payload));
    }
    return out_pay;
  }

  /**
   * @brief Delete a record using a given position.
   *
   * @param[in] pos A position where a target record is.
   * @param[out] out_pay An instance for storing an output payload.
   * @retval true if this node should be merged.
   * @retval false otherwise.
   */
  auto
  Delete(  //
      const int32_t pos,
      void* out_pay) noexcept  //
      -> bool
  {
    auto& meta = meta_arr_[pos];
    std::memcpy(out_pay, ShiftAddr(this, meta.GetPayOff()), meta.GetPayLen());
    if (level_ > 0) {
      assert(pos > 0);
      const auto rec_len = meta.rec_len;
      std::memmove(&meta, ShiftAddr(&meta, kMetaSize), kMetaSize * (rec_num_ - pos - 1UL));
      --rec_num_;
      usage_ -= rec_len + kMetaSize;
    } else {
      meta.deleted = true;
      usage_ -= pos > 0 ? meta.rec_len + kMetaSize : 0;
    }
    return kHeaderSize + usage_ < page_size_ / k8;
  }

  /*##########################################################################*
   * SMO APIs
   *##########################################################################*/

  /**
   * @brief Merge a given node into this node.
   *
   * @param grd The SIX-lock guard of this node.
   * @param r_node A right-merged (to be removed) node.
   * @param r_grd The SIX-lock guard of a right node.
   */
  void
  Merge(  //
      SIXGuard grd,
      Node* r_node,
      SIXGuard r_grd)
  {
    auto* tmp = AllocNode();
    tmp->offset_ = page_size_ - r_node->hk_len_;

    int32_t pos = 0;
    tmp->CopyHighKey(ShiftAddr(r_node, tmp->offset_), r_node->hk_len_);
    tmp->CopyLowRecord(this, pos);
    tmp->CopyRecords(this, pos, rec_num_);
    pos = 0;
    tmp->CopyRecords(r_node, pos, r_node->rec_num_);
    {
      auto&& x_grd = grd.UpgradeToX();
      CopyFrom(tmp);
      sib_node_ = r_node->sib_node_;
      VerIncrement<kSMOMask>(x_grd);
    }
    r_node->Remove(std::move(r_grd), this);

    FreeNode(tmp);
  }

  /**
   * @brief Set a remove flag to this node.
   *
   * @param grd The SIX-lock guard of this node.
   * @param next A side link for traversing to the next target node.
   */
  void
  Remove(  //
      SIXGuard grd,
      Node* next = nullptr)
  {
    auto&& x_grd = grd.UpgradeToX();
    removed_ = true;
    sib_node_ = next;
    VerIncrement<kSMOMask>(x_grd);
  }

  /*##########################################################################*
   * Public bulkload API
   *##########################################################################*/

  /**
   * @brief Bulkload records into this node.
   *
   * @tparam Entry A container for a key/payload pair.
   * @param iter The begin position of target records.
   * @param iter_end The end position of target records.
   */
  template <class BulkIter>
  void
  Bulkload(  //
      BulkIter& iter,
      const BulkIter& iter_end) noexcept
  {
    const uint32_t leaf_cap = page_size_ * 3 / 4;
    const uint32_t inner_cap = page_size_ - std::max(kMaxKeySize, page_size_ / 8U);

    offset_ = page_size_ - kMaxKeySize;
    const uint32_t max_usage = level_ > 0 ? inner_cap : leaf_cap;
    for (; iter < iter_end; ++iter) {
      const auto& [key, payload, key_len, pay_len] = ParseEntry(*iter);
      const uint32_t rec_len = key_len + pay_len;
      const uint32_t total_len = rec_len + kMetaSize;
      if (usage_ + total_len > max_usage) break;

      offset_ -= rec_len;
      std::memcpy(ShiftAddr(this, offset_), GetSrcAddr(key), key_len);
      std::memcpy(ShiftAddr(this, offset_ + key_len), &payload, pay_len);
      meta_arr_[rec_num_++] = Metadata{offset_, key_len, rec_len};
      usage_ += total_len;
    }
    leftmost_ = false;
  }

  /**
   * @brief Link this node with a right sibling node.
   *
   * @param sib_node A right sibling node.
   * @param key A highest key.
   * @param key_len The length of a highest key.
   */
  void
  LinkSiblingNode(  //
      Node* sib_node,
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
      Node* l_node,
      Node* r_node) noexcept
  {
    if (l_node == nullptr) return;

    while (true) {
      // link the border nodes
      const auto meta = r_node->meta_arr_[0];
      const auto key_len = meta.key_len;
      l_node->sib_node_ = r_node;
      l_node->hk_len_ = key_len;
      l_node->usage_ += key_len;
      auto* hk_addr = ShiftAddr(l_node, l_node->page_size_ - key_len);
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
      Node* node) noexcept
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
   * Public utility
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
      const auto* src_addr = ShiftAddr(this, page_size_ - hk_len_);
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

 private:
  /*##########################################################################*
   * Internal types
   *##########################################################################*/

  /// @brief A type for allocating variable-length data.
  using KeyWOPtr = std::remove_pointer_t<Key>;

  /*##########################################################################*
   * Internal constants
   *##########################################################################*/

  /// @brief A comparator instance.
  static constexpr Comp kLess{};

  /// @brief A flag for indicating whether the current SMO is splitting.
  static constexpr bool kSplit = true;

  /// @brief The maximum size of keys.
  static constexpr uint32_t kMaxKeySize = MaxSize<Key>();

  /*##########################################################################*
   * Internal APIs
   *##########################################################################*/

  /**
   * @return A temporary node
   */
  [[nodiscard]]
  auto
  AllocNode()  //
      -> Node*
  {
    auto* page = ::operator new(page_size_, static_cast<std::align_val_t>(align_val_));
    return new (page) Node{page_size_, align_val_};
  }

  /**
   * @param page A page for a temporary node.
   */
  void
  FreeNode(  //
      void* page)
  {
    ::operator delete(page, static_cast<std::align_val_t>(align_val_));
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
    const auto* src_addr = ShiftAddr(this, meta.offset);
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
   * @brief Clean up this node.
   *
   */
  void
  CleanUp()
  {
    auto* tmp = AllocNode();
    tmp->offset_ = page_size_ - hk_len_;

    int32_t pos = 0;
    tmp->CopyHighKey(ShiftAddr(this, tmp->offset_), hk_len_);
    tmp->CopyLowRecord(this, pos);
    tmp->CopyRecords(this, pos, rec_num_);
    CopyFrom(tmp);

    FreeNode(tmp);
  }

  /**
   * @brief Copy records and their relevant data from a TL temporary node.
   *
   */
  void
  CopyFrom(  //
      const Node* src) noexcept
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
      const void* src,
      const size_t key_len) noexcept
  {
    std::memcpy(ShiftAddr(this, page_size_ - key_len), src, key_len);
    usage_ += key_len;
    hk_len_ = key_len;
  }

  /**
   * @brief Copy a record that includes the lowest key in a node.
   *
   * @param src A source node page.
   */
  void
  CopyLowRecord(  //
      const Node* src,
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
      const Node* src,
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
  Node* sib_node_{};

  /// @brief An actual data block (it starts with record metadata).
  Metadata meta_arr_[0];
};

}  // namespace dbgroup::index::b_tree::component

#endif  // B_TREE_DBGROUP_B_TREE_COMPONENT_NODE_HPP_
