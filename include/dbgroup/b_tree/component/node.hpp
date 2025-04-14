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
#include "dbgroup/constants.hpp"
#include "dbgroup/index/utility.hpp"
#include "dbgroup/lock/utility.hpp"
#include "dbgroup/memory/utility.hpp"

// local sources
#include "dbgroup/b_tree/component/common.hpp"
#include "dbgroup/b_tree/component/metadata.hpp"
#include "dbgroup/b_tree/utility.hpp"

// we intentionally use a zero-length array for record metadata
#pragma GCC diagnostic ignored "-Warray-bounds"

namespace dbgroup::index::b_tree::component
{
/**
 * @brief A class for representing nodes in B+trees.
 *
 * @tparam Key A target key class.
 * @tparam Comp A comparator class for keys.
 * @tparam Lock A lock class for concurrency controls.
 * @tparam kUseOptimisticCC A flag for using optimistic concurrency controls.
 */
template <class Key,  //
          class Comp,
          ::dbgroup::lock::Lockable Lock,
          bool kUseOptimisticCC>
class Node
{
  /*##########################################################################*
   * Type aliases
   *##########################################################################*/

  using ScanKey = std::optional<std::tuple<const Key &, size_t, bool>>;

 public:
  /*##########################################################################*
   * Public types
   *##########################################################################*/

  using OptGuard = OptGuardExtractor<Lock, kUseOptimisticCC>::type;
  using SGuard = SGuardExtractor<Lock, !kUseOptimisticCC>::type;
  using SIXGuard = SIXGuardExtractor<Lock, !kUseOptimisticCC>::type;
  using XGuard = typename Lock::XGuard;
  using ReadGuard = std::conditional_t<kUseOptimisticCC, OptGuard, SGuard>;
  using CheckGuard = std::conditional_t<kUseOptimisticCC, OptGuard, SIXGuard>;

  /*##########################################################################*
   * Public constants
   *##########################################################################*/

  /// @brief A flag for indicating this node uses optimistic CC.
  static constexpr bool kUseOptCC = kUseOptimisticCC;

  /*##########################################################################*
   * Public constructors and assignment operators
   *##########################################################################*/

  /**
   * @brief Construct an empty node.
   *
   * @param level The level where this node is.
   */
  constexpr explicit Node(  //
      const size_t level = 0)
      : level_{static_cast<uint8_t>(level)}
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
      const size_t level,
      const Key &key,
      const size_t key_len,
      const void *l_child,
      const void *r_child)  //
      : rec_num_{2},
        offset_{static_cast<uint16_t>(kPageSize - (2 * kPtrSize + key_len))},
        usage_{static_cast<uint16_t>(kHeaderSize + 2 * (kMetaSize + kPtrSize) + key_len)},
        level_{static_cast<uint8_t>(level)}
  {
    constexpr size_t kLChildOffset = kPageSize - kPtrSize;
    constexpr size_t kRChildOffset = kPageSize - 2 * kPtrSize;

    std::memcpy(ShiftAddr(this, kLChildOffset), &l_child, kPtrSize);
    meta_arr_[0] = Metadata{kLChildOffset, 0, kPtrSize};

    const auto offset = kRChildOffset - key_len;
    std::memcpy(ShiftAddr(this, kRChildOffset), &r_child, kPtrSize);
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
      Node *l_node)
      : offset_{static_cast<uint16_t>(kPageSize - l_node->hk_len_)},
        level_{l_node->level_},
        leftmost_{false},
        sib_node_{l_node->sib_node_}
  {
    auto *tmp = new (&_tls_page) Node{};
    tmp->offset_ = kPageSize - MaxSize<Key>();

    // copy split-left records to a temporary region
    const double sep_ratio = l_node->leftmost_ && sib_node_     ? 0.25  // NOLINTBEGIN
                             : !l_node->leftmost_ && !sib_node_ ? 0.75
                                                                : 0.5;  // NOLINTEND
    const size_t sep_size = (l_node->usage_ - l_node->hk_len_) * sep_ratio;
    size_t pos = 0;
    tmp->template CopyRecords<kSplit>(l_node, pos, l_node->rec_num_, sep_size);
    const auto sep_meta = l_node->meta_arr_[pos];
    tmp->CopyHighKey(ShiftAddr(l_node, sep_meta.offset), sep_meta.key_len);

    // copy split-right records to this node
    CopyHighKey(ShiftAddr(l_node, offset_), l_node->hk_len_);
    CopyRecords(l_node, pos, l_node->rec_num_);

    // copy the split-left records to the left node
    l_node->CopyFromTmpNode();
    l_node->sib_node_ = this;
  }

  Node(const Node &) = delete;
  Node(Node &&) = delete;

  auto operator=(const Node &) -> Node & = delete;
  auto operator=(Node &&) -> Node & = delete;

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
  [[nodiscard]] constexpr auto
  GetLevel() const  //
      -> size_t
  {
    return level_;
  }

  /**
   * @return The number of records in this node.
   */
  [[nodiscard]] constexpr auto
  GetRecNum() const  //
      -> size_t
  {
    return rec_num_;
  }

  /**
   * @retval true if this node has been removed.
   * @retval false otherwise.
   */
  [[nodiscard]] constexpr auto
  Removed() const  //
      -> bool
  {
    return removed_;
  }

  /**
   * @return A sibling node.
   */
  [[nodiscard]] constexpr auto
  GetSibNode() const  //
      -> Node *
  {
    return sib_node_;
  }

  /**
   * @param pos A target record position.
   * @return A child node.
   */
  [[nodiscard]] auto
  GetChild(                    //
      const size_t pos) const  //
      -> Node *
  {
    Node *child;
    const auto meta = meta_arr_[pos];
    std::memcpy(&child, ShiftAddr(this, meta.offset + meta.key_len), kPtrSize);
    return child;
  }

  /**
   * @retval 1st: A separator key.
   * @retval 2nd: The length of a separator key.
   */
  [[nodiscard]] auto
  GetSeparatorKey() const  //
      -> std::pair<Key, size_t>
  {
    const auto *src_addr = ShiftAddr(this, kPageSize - hk_len_);
    if constexpr (IsVarLenData<Key>()) {
      auto *key = ::dbgroup::memory::Allocate<KeyWOPtr>(hk_len_);
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
   * @retval true if this node can include a given key.
   * @retval false otherwise.
   */
  [[nodiscard]] auto
  Include(                   //
      const Key &key) const  //
      -> bool
  {
    if (hk_len_ == 0) return true;

    Key h_key;
    const auto *src_addr = ShiftAddr(this, kPageSize - hk_len_);
    if constexpr (IsVarLenData<Key>()) {
      alignas(alignof(KeyWOPtr)) thread_local std::byte tls_key[kMaxVarDataSize];
      h_key = std::bit_cast<Key>(&tls_key);
      std::memcpy(h_key, src_addr, hk_len_);
    } else {
      std::memcpy(&h_key, src_addr, sizeof(Key));
    }
    return kLess(key, h_key);
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
  [[nodiscard]] auto
  CheckUniqueness(           //
      const Key &key) const  //
      -> std::tuple<bool, bool, size_t>
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
  [[nodiscard]] auto
  SearchChild(               //
      const Key &key) const  //
      -> Node *
  {
    auto [found, pos] = SearchRecord(key);
    pos -= static_cast<int64_t>(!found);

    Node *child;
    const auto meta = meta_arr_[pos];
    std::memcpy(&child, ShiftAddr(this, meta.offset + meta.key_len), kPtrSize);
    return child;
  }

  /**
   * @tparam Guard A desired guard type.
   * @return A guard instance for this node.
   */
  template <class Guard>
  [[nodiscard]] auto
  GetGuard()
  {
    if constexpr (kUseOptimisticCC) {
      if constexpr (std::is_same_v<Guard, typename Lock::OptGuard>) {
        return lock_.GetVersion();
      } else {
        return lock_.LockX();
      }
    } else {
      if constexpr (std::is_same_v<Guard, typename Lock::SGuard>) {
        return lock_.LockS();
      } else if constexpr (std::is_same_v<Guard, typename Lock::SIXGuard>) {
        return lock_.LockSIX();
      } else {
        return lock_.LockX();
      }
    }
  }

  /*##########################################################################*
   * Read APIs
   *##########################################################################*/

  /**
   * @brief Read the payload associated with a given key.
   *
   * @tparam Payload A class of stored payloads.
   * @param key A search key.
   * @retval A found payload if exist.
   * @retval std::nullopt otherwise.
   */
  template <class Payload>
  [[nodiscard]] auto
  Read(                      //
      const Key &key) const  //
      -> std::optional<Payload>
  {
    const auto [found, pos] = SearchRecord(key);
    const auto meta = meta_arr_[pos];
    if (!found || meta.deleted) return std::nullopt;

    Payload payload;
    const auto *addr = ShiftAddr(this, meta.GetPayOff());
    if constexpr (IsVarLenData<Payload>()) {
      using PayWOPtr = std::remove_pointer_t<Payload>;
      alignas(alignof(PayWOPtr)) thread_local std::byte tls_payload[kMaxVarDataSize];
      payload = std::bit_cast<Payload>(&tls_payload);
      std::memcpy(payload, addr, meta.GetPayLen());
    } else {
      std::memcpy(&payload, addr, sizeof(Payload));
    }
    return payload;
  }

  /*##########################################################################*
   * Write APIs
   *##########################################################################*/

  /**
   * @brief Write a given kay/payload pair.
   *
   * @param x_guard A lock guard to set version.
   * @param key A target key to be written.
   * @param key_len The length of a target key.
   * @param payload A target payload to be written.
   * @param pay_len The length of a target payload.
   * @retval kCompleted if a record is written.
   * @retval kNeedSplit if this node should be split before inserting a record.
   */
  auto
  Write(  //
      Lock::XGuard &x_guard,
      const Key &key,
      const size_t key_len,
      const void *payload,
      const size_t pay_len)  //
      -> NodeRC
  {
    const auto rec_len = key_len + pay_len;
    const auto total_len = rec_len + kMetaSize;
    if (usage_ + total_len > kPageSize) return kNeedSplit;
    if (kHeaderSize + kMetaSize * rec_num_ + total_len > offset_) {
      CleanUp();
      if constexpr (kUseOptimisticCC) {
        VerIncrement<kInsDelMask>(x_guard);
      }
    }

    const auto [found, pos] = SearchRecord(key);
    auto &meta = meta_arr_[pos];
    auto offset = offset_ - rec_len;
    if (found) {  // try to update a record
      if (meta.deleted) {
        usage_ += total_len;
        if constexpr (kUseOptimisticCC) {
          VerIncrement<kInsDelMask>(x_guard);
        }
      } else {
        usage_ += static_cast<int64_t>(rec_len - meta.rec_len);
      }
      if (rec_len <= meta.rec_len) {  // reuse a record
        offset = meta.offset;
      } else {  // insert as a new record
        std::memcpy(ShiftAddr(this, offset), ShiftAddr(this, meta.offset), key_len);
        offset_ -= rec_len;
      }
      std::memcpy(ShiftAddr(this, offset + key_len), payload, pay_len);
    } else {
      // insert a new record
      std::memcpy(ShiftAddr(this, offset), GetSrcAddr(key), key_len);
      std::memcpy(ShiftAddr(this, offset + key_len), payload, pay_len);
      std::memmove(ShiftAddr(&meta, kMetaSize), &meta, kMetaSize * (rec_num_ - pos));
      meta = Metadata{offset, key_len, rec_len};
      ++rec_num_;
      offset_ -= rec_len;
      usage_ += total_len;
      if constexpr (kUseOptimisticCC) {
        VerIncrement<kInsDelMask>(x_guard);
      }
    }
    meta = Metadata{offset, key_len, rec_len};

    return kCompleted;
  }

  /**
   * @brief Insert a given kay/payload pair.
   *
   * @param x_guard A lock guard to set version.
   * @param pos A position where a new record is inserted.
   * @param found A flag for indicating the same (deleted) key is found.
   * @param key A target key to be written.
   * @param key_len The length of a target key.
   * @param payload A target payload to be written.
   * @param pay_len The length of a target payload.
   * @retval kCompleted if a record is written.
   * @retval kNeedSplit if this node should be split before inserting a record.
   */
  auto
  Insert(  //
      Lock::XGuard &x_guard,
      size_t pos,
      bool found,
      const Key &key,
      const size_t key_len,
      const void *payload,
      const size_t pay_len)  //
      -> NodeRC
  {
    const auto rec_len = key_len + pay_len;
    const auto total_len = rec_len + kMetaSize;
    if (usage_ + total_len > kPageSize) return kNeedSplit;
    if (kHeaderSize + kMetaSize * rec_num_ + total_len > offset_) {
      CleanUp();
      std::tie(found, pos) = SearchRecord(key);
    }

    // insert a new record
    auto &meta = meta_arr_[pos];
    auto offset = offset_ - rec_len;
    usage_ += total_len;
    if (found) {                      // found a deleted record
      if (rec_len <= meta.rec_len) {  // reuse a record
        offset = meta.offset;
      } else {  // insert as a new record
        std::memcpy(ShiftAddr(this, offset), ShiftAddr(this, meta.offset), key_len);
        offset_ -= rec_len;
      }
      std::memcpy(ShiftAddr(this, offset + key_len), payload, pay_len);
    } else {
      std::memcpy(ShiftAddr(this, offset), GetSrcAddr(key), key_len);
      std::memcpy(ShiftAddr(this, offset + key_len), payload, pay_len);
      std::memmove(ShiftAddr(&meta, kMetaSize), &meta, kMetaSize * (rec_num_ - pos));
      ++rec_num_;
      offset_ -= rec_len;
    }
    meta = Metadata{offset, key_len, rec_len};

    if constexpr (kUseOptimisticCC) {
      VerIncrement<kInsDelMask>(x_guard);
    }
    return kCompleted;
  }

  /**
   * @brief Update a record using a given kay/payload pair.
   *
   * @param x_guard A lock guard to set version.
   * @param pos A position where a target record is.
   * @param key A target key to be written.
   * @param key_len The length of a target key.
   * @param payload A target payload to be written.
   * @param pay_len The length of a target payload.
   * @retval kCompleted if a record is written.
   * @retval kNeedSplit if this node should be split before inserting a record.
   */
  auto
  Update(  //
      Lock::XGuard &x_guard,
      size_t pos,
      const Key &key,
      const size_t key_len,
      const void *payload,
      const size_t pay_len)  //
      -> NodeRC
  {
    const auto rec_len = key_len + pay_len;
    const auto total_len = rec_len + kMetaSize;
    if (usage_ + total_len > kPageSize) return kNeedSplit;
    if (kHeaderSize + kMetaSize * rec_num_ + total_len > offset_) {
      CleanUp();
      pos = SearchRecord(key).second;
      if constexpr (kUseOptimisticCC) {
        VerIncrement<kInsDelMask>(x_guard);
      }
    }

    auto &meta = meta_arr_[pos];
    auto offset = offset_ - rec_len;
    if (rec_len <= meta.rec_len) {  // reuse a record
      offset = meta.offset;
    } else {  // insert as a new record
      std::memcpy(ShiftAddr(this, offset), ShiftAddr(this, meta.offset), key_len);
      offset_ -= rec_len;
    }
    std::memcpy(ShiftAddr(this, offset + key_len), payload, pay_len);
    usage_ += static_cast<int64_t>(rec_len - meta.rec_len);
    meta = Metadata{offset, key_len, rec_len};
    return kCompleted;
  }

  /*##########################################################################*
   * Public utility
   *##########################################################################*/

  /**
   * @return The explanation of this instance.
   */
  [[nodiscard]] auto
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
      const auto *src_addr = ShiftAddr(this, kPageSize - hk_len_);
      if constexpr (IsVarLenData<Key>()) {
        alignas(alignof(KeyWOPtr)) thread_local std::byte tls_key[kMaxVarDataSize];
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

  /*##########################################################################*
   * Internal getter for header information
   *##########################################################################*/

  /**
   * @param pos The position of a target record.
   * @return A key in a target record.
   */
  [[nodiscard]] auto
  GetKey(                      //
      const size_t pos) const  //
      -> Key
  {
    Key key;
    const auto meta = meta_arr_[pos];
    const auto *src_addr = ShiftAddr(this, meta.offset);
    if constexpr (IsVarLenData<Key>()) {
      alignas(alignof(KeyWOPtr)) thread_local std::byte tls_key[kMaxVarDataSize];
      key = std::bit_cast<Key>(&tls_key);
      std::memcpy(key, src_addr, meta.key_len);
    } else {
      std::memcpy(&key, src_addr, sizeof(Key));
    }
    return key;
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
  [[nodiscard]] auto
  SearchRecord(              //
      const Key &key) const  //
      -> std::pair<bool, size_t>
  {
    auto found = false;
    auto pos = static_cast<int64_t>(level_ > 0);
    int64_t end_pos = rec_num_ - 1;
    while (pos <= end_pos) {
      const int64_t cur_pos = (pos + end_pos) / 2;
      const auto &index_key = GetKey(cur_pos);
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
   * @brief Clean up this node.
   *
   */
  void
  CleanUp()
  {
    auto *tmp = new (&_tls_page) Node{};
    tmp->offset_ = kPageSize - hk_len_;

    size_t pos = 0;
    tmp->CopyHighKey(ShiftAddr(this, tmp->offset_), hk_len_);
    tmp->CopyRecords(this, pos, rec_num_);
    CopyFromTmpNode();
  }

  /**
   * @brief Copy records and their relevant data from a TL temporary node.
   *
   */
  void
  CopyFromTmpNode()
  {
    const auto *tmp = std::bit_cast<Node *>(&_tls_page);
    rec_num_ = tmp->rec_num_;
    offset_ = tmp->offset_;
    usage_ = tmp->usage_;
    hk_len_ = tmp->hk_len_;
    std::memcpy(meta_arr_, tmp->meta_arr_, kMetaSize * rec_num_);
    std::memcpy(ShiftAddr(this, offset_), ShiftAddr(tmp, offset_), kPageSize - offset_);
  }

  /**
   * @brief Copy a highest key from a given address.
   *
   * @param[in] src The address of a source key.
   * @param[in] key_len The length of a key.
   */
  void
  CopyHighKey(  //
      const void *src,
      const size_t key_len)
  {
    std::memcpy(ShiftAddr(this, kPageSize - key_len), src, key_len);
    usage_ += key_len;
    hk_len_ = key_len;
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
      const Node *src,
      size_t &pos,
      const size_t end_pos,
      [[maybe_unused]] const size_t sep_size = kPageSize)
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

  static_assert(  //
      !kUseOptimisticCC || ::dbgroup::lock::OptimisticallyLockable<Lock>,
      "A lock class must have optimistic lock APIs when using Optimistic CC.");

  /*##########################################################################*
   * Internal member variables
   *##########################################################################*/

  /// @brief The number of records in this node.
  uint16_t rec_num_{};

  /// @brief The total byte length of a record block.
  uint16_t offset_{};

  /// @brief The total usage of this node.
  uint16_t usage_{kHeaderSize};

  /// @brief The length of a highest key.
  uint16_t hk_len_{};

  /// @brief A level where this node is.
  uint8_t level_{};

  /// @brief A flag for indicating this node is removed from a tree.
  bool removed_{};

  /// @brief A flag for indicating the leftmost node.
  bool leftmost_{true};

  /// @brief A blank block.
  std::byte blank_[5]{};  // NOLINT

  /// @brief A lock for concurrency controls.
  Lock lock_{};

  /// @brief A right sibling node.
  Node *sib_node_{};

  /// @brief An actual data block (it starts with record metadata).
  Metadata meta_arr_[0];

  /*##########################################################################*
   * Thread local class variables
   *##########################################################################*/
  // NOLINTBEGIN

  /// @brief A temporary node page for SMOs.
  static thread_local inline Page _tls_page{};

  // NOLINTEND
};

}  // namespace dbgroup::index::b_tree::component

#endif  // B_TREE_DBGROUP_B_TREE_COMPONENT_NODE_HPP_
