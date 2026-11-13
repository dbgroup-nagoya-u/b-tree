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

#ifndef B_TREE_DBGROUP_B_TREE_COMPONENT_LEAF_NODE_HPP_
#define B_TREE_DBGROUP_B_TREE_COMPONENT_LEAF_NODE_HPP_

// C++ standard libraries
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <optional>
#include <tuple>
#include <utility>

// external libraries
#include <dbgroup/index/utility.hpp>
#include <dbgroup/lock/utility.hpp>
#include <dbgroup/memory/utility.hpp>

// local sources
#include "dbgroup/b_tree/component/common.hpp"
#include "dbgroup/b_tree/component/node_header.hpp"

namespace dbgroup::index::b_tree::component
{
/**
 * @brief A class for representing leaf nodes in B+trees.
 *
 * @tparam Key A target key class.
 * @tparam Comp A comparator class for keys.
 * @tparam Lock A lock class for concurrency controls.
 */
template <class Key,  //
          class Comp,
          ::dbgroup::lock::OptimisticallyLockable Lock>
class LeafNode : public NodeHeader<Key, Comp, Lock>
{
  /*##########################################################################*
   * Type aliases
   *##########################################################################*/

  using Header = NodeHeader<Key, Comp, Lock>;
  using ScanKey = std::optional<std::tuple<Key, size_t, bool>>;

  // header members used in this class
  using Header::GetKey;
  using Header::kLess;
  using Header::leftmost_;
  using Header::meta_arr_;
  using Header::page_size_;
  using Header::rec_num_;
  using Header::sib_node_;
  using Header::usage_;
  using typename Header::KeyWOPtr;

 public:
  /*##########################################################################*
   * Public constructors and assignment operators
   *##########################################################################*/

  /**
   * @brief Construct an empty node.
   *
   * @param page_size The page size for each node.
   * @param align_val The alignment size for each node.
   */
  constexpr LeafNode(  //
      const uint32_t page_size,
      const uint32_t align_val) noexcept
      : Header{page_size, align_val, 0}
  {
    static_assert(sizeof(LeafNode) == sizeof(Header));
  }

  /**
   * @brief Construct split nodes.
   *
   * @param[in,out] l_node A left (source) node.
   * @note This is a split right node.
   */
  explicit LeafNode(  //
      LeafNode* l_node)
      : Header{l_node}
  {
  }

  LeafNode(const LeafNode&) = delete;
  LeafNode(LeafNode&&) = delete;

  auto operator=(const LeafNode&) -> LeafNode& = delete;
  auto operator=(LeafNode&&) -> LeafNode& = delete;

  /*##########################################################################*
   * Public destructors
   *##########################################################################*/

  /**
   * @brief Destroy the node object.
   *
   */
  ~LeafNode() = default;

  /*##########################################################################*
   * Public getters
   *##########################################################################*/

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
    return this->SearchRecordImpl(key, 0);
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
    this->CopyPayloadTo(pos, &out_pay, sizeof(Payload));
    return out_pay;
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
      is_end = !sib_node_ || this->Include(key);
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
      const bool found,
      const Key& key,
      const size_t key_len,
      const void* payload,
      const size_t pay_len) noexcept  //
      -> bool
  {
    const auto total_len = key_len + pay_len + kMetaSize;
    if (found) {  // reuse the deleted record
      auto& meta = meta_arr_[pos];
      std::memcpy(ShiftAddr(this, meta.GetPayOff()), payload, pay_len);
      meta.deleted = false;
      usage_ += pos > 0 ? total_len : 0;
      return true;
    }

    if (!this->HasSpace(total_len)) return false;
    if (this->CleanUpIfNeeded(total_len)) {
      pos = SearchRecord(key).second;  // deleted records have been removed
    }
    this->InsertRecord(pos, key, key_len, payload, pay_len);
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
      const int32_t pos,
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
   * @brief Delete a record logically using a given position.
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
    meta.deleted = true;
    usage_ -= pos > 0 ? meta.rec_len + kMetaSize : 0;
    return this->NeedMerge();
  }

  /*##########################################################################*
   * Public bulkload APIs
   *##########################################################################*/

  /**
   * @brief Bulkload records into this node.
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
    this->BulkloadImpl(iter, iter_end, page_size_ * 3 / 4);
  }
};

}  // namespace dbgroup::index::b_tree::component

#endif  // B_TREE_DBGROUP_B_TREE_COMPONENT_LEAF_NODE_HPP_
