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
#include <bit>
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
#include "dbgroup/b_tree/component/metadata.hpp"
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
  using Header::offset_;
  using Header::page_size_;
  using Header::rec_num_;
  using Header::sib_node_;
  using Header::usage_;
  using typename Header::KeyWOPtr;
  using typename Header::SIXGuard;
  using typename Header::XGuard;

 public:
  /*##########################################################################*
   * Public constructors and assignment operators
   *##########################################################################*/

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
      pos = rec_num - 1;
    }
    return this->LinearSearchImpl(key, pos, found, 0);
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
   * @brief Copy the payload of a given position.
   *
   * @param[in] pos A target record position.
   * @param[out] out_pay A destination buffer of at least `max_len` bytes.
   * @param[in] max_len The maximum length of payloads.
   * @return The length of the copied payload.
   * @note This function may read corrupted metadata due to OCC, so it checks
   * the length before copying.
   */
  auto
  ReadPayload(  //
      const int32_t pos,
      void* out_pay,
      const uint32_t max_len) const noexcept  //
      -> uint32_t
  {
    const auto meta = meta_arr_[pos];
    const auto pay_len = meta.GetPayLen();
    const auto offset = meta.GetPayOff();
    if (pay_len > max_len || offset + pay_len > page_size_) [[unlikely]] {
      return 0;  // read corrupted data due to OCC
    }
    std::memcpy(out_pay, ShiftAddr(this, offset), pay_len);
    return pay_len;
  }

  /**
   * @brief Copy the payload of a given position if requested.
   *
   * @tparam Buf A payload buffer class (i.e., `BTreeVarLen::PayloadBuf`).
   * @param[in] pos A target record position.
   * @param[out] out A buffer for storing a payload (ignored if nullptr).
   */
  template <class Buf>
  void
  CopyOldPayload(  //
      const int32_t pos,
      Buf* out) const noexcept
  {
    if (!out) return;
    out->len = ReadPayload(pos, &(out->data[0]), sizeof(out->data));
  }

  /**
   * @brief Copy a record if exist.
   *
   * @param[in] pos A position to be copied.
   * @param[out] key A destination address for a key.
   * @param[out] key_len The length of a copied key.
   * @param[out] payload A destination address for a payload.
   * @param[out] pay_len The length of a copied payload.
   * @retval true if a target record is active (not deleted).
   * @retval false otherwise.
   */
  [[nodiscard]]
  auto
  CopyRecordTo(  //
      const int32_t pos,
      void* key,
      uint32_t& key_len,
      void* payload,
      uint32_t& pay_len) const noexcept  //
      -> bool
  {
    const auto meta = meta_arr_[pos];
    if (meta.deleted) return false;

    key_len = meta.key_len;
    pay_len = meta.GetPayLen();
    std::memcpy(key, ShiftAddr(this, meta.offset), key_len);
    std::memcpy(payload, ShiftAddr(this, meta.GetPayOff()), pay_len);
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
   * @note This node is not modified if this function returns false.
   */
  void
  Insert(  //
      int32_t pos,
      bool found,
      const Key& key,
      const size_t key_len,
      const void* payload,
      const size_t pay_len,
      SIXGuard& six_grd) noexcept
  {
    if (found && ReuseRecord(pos, key, key_len, payload, pay_len, six_grd)) return;

    const auto total_len = key_len + pay_len + kMetaSize;
    if (!this->HasSpace(total_len)) return;
    if (this->CleanUpIfNeeded(total_len, six_grd)) {
      this->LinearSearch(key, pos, found);  // deleted records have been removed
      if (found) {
        // the lowest record is retained even if deleted
        ReuseRecord(pos, key, key_len, payload, pay_len, six_grd);
        assert(!six_grd);
        return;
      }
    }
    assert(!found);  // reuse failure implies no contiguous space, i.e., clean-up
    this->InsertRecord(pos, key, key_len, payload, pay_len, six_grd);
  }

  /**
   * @brief Update the payload of a given position.
   *
   * @tparam Payload A payload class.
   * @tparam T A class of a given payload or delta.
   * @tparam Buf A payload buffer class (i.e., `BTreeVarLen::PayloadBuf`).
   * @param[in] pos A position where a target record is.
   * @param[in] key A target key.
   * @param[in] key_len The length of a target key.
   * @param[in] payload A target payload or delta.
   * @param[in] pay_len The length of a target payload or delta.
   * @param[out] old A buffer for storing a previous payload (ignored if nullptr).
   * @param[in] merger A function for merging payloads.
   * @retval true if the payload is updated.
   * @retval false if this node should be split before updating.
   * @note This node is not modified if this function returns false.
   */
  template <class Payload, class T, class Buf>
  void
  Update(  //
      const int32_t pos,
      [[maybe_unused]] const Key& key,
      [[maybe_unused]] const size_t key_len,
      const T& payload,
      const size_t pay_len,
      Buf* old,
      Merger<Payload, T> merger,
      SIXGuard& six_grd) noexcept
  {
    constexpr uint32_t kMaxPaySize = MaxPaySize<Payload>();

    const void* src = GetSrcAddr(payload);
    size_t len = pay_len;
    if (merger) {
      alignas(alignof(WOPtr<Payload>)) thread_local std::byte tls_old[kMaxPaySize];
      alignas(alignof(WOPtr<Payload>)) thread_local std::byte tls_new[kMaxPaySize];
      void* cur;
      uint32_t cur_len;
      uint32_t* len_ptr;
      if (old) {
        cur = &(old->data[0]);
        len_ptr = &old->len;
      } else {
        cur = &tls_old[0];
        len_ptr = &cur_len;
      }
      *len_ptr = ReadPayload(pos, cur, kMaxPaySize);
      len = merger(std::bit_cast<const WOPtr<Payload>*>(cur), *len_ptr,
                   static_cast<const WOPtr<T>*>(src), pay_len,
                   std::bit_cast<WOPtr<Payload>*>(&tls_new[0]));
      assert(len <= kMaxPaySize);
      src = &tls_new[0];
    } else {
      assert((std::is_same_v<T, Payload>));  // a delta requires a merger
      CopyOldPayload(pos, old);
    }

    if constexpr (IsVarLenData<Payload>()) {
      Replace(pos, key, key_len, src, len, six_grd);
    } else {
      assert(len == sizeof(Payload));
      auto&& x_grd = six_grd.UpgradeToX();
      std::memcpy(ShiftAddr(this, meta_arr_[pos].GetPayOff()), src, len);
    }
  }

  /**
   * @brief Delete a record logically using a given position.
   *
   * @param pos A position where a target record is.
   * @retval true if this node should be merged.
   * @retval false otherwise.
   */
  void
  Delete(  //
      const int32_t pos,
      SIXGuard& six_grd) noexcept
  {
    auto& meta = meta_arr_[pos];
    auto&& x_grd = six_grd.UpgradeToX();
    meta.deleted = true;
    usage_ -= pos > 0 ? meta.rec_len + kMetaSize : 0;
    VerIncrement<kInsDelMask>(x_grd);
    if (this->NeedMerge()) {
      six_grd = x_grd.DowngradeToSIX();
    }
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

 private:
  /*##########################################################################*
   * Internal utilities
   *##########################################################################*/

  /**
   * @brief Replace a payload with another one that may have a different length.
   *
   * @param pos A position where a target record is.
   * @param key A target key.
   * @param key_len The length of a target key.
   * @param payload A target payload to be written.
   * @param pay_len The length of a target payload.
   * @retval true if the payload is replaced.
   * @retval false if this node should be split before replacing.
   * @note This node is not modified if this function returns false.
   */
  void
  Replace(  //
      int32_t pos,
      const Key& key,
      const size_t key_len,
      const void* payload,
      const size_t pay_len,
      SIXGuard& six_grd) noexcept
  {
    if (ReuseRecord(pos, key, key_len, payload, pay_len, six_grd)) return;
    if (!this->HasSpace(key_len + pay_len + kMetaSize)) return;

    bool found;
    this->CleanUp(six_grd);
    this->LinearSearch(key, pos, found);  // deleted records have been removed
    if (found) {
      // the lowest record is retained even if deleted
      ReuseRecord(pos, key, key_len, payload, pay_len, six_grd);
    } else {
      this->InsertRecord(pos, key, key_len, payload, pay_len, six_grd);
    }
    assert(!six_grd);
  }

  /**
   * @brief Reuse a deleted record for writing a given payload.
   *
   * @param pos A position where a deleted record is.
   * @param key A target key.
   * @param key_len The length of a target key.
   * @param payload A target payload to be written.
   * @param pay_len The length of a target payload.
   * @retval true if a payload is written.
   * @retval false if this node does not have enough contiguous space.
   */
  auto
  ReuseRecord(  //
      const int32_t pos,
      const Key& key,
      const size_t key_len,
      const void* payload,
      const size_t pay_len,
      SIXGuard& six_grd) noexcept  //
      -> bool
  {
    auto& meta = meta_arr_[pos];
    const bool deleted = meta.deleted;
    const uint32_t old_len = meta.rec_len;
    const uint32_t rec_len = key_len + pay_len;
    XGuard x_grd;
    if (pay_len <= meta.GetPayLen()) {  // overwrite in place
      x_grd = six_grd.UpgradeToX();
      std::memcpy(ShiftAddr(this, meta.GetPayOff()), payload, pay_len);
      meta = Metadata{meta.offset, key_len, rec_len};
    } else if (kHeaderSize + kMetaSize * rec_num_ + rec_len <= offset_) {
      x_grd = six_grd.UpgradeToX();
      offset_ -= rec_len;  // move the record to the free space
      std::memcpy(ShiftAddr(this, offset_), GetSrcAddr(key), key_len);
      std::memcpy(ShiftAddr(this, offset_ + key_len), payload, pay_len);
      meta = Metadata{offset_, key_len, rec_len};
    } else {
      return false;
    }

    usage_ += rec_len;
    if (deleted && pos > 0) {
      usage_ += kMetaSize;
    } else {
      usage_ -= old_len;
    }
    VerIncrement<kInsDelMask>(x_grd);
    return true;
  }
};

}  // namespace dbgroup::index::b_tree::component

#endif  // B_TREE_DBGROUP_B_TREE_COMPONENT_LEAF_NODE_HPP_
