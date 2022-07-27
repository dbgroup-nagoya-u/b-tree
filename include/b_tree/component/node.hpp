
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

#ifndef B_TREE_COMPONENT_NODE_HPP
#define B_TREE_COMPONENT_NODE_HPP

#include <atomic>
#include <optional>
#include <utility>
#include <vector>

#include "common.hpp"
#include "lock/pessimistic_lock.hpp"
#include "metadata.hpp"

namespace dbgroup::index::b_tree::component
{

/**
 * @brief A class to represent nodes in Btree.
 *
 * @tparam Key a target key class.
 * @tparam Comp a comparetor class for keys.
 */
template <class Key, class Comp>
class PessimisticNode
{
 public:
  /*####################################################################################
   * Type aliases
   *##################################################################################*/

  using Lock = ::dbgroup::lock::PessimisticLock;

  /*####################################################################################
   * Public constructors and assignment operators
   *##################################################################################*/

  /**
   * @brief Construct a new node object
   *
   * @param is_leaf a flag to indicate whether a leaf node is constructed
   */
  explicit PessimisticNode(const bool is_leaf)
      : is_leaf_{static_cast<uint64_t>(is_leaf)}, record_count_{0}, block_size_{0}, deleted_size_{0}
  {
  }

  /**
   * @brief Construct a new node object when root split
   *
   * @param l_node a left child node which is previous root
   * @param r_node a right child node which is split into
   */
  PessimisticNode(  //
      PessimisticNode *l_node,
      const PessimisticNode *r_node)  //
      : is_leaf_{0}, record_count_{2}, block_size_{0}, deleted_size_{0}
  {
    // insert l_node
    const auto l_high_meta = l_node->high_meta_;
    const auto l_key_len = l_high_meta.key_len;
    auto offset = SetPayload(kPageSize, l_node);
    offset = CopyKeyFrom(l_node, l_high_meta, offset);
    meta_array_[0] = Metadata{offset, l_key_len, l_key_len + kPtrLen};

    // insert r_node
    offset = SetPayload(offset, r_node);
    meta_array_[1] = Metadata{offset, 0, kPtrLen};

    // update header information
    block_size_ += kPageSize - offset;
  }

  PessimisticNode(const PessimisticNode &) = delete;
  PessimisticNode(PessimisticNode &&) = delete;

  auto operator=(const PessimisticNode &) -> PessimisticNode & = delete;
  auto operator=(PessimisticNode &&) -> PessimisticNode & = delete;

  /*####################################################################################
   * Public destructors
   *##################################################################################*/

  /**
   * @brief Destroy the node object.
   *
   */
  ~PessimisticNode() = default;

  /*####################################################################################
   * new/delete operators
   *##################################################################################*/

  static auto
  operator new([[maybe_unused]] std::size_t n)  //
      -> void *
  {
    return ::operator new(kPageSize);
  }

  static void
  operator delete(void *p) noexcept
  {
    ::operator delete(p);
  }

  /*####################################################################################
   * Public getters for header information
   *##################################################################################*/

  /**
   * @return true if this is a leaf node
   * @return false otherwise
   */
  [[nodiscard]] constexpr auto
  IsLeaf() const  //
      -> bool
  {
    return is_leaf_;
  }

  /**
   * @param new_rec_len the length of a new record.
   * @retval true if this node requires splitting before inserting a new record.
   * @retval false otherwise.
   */
  [[nodiscard]] auto
  NeedSplit(const size_t new_rec_len)  //
      -> bool
  {
    // check whether the node has space for a new record
    const auto total_size = kMetaLen * (record_count_ + 1) + block_size_ + new_rec_len;
    if (total_size <= kPageSize - kHeaderLen) return false;
    if (total_size - deleted_size_ > kMaxUsedSpaceSize) return true;

    // this node has enough space but cleaning up is required
    CleanUp();
    return false;
  }

  /**
   * @retval true if this node requires merging.
   * @retval false otherwise.
   */
  [[nodiscard]] auto
  NeedMerge()  //
      -> bool
  {
    // check this node uses enough space
    if (GetUsedSize() < kMinUsedSpaceSize) return true;
    if (deleted_size_ <= kMaxDeletedSpaceSize) return false;

    // this node has a lot of dead space
    CleanUp();
    return false;
  }

  /**
   * @return true r_node can be merged into this node
   * @return false otherwise
   */
  auto
  CanMerge(PessimisticNode *r_node) const  //
      -> bool
  {
    return GetUsedSize() + r_node->GetUsedSize() < kMaxUsedSpaceSize;
  }

  /**
   * @return the number of records in this node.
   */
  [[nodiscard]] constexpr auto
  GetRecordCount() const  //
      -> size_t
  {
    return record_count_;
  }

  /**
   * @return pointer of next node with shared lock
   */
  [[nodiscard]] auto
  GetNextNodeForRead()  //
      -> PessimisticNode *
  {
    next_->mutex_.LockShared();
    mutex_.UnlockShared();
    return next_;
  }

  /**
   * @return pointer of next node with exclusive lock
   */
  [[nodiscard]] auto
  GetValidSplitNode(const Key &key)  //
      -> PessimisticNode *
  {
    auto *node = this;
    const auto &high_key = GetHighKey();
    if (!high_key || Comp{}(*high_key, key)) {
      node = next_;
      node->mutex_.Lock();
      mutex_.Unlock();
    }

    return node;
  }

  /**
   * @return a high key in this node
   */
  [[nodiscard]] auto
  GetHighKey() const  //
      -> std::optional<Key>
  {
    if (high_meta_.key_len == 0) return std::nullopt;
    return GetKey(high_meta_);
  }

  /*####################################################################################
   * Public getters for records
   *##################################################################################*/

  /**
   * @param pos the position of a target record.
   * @return a key in a target record.
   */
  [[nodiscard]] auto
  GetKey(const size_t pos) const  //
      -> Key
  {
    return GetKey(meta_array_[pos]);
  }

  /**
   * @tparam Payload a class of payload.
   * @param pos the position of a target record.
   * @return a payload in a target record.
   */
  template <class Payload>
  [[nodiscard]] auto
  GetPayload(const size_t pos) const  //
      -> Payload
  {
    return GetPayload<Payload>(meta_array_[pos]);
  }

  [[nodiscard]] constexpr auto
  GetChildForRead(const size_t pos)  //
      -> PessimisticNode *
  {
    auto *child = GetPayload<PessimisticNode *>(pos);
    child->mutex_.LockShared();
    mutex_.UnlockShared();
    return child;
  }

  [[nodiscard]] constexpr auto
  GetChildForWrite(const size_t pos)  //
      -> PessimisticNode *
  {
    auto *child = GetPayload<PessimisticNode *>(pos);
    child->mutex_.Lock();
    return child;
  }

  [[nodiscard]] auto
  GetMergeableRightChild(  //
      const PessimisticNode *l_node,
      const size_t l_pos)  //
      -> PessimisticNode *
  {
    // check there is a right-sibling node
    const auto r_pos = l_pos + 1;
    if (r_pos == record_count_) {
      // a rightmost node is cannot be merged
      mutex_.Unlock();
      return nullptr;
    }

    // check the right-sibling node has enough capacity for merging
    auto *r_node = GetChildForWrite(r_pos);
    if (!l_node->CanMerge(r_node)) {
      mutex_.Unlock();
      r_node->mutex_.Unlock();
      return nullptr;
    }

    return r_node;
  }

  [[nodiscard]] auto
  RecordIsDeleted(const size_t pos) const  //
      -> bool
  {
    return meta_array_[pos].is_deleted;
  }

  template <class Payload>
  [[nodiscard]] auto
  GetRecord(const size_t pos) const  //
      -> std::pair<Key, Payload>
  {
    return {GetKey(pos), GetPayload<Payload>(pos)};
  }

  /*####################################################################################
   * Public lock management APIs
   *##################################################################################*/

  void
  AcquireSharedLock()
  {
    mutex_.LockShared();
  }

  auto
  ReleaseSharedLock()  //
      -> void
  {
    mutex_.UnlockShared();
  }

  void
  AcquireExclusiveLock()
  {
    mutex_.Lock();
  }

  void
  ReleaseExclusiveLock()
  {
    mutex_.Unlock();
  }

  /*####################################################################################
   * Public utilities
   *##################################################################################*/

  /**
   * @brief Get the position of a specified key by using binary search. If there is no
   * specified key, this returns the minimum metadata position that is greater than the
   * specified key
   *
   * @param key a target key.
   * @return the pair of record's existence and a computed position.
   */
  [[nodiscard]] auto
  SearchRecord(const Key &key) const  //
      -> std::pair<NodeRC, size_t>
  {
    int64_t begin_pos = 0;
    int64_t end_pos = record_count_ - 1;
    while (begin_pos <= end_pos) {
      size_t pos = (begin_pos + end_pos) >> 1UL;  // NOLINT
      const auto &index_key = GetKey(meta_array_[pos]);

      if (Comp{}(key, index_key)) {  // a target key is in a left side
        end_pos = pos - 1;
      } else if (Comp{}(index_key, key)) {  // a target key is in a right side
        begin_pos = pos + 1;
      } else {  // find an equivalent key
        if (meta_array_[pos].is_deleted) return {kKeyAlreadyDeleted, pos};
        return {kKeyAlreadyInserted, pos};
      }
    }

    return {kKeyNotInserted, begin_pos};
  }

  /**
   * @brief Get the position of a specified key by using binary search. If there is no
   * specified key, this returns the minimum metadata index that is greater than the
   * specified key
   *
   * @param key a target key.
   * @param range_is_closed a flag to indicate that a target key is included.
   * @return the position of a specified key.
   */
  [[nodiscard]] auto
  SearchChild(  //
      const Key &key,
      const bool range_is_closed)  //
      -> size_t
  {
    int64_t begin_pos = 0;
    int64_t end_pos = record_count_ - 2;
    while (begin_pos <= end_pos) {
      size_t pos = (begin_pos + end_pos) >> 1UL;  // NOLINT

      const auto &index_key = GetKey(meta_array_[pos]);

      if (Comp{}(key, index_key)) {  // a target key is in a left side
        end_pos = pos - 1;
      } else if (Comp{}(index_key, key)) {  // a target key is in a right side
        begin_pos = pos + 1;
      } else {  // find an equivalent key
        if (!range_is_closed) ++pos;
        begin_pos = pos;
        break;
      }
    }

    return begin_pos;
  }

  /**
   * @brief Get the end position of records for scanning and check it has been finished.
   *
   * @param end_key a pair of a target key and its closed/open-interval flag.
   * @retval 1st: true if this node is end of scanning.
   * @retval 2nd: the end position for scanning.
   */
  [[nodiscard]] auto
  SearchEndPositionFor(const std::optional<std::pair<const Key &, bool>> &end_key)  //
      -> std::pair<bool, size_t>
  {
    const auto is_end = IsRightmostOf(end_key);
    size_t end_pos{};
    if (is_end && end_key) {
      const auto &[e_key, e_closed] = *end_key;
      const auto [rc, pos] = SearchRecord(e_key);
      end_pos = (rc == kKeyAlreadyInserted && e_closed) ? pos + 1 : pos;
    } else {
      end_pos = record_count_;
    }

    return {is_end, end_pos};
  }

  /*####################################################################################
   * Leaf read operations
   *##################################################################################*/

  /**
   * @brief Read a payload of a specified key if it exists.
   *
   * @tparam Payload a class of payload.
   * @param key a target key.
   * @param out_payload a reference to be stored a target payload.
   * @retval kSuccess if a target record is read.
   * @retval kKeyNotExist otherwise.
   */
  template <class Payload>
  auto
  Read(  //
      const Key &key,
      Payload &out_payload)  //
      -> ReturnCode
  {
    const auto [existence, pos] = SearchRecord(key);
    auto rc = kKeyNotExist;
    if (existence == kKeyAlreadyInserted) {
      const auto meta = meta_array_[pos];
      memcpy(&out_payload, GetPayloadAddr(meta), sizeof(Payload));
      rc = kSuccess;
    }

    mutex_.UnlockShared();
    return rc;
  }

  /*####################################################################################
   * Leaf write operations
   *##################################################################################*/

  /**
   * @brief Write (i.e., upsert) a specified kay/payload pair.
   *
   * If a specified key does not exist in a leaf node, this function performs an insert
   * operation. If a specified key has been already inserted, this function perfroms an
   * update operation.
   *
   * Note that if a target key/payload is binary data, it is required to specify its
   * length in bytes.
   *
   * @tparam Payload a class of payload.
   * @param key a target key to be written.
   * @param key_len the length of a target key.
   * @param payload a target payload to be written.
   */
  template <class Payload>
  void
  Write(  //
      const Key &key,
      const size_t key_len,
      const Payload &payload)
  {
    // search position where this key has to be set
    const auto [rc, pos] = SearchRecord(key);

    // perform insert or update operation
    if (rc == kKeyNotInserted) {
      InsertRecord(key, key_len, payload, pos);
    } else if (rc == kKeyAlreadyDeleted) {
      ReuseRecord(payload, pos);
    } else {  // update operation
      memcpy(GetPayloadAddr(meta_array_[pos]), &payload, sizeof(Payload));
    }

    mutex_.Unlock();
  }

  /**
   * @brief Insert a specified kay/payload pair.
   *
   * This function performs a uniqueness check in its processing. If a specified key does
   * not exist, this function insert a target payload into a target leaf node. If a
   * specified key exists in a target leaf node, this function does nothing and returns
   * kKeyAlreadyInserted as a return code.
   *
   * Note that if a target key/payload is binary data, it is required to specify its
   * length in bytes.
   *
   * @tparam Payload a class of payload.
   * @param key a target key to be written.
   * @param key_len the length of a target key.
   * @param payload a target payload to be written.
   * @retval kSuccess if a key/payload pair is written.
   * @retval kKeyExist otherwise.
   */
  template <class Payload>
  auto
  Insert(  //
      const Key &key,
      const size_t key_len,
      const Payload &payload)  //
      -> ReturnCode
  {
    // search position where this key has to be set
    const auto [existence, pos] = SearchRecord(key);

    // perform insert operation if possible
    auto rc = kSuccess;
    if (existence == kKeyNotInserted) {
      InsertRecord(key, key_len, payload, pos);
    } else if (existence == kKeyAlreadyDeleted) {
      ReuseRecord(payload, pos);
    } else {  // a target key has been inserted
      rc = kKeyExist;
    }

    mutex_.Unlock();
    return rc;
  }

  /**
   * @brief Update a target kay with a specified payload.
   *
   * This function performs a uniqueness check in its processing. If a specified key
   * exist, this function update a target payload. If a specified key does not exist in
   * a target leaf node, this function does nothing and returns kKeyNotInserted as a return
   * code.
   *
   * Note that if a target key/payload is binary data, it is required to specify its
   * length in bytes.
   *
   * @tparam Payload a class of payload.
   * @param key a target key to be written.
   * @param key_length the length of a target key.
   * @param payload a target payload to be written.
   * @retval kSuccess if a key/payload pair is written.
   * @retval kKeyNotExist otherwise.
   */
  template <class Payload>
  auto
  Update(  //
      const Key &key,
      const Payload &payload)  //
      -> ReturnCode
  {
    // check this node has a target record
    const auto [existence, pos] = SearchRecord(key);

    // perform update operation if possible
    auto rc = kKeyNotExist;
    if (existence == kKeyAlreadyInserted) {
      memcpy(GetPayloadAddr(meta_array_[pos]), &payload, sizeof(Payload));
      rc = kSuccess;
    }

    mutex_.Unlock();
    return rc;
  }

  /**
   * @brief Delete a target key from the index.
   *
   * This function performs a uniqueness check in its processing. If a specified key
   * exist, this function deletes it. If a specified key does not exist in a leaf node,
   * this function does nothing and returns kKeyNotInserted as a return code.
   *
   * @tparam Payload a class of payload.
   * @param key a target key to be written.
   * @retval kSuccess if a record is deleted.
   * @retval kKeyNotExist otherwise.
   */
  auto
  Delete(const Key &key)  //
      -> ReturnCode
  {
    // check this node has a target record
    const auto [existence, pos] = SearchRecord(key);

    // perform update operation if possible
    auto rc = kKeyNotExist;
    if (existence == kKeyAlreadyInserted) {
      meta_array_[pos].is_deleted = 1;
      deleted_size_ += meta_array_[pos].rec_len + kMetaLen;
      rc = kSuccess;
    }

    mutex_.Unlock();
    return rc;
  }

  /**
   * @brief Insert a new child node to this node.
   *
   * @param l_node a left child node.
   * @param r_node a right child (i.e., new) node.
   * @param pos the position of a left child node.
   */
  void
  InsertChild(  //
      const PessimisticNode *l_node,
      const PessimisticNode *r_node,
      const size_t pos)  //
  {
    // update a current pointer and prepare space for a left child
    memcpy(GetPayloadAddr(meta_array_[pos]), &r_node, kPtrLen);
    memmove(&(meta_array_[pos + 1]), &(meta_array_[pos]), kMetaLen * (record_count_ - pos));

    // insert a left child
    const auto l_high_meta = l_node->high_meta_;
    const auto key_len = l_high_meta.key_len;
    const auto rec_len = key_len + kPtrLen;
    auto offset = SetPayload(kPageSize - block_size_, l_node);
    offset = CopyKeyFrom(l_node, l_high_meta, offset);
    meta_array_[pos] = Metadata{offset, key_len, rec_len};

    // update header information
    ++record_count_;
    block_size_ += rec_len;

    mutex_.Unlock();
  }

  /**
   * @brief Delete child node in this node.
   *
   * @param l_node a left child node.
   * @param pos the position of a right child node.
   */
  void
  DeleteChild(  //
      const PessimisticNode *l_node,
      const size_t pos)  //
  {
    const auto del_rec_len = meta_array_[pos].rec_len;  // keep a record length to be deleted

    // delete a right child by shifting metadata
    memmove(&(meta_array_[pos]), &(meta_array_[pos + 1]), kMetaLen * (record_count_ - pos));
    memcpy(GetPayloadAddr(meta_array_[pos]), &l_node, kPtrLen);

    // update this header
    --record_count_;
    deleted_size_ += del_rec_len;

    mutex_.Unlock();
  }

  /*####################################################################################
   * Public structure modification operations
   *##################################################################################*/

  /**
   * @brief Split this node into right nodes.
   *
   * @param r_node a split right node.
   */
  void
  Split(PessimisticNode *r_node)
  {
    const auto half_size = (kMetaLen * record_count_ + block_size_ - deleted_size_) / 2;

    // copy left half records to a temporal node
    size_t offset = kPageSize;
    size_t pos = 0;
    size_t used_size = 0;
    while (used_size < half_size) {
      const auto meta = meta_array_[pos++];
      if (!meta.is_deleted) {
        offset = temp_node_->CopyRecordFrom(this, meta, offset);
        used_size += meta.rec_len + kMetaLen;
      }
    }
    const auto sep_key_len = meta_array_[pos - 1].key_len;
    temp_node_->high_meta_ = Metadata{offset, sep_key_len, sep_key_len};
    if (!is_leaf_) {
      const auto rightmost_pos = temp_node_->record_count_ - 1;
      temp_node_->meta_array_[rightmost_pos] = Metadata{offset + sep_key_len, 0, kPtrLen};
    }

    // copy right half records to a right node
    auto r_offset = r_node->CopyHighKeyFrom(this);
    r_offset = r_node->CopyRecordsFrom(this, pos, record_count_, r_offset);

    // update a right header
    r_node->block_size_ = kPageSize - r_offset;
    r_node->next_ = next_;

    // update a header
    block_size_ = kPageSize - offset;
    deleted_size_ = 0;
    record_count_ = temp_node_->record_count_;
    next_ = r_node;

    // copy temporal node to this node
    memcpy(&high_meta_, &(temp_node_->high_meta_), kMetaLen * (record_count_ + 1));
    memcpy(ShiftAddr(this, offset), ShiftAddr(temp_node_.get(), offset), block_size_);

    // reset a temp node
    temp_node_->record_count_ = 0;
  }

  /**
   * @brief Merge given node into this node.
   *
   * @param r_node a right node to be merged.
   */
  void
  Merge(const PessimisticNode *r_node)
  {
    // copy a highest key of a merged node to a temporal node
    auto offset = temp_node_->CopyHighKeyFrom(r_node);

    // copy consolidated records to the original node
    offset = temp_node_->CopyRecordsFrom(this, 0, record_count_, offset);
    record_count_ = temp_node_->record_count_;
    if (!is_leaf_) {
      offset = temp_node_->CopyKeyFrom(this, high_meta_, offset);
      const auto key_len = high_meta_.key_len;
      temp_node_->meta_array_[record_count_ - 1] = Metadata{offset, key_len, key_len + kPtrLen};
    }
    memcpy(&high_meta_, &(temp_node_->high_meta_), kMetaLen * (record_count_ + 1));
    memcpy(ShiftAddr(this, offset), ShiftAddr(temp_node_.get(), offset), kPageSize - offset);

    // copy right records to this nodes
    offset = CopyRecordsFrom(r_node, 0, r_node->record_count_, offset);

    // update a header
    block_size_ = kPageSize - offset;
    deleted_size_ = 0;
    next_ = r_node->next_;

    // reset a temp node
    temp_node_->record_count_ = 0;
  }

  /*####################################################################################
   * Public bulkload APIs
   *##################################################################################*/

  /**
   * @brief Create a leaf node with the maximum number of records for bulkloading.
   *
   * @tparam Payload a target payload class.
   * @param iter the begin position of target records.
   * @param iter_end the end position of target records.
   * @param l_node the left node of this node.
   */
  template <class Entry, class Payload>
  void
  Bulkload(  //
      typename std::vector<Entry>::const_iterator &iter,
      const typename std::vector<Entry>::const_iterator &iter_end,
      const bool is_rightmost,
      PessimisticNode *l_node = nullptr)
  {
    // extract and insert entries for the leaf node
    size_t node_size = kHeaderLen;
    auto offset = kPageSize;
    while (iter < iter_end) {
      const auto &[key, payload, key_len] = ParseEntry<Payload>(*iter);
      // const auto key_len = iter->GetKeyLength();
      const auto rec_len = key_len + sizeof(Payload);

      // check whether the node has sufficient space
      node_size += rec_len + kMetaLen;
      if (node_size > kPageSize - kMinFreeSpaceSize) break;

      // insert an entry to the leaf node
      auto tmp_offset = SetPayload(offset, payload);
      tmp_offset = SetKey(tmp_offset, key, key_len);
      meta_array_[record_count_] = Metadata{tmp_offset, key_len, rec_len};
      offset -= rec_len;

      ++record_count_;
      ++iter;
    }

    // set a highest key if needed
    if (iter < iter_end || !is_rightmost) {
      const auto high_meta = meta_array_[record_count_ - 1];
      const auto high_key_len = high_meta.key_len;
      high_meta_ = Metadata{high_meta.offset, high_key_len, high_key_len};
      offset -= high_key_len;
    }

    block_size_ = kPageSize - offset;
    // set this node to next node of l_node
    if (l_node) l_node->next_ = this;
  }

  /**
   * @brief Create a inner node with the maximum number of records for bulkloading.
   *
   * @param iter the begin position of target records.
   * @param iter_end the end position of target records.
   * @param l_node the left node of this node.
   */
  void
  Bulkload(  //
      typename std::vector<PessimisticNode *>::const_iterator &iter,
      const typename std::vector<PessimisticNode *>::const_iterator &iter_end,
      PessimisticNode *l_node = nullptr)
  {
    size_t node_size = kHeaderLen;
    auto offset = kPageSize;
    while (iter < iter_end) {
      const auto key_len = (*iter)->high_meta_.key_len;
      const auto rec_len = key_len + sizeof(PessimisticNode *);

      // check whether the node has sufficient space
      node_size += rec_len + kMetaLen;
      if (node_size > kPageSize - kMinFreeSpaceSize) break;

      // insert an entry to the inner node
      auto tmp_offset = SetPayload(offset, *iter);
      if (key_len) tmp_offset = SetKey(tmp_offset, (*iter)->GetHighKey().value(), key_len);
      meta_array_[record_count_] = Metadata{tmp_offset, key_len, rec_len};
      offset -= rec_len;

      ++record_count_;
      ++iter;
    }

    // set a highest key
    const auto high_meta = meta_array_[record_count_ - 1];
    const auto high_key_len = high_meta.key_len;
    high_meta_ = Metadata{high_meta.offset, high_key_len, high_key_len};
    offset -= high_key_len;
    block_size_ = kPageSize - offset;

    // set this node to next node of l_node
    if (l_node) l_node->next_ = this;
  }

 private:
  /*####################################################################################
   * Internal constants
   *##################################################################################*/

  /// the length of child pointers.
  static constexpr size_t kPtrLen = sizeof(PessimisticNode *);

  /// the length of metadata.
  static constexpr size_t kMetaLen = sizeof(Metadata);

  /// the length of a header in each node page.
  static constexpr size_t kHeaderLen = sizeof(PessimisticNode);

  /// the maximum usage of record space.
  static constexpr size_t kMaxUsedSpaceSize = kPageSize - (kHeaderLen + kMinFreeSpaceSize);

  /*####################################################################################
   * Internal getter for header information
   *##################################################################################*/

  /**
   * @param end_key a pair of a target key and its closed/open-interval flag.
   * @retval true if this node is a rightmost node for the given key.
   * @retval false otherwise.
   */
  [[nodiscard]] auto
  IsRightmostOf(const std::optional<std::pair<const Key &, bool>> &end_key) const  //
      -> bool
  {
    if (!next_) return true;     // the rightmost node
    if (!end_key) return false;  // perform full scan
    const auto &high_key = GetKey(high_meta_);
    return !Comp{}(high_key, end_key->first);
  }

  [[nodiscard]] auto
  GetUsedSize() const  //
      -> size_t
  {
    return kMetaLen * record_count_ + block_size_ - deleted_size_;
  }

  /*####################################################################################
   * Internal getter/setters for records
   *##################################################################################*/

  /**
   * @param meta metadata of a corresponding record.
   * @return an address of a target key.
   */
  [[nodiscard]] constexpr auto
  GetKeyAddr(const Metadata meta) const  //
      -> void *
  {
    return ShiftAddr(this, meta.offset);
  }

  /**
   * @param meta metadata of a corresponding record.
   * @return a key in a target record.
   */
  [[nodiscard]] auto
  GetKey(const Metadata meta) const  //
      -> Key
  {
    if constexpr (IsVarLenData<Key>()) {
      return reinterpret_cast<Key>(GetKeyAddr(meta));
    } else {
      Key key{};
      memcpy(&key, GetKeyAddr(meta), sizeof(Key));
      return key;
    }
  }

  /**
   * @param meta metadata of a corresponding record.
   * @return an address of a target payload.
   */
  [[nodiscard]] constexpr auto
  GetPayloadAddr(const Metadata meta) const  //
      -> void *
  {
    return ShiftAddr(this, meta.offset + meta.key_len);
  }

  /**
   * @tparam Payload a class of payload.
   * @param meta metadata of a corresponding record.
   * @return a payload in a target record.
   */
  template <class Payload>
  [[nodiscard]] auto
  GetPayload(const Metadata meta) const  //
      -> Payload
  {
    Payload payload{};
    memcpy(&payload, GetPayloadAddr(meta), sizeof(Payload));
    return payload;
  }

  /**
   * @brief Set a target key directly
   *
   * @param offset an offset to set a target key
   * @param key a target key to be set
   * @param key_len a target key length to be set
   */
  auto
  SetKey(  //
      size_t offset,
      const Key &key,
      [[maybe_unused]] const size_t key_len)  //
      -> size_t
  {
    if constexpr (IsVarLenData<Key>()) {
      offset -= key_len;
      memcpy(ShiftAddr(this, offset), key, key_len);
    } else {
      offset -= sizeof(Key);
      memcpy(ShiftAddr(this, offset), &key, sizeof(Key));
    }

    return offset;
  }

  /**
   * @brief Set a target payload directly
   *
   * @tparam Payload a class of payload
   * @param offset an offset to set a target payload
   * @param payload a target payload to be set
   * @param pay_len a target payload length to be set
   */
  template <class Payload>
  auto
  SetPayload(  //
      size_t offset,
      const Payload &payload)  //
      -> size_t
  {
    offset -= sizeof(Payload);
    memcpy(ShiftAddr(this, offset), &payload, sizeof(Payload));
    return offset;
  }

  /*####################################################################################
   * Internal utility functions
   *##################################################################################*/

  template <class Payload>
  void
  InsertRecord(  //
      const Key &key,
      const size_t key_len,
      const Payload &payload,
      const size_t pos)
  {
    const auto rec_len = key_len + sizeof(Payload);

    // insert a new record
    auto offset = kPageSize - block_size_;
    offset = SetPayload(offset, payload);
    offset = SetKey(offset, key, key_len);
    memmove(&(meta_array_[pos + 1]), &(meta_array_[pos]), kMetaLen * (record_count_ - pos));
    meta_array_[pos] = Metadata{offset, key_len, rec_len};

    // update header information
    ++record_count_;
    block_size_ += rec_len;
  }

  template <class Payload>
  void
  ReuseRecord(  //
      const Payload &payload,
      const size_t pos)
  {
    // reuse a deleted record
    const auto meta = meta_array_[pos];
    meta_array_[pos].is_deleted = 0;
    memcpy(GetPayloadAddr(meta), &payload, sizeof(Payload));

    // update header information
    deleted_size_ -= meta.rec_len + kMetaLen;
  }

  /**
   * @brief Clean up this node.
   *
   * @tparam Payload a class of payload.
   */
  void
  CleanUp()
  {
    // copy records to a temporal node
    auto offset = temp_node_->CopyHighKeyFrom(this);
    offset = temp_node_->CopyRecordsFrom(this, 0, record_count_, offset);

    // update a header
    block_size_ = kPageSize - offset;
    deleted_size_ = 0;
    record_count_ = temp_node_->record_count_;

    // copy cleaned up records to the original node
    memcpy(&high_meta_, &(temp_node_->high_meta_), kMetaLen * (record_count_ + 1));
    memcpy(ShiftAddr(this, offset), ShiftAddr(temp_node_.get(), offset), block_size_);

    // reset a temp node
    temp_node_->record_count_ = 0;
  }

  /**
   * @brief Copy a key from a given node.
   *
   * @param node an original node that has a target key.
   * @param meta the corresponding metadata of a target key.
   * @param offset the current offset of this node.
   * @return the updated offset value.
   */
  auto
  CopyKeyFrom(  //
      const PessimisticNode *node,
      const Metadata meta,
      size_t offset)  //
      -> size_t
  {
    // copy a record from the given node
    if constexpr (IsVarLenData<Key>()) {
      const auto key_len = meta.key_len;
      offset -= key_len;
      memcpy(ShiftAddr(this, offset), node->GetKeyAddr(meta), key_len);
    } else {
      offset -= sizeof(Key);
      memcpy(ShiftAddr(this, offset), node->GetKeyAddr(meta), sizeof(Key));
    }

    return offset;
  }

  /**
   * @brief Copy a high key from a given node.
   *
   * @param node an original node that has a high key.
   * @return the updated offset value.
   */
  auto
  CopyHighKeyFrom(const PessimisticNode *node)  //
      -> size_t
  {
    const auto meta = node->high_meta_;
    const auto high_key_len = meta.key_len;
    const auto offset = (high_key_len > 0) ? CopyKeyFrom(node, meta, kPageSize) : kPageSize;
    high_meta_ = Metadata{offset, high_key_len, high_key_len};

    return offset;
  }

  /**
   * @brief Copy a record from a given node.
   *
   * @param orig_node an original node that has a target record.
   * @param meta the corresponding metadata of a target record.
   * @param offset the current offset of this node.
   * @return the updated offset value.
   */
  auto
  CopyRecordFrom(  //
      const PessimisticNode *node,
      const Metadata meta,
      size_t offset)  //
      -> size_t
  {
    // copy a record from the given node
    const auto rec_len = meta.rec_len;
    offset -= rec_len;
    memcpy(ShiftAddr(this, offset), node->GetKeyAddr(meta), rec_len);

    // set new metadata
    meta_array_[record_count_++] = Metadata{offset, meta.key_len, rec_len};

    return offset;
  }

  /**
   * @brief Copy records from a given node.
   *
   * @param orig_node an original node that has target records.
   * @param begin_pos the begin position of target records.
   * @param end_pos the end position of target records.
   * @param rec_count the current number of records in this node.
   * @param offset the current offset of this node.
   * @return the updated offset value.
   */
  auto
  CopyRecordsFrom(  //
      const PessimisticNode *orig_node,
      const size_t begin_pos,
      const size_t end_pos,
      size_t offset)  //
      -> size_t
  {
    // copy records from the given node
    for (size_t i = begin_pos; i < end_pos; ++i) {
      const auto meta = orig_node->meta_array_[i];
      if (!meta.is_deleted) {
        offset = CopyRecordFrom(orig_node, meta, offset);
      }
    }

    return offset;
  }

  /**
   * @brief Parse an entry of bulkload according to key's type.
   *
   * @tparam Payload a payload type.
   * @tparam Entry std::pair or std::tuple for containing entries.
   * @param entry a bulkload entry.
   * @retval 1st: a target key.
   * @retval 2nd: a target payload.
   * @retval 3rd: the length of a target key.
   */
  template <class Payload, class Entry>
  constexpr auto
  ParseEntry(const Entry &entry)  //
      -> std::tuple<Key, Payload, size_t>
  {
    if constexpr (IsVarLenData<Key>()) {
      return entry;
    } else {
      const auto &[key, payload] = entry;
      return {key, payload, sizeof(Key)};
    }
  }

  /*####################################################################################
   * Internal member variables
   *##################################################################################*/

  /// a flag to indicate whether this node is a leaf or internal node.
  uint64_t is_leaf_ : 1;

  /// the number of records in this node.
  uint64_t record_count_ : 16;

  /// the total byte length of records in a node.
  uint64_t block_size_ : 16;

  /// the total byte length of deleted records in a node.
  uint64_t deleted_size_ : 16;

  /// a black block for alignment.
  uint64_t : 0;

  /// the latch this node.
  Lock mutex_{};

  /// the pointer to the next node.
  PessimisticNode *next_{nullptr};

  /// the metadata of a highest key.
  Metadata high_meta_{kPageSize, 0, 0};

  /// an actual data block (it starts with record metadata).
  Metadata meta_array_[0];

  // temporary node for SMO
  static thread_local inline std::unique_ptr<PessimisticNode> temp_node_ =  // NOLINT
      std::make_unique<PessimisticNode>(0);
};

}  // namespace dbgroup::index::b_tree::component

#endif  // B_TREE_COMPONENT_NODE_HPP
