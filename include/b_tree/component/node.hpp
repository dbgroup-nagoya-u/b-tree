
/*
 * Copyright 2021 Database Group, Nagoya University
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
  /*####################################################################################
   * Type aliases
   *##################################################################################*/

  using Lock = ::dbgroup::lock::PessimisticLock;

 public:
  /*################################################################################################
   * Public constructors/destructors
   *##############################################################################################*/

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
    constexpr auto kPayLen = sizeof(PessimisticNode *);

    // insert l_node
    const auto l_high_meta = l_node->high_meta_;
    const auto l_key_len = l_high_meta.key_length;
    const auto l_rec_len = l_key_len + kPayLen;
    auto offset = SetPayload(kPageSize, l_node);
    offset = CopyKeyFrom(l_node, l_high_meta, offset);
    meta_array_[0] = Metadata{offset, l_key_len, l_rec_len};

    // insert r_node
    offset = SetPayload(offset, r_node);
    meta_array_[1] = Metadata{offset, 0, kPayLen};

    block_size_ += l_rec_len + kPayLen;
  }

  PessimisticNode(const PessimisticNode &) = delete;
  PessimisticNode(PessimisticNode &&) = delete;

  auto operator=(const PessimisticNode &) -> PessimisticNode & = delete;
  auto operator=(PessimisticNode &&) -> PessimisticNode & = delete;

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

  /**
   * @return true r_node can be merged into this node
   * @return false otherwise
   */
  auto
  CanMerge(PessimisticNode *r_node)  //
      -> bool
  {
    constexpr auto kMetaLen = sizeof(Metadata);

    const auto l_size = kMetaLen * record_count_ + block_size_ - deleted_size_;
    const auto r_size =
        kMetaLen * r_node->record_count_ + r_node->block_size_ - r_node->deleted_size_;

    const auto can_merge = kPageSize - (kHeaderLength + kMinFreeSpaceSize) > l_size + r_size;

    if (!can_merge) r_node->mutex_.Unlock();
    return can_merge;
  }

  /*################################################################################################
   * Public getters
   *##############################################################################################*/

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
   * @return pointer of next node with exclusive lock
   */
  [[nodiscard]] auto
  GetValidSplitNode(const Key &key)  //
      -> PessimisticNode *
  {
    const auto &high_key = GetHighKey();
    if (!high_key || !Comp{}(high_key.value(), key)) {
      next_->mutex_.Unlock();
      return this;
    }
    mutex_.Unlock();
    return next_;
  }

  /**
   * @return pointer of next node with shared lock
   */
  [[nodiscard]] constexpr auto
  GetNextNodeForRead()  //
      -> PessimisticNode *
  {
    next_->mutex_.LockShared();
    mutex_.UnlockShared();
    return next_;
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
  HasSufficientSpace(const bool ops_is_del)  //
      -> bool
  {
    constexpr auto kMetaLen = sizeof(Metadata);
    constexpr auto kKeyLen = (IsVariableLengthData<Key>()) ? kMaxVarDataSize : sizeof(Key);
    constexpr auto kRecLen = kKeyLen + sizeof(PessimisticNode *) + kMetaLen;

    // check if the node has sufficient space
    const auto size = kMetaLen * record_count_ + block_size_;
    const auto has_sufficient_space =
        (!ops_is_del && (size <= kPageSize - (kHeaderLength + kRecLen)))
        || (ops_is_del && (size - deleted_size_ >= kRecLen + kMinUsedSpaceSize));
    return has_sufficient_space;
  }

  [[nodiscard]] constexpr auto
  GetChildForWrite(  //
      const size_t pos,
      const bool ops_is_del,
      const bool unlock_parent = true)  //
      -> std::pair<PessimisticNode *, bool>
  {
    auto *child = GetPayload<PessimisticNode *>(pos);
    child->mutex_.Lock();

    // check if the node has sufficient space
    const auto should_smo = !child->HasSufficientSpace(ops_is_del);
    if (!should_smo && unlock_parent) mutex_.Unlock();
    return {child, should_smo};
  }

  /**
   * @param position the position of record metadata to be get.
   * @return Metadata: record metadata.
   */
  [[nodiscard]] constexpr auto
  GetMetadata(const size_t position) const  //
      -> Metadata
  {
    return meta_array_[position];
  }

  /**
   * @param meta metadata of a corresponding record.
   * @return a key in a target record.
   */
  [[nodiscard]] auto
  GetKey(const Metadata meta) const  //
      -> Key
  {
    if constexpr (IsVariableLengthData<Key>()) {
      return reinterpret_cast<Key>(GetKeyAddr(meta));
    } else {
      Key key{};
      memcpy(&key, GetKeyAddr(meta), sizeof(Key));
      return key;
    }
  }

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
   * @return a high key in this node
   */
  [[nodiscard]] auto
  GetHighKey() const  //
      -> std::optional<Key>
  {
    if (high_meta_.key_length == 0) return std::nullopt;
    return GetKey(high_meta_);
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
   * @tparam Payload a class of payload.
   * @param pos the position of a target record.
   * @return a payload in a target record.
   */
  template <class Payload>
  [[nodiscard]] auto
  GetPayload(const size_t pos) const  //
      -> Payload
  {
    Payload payload{};
    memcpy(&payload, GetPayloadAddr(meta_array_[pos]), sizeof(Payload));
    return payload;
  }

  template <class Payload>
  [[nodiscard]] auto
  GetRecord(const size_t pos) const  //
      -> std::pair<Key, Payload>
  {
    return {GetKey(pos), GetPayload<Payload>(pos)};
  }

  /*################################################################################################
   * Public APIs
   *##############################################################################################*/

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
   * @retval kCompleted if a key exists.
   * @retval kKeyNotInserted if a key does not exist.
   */
  template <class Payload>
  auto
  Read(  //
      const Key &key,
      Payload &out_payload)  //
      -> NodeRC
  {
    auto [rc, pos] = SearchRecord(key);
    if (rc == kKeyAlreadyInserted) {
      const auto meta = meta_array_[pos];
      memcpy(&out_payload, GetPayloadAddr(meta), sizeof(Payload));
      rc = kCompleted;
    } else {
      rc = kKeyNotInserted;
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
   * @param key_length the length of a target key.
   * @param payload a target payload to be written.
   * @retval kCompleted if a key/payload pair is written.
   * @retval kNeedConsolidation if a target node requires split.
   */
  template <class Payload>
  auto
  Write(  //
      const Key &key,
      const size_t key_length,
      const Payload &payload)  //
      -> NodeRC
  {
    constexpr auto kMetaLen = sizeof(Metadata);

    const auto total_length = key_length + sizeof(Payload);

    if (auto rc = GetSpaceStatus(total_length); rc != kHasSpace) {
      mutex_.Unlock();
      return rc;
    }

    // search position where this key has to be set
    const auto [rc, pos] = SearchRecord(key);

    if (rc == kKeyNotInserted) {
      // insert
      auto offset = kPageSize - block_size_;
      offset = SetPayload(offset, payload);
      offset = SetKey(offset, key, key_length);
      block_size_ += total_length;

      memmove(&(meta_array_[pos + 1]), &(meta_array_[pos]), kMetaLen * (record_count_ - pos));
      meta_array_[pos] = Metadata{offset, key_length, total_length};
      record_count_ += 1;
    } else {
      // update
      memcpy(GetPayloadAddr(meta_array_[pos]), &payload, sizeof(Payload));
      if (rc == kKeyAlreadyDeleted) {
        meta_array_[pos].is_deleted = 0;
        deleted_size_ -= meta_array_[pos].total_length + kMetaLen;
      }
    }

    mutex_.Unlock();
    return kCompleted;
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
   * @param key_length the length of a target key.
   * @param payload a target payload to be written.
   * @retval kCompleted if a key/payload pair is written.
   * @retval kKeyAlreadyInserted if a specified key exists.
   * @retval kNeedConsolidation if a target node requires split.
   */
  template <class Payload>
  auto
  Insert(  //
      const Key &key,
      const size_t key_length,
      const Payload &payload)  //
      -> NodeRC
  {
    constexpr auto kMetaLen = sizeof(Metadata);

    const auto total_length = key_length + sizeof(Payload);

    if (auto rc = GetSpaceStatus(total_length); rc != kHasSpace) {
      mutex_.Unlock();
      return rc;
    }
    // search position where this key has to be set
    const auto [rc, pos] = SearchRecord(key);

    if (rc == kKeyAlreadyInserted) {
      mutex_.Unlock();
      return kKeyAlreadyInserted;
    }

    // reuse dead record
    if (rc == kKeyAlreadyDeleted) {
      memcpy(GetPayloadAddr(meta_array_[pos]), &payload, sizeof(Payload));
      meta_array_[pos].is_deleted = 0;
      deleted_size_ -= meta_array_[pos].total_length + kMetaLen;
      mutex_.Unlock();
      return kCompleted;
    }

    // insert record
    auto offset = kPageSize - block_size_;
    offset = SetPayload(offset, payload);
    offset = SetKey(offset, key, key_length);

    block_size_ += total_length;

    // insert metadata
    memmove(&(meta_array_[pos + 1]), &(meta_array_[pos]), kMetaLen * (record_count_ - pos));
    meta_array_[pos] = Metadata{offset, key_length, total_length};
    record_count_ += 1;

    mutex_.Unlock();
    return kCompleted;
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
   * @retval kCompleted if a key/payload pair is written.
   * @retval kKeyNotInserted if a specified key does not exist.
   * @retval kNeedConsolidation if a target node requires split.
   */
  template <class Payload>
  auto
  Update(  //
      const Key &key,
      const Payload &payload)  //
      -> NodeRC
  {
    // search position where this key has to be set
    const auto [rc, pos] = SearchRecord(key);

    if (rc == kKeyNotInserted || rc == kKeyAlreadyDeleted) {
      mutex_.Unlock();
      return kKeyNotInserted;
    }

    // update payload
    memcpy(GetPayloadAddr(meta_array_[pos]), &payload, sizeof(Payload));
    mutex_.Unlock();
    return kCompleted;
  }

  /**
   * @brief Delete a target key from the index.
   *
   * This function performs a uniqueness check in its processing. If a specified key
   * exist, this function deletes it. If a specified key does not exist in a leaf node,
   * this function does nothing and returns kKeyNotInserted as a return code.
   *
   * Note that if a target key is binary data, it is required to specify its length in
   * bytes.
   *
   * @tparam Payload a class of payload.
   * @param key a target key to be written.
   * @retval kCompleted if a key/payload pair is written.
   * @retval kKeyNotInserted if a specified key does not exist.
   */
  auto
  Delete(const Key &key)  //
      -> NodeRC
  {
    constexpr auto kMetaLen = sizeof(Metadata);

    const auto [rc, pos] = SearchRecord(key);
    if (rc == kKeyNotInserted || rc == kKeyAlreadyDeleted) {
      mutex_.Unlock();
      return kKeyNotInserted;
    }

    // update a header
    meta_array_[pos].is_deleted = 1;
    deleted_size_ += meta_array_[pos].total_length + kMetaLen;

    const auto used_size = kMetaLen * record_count_ + block_size_ - deleted_size_;

    mutex_.Unlock();
    return (used_size < kMinUsedSpaceSize) ? kNeedMerge : kCompleted;
  }

  /**
   * @brief Insert a new child node to this node.
   *
   * @param l_node a left child node.
   * @param r_node a right child (i.e., new) node.
   * @param pos the position of a left child node.
   * @retval kCompleted if a child node is inserted.
   * @retval kNeedConsolidation if this node needs consolidation for future insertion.
   * @retval kNeedSplit if this node needs splitting for future insertion.
   */
  void
  InsertChild(  //
      PessimisticNode *l_node,
      const PessimisticNode *r_node,
      const size_t pos)  //
  {
    constexpr auto kPayLen = sizeof(PessimisticNode *);
    const auto l_high_meta = l_node->high_meta_;
    const auto key_len = l_high_meta.key_length;
    const auto rec_len = key_len + kPayLen;

    // insert a right child by updating an original record
    memcpy(GetPayloadAddr(meta_array_[pos]), &r_node, kPayLen);

    // insert a left child
    auto offset = SetPayload(kPageSize - block_size_, l_node);
    offset = CopyKeyFrom(l_node, l_high_meta, offset);

    // add metadata for a left child
    memmove(&(meta_array_[pos + 1]), &(meta_array_[pos]), sizeof(Metadata) * (record_count_ - pos));
    meta_array_[pos] = Metadata{offset, key_len, rec_len};

    // update this header
    block_size_ += rec_len;
    ++record_count_;

    mutex_.Unlock();
  }

  /**
   * @brief Delete child node in this node.
   *
   * @param l_node a left child node.
   * @param r_node a right child (i.e., to be deleted) node.
   * @param pos the position of a right child node.
   * @retval kHasSpace if this node does not need any SMOs.
   * @retval kNeedConsolidation if this node needs consolidation for future insertion.
   * @retval kNeedSplit if this node needs splitting for future insertion.
   */
  void
  DeleteChild(  //
      PessimisticNode *l_node,
      const size_t pos)  //
  {
    constexpr auto kPayLen = sizeof(PessimisticNode *);

    // update payload
    memcpy(GetPayloadAddr(meta_array_[pos]), &l_node, kPayLen);

    const auto key_len = meta_array_[pos - 1].key_length;

    // delete metadata
    memmove(&(meta_array_[pos - 1]), &(meta_array_[pos]), sizeof(Metadata) * (record_count_ - pos));

    // update this header
    deleted_size_ += key_len + kPayLen;
    --record_count_;

    mutex_.Unlock();
  }

  /*####################################################################################
   * Public structure modification operations
   *##################################################################################*/

  /**
   * @brief Consolidate this node.
   *
   * @tparam Payload a class of payload.
   */
  template <class Payload>
  void
  Consolidate()
  {
    constexpr auto kMetaLen = sizeof(Metadata);

    // copy records to a temporal node
    auto offset = temp_node_->CopyHighKeyFrom(this, high_meta_);
    offset = temp_node_->template CopyRecordsFrom<Payload>(this, 0, record_count_, 0, offset);

    // update a header
    block_size_ = kPageSize - offset;
    deleted_size_ = 0;
    record_count_ = temp_node_->record_count_;

    // copy consolidated records to the original node
    memcpy(meta_array_, temp_node_->meta_array_, kMetaLen * (record_count_));
    memcpy(ShiftAddr(this, offset), ShiftAddr(temp_node_.get(), offset), block_size_);
  }

  /**
   * @brief Split this node into right nodes.
   *
   * @tparam Payload a class of payload.
   * @param r_node a split right node.
   */
  template <class Payload>
  void
  Split(PessimisticNode *r_node)
  {
    r_node->mutex_.Lock();

    constexpr auto kMetaLen = sizeof(Metadata);
    // copy left half records to a temporal node
    size_t offset = kPageSize;
    size_t l_count = 0;
    size_t pos = 0;
    temp_node_->record_count_ = 0;
    const auto half_size = (kMetaLen * record_count_ + block_size_ - deleted_size_) / 2;
    size_t used_size = 0;
    while (half_size > used_size) {
      const auto target_meta = meta_array_[pos++];
      if (!target_meta.is_deleted) {
        offset = temp_node_->template CopyRecordFrom<Payload>(this, target_meta, l_count++, offset);
        used_size += target_meta.total_length + kMetaLen;
      }
    }
    const auto sep_key_len = meta_array_[pos - 1].key_length;
    temp_node_->high_meta_ = Metadata{offset, sep_key_len, sep_key_len};

    // copy right half records to a right node
    auto r_offset = r_node->CopyHighKeyFrom(this, high_meta_);
    r_offset = r_node->template CopyRecordsFrom<Payload>(this, pos, record_count_, 0, r_offset);

    // update a right header
    r_node->block_size_ = kPageSize - r_offset;
    r_node->next_ = next_;

    // update a header
    block_size_ = kPageSize - offset;
    deleted_size_ = 0;
    record_count_ = l_count;
    next_ = r_node;

    // copy temporal node to this node
    memcpy(&high_meta_, &temp_node_->high_meta_, kMetaLen * (record_count_ + 1));
    memcpy(ShiftAddr(this, offset), ShiftAddr(temp_node_.get(), offset), block_size_);
  }

  /**
   * @brief Merge given node into this node.
   *
   * @tparam Payload a class of payload.
   * @param r_node a right node to be merged.
   */
  template <class Payload>
  void
  Merge(const PessimisticNode *r_node)
  {
    constexpr auto kMetaLen = sizeof(Metadata);

    // copy a highest key of a merged node to a temporal node
    auto offset = temp_node_->CopyHighKeyFrom(r_node, r_node->high_meta_);

    // copy consolidated records to the original node
    offset = temp_node_->template CopyRecordsFrom<Payload>(this, 0, record_count_, 0, offset);
    record_count_ = temp_node_->record_count_;
    memcpy(&high_meta_, &temp_node_->high_meta_, kMetaLen * (record_count_ + 1));
    memcpy(ShiftAddr(this, offset), ShiftAddr(temp_node_.get(), offset), kPageSize - offset);

    // copy right records to this nodes
    const auto r_count = r_node->GetRecordCount();
    offset = CopyRecordsFrom<Payload>(r_node, 0, r_count, record_count_, offset);

    // update a header
    block_size_ = kPageSize - offset;
    deleted_size_ = 0;
    next_ = r_node->next_;
  }

  /*####################################################################################
   * Public bulkload API
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
    size_t node_size = kHeaderLength;
    auto offset = kPageSize;
    while (iter < iter_end) {
      const auto &[key, payload, key_len] = ParseEntry<Payload>(*iter);
      // const auto key_len = iter->GetKeyLength();
      const auto rec_len = key_len + sizeof(Payload);

      // check whether the node has sufficient space
      node_size += rec_len + sizeof(Metadata);
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
      const auto high_key_len = high_meta.key_length;
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
    size_t node_size = kHeaderLength;
    auto offset = kPageSize;
    while (iter < iter_end) {
      const auto key_len = (*iter)->high_meta_.key_length;
      const auto rec_len = key_len + sizeof(PessimisticNode *);

      // check whether the node has sufficient space
      node_size += rec_len + sizeof(Metadata);
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
    const auto high_key_len = high_meta.key_length;
    high_meta_ = Metadata{high_meta.offset, high_key_len, high_key_len};
    offset -= high_key_len;
    block_size_ = kPageSize - offset;

    // set this node to next node of l_node
    if (l_node) l_node->next_ = this;
  }

 private:
  /*################################################################################################
   * Internal getter/setters
   *##############################################################################################*/

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
   * @return an address of a target payload.
   */
  [[nodiscard]] constexpr auto
  GetPayloadAddr(const Metadata meta) const  //
      -> void *
  {
    return ShiftAddr(this, meta.offset + meta.key_length);
  }

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
    return !Comp{}(GetKey(high_meta_), end_key->first);
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
    if constexpr (IsVariableLengthData<Key>()) {
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

  /*################################################################################################
   * Internal utility functions
   *##############################################################################################*/

  /**
   * @retval true if this node does not have sufficient free space.
   * @retval false otherwise.
   */
  [[nodiscard]] constexpr auto
  GetSpaceStatus(const size_t new_record_size) const  //
      -> NodeRC
  {
    constexpr auto kMetaLen = sizeof(Metadata);

    // check whether the node has space for a new record
    auto total_size = kMetaLen * (record_count_ + 1) + block_size_ + new_record_size;
    if (total_size <= kPageSize - kHeaderLength) return kHasSpace;

    // the node needs any SMO
    return (total_size - deleted_size_ > kPageSize - (kHeaderLength + kMinFreeSpaceSize))
               ? kNeedSplit
               : kNeedConsolidation;
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
    if constexpr (IsVariableLengthData<Key>()) {
      const auto key_len = meta.key_length;
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
   * @param meta the corresponding metadata of a high key.
   * @param offset the current offset of this node.
   * @return the updated offset value.
   */
  auto
  CopyHighKeyFrom(  //
      const PessimisticNode *orig_node,
      const Metadata target_meta)  //
      -> size_t
  {
    const auto offset = CopyKeyFrom(orig_node, target_meta, kPageSize);
    high_meta_ = target_meta;
    high_meta_.offset = offset;

    return offset;
  }

  /**
   * @brief Copy a record from a given node.
   *
   * @tparam Payload a class of payload.
   * @param orig_node an original node that has a target record.
   * @param meta the corresponding metadata of a target record.
   * @param rec_count the current number of records in this node.
   * @param offset the current offset of this node.
   * @return the updated offset value.
   */
  template <class Payload>
  auto
  CopyRecordFrom(  //
      const PessimisticNode *node,
      const Metadata meta,
      const size_t rec_count,
      size_t offset)  //
      -> size_t
  {
    // copy a record from the given node
    if constexpr (IsVariableLengthData<Key>()) {
      const auto rec_len = meta.total_length;
      offset -= rec_len;
      memcpy(ShiftAddr(this, offset), node->GetKeyAddr(meta), rec_len);
    } else {
      constexpr auto kRecLen = sizeof(Key) + sizeof(Payload);
      offset -= kRecLen;
      memcpy(ShiftAddr(this, offset), node->GetKeyAddr(meta), kRecLen);
    }

    // set new metadata
    meta_array_[rec_count] = meta;
    meta_array_[rec_count].offset = offset;
    record_count_ = rec_count + 1;

    return offset;
  }

  /**
   * @brief Copy records from a given node.
   *
   * @tparam Payload a class of payload.
   * @param orig_node an original node that has target records.
   * @param begin_pos the begin position of target records.
   * @param end_pos the end position of target records.
   * @param rec_count the current number of records in this node.
   * @param offset the current offset of this node.
   * @return the updated offset value.
   */
  template <class Payload>
  auto
  CopyRecordsFrom(  //
      const PessimisticNode *orig_node,
      const size_t begin_pos,
      const size_t end_pos,
      size_t rec_count,
      size_t offset)  //
      -> size_t
  {
    record_count_ = rec_count;
    // copy records from the given node
    for (size_t i = begin_pos; i < end_pos; ++i) {
      const auto target_meta = orig_node->meta_array_[i];
      if (!target_meta.is_deleted) {
        offset = CopyRecordFrom<Payload>(orig_node, target_meta, rec_count++, offset);
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
    if constexpr (IsVariableLengthData<Key>()) {
      return entry;
    } else {
      const auto &[key, payload] = entry;
      return {key, payload, sizeof(Key)};
    }
  }

  /*################################################################################################
   * Internal member variables
   *##############################################################################################*/

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
