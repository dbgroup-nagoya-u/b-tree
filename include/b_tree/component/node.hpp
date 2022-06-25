
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

#include <optional>
#include <shared_mutex>
#include <utility>

#include "common.hpp"
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
class Node
{
 public:
  /*################################################################################################
   * Public constructors/destructors
   *##############################################################################################*/

  /**
   * @brief Construct a new node object
   *
   * @param is_leaf a flag to indicate whether a leaf node is constructed
   */
  explicit Node(const bool is_leaf)
      : is_leaf_{static_cast<uint64_t>(is_leaf)},
        record_count_{0},
        deleted_count_{0},
        block_size_{0},
        deleted_size_{0}
  {
  }

  /**
   * @brief Construct a new node object when root split
   *
   * @param l_node a left child node which is previous root
   * @param r_node a right child node which is split into
   */
  Node(  //
      Node *l_node,
      const Node *r_node)  //
      : is_leaf_{0}, record_count_{2}, deleted_count_{0}, block_size_{0}, deleted_size_{0}
  {
    constexpr auto kPayLen = sizeof(Node *);

    // insert l_node
    const auto l_high_meta = l_node->high_meta_;
    const auto l_key_len = l_high_meta.key_length;
    const auto l_rec_len = l_key_len + kPayLen;
    auto offset = SetPayload(kPageSize, l_node);
    offset = CopyKeyFrom(l_node, l_high_meta, offset);
    meta_array_[0] = Metadata{0, offset, l_key_len, l_rec_len};

    // insert r_node
    offset = SetPayload(offset, r_node);
    meta_array_[1] = Metadata{0, offset, 0, kPayLen};

    block_size_ += l_rec_len + kPayLen;
  }

  Node(const Node &) = delete;
  Node(Node &&) = delete;

  auto operator=(const Node &) -> Node & = delete;
  auto operator=(Node &&) -> Node & = delete;

  /**
   * @brief Destroy the node object.
   *
   */
  ~Node() = default;

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
  CanMerge(Node *r_node)  //
      -> bool
  {
    const auto l_size =
        (sizeof(Metadata) * (record_count_ - deleted_count_)) + block_size_ - deleted_size_;
    const auto r_size = (sizeof(Metadata) * (r_node->record_count_ - r_node->deleted_count_))
                        + r_node->block_size_ - r_node->deleted_size_;

    const auto can_merge = kPageSize - (kHeaderLength + kMinFreeSpaceSize) > l_size + r_size;
    if (!can_merge) {
      r_node->mutex_.unlock();
      mutex_.unlock();
    }

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
    mutex_.lock_shared();
  }

  auto
  ReleaseSharedLock()  //
      -> void
  {
    mutex_.unlock_shared();
  }

  void
  AcquireExclusiveLock()
  {
    mutex_.lock();
  }

  void
  ReleaseExclusiveLock()
  {
    mutex_.unlock();
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
   * @return the number of valid records in this node.
   */
  [[nodiscard]] constexpr auto
  GetValidRecordCount() const  //
      -> size_t
  {
    return record_count_ - deleted_count_;
  }

  /**
   * @return nearest valid position to the right of pos
   */
  [[nodiscard]] constexpr auto
  GetRightValidRecordPos(const size_t pos) const  //
      -> size_t
  {
    for (size_t i = pos; i < record_count_; i++) {
      if (!meta_array_[i].is_deleted) return i;
    }
  }

  /**
   * @return pointer of next node with exclusive lock
   */
  [[nodiscard]] auto
  GetValidSplitNode(const Key &key)  //
      -> Node *
  {
    const auto &high_key = GetHighKey();
    if (!high_key || !Comp{}(high_key.value(), key)) {
      next_->mutex_.unlock();
      return this;
    }
    mutex_.unlock();
    return next_;
  }

  /**
   * @return pointer of next node with shared lock
   */
  [[nodiscard]] constexpr auto
  GetNextNodeForRead()  //
      -> Node *
  {
    next_->mutex_.lock_shared();
    mutex_.unlock_shared();
    return next_;
  }

  [[nodiscard]] constexpr auto
  GetChildForRead(const size_t pos)  //
      -> Node *
  {
    auto *child = GetPayload<Node *>(pos);
    child->mutex_.lock_shared();
    mutex_.unlock_shared();
    return child;
  }

  [[nodiscard]] constexpr auto
  GetChildForWrite(  //
      const size_t pos,
      const bool ops_is_del)  //
      -> std::pair<Node *, bool>
  {
    constexpr auto kKeyLen = (IsVariableLengthData<Key>()) ? kMaxVarDataSize : sizeof(Key);
    constexpr auto kRecLen = kKeyLen + sizeof(Node *) + sizeof(Metadata);

    auto *child = GetPayload<Node *>(pos);
    child->mutex_.lock();

    // check if the node has sufficient space
    const auto size = sizeof(Metadata) * child->record_count_ + child->block_size_;
    const auto keep_lock =
        (!ops_is_del && (size > kPageSize - (kHeaderLength + kRecLen)))
        || (ops_is_del
            && (size - child->deleted_size_ - sizeof(Metadata) * child->deleted_count_
                < kRecLen + kMinUsedSpaceSize));

    return {child, keep_lock};
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

    mutex_.unlock_shared();
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
    const auto total_length = key_length + sizeof(Payload);

    if (auto rc = GetSpaceStatus(total_length); rc != kHasSpace) return rc;

    // search position where this key has to be set
    const auto [rc, pos] = SearchRecord(key);

    if (rc == kKeyNotInserted) {
      // insert
      auto offset = kPageSize - block_size_;
      offset = SetPayload(offset, payload);
      offset = SetKey(offset, key, key_length);
      block_size_ += total_length;

      memmove(&(meta_array_[pos + 1]), &(meta_array_[pos]),
              sizeof(Metadata) * (record_count_ - pos));
      meta_array_[pos] = Metadata{0, offset, key_length, total_length};
      record_count_ += 1;
    } else {
      // update
      memcpy(GetPayloadAddr(meta_array_[pos]), &payload, sizeof(Payload));
      if (rc == kKeyAlreadyDeleted) {
        meta_array_[pos].is_deleted = 0;
        deleted_count_--;
      }
    }

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
    const auto total_length = key_length + sizeof(Payload);

    if (auto rc = GetSpaceStatus(total_length); rc != kHasSpace) return rc;

    // search position where this key has to be set
    const auto [rc, pos] = SearchRecord(key);

    if (rc == kKeyAlreadyInserted) return kKeyAlreadyInserted;

    // reuse dead record
    if (rc == kKeyAlreadyDeleted) {
      memcpy(GetPayloadAddr(meta_array_[pos]), &payload, sizeof(Payload));
      meta_array_[pos].is_deleted = 0;
      deleted_count_--;
      return kCompleted;
    }

    // insert record
    auto offset = kPageSize - block_size_;
    offset = SetPayload(offset, payload);
    offset = SetKey(offset, key, key_length);

    block_size_ += total_length;

    // insert metadata
    memmove(&(meta_array_[pos + 1]), &(meta_array_[pos]), sizeof(Metadata) * (record_count_ - pos));
    meta_array_[pos] = Metadata{0, offset, key_length, total_length};
    record_count_ += 1;

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

    if (rc == kKeyNotInserted || rc == kKeyAlreadyDeleted) return kKeyNotInserted;

    // update payload
    memcpy(GetPayloadAddr(meta_array_[pos]), &payload, sizeof(Payload));

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
    const auto [rc, pos] = SearchRecord(key);
    if (rc == kKeyNotInserted || rc == kKeyAlreadyDeleted) return kKeyNotInserted;

    // update a header
    meta_array_[pos].is_deleted = 1;
    deleted_size_ += meta_array_[pos].total_length;
    deleted_count_++;

    const auto used_size =
        sizeof(Metadata) * (record_count_ - deleted_count_) + block_size_ - deleted_size_;
    if (used_size < kMinUsedSpaceSize) return kNeedMerge;

    return kCompleted;
  }

  /**
   * @brief Insert a new child node to this node.
   *
   * @param l_node a left child node.
   * @param r_node a right child (i.e., new) node.
   * @param l_node_pos the position of a left child node.
   * @param r_node_pos the position of a right child node.
   * @retval kCompleted if a child node is inserted.
   * @retval kNeedConsolidation if this node needs consolidation for future insertion.
   * @retval kNeedSplit if this node needs splitting for future insertion.
   */
  auto
  InsertChild(  //
      Node *l_node,
      const Node *r_node,
      const size_t l_node_pos,
      const size_t r_node_pos)  //
      -> NodeRC
  {
    constexpr auto kPayLen = sizeof(Node *);
    const auto l_high_meta = l_node->high_meta_;
    const auto key_len = l_high_meta.key_length;
    const auto rec_len = key_len + kPayLen;

    if (auto rc = GetSpaceStatus(rec_len); rc != kHasSpace) return rc;

    // insert a right child by updating an original record
    memcpy(GetPayloadAddr(meta_array_[r_node_pos]), &r_node, kPayLen);

    // insert a left child
    auto offset = SetPayload(kPageSize - block_size_, l_node);
    offset = CopyKeyFrom(l_node, l_high_meta, offset);

    // add metadata for a left child
    memmove(&(meta_array_[l_node_pos + 1]), &(meta_array_[l_node_pos]),
            sizeof(Metadata) * (record_count_ - l_node_pos));
    meta_array_[l_node_pos] = Metadata{0, offset, key_len, rec_len};

    // update this header
    block_size_ += rec_len;
    ++record_count_;
    return kCompleted;
  }

  /**
   * @brief Delete child node in this node.
   *
   * @param l_node a left child node.
   * @param l_node_pos the position of a left child node.
   * @param r_node_pos the position of a right child node.
   * @retval kHasSpace if this node does not need any SMOs.
   * @retval kNeedConsolidation if this node needs consolidation for future insertion.
   * @retval kNeedSplit if this node needs splitting for future insertion.
   */
  auto
  DeleteChild(  //
      Node *l_node,
      const size_t l_node_pos,
      const size_t r_node_pos)  //
      -> NodeRC
  {
    constexpr auto kPayLen = sizeof(Node *);
    // update payload
    memcpy(GetPayloadAddr(meta_array_[r_node_pos]), &l_node, kPayLen);

    const auto key_len = meta_array_[l_node_pos].key_length;

    meta_array_[l_node_pos].is_deleted = 1;
    deleted_count_++;

    // update this header
    deleted_size_ += key_len + kPayLen;

    // check if this node should be merged
    const auto used_size =
        sizeof(Metadata) * (record_count_ - deleted_count_) + block_size_ - deleted_size_;
    return (used_size < kMinUsedSpaceSize) ? kNeedMerge : kCompleted;
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
    // copy records to a temporal node
    auto offset = temp_node_->CopyHighKeyFrom(this, high_meta_);
    offset = temp_node_->template CopyRecordsFrom<Payload>(this, 0, record_count_, 0, offset);

    // update a header
    block_size_ = kPageSize - offset;
    deleted_size_ = 0;
    record_count_ -= deleted_count_;
    deleted_count_ = 0;

    // copy consolidated records to the original node
    memcpy(meta_array_, temp_node_->meta_array_, sizeof(Metadata) * (record_count_));
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
  Split(Node *r_node)
  {
    // copy left half records to a temporal node
    size_t offset = kPageSize;
    size_t l_count = 0;
    size_t count = 0;
    const auto target_offset = offset - (block_size_ - deleted_size_) / 2;
    while (offset > target_offset) {
      const auto target_meta = meta_array_[count++];
      if (!target_meta.is_deleted)
        offset = temp_node_->template CopyRecordFrom<Payload>(this, target_meta, l_count++, offset);
    }
    const auto sep_key_len = meta_array_[count - 1].key_length;
    temp_node_->high_meta_ = Metadata{0, offset, sep_key_len, sep_key_len};

    // copy right half records to a right node
    auto r_offset = r_node->CopyHighKeyFrom(this, high_meta_);
    r_offset = r_node->template CopyRecordsFrom<Payload>(this, count, record_count_, 0, r_offset);

    // update a right header
    r_node->block_size_ = kPageSize - r_offset;
    r_node->record_count_ = record_count_ - deleted_count_ - l_count;
    r_node->next_ = next_;

    // update a header
    block_size_ = kPageSize - offset;
    deleted_size_ = 0;
    deleted_count_ = 0;
    record_count_ = l_count;
    next_ = r_node;

    // copy temporal node to this node
    memcpy(&high_meta_, &temp_node_->high_meta_, sizeof(Metadata) * (record_count_ + 1));
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
  Merge(const Node *r_node)
  {
    // copy a highest key of a merged node to a temporal node
    auto offset = temp_node_->CopyHighKeyFrom(r_node, r_node->high_meta_);

    // copy consolidated records to the original node
    offset = temp_node_->template CopyRecordsFrom<Payload>(this, 0, record_count_, 0, offset);
    record_count_ -= deleted_count_;
    memcpy(&high_meta_, &temp_node_->high_meta_, sizeof(Metadata) * (record_count_ + 1));
    memcpy(ShiftAddr(this, offset), ShiftAddr(temp_node_.get(), offset), kPageSize - offset);

    // copy right records to this nodes
    const auto r_count = r_node->GetRecordCount();
    offset = CopyRecordsFrom<Payload>(r_node, 0, r_count, record_count_, offset);

    // update a header
    record_count_ += r_node->GetValidRecordCount();
    block_size_ = kPageSize - offset;
    deleted_size_ = 0;
    deleted_count_ = 0;
    next_ = r_node->next_;

    mutex_.unlock();
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
    if (high_meta_.key_length == 0) return true;  // the rightmost node
    if (!end_key) return false;                   // perform full scan
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
    // check whether the node has space for a new record
    auto total_size = sizeof(Metadata) * (record_count_ + 1) + block_size_ + new_record_size;
    if (total_size <= kPageSize - kHeaderLength) return kHasSpace;

    // the node needs any SMO
    return (total_size - deleted_size_ - sizeof(Metadata) * deleted_count_
            > kPageSize - (kHeaderLength + kMinFreeSpaceSize))
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
      const Node *node,
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
      const Node *orig_node,
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
      const Node *node,
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
      const Node *orig_node,
      const size_t begin_pos,
      const size_t end_pos,
      size_t rec_count,
      size_t offset)  //
      -> size_t
  {
    // copy records from the given node
    for (size_t i = begin_pos; i < end_pos; ++i) {
      const auto target_meta = orig_node->meta_array_[i];
      if (!target_meta.is_deleted)
        offset = CopyRecordFrom<Payload>(orig_node, target_meta, rec_count++, offset);
    }

    return offset;
  }

  /*################################################################################################
   * Internal member variables
   *##############################################################################################*/

  /// a flag to indicate whether this node is a leaf or internal node.
  uint64_t is_leaf_ : 1;

  /// the number of records in this node.
  uint64_t record_count_ : 16;

  /// the number of deleted records in this node.
  uint64_t deleted_count_ : 16;

  /// the total byte length of records in a node.
  uint64_t block_size_ : 16;

  /// the total byte length of deleted records in a node.
  uint64_t deleted_size_ : 16;

  /// a black block for alignment.
  uint64_t : 0;

  /// the latch this node.
  std::shared_mutex mutex_{};

  /// the pointer to the next node.
  Node *next_{nullptr};

  /// the metadata of a highest key.
  Metadata high_meta_{0, kPageSize, 0, 0};

  /// an actual data block (it starts with record metadata).
  Metadata meta_array_[0];

  // temporary node for SMO
  static thread_local inline std::unique_ptr<Node> temp_node_ =  // NOLINT
      std::make_unique<Node>(0);
};

}  // namespace dbgroup::index::b_tree::component

#endif  // B_TREE_COMPONENT_NODE_HPP
