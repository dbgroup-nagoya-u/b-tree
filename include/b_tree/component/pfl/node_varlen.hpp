
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

#ifndef B_TREE_COMPONENT_PFL_NODE_VARLEN_HPP
#define B_TREE_COMPONENT_PFL_NODE_VARLEN_HPP

#include <atomic>
#include <optional>
#include <utility>
#include <vector>

// organization libraries
#include "lock/pessimistic_lock.hpp"

// local sources
#include "b_tree/component/common.hpp"
#include "b_tree/component/metadata.hpp"

namespace dbgroup::index::b_tree::component::pfl
{

/**
 * @brief A class for representing nodes in B+trees.
 *
 * This class uses pessimistic fine-grained locking for concurrency controls and can
 * contain variable-length data.
 *
 * @tparam Key a target key class.
 * @tparam Comp a comparetor class for keys.
 */
template <class Key, class Comp>
class NodeVarLen
{
 public:
  /*####################################################################################
   * Type aliases
   *##################################################################################*/

  using Node = NodeVarLen;

  /*####################################################################################
   * Public constructors and assignment operators
   *##################################################################################*/

  /**
   * @brief Construct an empty node object.
   *
   * @param is_leaf a flag to indicate whether a leaf node is constructed.
   */
  constexpr explicit NodeVarLen(const uint32_t is_leaf)
      : is_leaf_{is_leaf}, is_removed_{0}, block_size_{0}
  {
  }

  /**
   * @brief Construct a new root node.
   *
   * @param l_node a left child node which is the previous root node.
   * @param r_node a right child node.
   */
  NodeVarLen(  //
      const Key &l_key,
      const size_t l_key_len,
      const NodeVarLen *l_node,
      const NodeVarLen *r_node)  //
      : is_leaf_{0}, is_removed_{0}, record_count_{2}
  {
    // insert l_node
    auto offset = SetPayload(kPageSize, &l_node, kPtrLen);
    offset = SetKey(offset, l_key, l_key_len);
    meta_array_[0] = Metadata{offset, l_key_len, l_key_len + kPtrLen};

    // insert r_node
    offset = SetPayload(offset, &r_node, kPtrLen);
    meta_array_[1] = Metadata{offset, 0, kPtrLen};

    // update header information
    block_size_ = kPageSize - offset;
  }

  NodeVarLen(const NodeVarLen &) = delete;
  NodeVarLen(NodeVarLen &&) = delete;

  auto operator=(const NodeVarLen &) -> NodeVarLen & = delete;
  auto operator=(NodeVarLen &&) -> NodeVarLen & = delete;

  /*####################################################################################
   * Public destructors
   *##################################################################################*/

  /**
   * @brief Destroy the node object.
   *
   */
  ~NodeVarLen() = default;

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
   * @return true if this is a leaf node.
   * @return false otherwise.
   */
  [[nodiscard]] constexpr auto
  IsLeaf() const  //
      -> bool
  {
    return is_leaf_;
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
   * @return the next node with a shared lock.
   */
  [[nodiscard]] auto
  GetNextNodeForRead()  //
      -> Node *
  {
    next_->mutex_.LockS();
    mutex_.UnlockS();
    return next_;
  }

  /**
   * @brief Get a split node that includes a target key.
   *
   * The returned node is locked with an SIX lock and the other is unlocked.
   *
   * @param key a search key.
   * @return this node or a right sibling one.
   */
  [[nodiscard]] auto
  GetValidSplitNode(const Key &key)  //
      -> std::tuple<Node *, Key, size_t>
  {
    const auto &[sep_key, sep_key_len] = GetHighKeyForSMOs();

    auto *node = this;
    if (Comp{}(sep_key, key)) {
      node = next_;
      node->mutex_.LockX();
      mutex_.UnlockX();
    }

    return {node, sep_key, sep_key_len};
  }

  /**
   * @brief Get a mergeable-sibling node if exist.
   *
   * The returned node is locked with an exclusive lock. If there is no mergeable node,
   * this node is unlocked.
   *
   * @retval a mergeable-sibling node if exist.
   * @retval nullptr otherwise.
   */
  [[nodiscard]] auto
  GetMergeableSiblingNode()  //
      -> Node *
  {
    // check there is a right-sibling node
    if (next_ == nullptr) {
      mutex_.UnlockX();
      return nullptr;
    }

    // check the right-sibling node has enough capacity for merging
    next_->mutex_.LockX();
    if (GetUsedSize() + next_->GetUsedSize() >= kMaxUsedSpaceSize) {
      mutex_.UnlockX();
      next_->mutex_.UnlockX();
      return nullptr;
    }

    return next_;
  }

  /**
   * @retval 1st: a highest key.
   * @retval 2nd: the length of the highest key.
   */
  [[nodiscard]] auto
  GetHighKeyForSMOs() const  //
      -> std::pair<Key, size_t>
  {
    if constexpr (IsVarLenData<Key>()) {
      // allocate space dynamically to keep a copied key
      auto *h_key = reinterpret_cast<Key>(::operator new(h_key_len_));
      memcpy(h_key, GetHighKeyAddr(), h_key_len_);
      return {h_key, h_key_len_};
    } else {
      return {GetHighKey(), h_key_len_};
    }
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

  /**
   * @brief Get a child node in a given position.
   *
   * The node is unlocked after this function.
   *
   * @param pos the position of a child node.
   * @return the child node.
   */
  [[nodiscard]] auto
  GetChild(const size_t pos)  //
      -> Node *
  {
    auto *child = GetPayload<Node *>(pos);
    mutex_.UnlockS();
    return child;
  }

  /**
   * @brief Get a leftmost child node.
   *
   * The node is unlocked after this function.
   *
   * @return the child node.
   */
  [[nodiscard]] auto
  GetLeftmostChild()  //
      -> Node *
  {
    // check this node is not removed
    auto *node = this;
    while (true) {
      node->mutex_.LockS();
      if (is_removed_ == 0) break;

      // the node was removed, so go to the next node
      auto *next = node->next_;
      node->mutex_.UnlockS();
      node = next;
    }

    auto *child = GetPayload<Node *>(0);
    node->mutex_.UnlockS();
    return child;
  }

  /**
   * @param pos the position of a record.
   * @retval true if the record is deleted.
   * @retval false otherwise.
   */
  [[nodiscard]] auto
  RecordIsDeleted(const size_t pos) const  //
      -> bool
  {
    return meta_array_[pos].is_deleted;
  }

  /**
   * @tparam Payload a class of payload.
   * @param pos the position of a target record.
   * @retval 1st: a key in a target record.
   * @retval 2nd: a payload in a target record.
   */
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

  /**
   * @brief Acquire a shared lock for this node.
   *
   */
  void
  LockS()
  {
    mutex_.LockS();
  }

  /**
   * @brief Release the shared lock for this node.
   *
   */
  auto
  UnlockS()  //
      -> void
  {
    mutex_.UnlockS();
  }

  /**
   * @brief Release the shared lock for this node.
   *
   */
  auto
  UnlockX()  //
      -> void
  {
    mutex_.UnlockX();
  }

  /**
   * @brief Acquire a shared lock with intent-exclusive locking for this node.
   *
   */
  void
  LockSIX()
  {
    mutex_.LockSIX();
  }

  /**
   * @brief Release the shared lock with intent-exclusive locking for this node.
   *
   */
  void
  UnlockSIX()
  {
    mutex_.UnlockSIX();
  }

  /**
   * @brief Upgrade the SIX lock to an exclusive lock for this node.
   *
   */
  void
  UpgradeToX()
  {
    mutex_.UpgradeToX();
  }

  /*####################################################################################
   * Public utilities
   *##################################################################################*/

  /**
   * @brief Get the position of a specified key by using binary search.
   *
   * If there is no specified key in this node, this returns the minimum position that
   * is greater than the specified key.
   *
   * @param key a search key.
   * @retval 1st: record's existence.
   * @retval 2nd: the searched position.
   */
  [[nodiscard]] auto
  SearchRecord(const Key &key) const  //
      -> std::pair<NodeRC, size_t>
  {
    const auto inner_diff = static_cast<size_t>(!static_cast<bool>(is_leaf_));

    int64_t begin_pos = 0;
    int64_t end_pos = record_count_ - 1 - inner_diff;
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
   * @brief Get a child node of a specified key by using binary search.
   *
   * If there is no specified key in this node, this returns the minimum position that
   * is greater than the specified key.
   *
   * @param key a search key.
   * @param is_closed a flag for indicating closed-interval.
   * @return the child node that includes the given key.
   */
  [[nodiscard]] auto
  SearchChild(  //
      const Key &key,
      const bool is_closed)  //
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
        if (!is_closed) ++pos;
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

  /**
   * @brief Check the key range of a given node and traverse side links if needed.
   *
   * This function acquires a shared lock for a returned node.
   *
   * @param node a current node to be locked.
   * @param key a search key.
   * @param is_closed a flag for indicating closed-interval.
   * @return a node whose key range includes the search key.
   */
  [[nodiscard]] static auto
  CheckKeyRangeAndLockForRead(  //
      Node *node,
      const Key &key,
      const bool is_closed)  //
      -> Node *
  {
    while (true) {
      node->mutex_.LockS();

      // check the node is not removed
      if (node->is_removed_ == 0) {
        // check the node includes a target key
        if (node->h_key_len_ == 0) return node;
        const auto &high_key = node->GetHighKey();
        if (Comp{}(key, high_key) || (is_closed && !Comp{}(high_key, key))) return node;
      }

      // go to the next node
      auto *next = node->next_;
      node->mutex_.UnlockS();
      node = next;
      if (node == nullptr) return node;
    }
  }

  /**
   * @brief Check the key range of a given node and traverse side links if needed.
   *
   * This function acquires an exclusive lock for a returned node.
   *
   * @param node a current node to be locked.
   * @param key a search key.
   * @return a node whose key range includes the search key.
   */
  [[nodiscard]] static auto
  CheckKeyRangeAndLockForWrite(  //
      Node *node,
      const Key &key)  //
      -> Node *
  {
    while (true) {
      node->mutex_.LockX();

      // check the node is not removed
      if (node->is_removed_ == 0) {
        // check the node includes a target key
        if (node->h_key_len_ == 0 || !Comp{}(node->GetHighKey(), key)) return node;
      }

      // go to the next node
      auto *next = node->next_;
      node->mutex_.UnlockX();
      node = next;
      if (node == nullptr) return node;
    }
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

    mutex_.UnlockS();
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
   * @param key a target key to be written.
   * @param key_len the length of a target key.
   * @param payload a target payload to be written.
   * @param pay_len the length of a target payload.
   */
  auto
  Write(  //
      const Key &key,
      const size_t key_len,
      const void *payload,
      const size_t pay_len)  //
      -> NodeRC
  {
    const auto rec_len = key_len + pay_len;
    CleanUpIfNeeded(rec_len);

    // search position where this key has to be set
    const auto [existence, pos] = SearchRecord(key);
    if (existence == kKeyNotInserted) return InsertRecord(key, key_len, payload, pay_len, pos);

    // there is a record with the same key, so reuse/update it
    if (existence == kKeyAlreadyDeleted) {
      ReuseRecord(payload, pay_len, pos);
    } else {
      memcpy(GetPayloadAddr(meta_array_[pos]), payload, pay_len);
    }

    mutex_.UnlockX();
    return kCompleted;
  }

  /**
   * @brief Insert a specified kay/payload pair.
   *
   * If a specified key does not exist, this function insert a target payload into a
   * target leaf node. If a specified key exists in a target leaf node, this function
   * does nothing and returns kKeyExist as a return code.
   *
   * @param key a target key to be written.
   * @param key_len the length of a target key.
   * @param payload a target payload to be written.
   * @param pay_len the length of a target payload.
   * @retval kSuccess if a key/payload pair is written.
   * @retval kKeyExist otherwise.
   */
  auto
  Insert(  //
      const Key &key,
      const size_t key_len,
      const void *payload,
      const size_t pay_len)  //
      -> NodeRC
  {
    const auto rec_len = key_len + pay_len;
    CleanUpIfNeeded(rec_len);

    // search position where this key has to be set
    const auto [existence, pos] = SearchRecord(key);
    if (existence == kKeyNotInserted) return InsertRecord(key, key_len, payload, pay_len, pos);

    auto rc = kKeyAlreadyInserted;
    if (existence == kKeyAlreadyDeleted) {
      // there is a deleted record with the same key, so reuse it
      ReuseRecord(payload, pay_len, pos);
      rc = kCompleted;
    }

    mutex_.UnlockX();
    return rc;
  }

  /**
   * @brief Update a target kay with a specified payload.
   *
   * If a specified key exist, this function update a target payload. If a specified key
   * does not exist in a target leaf node, this function does nothing and returns
   * kKeyNotExist as a return code.
   *
   * @param key a target key to be written.
   * @param payload a target payload to be written.
   * @param pay_len the length of a target payload.
   * @retval kSuccess if a key/payload pair is written.
   * @retval kKeyNotExist otherwise.
   */
  auto
  Update(  //
      const Key &key,
      const void *payload,
      const size_t pay_len)  //
      -> ReturnCode
  {
    // check this node has a target record
    const auto [existence, pos] = SearchRecord(key);

    // perform update operation if possible
    auto rc = kKeyNotExist;
    if (existence == kKeyAlreadyInserted) {
      memcpy(GetPayloadAddr(meta_array_[pos]), payload, pay_len);
      rc = kSuccess;
    }

    mutex_.UnlockX();
    return rc;
  }

  /**
   * @brief Delete a target key from the index.
   *
   * If a specified key exist, this function deletes it. If a specified key does not
   * exist in a leaf node, this function does nothing and returns kKeyNotExist as a
   * return code.
   *
   * @param key a target key to be written.
   * @retval kSuccess if a record is deleted.
   * @retval kKeyNotExist otherwise.
   */
  auto
  Delete(const Key &key)  //
      -> NodeRC
  {
    // check this node has a target record
    const auto [existence, pos] = SearchRecord(key);
    if (existence != kKeyAlreadyInserted) {
      mutex_.UnlockX();
      return kKeyNotInserted;
    }

    // perform delete operation
    meta_array_[pos].is_deleted = 1;
    deleted_size_ += meta_array_[pos].rec_len + kMetaLen;
    if (NeedMerge()) return kNeedMerge;

    mutex_.UnlockX();
    return kCompleted;
  }

  /**
   * @brief Insert a new child node into this node.
   *
   * @param l_node a left child node.
   * @param r_node a right child (i.e., new) node.
   * @param pos the position of the left child node.
   */
  auto
  InsertChild(  //
      const Node *r_node,
      const Key &sep_key,
      const size_t sep_key_len)  //
      -> NodeRC
  {
    // check free space in this node
    const auto rec_len = sep_key_len + kPtrLen;
    CleanUpIfNeeded(rec_len);

    // check an inserting position and concurrent SMOs
    const auto [existence, pos] = SearchRecord(sep_key);
    if (existence == kKeyAlreadyInserted) {
      // previous merging has not been applied, so unlock and retry
      mutex_.UnlockX();
      return kNeedWaitAndRetry;
    }

    // recheck free space in this node
    if (NeedSplit(rec_len)) return kNeedSplit;  // this node is full, so perform splitting

    // insert a split-right child node
    const auto *l_node = GetPayload<Node *>(pos);
    memcpy(GetPayloadAddr(meta_array_[pos]), &r_node, kPtrLen);           // update as a right child
    auto offset = SetPayload(kPageSize - block_size_, &l_node, kPtrLen);  // insert a left node
    offset = SetKey(offset, sep_key, sep_key_len);
    memmove(&(meta_array_[pos + 1]), &(meta_array_[pos]), kMetaLen * (record_count_ - pos));
    meta_array_[pos] = Metadata{offset, sep_key_len, rec_len};

    // update header information
    ++record_count_;
    block_size_ += rec_len;

    mutex_.UnlockX();
    return kCompleted;
  }

  /**
   * @brief Delete a child node from this node.
   *
   * @param l_node a left child node.
   * @param pos the position of the left child node.
   */
  auto
  DeleteChild(  //
      Node *r_child,
      const Key &del_key)  //
      -> NodeRC
  {
    // check a child node to be deleted is not rightmost
    if (h_key_len_ > 0 && !Comp{}(del_key, GetHighKey())) {
      mutex_.UnlockX();
      return kAbortMerge;
    }

    // check a position to be deleted and concurrent SMOs
    const auto [existence, l_pos] = SearchRecord(del_key);
    if (existence == kKeyNotInserted) {
      // previous splitting has not been applied, so unlock and retry
      mutex_.UnlockX();
      return kNeedWaitAndRetry;
    }

    // merging have succeeded, so unlock child nodes
    r_child->next_->mutex_.UnlockSIX();
    r_child->mutex_.UnlockSIX();

    // get a current left child node
    const auto l_meta = meta_array_[l_pos];
    const auto *l_child = GetPayload<Node *>(l_meta);

    // delete metadata of a deleted key
    const auto del_rec_len = meta_array_[l_pos].rec_len;  // keep a record length to be deleted
    const auto r_pos = l_pos + 1;
    memmove(&(meta_array_[l_pos]), &(meta_array_[r_pos]), kMetaLen * (record_count_ - r_pos));
    memcpy(GetPayloadAddr(meta_array_[l_pos]), &l_child, kPtrLen);

    // update this header
    --record_count_;
    deleted_size_ += del_rec_len;

    if (NeedMerge()) return kNeedMerge;

    mutex_.UnlockX();
    return kCompleted;
  }

  /*####################################################################################
   * Public structure modification operations
   *##################################################################################*/

  /**
   * @brief Split this node into a right node.
   *
   * @param r_node a split right node.
   */
  void
  Split(Node *r_node)
  {
    const auto half_size = (kMetaLen * record_count_ + block_size_ - deleted_size_) / 2;

    /*------------------------------
     * temporal left-split
     *----------------------------*/

    // set a lowest key for a split-left node
    auto offset = temp_node_->CopyLowKeyFrom(this);

    // copy left half records to a temporal node
    size_t pos = 0;
    for (size_t used_size = 0; used_size < half_size; ++pos) {
      const auto meta = meta_array_[pos];
      if (!meta.is_deleted) {
        offset = temp_node_->CopyRecordFrom(this, meta, offset);
        used_size += meta.rec_len + kMetaLen;
      }
    }
    const auto sep_key_len = meta_array_[pos - 1].key_len;
    temp_node_->h_key_offset_ = offset;
    temp_node_->h_key_len_ = sep_key_len;
    if (is_leaf_) {
      offset -= sep_key_len;  // reserve space for a highest key
    } else {
      const auto rightmost_pos = temp_node_->record_count_ - 1;
      temp_node_->meta_array_[rightmost_pos] = Metadata{offset + sep_key_len, 0, kPtrLen};
    }

    /*------------------------------
     * right-split
     *----------------------------*/

    // set lowest/highest keys for a split-right node
    auto r_offset = r_node->CopyHighKeyFrom(temp_node_.get(), kPageSize);
    r_node->l_key_offset_ = r_offset;
    r_node->l_key_len_ = sep_key_len;
    r_offset = r_node->CopyHighKeyFrom(this, r_offset);
    r_node->h_key_offset_ = r_offset;
    r_node->h_key_len_ = h_key_len_;

    // copy right half records to a right node
    r_offset = r_node->CopyRecordsFrom(this, pos, record_count_, r_offset);

    // update a right header
    r_node->block_size_ = kPageSize - r_offset;
    r_node->next_ = next_;

    /*------------------------------
     * reflect left-split
     *----------------------------*/

    // update a header
    block_size_ = kPageSize - offset;
    deleted_size_ = 0;
    record_count_ = temp_node_->record_count_;
    next_ = r_node;

    // copy temporal node to this node
    h_key_offset_ = temp_node_->h_key_offset_;
    h_key_len_ = sep_key_len;
    memcpy(meta_array_, temp_node_->meta_array_, kMetaLen * record_count_);
    memcpy(ShiftAddr(this, offset), ShiftAddr(temp_node_.get(), offset), block_size_);

    // reset a temp node
    temp_node_->record_count_ = 0;
  }

  /**
   * @brief Merge a given node into this node.
   *
   * @param r_node a right node to be merged.
   */
  void
  Merge(Node *r_node)
  {
    // copy lowest/highest keys of merged nodes to a temporal node
    auto offset = temp_node_->CopyLowKeyFrom(this);
    offset = temp_node_->CopyHighKeyFrom(r_node, offset);
    temp_node_->h_key_offset_ = offset;

    // copy consolidated records to the original node
    offset = temp_node_->CopyRecordsFrom(this, 0, record_count_, offset);
    record_count_ = temp_node_->record_count_;
    if (!is_leaf_) {
      offset = temp_node_->CopyHighKeyFrom(this, offset);
      temp_node_->meta_array_[record_count_ - 1] =
          Metadata{offset, h_key_len_, h_key_len_ + kPtrLen};
    }
    h_key_offset_ = temp_node_->h_key_offset_;
    h_key_len_ = r_node->h_key_len_;
    memcpy(meta_array_, temp_node_->meta_array_, kMetaLen * record_count_);
    memcpy(ShiftAddr(this, offset), ShiftAddr(temp_node_.get(), offset), kPageSize - offset);

    // copy right records to this nodes
    offset = CopyRecordsFrom(r_node, 0, r_node->record_count_, offset);

    // update a header
    block_size_ = kPageSize - offset;
    deleted_size_ = 0;
    next_ = r_node->next_;

    // update a header of a right node
    r_node->is_removed_ = 1;
    r_node->next_ = this;

    mutex_.DowngradeToSIX();
    r_node->mutex_.DowngradeToSIX();

    // reset a temp node
    temp_node_->record_count_ = 0;
  }

  void
  AbortMerge(  //
      Node *r_node,
      const Key &l_key,
      const size_t l_key_len)
  {
    // upgrade locks to abort merging
    mutex_.UpgradeToX();
    r_node->mutex_.UpgradeToX();

    // revert header information of a right node
    r_node->is_removed_ = 0;
    r_node->next_ = next_;
    r_node->mutex_.UnlockX();

    // check the top position of a record block
    size_t offset{};
    if (is_leaf_) {
      // set a highest key for leaf nodes
      auto [existence, pos] = SearchRecord(l_key);
      record_count_ = (existence == kKeyAlreadyInserted) ? pos + 1 : (((pos > 0)) ? pos : 0);
      offset = (record_count_ > 0) ? meta_array_[record_count_ - 1].offset : kPageSize;
      offset = SetKey(offset, l_key, l_key_len);
    } else {
      // remove a key of a rightmost record in internal nodes
      const auto pos = SearchRecord(l_key).second;
      offset = meta_array_[pos].offset;
      meta_array_[pos] = Metadata{offset + l_key_len, 0, kPtrLen};
      record_count_ = pos + 1;
    }

    // update header information
    block_size_ = kPageSize - offset;
    deleted_size_ = h_key_len_;
    next_ = r_node;
    h_key_offset_ = offset;
    h_key_len_ = l_key_len;

    mutex_.UnlockX();
  }

  auto
  RemoveRoot()  //
      -> Node *
  {
    mutex_.UpgradeToX();

    // get a child node as a new root node
    auto *child = GetPayload<Node *>(0);
    child->LockSIX();

    // remove this node from a tree
    is_removed_ = 1;
    next_ = nullptr;

    mutex_.UnlockX();
    return child;
  }

  /*####################################################################################
   * Public bulkload APIs
   *##################################################################################*/

  /**
   * @brief Create a leaf node with the maximum number of records for bulkloading.
   *
   * @tparam Entry a pair/tuple of keys and payloads.
   * @tparam Payload a target payload class.
   * @param iter the begin position of target records.
   * @param iter_end the end position of target records.
   * @param is_rightmost a flag for indicating this node can be rightmost.
   * @param l_node the left sibling node of this node.
   */
  template <class Entry, class Payload>
  void
  Bulkload(  //
      typename std::vector<Entry>::const_iterator &iter,
      const typename std::vector<Entry>::const_iterator &iter_end,
      const bool is_rightmost,
      Node *l_node)
  {
    // extract and insert entries for the leaf node
    auto offset = (l_node != nullptr) ? l_node->LinkNext(this) : kPageSize - kMaxKeyLen;
    auto node_size = kHeaderLen;
    for (; iter < iter_end; ++iter) {
      const auto &[key, payload, key_len] = ParseEntry<Payload>(*iter);
      const auto rec_len = key_len + sizeof(Payload);

      // check whether the node has sufficient space
      node_size += rec_len + kMetaLen;
      if (node_size + key_len > kPageSize - kMinFreeSpaceSize) break;

      // insert an entry to the leaf node
      offset = SetPayload(offset, &payload, sizeof(Payload));
      offset = SetKey(offset, key, key_len);
      meta_array_[record_count_++] = Metadata{offset, key_len, rec_len};
    }

    // set a highest key if this node is not rightmost
    if (iter < iter_end || !is_rightmost) {
      h_key_offset_ = offset;
      h_key_len_ = meta_array_[record_count_ - 1].key_len;
      offset -= h_key_len_;  // reserve space for a highest key
    }

    // update header information
    block_size_ = kPageSize - offset;
  }

  /**
   * @brief Create a inner node with the maximum number of records for bulkloading.
   *
   * @param iter the begin position of target records.
   * @param iter_end the end position of target records.
   */
  void
  Bulkload(  //
      typename std::vector<Node *>::const_iterator &iter,
      const typename std::vector<Node *>::const_iterator &iter_end,
      Node *l_node)
  {
    auto offset = (l_node != nullptr) ? l_node->LinkNext(this) : kPageSize - kMaxKeyLen;
    auto node_size = kHeaderLen;
    for (; iter < iter_end; ++iter) {
      const auto *child = *iter;
      const auto high_key_len = child->h_key_len_;
      const auto rec_len = high_key_len + kPtrLen;

      // check whether the node has sufficient space
      node_size += rec_len + kMetaLen;
      if (node_size > kPageSize - kMinFreeSpaceSize) break;

      // insert an entry to the inner node
      offset = SetPayload(offset, &child, kPtrLen);
      offset = CopyHighKeyFrom(child, offset);
      meta_array_[record_count_++] = Metadata{offset, high_key_len, rec_len};
    }

    // move a highest key to header
    const auto rightmost_pos = record_count_ - 1;
    h_key_offset_ = offset;
    h_key_len_ = meta_array_[rightmost_pos].key_len;
    meta_array_[rightmost_pos] = Metadata{offset + h_key_len_, 0, kPtrLen};

    // update header information
    block_size_ = kPageSize - offset;
  }

 private:
  /*####################################################################################
   * Internal constants
   *##################################################################################*/

  /// the maximum length of keys.
  static constexpr size_t kMaxKeyLen = (IsVarLenData<Key>()) ? kMaxVarLenDataSize : sizeof(Key);

  /// the length of child pointers.
  static constexpr size_t kPtrLen = sizeof(Node *);

  /// the length of metadata.
  static constexpr size_t kMetaLen = sizeof(Metadata);

  /// the length of a header in each node page.
  static constexpr size_t kHeaderLen = sizeof(Node);

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
    return !Comp{}(GetHighKey(), end_key->first);
  }

  void
  CleanUpIfNeeded(const size_t new_rec_len)
  {
    // check whether the node has space for a new record
    const auto total_size = kMetaLen * (record_count_ + 1) + block_size_ + new_rec_len;
    const auto used_size = total_size - deleted_size_;
    if (total_size <= kPageSize - kHeaderLen || used_size > kMaxUsedSpaceSize) return;

    // this node has a lot of dead space
    CleanUp();
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
    return total_size > kPageSize - kHeaderLen;
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
   * @return Current usage of the record block.
   */
  [[nodiscard]] auto
  GetUsedSize() const  //
      -> size_t
  {
    return kMetaLen * record_count_ + block_size_ - deleted_size_;
  }

  /**
   * @return an address of a lowest key.
   */
  [[nodiscard]] constexpr auto
  GetLowKeyAddr() const  //
      -> void *
  {
    return ShiftAddr(this, l_key_offset_);
  }

  /**
   * @retval a lowest key in this node.
   */
  [[nodiscard]] auto
  GetLowKey() const  //
      -> Key
  {
    if constexpr (IsVarLenData<Key>()) {
      return reinterpret_cast<Key>(GetLowKeyAddr());
    } else {
      Key key{};
      memcpy(&key, GetLowKeyAddr(), sizeof(Key));
      return key;
    }
  }

  /**
   * @return an address of a highest key.
   */
  [[nodiscard]] constexpr auto
  GetHighKeyAddr() const  //
      -> void *
  {
    return ShiftAddr(this, h_key_offset_);
  }

  /**
   * @retval a highest key in this node.
   */
  [[nodiscard]] auto
  GetHighKey() const  //
      -> Key
  {
    if constexpr (IsVarLenData<Key>()) {
      return reinterpret_cast<Key>(GetHighKeyAddr());
    } else {
      Key key{};
      memcpy(&key, GetHighKeyAddr(), sizeof(Key));
      return key;
    }
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
   * @brief Set a target key directly.
   *
   * @param offset an offset to the top of the record block.
   * @param key a target key to be set.
   * @param key_len the length of the key.
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
   * @brief Set a target payload directly.
   *
   * @param offset an offset to the top of the record block.
   * @param payload a target payload to be written.
   * @param pay_len the length of a target payload.
   */
  auto
  SetPayload(  //
      size_t offset,
      const void *payload,
      const size_t pay_len)  //
      -> size_t
  {
    offset -= pay_len;
    memcpy(ShiftAddr(this, offset), payload, pay_len);
    return offset;
  }

  /*####################################################################################
   * Internal utility functions
   *##################################################################################*/

  /**
   * @brief Insert a given record into this node.
   *
   * @param key a target key to be set.
   * @param key_len the length of the key.
   * @param payload a target payload to be written.
   * @param pay_len the length of a target payload.
   * @param pos an insertion position.
   */
  auto
  InsertRecord(  //
      const Key &key,
      const size_t key_len,
      const void *payload,
      const size_t pay_len,
      const size_t pos)  //
      -> NodeRC
  {
    const auto rec_len = key_len + pay_len;
    if (NeedSplit(rec_len)) return kNeedSplit;

    // insert a new record
    auto offset = kPageSize - block_size_;
    offset = SetPayload(offset, payload, pay_len);
    offset = SetKey(offset, key, key_len);
    memmove(&(meta_array_[pos + 1]), &(meta_array_[pos]), kMetaLen * (record_count_ - pos));
    meta_array_[pos] = Metadata{offset, key_len, rec_len};

    // update header information
    ++record_count_;
    block_size_ += rec_len;

    mutex_.UnlockX();
    return kCompleted;
  }

  /**
   * @brief Reuse the deleted record to insert a new record.
   *
   * @param payload a target payload to be written.
   * @param pay_len the length of a target payload.
   * @param pos the position of the deleted record.
   */
  void
  ReuseRecord(  //
      const void *payload,
      const size_t pay_len,
      const size_t pos)
  {
    // reuse a deleted record
    const auto meta = meta_array_[pos];
    meta_array_[pos].is_deleted = 0;
    memcpy(GetPayloadAddr(meta), payload, pay_len);

    // update header information
    deleted_size_ -= meta.rec_len + kMetaLen;
  }

  /**
   * @brief Clean up this node.
   *
   */
  void
  CleanUp()
  {
    // copy records to a temporal node
    auto offset = temp_node_->CopyHighKeyFrom(this, kPageSize);
    temp_node_->h_key_offset_ = offset;
    offset = temp_node_->CopyRecordsFrom(this, 0, record_count_, offset);

    // update a header
    block_size_ = kPageSize - offset;
    deleted_size_ = 0;
    record_count_ = temp_node_->record_count_;
    h_key_offset_ = temp_node_->h_key_offset_;

    // copy cleaned up records to the original node
    memcpy(meta_array_, temp_node_->meta_array_, kMetaLen * record_count_);
    memcpy(ShiftAddr(this, offset), ShiftAddr(temp_node_.get(), offset), block_size_);

    // reset a temp node
    temp_node_->record_count_ = 0;
  }

  /**
   * @brief Copy a key from a given node.
   *
   * @param node an original node that has a target key.
   * @param meta the corresponding metadata of a target key.
   * @param offset an offset to the top of the record block.
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
    const auto key_len = meta.key_len;
    offset -= key_len;
    memcpy(ShiftAddr(this, offset), node->GetKeyAddr(meta), key_len);

    return offset;
  }

  /**
   * @brief Copy a lowest key from a given node.
   *
   * @param node an original node that has a lowest key.
   * @return the updated offset value.
   */
  auto
  CopyLowKeyFrom(const Node *node)  //
      -> size_t
  {
    const auto key_len = node->l_key_len_;
    auto offset = kPageSize;
    if (key_len > 0) {
      offset -= key_len;
      memcpy(ShiftAddr(this, offset), node->GetLowKeyAddr(), key_len);
    }
    l_key_offset_ = offset;
    l_key_len_ = key_len;

    return offset;
  }

  /**
   * @brief Copy a highest key from a given node.
   *
   * @param node an original node that has a highest key.
   * @return the updated offset value.
   */
  auto
  CopyHighKeyFrom(  //
      const Node *node,
      size_t offset)  //
      -> size_t
  {
    const auto key_len = node->h_key_len_;
    if (key_len > 0) {
      offset -= key_len;
      memcpy(ShiftAddr(this, offset), node->GetHighKeyAddr(), key_len);
    }

    return offset;
  }

  /**
   * @brief Copy a record from a given node.
   *
   * @param node an original node that has a target record.
   * @param meta the corresponding metadata of a target record.
   * @param offset an offset to the top of the record block.
   * @return the updated offset value.
   */
  auto
  CopyRecordFrom(  //
      const Node *node,
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
   * @param node an original node that has target records.
   * @param begin_pos the begin position of target records.
   * @param end_pos the end position of target records.
   * @param offset an offset to the top of the record block.
   * @return the updated offset value.
   */
  auto
  CopyRecordsFrom(  //
      const Node *node,
      const size_t begin_pos,
      const size_t end_pos,
      size_t offset)  //
      -> size_t
  {
    // copy records from the given node
    for (size_t i = begin_pos; i < end_pos; ++i) {
      const auto meta = node->meta_array_[i];
      if (!meta.is_deleted) {
        offset = CopyRecordFrom(node, meta, offset);
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

  /**
   * @brief Link this node to a right sibling node.
   *
   * @param r_node a right sibling node.
   * @return an offset to a lowest key in a right sibling node.
   */
  auto
  LinkNext(Node *r_node)  //
      -> size_t
  {
    // set a sibling link in a left node
    next_ = r_node;

    // copy a highest key in a left node as a lowest key in a right node
    const auto offset = kPageSize - h_key_len_;
    memcpy(ShiftAddr(r_node, offset), GetHighKeyAddr(), h_key_len_);
    r_node->l_key_offset_ = offset;
    r_node->l_key_len_ = h_key_len_;

    return offset;
  }

  /*####################################################################################
   * Internal member variables
   *##################################################################################*/

  /// a flag for indicating this node is a leaf or internal node.
  uint32_t is_leaf_ : 1;

  /// a flag for indicating this node is removed from a tree.
  uint32_t is_removed_ : 1;

  /// the total byte length of records in a node.
  uint32_t block_size_ : 30;

  /// the total byte length of deleted records in a node.
  uint16_t deleted_size_{0};

  /// the number of records in this node.
  uint16_t record_count_{0};

  /// a lock for concurrency controls.
  ::dbgroup::lock::PessimisticLock mutex_{};

  /// the pointer to the next node.
  Node *next_{nullptr};

  uint16_t l_key_offset_{kPageSize};

  uint16_t l_key_len_{0};

  uint16_t h_key_offset_{kPageSize};

  uint16_t h_key_len_{0};

  /// an actual data block (it starts with record metadata).
  Metadata meta_array_[0];

  // a temporary node for SMOs.
  static thread_local inline std::unique_ptr<Node> temp_node_ =  // NOLINT
      std::make_unique<Node>(0);
};

}  // namespace dbgroup::index::b_tree::component::pfl

#endif  // B_TREE_COMPONENT_PFL_NODE_VARLEN_HPP
