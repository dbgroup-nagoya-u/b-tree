
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

#ifndef B_TREE_COMPONENT_OML_NODE_FIXLEN_HPP
#define B_TREE_COMPONENT_OML_NODE_FIXLEN_HPP

#include <atomic>
#include <optional>
#include <utility>
#include <vector>

// organization libraries
#include "lock/optimistic_lock.hpp"

// local sources
#include "b_tree/component/common.hpp"
#include "b_tree/component/metadata.hpp"

namespace dbgroup::index::b_tree::component::oml
{

/**
 * @brief A class for representing nodes in B+trees.
 *
 * This class uses optimistic multi-layer locking for concurrency controls and
 * optimizes a page layout for fixed-length data.
 *
 * @tparam Key a target key class.
 * @tparam Comp a comparetor class for keys.
 */
template <class Key, class Comp>
class NodeFixLen
{
 public:
  /*####################################################################################
   * Type aliases
   *##################################################################################*/

  using Node = NodeFixLen;

  // aliases for bulkloading
  template <class Entry>
  using BulkIter = typename std::vector<Entry>::const_iterator;

  /*####################################################################################
   * Public constructors and assignment operators
   *##################################################################################*/

  /**
   * @brief Construct an empty node object.
   *
   * @param is_leaf a flag to indicate whether a leaf node is constructed.
   */
  constexpr explicit NodeFixLen(const uint32_t is_leaf)
      : is_leaf_{is_leaf}, is_removed_{0}, block_size_{0}, has_low_key_{0}, has_high_key_{0}
  {
  }

  /**
   * @brief Construct a new root node.
   *
   * @param l_node a left child node which is the previous root node.
   * @param r_node a right child node.
   */
  NodeFixLen(  //
      const Key &l_key,
      [[maybe_unused]] const size_t l_key_len,
      const NodeFixLen *l_node,
      const NodeFixLen *r_node)  //
      : is_leaf_{0},
        is_removed_{0},
        block_size_{2 * kPtrLen},
        record_count_{2},
        has_low_key_{0},
        has_high_key_{0}
  {
    keys_[0] = l_key;
    SetPayload(kPageSize, &l_node);
    SetPayload(kPageSize - kPtrLen, &r_node);
  }

  NodeFixLen(const NodeFixLen &) = delete;
  NodeFixLen(NodeFixLen &&) = delete;

  auto operator=(const NodeFixLen &) -> NodeFixLen & = delete;
  auto operator=(NodeFixLen &&) -> NodeFixLen & = delete;

  /*####################################################################################
   * Public destructors
   *##################################################################################*/

  /**
   * @brief Destroy the node object.
   *
   */
  ~NodeFixLen() = default;

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

  [[nodiscard]] constexpr auto
  IsRemoved() const  //
      -> bool
  {
    return is_removed_;
  }

  [[nodiscard]] auto
  NeedRetryInner(const Key key,             //
                 const uint64_t ver) const  //
      -> NodeRC
  {
    const auto &high_key = GetHighKey();
    if (!mutex_.HasSameVersion(ver)) {
      return kNeedRetry;
    } else if (is_removed_ || (high_key && Comp{}(*high_key, key))) {
      return kNeedRootRetry;
    }
    return kCompleted;
  }

  [[nodiscard]] auto
  NeedRetryLeaf(const Key key,             //
                const uint64_t ver) const  //
      -> NodeRC
  {
    const auto &high_key = GetHighKey();
    if (!mutex_.HasSameVersion(ver)) {
      return kNeedRetry;
    } else if (is_removed_ || (high_key && Comp{}(*high_key, key))) {
      return kNeedNextRetry;
    }
    return kCompleted;
  }

  /**
   * @param new_rec_len the length of a new record.
   * @retval true if this node requires splitting before inserting a new record.
   * @retval false otherwise.
   */
  [[nodiscard]] auto
  NeedSplit(const size_t new_rec_len)  //
      -> std::pair<bool, uint64_t>
  {
    while (true) {
      const auto ver = mutex_.GetVersion();
      auto need_split = false;

      // check whether the node has space for a new record
      const auto total_size = kMetaLen * (record_count_ + 1) + block_size_ + new_rec_len;
      if (total_size <= kPageSize - kHeaderLen) {
        need_split = false;
      } else if (total_size - deleted_size_ > kMaxUsedSpaceSize) {
        need_split = true;
      } else {
        // this node has enough space but cleaning up is required
        if (mutex_.TryLockX(ver)) {
          CleanUp();
          mutex_.UnlockX();
          continue;
        }
        need_split = false;
      }
      if (mutex_.HasSameVersion(ver)) return {need_split, ver};
    }
  }

  /**
   * @retval true if this node requires merging.
   * @retval false otherwise.
   */
  [[nodiscard]] auto
  NeedMerge()  //
      -> std::pair<bool, uint64_t>
  {
    while (true) {
      const auto ver = mutex_.GetVersion();
      auto need_merge = false;

      // check this node uses enough space
      if (GetUsedSize() < kMinUsedSpaceSize) {
        need_merge = false;
      } else {
        // this node has a lot of dead space
        if (mutex_.TryLockX(ver)) {
          CleanUp();
          mutex_.UnlockX();
        }
        need_merge = false;
      }
      if (mutex_.HasSameVersion(ver)) return {need_merge, ver};
    }
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
      -> Node *
  {
    auto *node = this;
    const auto &high_key = GetHighKey();
    if (node->is_removed_) {
      node = next_;
    } else if (!high_key_ || Comp{}(GetHighKey(), key)) {
      node = next_;
    }

    return node;
  }

  /**
   * @brief Set the length of payloads for leaf read/write operations.
   *
   * @param pay_len the length of payloads.
   */
  void
  SetPayloadLength(const size_t pay_len)
  {
    pay_len_ = pay_len;
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
    return keys_[pos];
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
    memcpy(&payload, GetPayloadAddr(pos), sizeof(Payload));
    return payload;
  }

  /**
   * @brief Get a child node in a given position.
   *
   * The returned child node is locked with a shared lock and this node is unlocked.
   *
   * @param pos the position of a child node.
   * @return the child node with a shared lock.
   */
  [[nodiscard]] constexpr auto
  GetChildForRead(const size_t pos)  //
      -> Node *
  {
    auto *child = GetPayload<Node *>(pos);
    child->mutex_.LockS();
    mutex_.UnlockS();
    return child;
  }

  /**
   * @brief Get a child node in a given position.
   *
   * @param pos the position of a child node.
   * @return the child node with an SIX lock.
   */
  [[nodiscard]] constexpr auto
  GetChildForWrite(const size_t pos)  //
      -> Node *
  {
    auto *child = GetPayload<Node *>(pos);
    child->mutex_.LockSIX();
    return child;
  }

  /**
   * @brief Get a mergeable child node if exist.
   *
   * The returned child node is locked with an SIX lock. If there is no mergeable child
   * node, this node is unlocked.
   *
   * @param l_node a left child node to be merged.
   * @param l_pos the position of the left child node.
   * @retval a mergeable child node if exist.
   * @retval nullptr otherwise.
   */
  [[nodiscard]] auto
  GetMergeableRightChild(  //
      const Node *l_node,
      const size_t l_pos)  //
      -> Node *
  {
    // check there is a right-sibling node
    const auto r_pos = l_pos + 1;
    if (r_pos == record_count_) {
      // a rightmost node is cannot be merged
      mutex_.UnlockSIX();
      return nullptr;
    }

    // check the right-sibling node has enough capacity for merging
    auto *r_node = GetChildForWrite(r_pos);
    if (!l_node->CanMerge(r_node)) {
      mutex_.UnlockSIX();
      r_node->mutex_.UnlockSIX();
      return nullptr;
    }

    UpgradeToX();
    r_node->UpgradeToX();
    return r_node;
  }

  /**
   * @param pos the position of a record.
   * @return false.
   */
  [[nodiscard]] constexpr auto
  RecordIsDeleted([[maybe_unused]] const size_t pos) const  //
      -> bool
  {
    return false;
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

  auto
  LockX()  //
      -> void
  {
    mutex_.LockX();
  }

  auto
  TryLockX(const uint64_t ver)  //
      -> bool
  {
    return mutex_.TryLockX(ver);
  }

  void
  DowngradeToSIX()
  {
    mutex_.DowngradeToSIX();
  }

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

  auto
  TryLockSIX(const uint64_t ver)  //
      -> bool
  {
    return mutex_.TryLockSIX(ver);
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

  /**
   * @brief Release the shared lock with intent-exclusive locking for this node.
   *
   */
  void
  UnlockSIX()
  {
    mutex_.UnlockSIX();
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
    int64_t begin_pos = 0;
    int64_t end_pos = record_count_ - 1;
    while (begin_pos <= end_pos) {
      size_t pos = (begin_pos + end_pos) >> 1UL;  // NOLINT
      const auto &index_key = keys_[pos];

      if (Comp{}(key, index_key)) {  // a target key is in a left side
        end_pos = pos - 1;
      } else if (Comp{}(index_key, key)) {  // a target key is in a right side
        begin_pos = pos + 1;
      } else {  // find an equivalent key
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
      const auto &index_key = keys_[pos];

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
      memcpy(&out_payload, GetPayloadAddr(pos), sizeof(Payload));
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
  void
  Write(  //
      const Key &key,
      [[maybe_unused]] const size_t key_len,
      const void *payload,
      [[maybe_unused]] const size_t pay_len)
  {
    // search position where this key has to be set
    const auto [rc, pos] = SearchRecord(key);

    // perform insert or update operation
    if (rc == kKeyNotInserted) {
      InsertRecord(key, payload, pos);
    } else {  // update operation
      memcpy(GetPayloadAddr(pos), payload, pay_len_);
    }

    mutex_.UnlockX();
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
      [[maybe_unused]] const size_t key_len,
      const void *payload,
      [[maybe_unused]] const size_t pay_len)  //
      -> ReturnCode
  {
    // search position where this key has to be set
    const auto [existence, pos] = SearchRecord(key);

    // perform insert operation if possible
    auto rc = kKeyExist;
    if (existence == kKeyNotInserted) {
      InsertRecord(key, payload, pos);
      rc = kSuccess;
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
      [[maybe_unused]] const size_t pay_len)  //
      -> ReturnCode
  {
    // check this node has a target record
    const auto [existence, pos] = SearchRecord(key);

    // perform update operation if possible
    auto rc = kKeyNotExist;
    if (existence == kKeyAlreadyInserted) {
      memcpy(GetPayloadAddr(pos), payload, pay_len_);
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
      -> ReturnCode
  {
    // check this node has a target record
    const auto [existence, pos] = SearchRecord(key);

    auto rc = kKeyNotExist;
    if (existence == kKeyAlreadyInserted) {
      // remove a key and a payload
      const auto move_num = record_count_ - 1 - pos;
      memmove(&(keys_[pos]), &(keys_[pos + 1]), kKeyLen * (move_num + has_high_key_));
      auto *top_addr = ShiftAddr(this, kPageSize - block_size_);
      memmove(ShiftAddr(top_addr, pay_len_), top_addr, pay_len_ * move_num);

      // update header information
      --record_count_;
      block_size_ -= pay_len_;

      rc = kSuccess;
    }

    mutex_.UnlockX();
    return rc;
  }

  /**
   * @brief Insert a new child node into this node.
   *
   * @param l_node a left child node.
   * @param r_node a right child (i.e., new) node.
   * @param pos the position of the left child node.
   */
  void
  InsertChild(  //
      const Node *l_node,
      const Node *r_node,
      const size_t pos)  //
  {
    // insert a right child
    const auto move_num = record_count_ - 1 - pos;
    const auto move_size = kPtrLen * move_num;
    const auto top_offset = kPageSize - block_size_;
    memmove(ShiftAddr(this, top_offset - kPtrLen), ShiftAddr(this, top_offset), move_size);
    SetPayload(top_offset + move_size, &r_node);

    // insert a separator key
    memmove(&(keys_[pos + 1]), &(keys_[pos]), kKeyLen * (move_num + has_high_key_));
    keys_[pos] = l_node->GetHighKey();

    // update header information
    ++record_count_;
    block_size_ += kPtrLen;

    mutex_.UnlockX();
  }

  /**
   * @brief Delete a child node from this node.
   *
   * @param l_node a left child node.
   * @param pos the position of the left child node.
   */
  void
  DeleteChild(  //
      [[maybe_unused]] const Node *l_node,
      const size_t pos)  //
  {
    // delete a separator key
    const auto move_num = record_count_ - 2 - pos;
    memmove(&(keys_[pos]), &(keys_[pos + 1]), kKeyLen * (move_num + has_high_key_));

    // delete a right child
    auto *top_addr = ShiftAddr(this, kPageSize - block_size_);
    memmove(ShiftAddr(top_addr, kPtrLen), top_addr, kPtrLen * move_num);

    // update this header
    --record_count_;
    block_size_ -= kPtrLen;

    mutex_.UnlockX();
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
    const auto l_count = record_count_ / 2;
    const auto r_count = record_count_ - l_count;
    const auto is_inner = static_cast<size_t>(!static_cast<bool>(is_leaf_));

    // copy right half records to a right node
    r_node->pay_len_ = pay_len_;
    auto r_offset = r_node->CopyRecordsFrom(this, l_count, record_count_, kPageSize);
    r_node->keys_[r_count - is_inner] = keys_[record_count_ - is_inner];

    // update a right header
    r_node->block_size_ = kPageSize - r_offset;
    r_node->next_ = next_;
    r_node->has_high_key_ = has_high_key_;

    mutex_.UpgradeToX();  // upgrade the lock to modify the left node

    // update a header
    block_size_ -= r_node->block_size_;
    record_count_ = l_count;
    next_ = r_node;
    has_high_key_ = 1;
    keys_[l_count - is_inner] = keys_[l_count - 1];

    mutex_.DowngradeToSIX();
  }

  /**
   * @brief Merge a given node into this node.
   *
   * @param r_node a right node to be merged.
   */
  void
  Merge(const Node *r_node)
  {
    mutex_.UpgradeToX();

    // copy right records to this nodes
    auto offset = CopyRecordsFrom(r_node, 0, r_node->record_count_, kPageSize - block_size_);

    // update a header
    block_size_ = kPageSize - offset;
    next_ = r_node->next_;
    has_high_key_ = r_node->has_high_key_;

    mutex_.DowngradeToSIX();
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
      BulkIter<Entry> &iter,
      const BulkIter<Entry> &iter_end,
      const bool is_rightmost,
      Node *l_node)
  {
    constexpr auto kPayLen = sizeof(Payload);
    constexpr auto kRecLen = kKeyLen + kPayLen;
    constexpr auto kNodeCap = (kPageSize - kHeaderLen - kKeyLen - kMinFreeSpaceSize) / kRecLen;

    // set header information for treating payloads
    pay_len_ = kPayLen;

    // extract and insert entries for the leaf node
    const size_t n = std::distance(iter, iter_end);
    const auto has_last_record = n < kNodeCap;
    const auto rec_count = (has_last_record) ? n : kNodeCap;
    auto offset = kPageSize;
    for (size_t i = 0; i < rec_count; ++i, ++iter) {
      // insert an entry to the leaf node
      const auto &[key, payload] = *iter;
      keys_[record_count_++] = key;
      offset = SetPayload(offset, &payload);
    }

    // set a highest key if this node is not rightmost
    if (!has_last_record || !is_rightmost) {
      has_high_key_ = 1;
      keys_[rec_count] = keys_[rec_count - 1];
    }

    // update header information
    block_size_ = kPageSize - offset;
    if (l_node != nullptr) {
      l_node->next_ = this;
    }
  }

  /**
   * @brief Create a inner node with the maximum number of records for bulkloading.
   *
   * @param iter the begin position of target records.
   * @param iter_end the end position of target records.
   */
  void
  Bulkload(  //
      BulkIter<Node *> &iter,
      const BulkIter<Node *> &iter_end)
  {
    constexpr auto kRecLen = kKeyLen + kPtrLen;
    constexpr auto kNodeCap = (kPageSize - kHeaderLen - kKeyLen - kMinFreeSpaceSize) / kRecLen;

    // check the number of child nodes to be loaded
    const size_t n = std::distance(iter, iter_end);
    const auto rec_count = (n < kNodeCap) ? n : kNodeCap;
    const auto *last_child = *std::next(iter, rec_count - 1);
    has_high_key_ = last_child->has_high_key_;

    // extract and insert entries for the leaf node
    auto offset = kPageSize;
    for (size_t i = 0; i < rec_count; ++i, ++iter) {
      // insert an entry to the inner node
      const auto *child = *iter;
      keys_[record_count_++] = child->GetHighKey();
      offset = SetPayload(offset, &child);
    }

    // update header information
    block_size_ = kPageSize - offset;
  }

  /**
   * @brief Link border nodes between partial trees.
   *
   * @param l_node a highest border node in a left tree.
   * @param r_node a highest border node in a right tree.
   */
  static void
  LinkVerticalBorderNodes(  //
      Node *l_node,
      Node *r_node)
  {
    while (true) {
      if (l_node->is_leaf_) {
        l_node->next_ = r_node;
        return;
      }

      // go down to the lower level
      l_node = l_node->template GetPayload<Node *>(l_node->record_count_ - 1);
      r_node = r_node->template GetPayload<Node *>(0);
    }
  }

 private:
  /*####################################################################################
   * Internal constants
   *##################################################################################*/

  /// the length of keys.
  static constexpr size_t kKeyLen = sizeof(Key);

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

  /**
   * @return Current usage of the record block.
   */
  [[nodiscard]] auto
  GetUsedSize() const  //
      -> size_t
  {
    return kKeyLen * (record_count_ + (is_leaf_ & has_high_key_)) + block_size_;
  }

  /**
   * @retval a highest key in this node if exist.
   * @retval std::nullopt otherwise.
   */
  [[nodiscard]] auto
  GetHighKey() const  //
      -> const Key &
  {
    const auto is_inner = static_cast<size_t>(!static_cast<bool>(is_leaf_));
    return keys_[record_count_ - is_inner];
  }

  /*####################################################################################
   * Internal getter/setters for records
   *##################################################################################*/

  /**
   * @param pos the position of a target record.
   * @return an address of a target payload.
   */
  [[nodiscard]] constexpr auto
  GetPayloadAddr(const size_t pos) const  //
      -> void *
  {
    return ShiftAddr(this, kPageSize - (pos + 1) * pay_len_);
  }

  /**
   * @brief Set a target payload directly.
   *
   * @param offset an offset to the top of the record block.
   * @param payload a target payload to be written.
   */
  auto
  SetPayload(  //
      size_t offset,
      const void *payload)  //
      -> size_t
  {
    offset -= pay_len_;
    memcpy(ShiftAddr(this, offset), payload, pay_len_);
    return offset;
  }

  /*####################################################################################
   * Internal utility functions
   *##################################################################################*/

  /**
   * @brief Insert a given record into this node.
   *
   * @param key a target key to be set.
   * @param payload a target payload to be written.
   * @param pos an insertion position.
   */
  void
  InsertRecord(  //
      const Key &key,
      const void *payload,
      const size_t pos)  //
      ->NodeRc
  {
    const auto rec_len = kKeyLen + pay_len_;
    if (GetUsedSize() + rec_len >) const auto move_num = record_count_ - pos;

    // insert a new key
    memmove(&(keys_[pos + 1]), &(keys_[pos]), kKeyLen * (move_num + has_high_key_));
    keys_[pos] = key;

    // insert a new payload
    const auto top_offset = kPageSize - block_size_;
    const auto move_size = pay_len_ * move_num;
    memmove(ShiftAddr(this, top_offset - pay_len_), ShiftAddr(this, top_offset), move_size);
    SetPayload(top_offset + move_size, payload);

    // update header information
    ++record_count_;
    block_size_ += pay_len_;
  }

  /**
   * @brief Clean up this node.
   *
   */
  void
  CleanUp()
  {
    // copy records to a temporal node
    auto offset = temp_node_->CopyRecordsFrom(this, 0, record_count_, offset);

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
   * @brief Copy a record from a given node.
   *
   * @param node an original node that has a target record.
   * @param pos the position of a target record.
   * @param offset an offset to the top of the record block.
   * @return the updated offset value.
   */
  auto
  CopyRecordFrom(  //
      const Node *node,
      const size_t pos,
      size_t offset)  //
      -> size_t
  {
    // copy a record from the given node
    keys_[record_count_++] = node->keys_[pos];
    offset -= pay_len_;
    memcpy(ShiftAddr(this, offset), node->GetPayloadAddr(pos), pay_len_);

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
      offset = CopyRecordFrom(node, i, offset);
    }

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

  /// the total byte length of records in a node.
  uint16_t deleted_size_{0};

  /// the number of records in this node.
  uint16_t record_count_{0};

  /// the length of payloads in this node.
  uint16_t pay_len_{kPtrLen};

  /// a lock for concurrency controls.
  ::dbgroup::lock::OptimisticLock mutex_{};

  /// the pointer to the next node.
  Node *next_{nullptr};

  /// a flag for a lowest key.
  uint64_t has_low_key_ : 1;

  /// a flag for a highest key.
  uint64_t has_high_key_ : 1;

  /// an alignment block.
  uint64_t : 0;

  /// an actual data block (it starts with record keys).
  Key keys_[0];

  /// the metadata of a highest key.
  Metadata high_meta_{kPageSize, 0, 0};

  /// a temporary node for SMOs.
  static thread_local inline std::unique_ptr<Node> temp_node_ =  // NOLINT
      std::make_unique<Node>(0);
};

}  // namespace dbgroup::index::b_tree::component::oml

#endif  // B_TREE_COMPONENT_OML_NODE_FIXLEN_HPP
