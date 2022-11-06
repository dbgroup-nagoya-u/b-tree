
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

#ifndef B_TREE_COMPONENT_PML_NODE_VARLEN_HPP
#define B_TREE_COMPONENT_PML_NODE_VARLEN_HPP

#include <atomic>
#include <optional>
#include <utility>
#include <vector>

// organization libraries
#include "lock/pessimistic_lock.hpp"

// local sources
#include "b_tree/component/common.hpp"
#include "b_tree/component/metadata.hpp"

namespace dbgroup::index::b_tree::component::pml
{

/**
 * @brief A class for representing nodes in B+trees.
 *
 * This class uses pessimistic multi-layer locking for concurrency controls and can
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
  using ScanKey = std::optional<std::tuple<const Key &, size_t, bool>>;
  // aliases for bulkloading
  template <class Entry>
  using BulkIter = typename std::vector<Entry>::const_iterator;
  using NodeEntry = std::tuple<Key, Node *, size_t>;

  /*####################################################################################
   * Public constructors and assignment operators
   *##################################################################################*/

  /**
   * @brief Construct an empty node object.
   *
   * @param is_inner a flag to indicate whether a leaf node is constructed.
   */
  constexpr explicit NodeVarLen(const uint32_t is_inner) : is_inner_{is_inner}, block_size_{0} {}

  /**
   * @brief Construct a new root node.
   *
   * @param l_node a left child node which is the previous root node.
   * @param r_node a right child node.
   */
  NodeVarLen(  //
      const NodeVarLen *l_node,
      const NodeVarLen *r_node)  //
      : is_inner_{1}, record_count_{2}
  {
    // insert l_node
    const auto l_high_meta = l_node->high_meta_;
    const auto l_key_len = l_high_meta.key_len;
    auto offset = SetPayload(kPageSize, &l_node, kPtrLen);
    meta_array_[0] = Metadata{offset, 0, kPtrLen};

    // insert r_node
    offset = SetPayload(offset, &r_node, kPtrLen);
    offset = CopyKeyFrom(l_node, l_high_meta, offset);
    meta_array_[1] = Metadata{offset, l_key_len, l_key_len + kPtrLen};

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
   * @return true if this is a inner node.
   * @return false otherwise.
   */
  [[nodiscard]] constexpr auto
  IsInner() const  //
      -> bool
  {
    return is_inner_;
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
    mutex_.UpgradeToX();
    CleanUp();
    mutex_.DowngradeToSIX();
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
    mutex_.UpgradeToX();
    CleanUp();
    mutex_.DowngradeToSIX();
    return false;
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
    auto *next = next_;
    next_->mutex_.LockS();
    mutex_.UnlockS();
    return next;
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
  GetValidSplitNode(  //
      const Key &key,
      Node *r_node)  //
      -> Node *
  {
    auto *node = this;
    const auto &high_key = GetHighKey();
    if (!high_key || !Comp{}(key, *high_key)) {
      node = r_node;
      mutex_.UnlockSIX();
    } else {
      r_node->mutex_.UnlockSIX();
    }

    return node;
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
    if (l_node->GetUsedSize() + r_node->GetUsedSize() >= kMaxUsedSpaceSize) {
      mutex_.UnlockSIX();
      r_node->mutex_.UnlockSIX();
      return nullptr;
    }

    UpgradeToX();
    return r_node;
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

  /**
   * @brief Downgrade the X lock to an SIX lock for this node.
   *
   */
  void
  DowngradeToSIX()
  {
    mutex_.DowngradeToSIX();
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
    int64_t begin_pos = is_inner_;
    int64_t end_pos = record_count_ - 1;
    while (begin_pos <= end_pos) {
      size_t pos = (begin_pos + end_pos) >> 1UL;  // NOLINT
      const auto &index_key = GetKey(meta_array_[pos]);

      if (Comp{}(key, index_key)) {  // a target key is in a left side
        end_pos = pos - 1;
      } else if (Comp{}(index_key, key)) {  // a target key is in a right side
        begin_pos = pos + 1;
      } else if (meta_array_[pos].is_deleted) {
        return {kKeyAlreadyDeleted, pos};
      } else {
        return {kKeyAlreadyInserted, pos};
      }
    }

    return {kKeyNotInserted, begin_pos - is_inner_};
  }

  /**
   * @brief Get the end position of records for scanning and check it has been finished.
   *
   * @param end_key a pair of a target key and its closed/open-interval flag.
   * @retval 1st: true if this node is end of scanning.
   * @retval 2nd: the end position for scanning.
   */
  [[nodiscard]] auto
  SearchEndPositionFor(const ScanKey &end_key)  //
      -> std::pair<bool, size_t>
  {
    const auto is_end = IsRightmostOf(end_key);
    size_t end_pos{};
    if (is_end && end_key) {
      const auto &[e_key, e_key_len, e_closed] = *end_key;
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
      const size_t key_len,
      const void *payload,
      const size_t pay_len)
  {
    // search position where this key has to be set
    const auto [rc, pos] = SearchRecord(key);

    // perform insert or update operation
    if (rc == kKeyNotInserted) {
      InsertRecord(key, key_len, payload, pay_len, pos);
    } else if (rc == kKeyAlreadyDeleted) {
      ReuseRecord(payload, pay_len, pos);
    } else {  // update operation
      memcpy(GetPayloadAddr(meta_array_[pos]), payload, pay_len);
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
      const size_t key_len,
      const void *payload,
      const size_t pay_len)  //
      -> ReturnCode
  {
    // search position where this key has to be set
    const auto [existence, pos] = SearchRecord(key);

    // perform insert operation if possible
    auto rc = kSuccess;
    if (existence == kKeyNotInserted) {
      InsertRecord(key, key_len, payload, pay_len, pos);
    } else if (existence == kKeyAlreadyDeleted) {
      ReuseRecord(payload, pay_len, pos);
    } else {  // a target key has been inserted
      rc = kKeyExist;
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
    // prepare space for a right child
    memmove(&(meta_array_[pos + 2]), &(meta_array_[pos + 1]), kMetaLen * (record_count_ - pos - 1));

    // insert a right child
    const auto l_high_meta = l_node->high_meta_;
    const auto key_len = l_high_meta.key_len;
    const auto rec_len = key_len + kPtrLen;
    auto offset = SetPayload(kPageSize - block_size_, &r_node, kPtrLen);
    offset = CopyKeyFrom(l_node, l_high_meta, offset);
    meta_array_[pos + 1] = Metadata{offset, key_len, rec_len};

    // update header information
    ++record_count_;
    block_size_ += rec_len;

    mutex_.UnlockX();
  }

  /**
   * @brief Delete a child node from this node.
   *
   * @param pos the position of the left child node.
   */
  void
  DeleteChild(const size_t pos)  //
  {
    const auto del_rec_len = meta_array_[pos + 1].rec_len;  // keep a record length to be deleted

    // delete a right child by shifting metadata
    memmove(&(meta_array_[pos + 1]), &(meta_array_[pos + 2]), kMetaLen * (record_count_ - 2 - pos));

    // update this header
    --record_count_;
    deleted_size_ += del_rec_len;

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
    offset = temp_node_->CopyKeyFrom(this, meta_array_[pos], offset);
    const auto sep_key_len = meta_array_[pos].key_len;
    temp_node_->high_meta_ = Metadata{offset, sep_key_len, sep_key_len};

    // copy right half records to a right node
    auto r_offset = r_node->CopyHighKeyFrom(this);
    r_offset = r_node->CopyRecordsFrom(this, pos, record_count_, r_offset);

    // update a right header
    r_node->block_size_ = kPageSize - r_offset;
    if (!is_inner_) {
      r_node->next_ = next_;
    }

    mutex_.UpgradeToX();  // upgrade the lock to modify the left node

    // update a header
    block_size_ = kPageSize - offset;
    deleted_size_ = 0;
    record_count_ = temp_node_->record_count_;
    if (!is_inner_) {
      next_ = r_node;
    }

    // copy temporal node to this node
    memcpy(&high_meta_, &(temp_node_->high_meta_), kMetaLen * (record_count_ + 1));
    memcpy(ShiftAddr(this, offset), ShiftAddr(temp_node_.get(), offset), block_size_);

    // reset a temp node
    temp_node_->record_count_ = 0;

    mutex_.DowngradeToSIX();
    r_node->mutex_.LockSIX();
  }

  /**
   * @brief Merge a given node into this node.
   *
   * @param r_node a right node to be merged.
   */
  void
  Merge(Node *r_node)
  {
    mutex_.UpgradeToX();
    r_node->UpgradeToX();

    // copy a highest key of a merged node to a temporal node
    auto offset = temp_node_->CopyHighKeyFrom(r_node);

    // copy consolidated records to the original node
    offset = temp_node_->CopyRecordsFrom(this, 0, record_count_, offset);
    record_count_ = temp_node_->record_count_;
    memcpy(&high_meta_, &(temp_node_->high_meta_), kMetaLen * (record_count_ + 1));
    memcpy(ShiftAddr(this, offset), ShiftAddr(temp_node_.get(), offset), kPageSize - offset);

    // copy right records to this nodes
    offset = CopyRecordsFrom(r_node, 0, r_node->record_count_, offset);

    // update a header
    block_size_ = kPageSize - offset;
    deleted_size_ = 0;
    if (!is_inner_) {
      next_ = r_node->next_;
    }

    // reset a temp node
    temp_node_->record_count_ = 0;

    mutex_.DowngradeToSIX();
  }

  /*####################################################################################
   * Public bulkload APIs
   *##################################################################################*/
  /**
   * @brief Create a node with the maximum number of records for bulkloading.
   *
   * @tparam Entry a container of a key/payload pair.
   * @param iter the begin position of target records.
   * @param iter_end the end position of target records.
   * @param prev_node a left sibling node.
   * @param nodes the container of construcred nodes.
   */
  template <class Entry>
  void
  Bulkload(  //
      BulkIter<Entry> &iter,
      const BulkIter<Entry> &iter_end,
      Node *prev_node,
      std::vector<NodeEntry> &nodes)
  {
    using Payload = std::tuple_element_t<1, Entry>;

    constexpr auto kMaxKeyLen = (IsVarLenData<Key>()) ? kMaxVarLenDataSize : sizeof(Key);
    const auto kPayLen = sizeof(Payload);

    // extract and insert entries into this node
    auto offset = kPageSize - kMaxKeyLen;  // reserve the space for a highest key
    auto node_size = kHeaderLen + kMaxKeyLen;
    for (; iter < iter_end; ++iter) {
      const auto &[key, payload, key_len] = ParseEntry(*iter);
      const auto rec_len = key_len + kPayLen;

      // check whether the node has sufficent space
      node_size += rec_len + sizeof(Metadata);
      if (node_size + kMaxKeyLen > kPageSize) break;

      // insert an entry into this node
      offset = SetPayload(offset, &payload, kPayLen);
      offset = SetKey(offset, key, key_len);
      meta_array_[record_count_++] = Metadata{offset, key_len, rec_len};
    }

    block_size_ = kPageSize - offset;

    // link the sibling nodes if exist
    if (prev_node != nullptr) {
      prev_node->LinkNext(this);
    }

    nodes.emplace_back(GetKey(0), this, meta_array_[0].key_len);
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
      l_node->LinkNext(r_node);
      if (!l_node->is_inner_) {
        return;
      }

      // go down to the lower level
      l_node = l_node->template GetPayload<Node *>(l_node->record_count_ - 1);
      r_node = r_node->template GetPayload<Node *>(0);
    }
  }

  /**
   * @brief Remove the leftmost keys from the leftmost nodes.
   *
   */
  static void
  RemoveLeftmostKeys(Node *node)
  {
    while (!node->IsInner()) {
      // remove the leftmost key in a record region of an inner node
      const auto meta = node->meta_array_[0];
      const auto offset = meta.offset;
      const auto key_len = meta.key_len;
      const size_t rec_len = meta.rec_len - key_len;
      node->meta_array_[0] = Metadata{offset + key_len, 0, rec_len};

      // go down to the lower level
      node = node->template GetPayload<Node *>(0);
    }
  }

 private:
  /*####################################################################################
   * Internal constants
   *##################################################################################*/

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
  IsRightmostOf(const ScanKey &end_key) const  //
      -> bool
  {
    if (!next_) return true;     // the rightmost node
    if (!end_key) return false;  // perform full scan
    const auto &high_key = GetKey(high_meta_);
    return Comp{}(std::get<0>(*end_key), high_key);
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
   * @retval a highest key in this node if exist.
   * @retval std::nullopt otherwise.
   */
  [[nodiscard]] auto
  GetHighKey() const  //
      -> std::optional<Key>
  {
    if (high_meta_.key_len == 0) return std::nullopt;
    return GetKey(high_meta_);
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
  void
  InsertRecord(  //
      const Key &key,
      const size_t key_len,
      const void *payload,
      const size_t pay_len,
      const size_t pos)
  {
    const auto rec_len = key_len + pay_len;

    // insert a new record
    auto offset = kPageSize - block_size_;
    offset = SetPayload(offset, payload, pay_len);
    offset = SetKey(offset, key, key_len);
    memmove(&(meta_array_[pos + 1]), &(meta_array_[pos]), kMetaLen * (record_count_ - pos));
    meta_array_[pos] = Metadata{offset, key_len, rec_len};

    // update header information
    ++record_count_;
    block_size_ += rec_len;
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
   * @brief Copy a high key from a given node.
   *
   * @param node an original node that has a high key.
   * @return the updated offset value.
   */
  auto
  CopyHighKeyFrom(const Node *node)  //
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
  template <class Entry>
  constexpr auto
  ParseEntry(const Entry &entry)  //
      -> std::tuple<Key, std::tuple_element_t<1, Entry>, size_t>
  {
    constexpr auto kTupleSize = std::tuple_size_v<Entry>;
    static_assert(2 <= kTupleSize && kTupleSize <= 3);

    if constexpr (kTupleSize == 3) {
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
    // set a highest key in a left node
    const auto high_key_len = r_node->meta_array_[0].key_len;
    auto offset = kPageSize - block_size_;
    offset = SetKey(offset, r_node->GetKey(0), high_key_len);
    block_size_ += high_key_len;
    high_meta_ = Metadata{offset, high_key_len, high_key_len};

    // set a sibling link in a left node
    if (IsInner()) return offset;

    next_ = r_node;

    return offset;
  }

  /*####################################################################################
   * Internal member variables
   *##################################################################################*/

  /// a flag for indicating this node is a leaf or internal node.
  uint32_t is_inner_ : 1;

  /// the total byte length of records in a node.
  uint32_t block_size_ : 31;

  /// the total byte length of deleted records in a node.
  uint16_t deleted_size_{0};

  /// the number of records in this node.
  uint16_t record_count_{0};

  /// a lock for concurrency controls.
  ::dbgroup::lock::PessimisticLock mutex_{};

  /// the pointer to the next node.
  Node *next_{nullptr};

  /// the metadata of a highest key.
  Metadata high_meta_{kPageSize, 0, 0};

  /// an actual data block (it starts with record metadata).
  Metadata meta_array_[0];

  // a temporary node for SMOs.
  static thread_local inline std::unique_ptr<Node> temp_node_ =  // NOLINT
      std::make_unique<Node>(0);
};

}  // namespace dbgroup::index::b_tree::component::pml

#endif  // B_TREE_COMPONENT_PML_NODE_VARLEN_HPP
