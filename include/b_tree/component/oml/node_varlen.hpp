
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

#ifndef B_TREE_COMPONENT_OML_NODE_VARLEN_HPP
#define B_TREE_COMPONENT_OML_NODE_VARLEN_HPP

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
 * This class uses optimistic multi-layer locking for concurrency controls and can
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
      const NodeVarLen *l_node,
      const NodeVarLen *r_node)  //
      : is_leaf_{0}, is_removed_{0}, record_count_{2}
  {
    // insert l_node
    const auto l_high_meta = l_node->high_meta_;
    const auto l_key_len = l_high_meta.key_len;
    auto offset = SetPayload(kPageSize, &l_node, kPtrLen);
    offset = CopyKeyFrom(l_node, l_high_meta, offset);
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
      -> std::tuple<Node *, Key, size_t, uint64_t>
  {
    Node *node{};
    uint64_t ver{};
    const auto &[sep_key, sep_key_len] = GetHighKeyForSMOs();

    if (Comp{}(sep_key, key)) {
      node = next_;
      ver = next_->mutex_.UnlockX();
      mutex_.UnlockX();
    } else {
      node = this;
      next_->mutex_.UnlockX();
      ver = mutex_.UnlockX();
    }

    return {node, sep_key, sep_key_len, ver};
  }

  /**
   * @brief
   *
   * @param new_rec_len the length of a new record.
   * @retval 1st: kNeedSplit if this node requires splitting before modification.
   * @retval 1st: kNeedMerge if this node requires splitting before modification.
   * @retval 1st: kNeedRetry if this node is removed.
   * @retval 1st: kCompleted otherwise.
   * @retval 2nd: a version value for optimistic concurrency controls.
   */
  [[nodiscard]] auto
  CheckNodeStatus(const size_t new_rec_len)  //
      -> std::pair<NodeRC, uint64_t>
  {
    NodeRC rc{};
    uint64_t ver{};
    while (true) {
      ver = mutex_.GetVersion();

      // check this node is not removed
      if (is_removed_) {
        if (!mutex_.HasSameVersion(ver)) continue;
        rc = kNeedRetry;
        break;
      }

      // check a certain SMO is needed
      rc = CheckSMOs(new_rec_len, ver);
      if (rc != kNeedRetry) break;
    }

    return {rc, ver};
  }

  /**
   * @retval 1st: a highest key.
   * @retval 2nd: the length of the highest key.
   */
  [[nodiscard]] auto
  GetHighKeyForSMOs() const  //
      -> std::pair<Key, size_t>
  {
    const auto h_key_len = high_meta_.key_len;
    if constexpr (IsVarLenData<Key>()) {
      // allocate space dynamically to keep a copied key
      auto *h_key = reinterpret_cast<Key>(::operator new(h_key_len));
      memcpy(h_key, GetKeyAddr(high_meta_), h_key_len);
      return {h_key, h_key_len};
    } else {
      return {*GetHighKey(), h_key_len};
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
   * @brief Get a leftmost child node.
   *
   * @return the child node.
   */
  [[nodiscard]] auto
  GetLeftmostChild() const  //
      -> Node *
  {
    while (true) {
      const auto ver = mutex_.GetVersion();

      // check this node is not removed
      Node *child = nullptr;
      if (is_removed_ == 0) {
        child = GetPayload<Node *>(0);
      }

      if (mutex_.HasSameVersion(ver)) return child;
    }
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
      Node *l_node,
      const size_t l_pos)  //
      -> Node *
  {
    // check there is a right-sibling node
    const auto r_pos = l_pos + 1;
    if (r_pos == record_count_) {
      // a rightmost node is cannot be merged
      mutex_.UnlockSIX();
      l_node->mutex_.UnlockSIX();
      return nullptr;
    }

    // check the right-sibling node has enough capacity for merging
    auto *r_node = GetPayload<Node *>(r_pos);
    r_node->mutex_.LockSIX();
    if (l_node->GetUsedSize() + r_node->GetUsedSize() >= kMaxUsedSpaceSize) {
      mutex_.UnlockSIX();
      l_node->mutex_.UnlockSIX();
      r_node->mutex_.UnlockSIX();
      return nullptr;
    }

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
  void
  UnlockS()
  {
    mutex_.UnlockS();
  }

  /**
   * @brief Release the exclusive lock for this node.
   *
   * @return a version value after the exclusive lock is released.
   */
  auto
  UnlockX()  //
      -> uint64_t
  {
    return mutex_.UnlockX();
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
   * @param ver an expected version value.
   * @retval true if a shared lock with intent-exclusive locking is acquired.
   * @retval false otherwise.
   */
  auto
  TryLockSIX(const uint64_t ver)  //
      -> bool
  {
    return mutex_.TryLockSIX(ver);
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
   * @retval 1st: kKeyAlreadyInserted if a record is found.
   * @retval 1st: kKeyAlreadyDeleted if a found record is deleted.
   * @retval 1st: kKeyNotInserted otherwise.
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
      } else if (!meta_array_[pos].is_deleted) {
        return {kKeyAlreadyInserted, pos};
      } else {
        return {kKeyAlreadyDeleted, pos};
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
      -> Node *
  {
    Node *child{};
    while (true) {
      const auto ver = mutex_.GetVersion();

      // check the current node is not removed
      if (is_removed_) {
        if (!mutex_.HasSameVersion(ver)) continue;
        child = nullptr;  // retry from a root node
        break;
      }

      // check the current node has a target key
      const auto &high_key = GetHighKey();
      if (high_key && (Comp{}(*high_key, key) || (!is_closed && !Comp{}(key, *high_key)))) {
        if (!mutex_.HasSameVersion(ver)) continue;
        child = nullptr;  // retry from a root node
        break;
      }

      // search a child node
      auto [rc, pos] = SearchRecord(key);
      if (!is_closed && rc == kKeyAlreadyInserted) {
        ++pos;
      }
      child = GetPayload<Node *>(pos);

      if (mutex_.HasSameVersion(ver)) break;
    }

    return child;
  }

  /**
   * @brief Get a child node of a specified key by using binary search.
   *
   * If there is no specified key in this node, this returns the minimum position that
   * is greater than the specified key.
   *
   * @param key a search key.
   * @param ver an expected version value.
   * @param rec_len the length of a record to be inserted.
   * @retval 1st: the position of a child node.
   * @retval 2nd: nullptr if the key is out of range of this node.
   * @retval 2nd: the child node that includes the key otherwise.
   */
  [[nodiscard]] auto
  SearchChild(  //
      const Key &key,
      uint64_t &ver,
      const size_t rec_len)  //
      -> std::pair<size_t, Node *>
  {
    while (true) {
      while (mutex_.HasSameVersion(ver)) {
        // check the current node has a target key
        const auto &high_key = GetHighKey();
        if (high_key && Comp{}(*high_key, key)) {
          if (!mutex_.HasSameVersion(ver)) break;
          return {0, nullptr};  // retry from a root node
        }

        // search a child node
        const auto pos = SearchRecord(key).second;
        auto *child = GetPayload<Node *>(pos);

        if (!mutex_.HasSameVersion(ver)) break;
        return {pos, child};
      }

      // check a current node status
      const auto [rc, new_ver] = CheckNodeStatus(rec_len);
      if (rc != kCompleted) return {0, nullptr};  // retry from a root node

      // this node is modified, but it may includes a target key
      ver = new_ver;
    }
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
   * @param node a current node to be checked.
   * @param key a search key.
   * @param is_closed a flag for indicating closed-interval.
   * @return a version value for optimistic concurrency controls.
   */
  [[nodiscard]] static auto
  CheckKeyRange(  //
      Node *&node,
      const Key &key,
      const bool is_closed)  //
      -> uint64_t
  {
    uint64_t ver{};
    while (true) {
      ver = node->mutex_.GetVersion();
      const auto &high_key = node->GetHighKey();

      // check the node is not removed
      if (node->is_removed_ == 0) {
        // check the node includes a target key
        if (!high_key) {
          if (node->mutex_.HasSameVersion(ver)) break;
          continue;
        }
        if (Comp{}(key, *high_key) || (is_closed && !Comp{}(*high_key, key))) {
          if (node->mutex_.HasSameVersion(ver)) break;
          continue;
        }
      }

      // go to the next node
      auto *next = node->next_;
      if (!node->mutex_.HasSameVersion(ver)) continue;

      node = next;
      if (node == nullptr) break;
    }

    return ver;
  }

  /**
   * @brief Check the key range of a given node and traverse side links if needed.
   *
   * This function acquires a shared lock for a returned node.
   *
   * @param node a current node to be locked.
   * @param key a search key.
   * @param is_closed a flag for indicating closed-interval.
   */
  static void
  CheckKeyRangeAndLockForRead(  //
      Node *&node,
      const Key &key,
      const bool is_closed)
  {
    while (true) {
      auto ver = CheckKeyRange(node, key, is_closed);
      if (node->mutex_.TryLockS(ver)) return;
    }
  }

  /**
   * @brief Check the key range of a given node and traverse side links if needed.
   *
   * This function acquires a shared lock with intent-exclusive locking for the node.
   *
   * @param node a current node to be locked.
   * @param key a search key.
   * @param new_rec_len the length of a new record.
   * @retval kNeedRetry if splitting is required.
   * @retval kCompleted otherwise.
   */
  [[nodiscard]] static auto
  CheckKeyRangeAndLockForWrite(  //
      Node *&node,
      const Key &key,
      const size_t new_rec_len)  //
      -> NodeRC
  {
    while (true) {
      auto ver = CheckKeyRange(node, key, kClosed);
      const auto rc = node->CheckSMOs(new_rec_len, ver);
      if (rc == kNeedRetry) continue;
      if (rc == kNeedSplit) return kNeedRetry;  // ignore merging in the leaf traversal

      // if there is a space for a new record, acquire lock
      if (!node->TryLockSIX(ver)) continue;  // retry on the leaf level

      return kCompleted;
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
   * @retval kCompleted if a target record is read.
   * @retval kKeyNotInserted otherwise.
   */
  template <class Payload>
  static auto
  Read(Node *&node,
       const Key &key,
       Payload &out_payload)  //
      -> NodeRC
  {
    while (true) {
      const auto ver = CheckKeyRange(node, key, kClosed);

      const auto [existence, pos] = node->SearchRecord(key);
      if (existence == kKeyAlreadyInserted) {
        const auto meta = node->meta_array_[pos];
        memcpy(&out_payload, node->GetPayloadAddr(meta), sizeof(Payload));
        if (node->mutex_.HasSameVersion(ver)) return kCompleted;
      }

      if (node->mutex_.HasSameVersion(ver)) return kKeyNotInserted;
    }
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

    mutex_.UpgradeToX();

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
   * @retval kCompleted if a key/payload pair is written.
   * @retval kNeedRetry if a leaf SMO is required.
   * @retval kKeyAlreadyInserted otherwise.
   */
  static auto
  Insert(  //
      Node *&node,
      const Key &key,
      const size_t key_len,
      const void *payload,
      const size_t pay_len)  //
      -> NodeRC
  {
    const auto rec_len = key_len + pay_len;

    while (true) {
      auto ver = CheckKeyRange(node, key, kClosed);
      const auto rc = node->CheckSMOs(rec_len, ver);
      if (rc == kNeedRetry) continue;           // retry on the leaf level
      if (rc == kNeedSplit) return kNeedRetry;  // ignore merging in the leaf traversal

      // search position where this key has to be set
      const auto [existence, pos] = node->SearchRecord(key);
      if (existence == kKeyAlreadyInserted) {
        if (node->mutex_.HasSameVersion(ver)) return kKeyAlreadyInserted;
        continue;
      }

      // a target record does not exist, so try to acquire an exclusive lock
      if (!node->mutex_.TryLockX(ver)) continue;

      // perform insert operation if possible
      if (existence == kKeyNotInserted) {
        node->InsertRecord(key, key_len, payload, pay_len, pos);
      } else if (existence == kKeyAlreadyDeleted) {
        node->ReuseRecord(payload, pay_len, pos);
      }

      node->mutex_.UnlockX();
      return kCompleted;
    }
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
  static auto
  Update(  //
      Node *&node,
      const Key &key,
      const void *payload,
      const size_t pay_len)  //
      -> ReturnCode
  {
    while (true) {
      const auto ver = CheckKeyRange(node, key, kClosed);

      // check this node has a target record
      const auto [existence, pos] = node->SearchRecord(key);
      if (existence != kKeyAlreadyInserted) {
        if (node->mutex_.HasSameVersion(ver)) return kKeyNotExist;
        continue;
      }

      // a target record exists, so try to acquire an exclusive lock
      if (!node->mutex_.TryLockX(ver)) continue;

      memcpy(node->GetPayloadAddr(node->meta_array_[pos]), payload, pay_len);
      node->mutex_.UnlockX();
      return kSuccess;
    }
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
  static auto
  Delete(  //
      Node *&node,
      const Key &key)  //
      -> ReturnCode
  {
    while (true) {
      const auto ver = CheckKeyRange(node, key, kClosed);

      // check this node has a target record
      const auto [existence, pos] = node->SearchRecord(key);
      if (existence != kKeyAlreadyInserted) {
        if (node->mutex_.HasSameVersion(ver)) return kKeyNotExist;
        continue;
      }

      // a target record exists, so try to acquire an exclusive lock
      if (!node->mutex_.TryLockX(ver)) continue;

      node->meta_array_[pos].is_deleted = 1;
      node->deleted_size_ += node->meta_array_[pos].rec_len + kMetaLen;
      node->mutex_.UnlockX();
      return kSuccess;
    }
  }

  /**
   * @brief Insert a new child node into this node.
   *
   * @param l_node a left child node.
   * @param r_node a right child (i.e., new) node.
   * @param sep_key a separator key.
   * @param sep_key_len the length of the separator key.
   * @param pos the position of the left child node.
   */
  void
  InsertChild(  //
      const Node *l_node,
      const Node *r_node,
      const Key &sep_key,
      const size_t sep_key_len,
      const size_t pos)  //
  {
    mutex_.UpgradeToX();

    // update a current pointer and prepare space for a left child
    memcpy(GetPayloadAddr(meta_array_[pos]), &r_node, kPtrLen);
    memmove(&(meta_array_[pos + 1]), &(meta_array_[pos]), kMetaLen * (record_count_ - pos));

    // insert a left child
    const auto rec_len = sep_key_len + kPtrLen;
    auto offset = SetPayload(kPageSize - block_size_, &l_node, kPtrLen);
    offset = SetKey(offset, sep_key, sep_key_len);
    meta_array_[pos] = Metadata{offset, sep_key_len, rec_len};

    // update header information
    ++record_count_;
    block_size_ += rec_len;

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
      const Node *l_node,
      const size_t pos)  //
  {
    const auto del_rec_len = meta_array_[pos].rec_len;  // keep a record length to be deleted

    mutex_.UpgradeToX();

    // delete a right child by shifting metadata
    memmove(&(meta_array_[pos]), &(meta_array_[pos + 1]), kMetaLen * (record_count_ - 1 - pos));
    memcpy(GetPayloadAddr(meta_array_[pos]), &l_node, kPtrLen);

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
    const auto sep_key_len = meta_array_[pos - 1].key_len;
    temp_node_->high_meta_ = Metadata{offset, sep_key_len, sep_key_len};
    if (!is_leaf_) {
      const auto rightmost_pos = temp_node_->record_count_ - 1;
      temp_node_->meta_array_[rightmost_pos] = Metadata{offset + sep_key_len, 0, kPtrLen};
    }

    // prevent deadlock
    mutex_.UpgradeToX();

    // copy right half records to a right node
    r_node->mutex_.LockX();
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
   * @brief Merge a given node into this node.
   *
   * @param r_node a right node to be merged.
   */
  auto
  Merge(Node *r_node)  //
      -> uint64_t
  {
    // copy a highest key of a merged node to a temporal node
    auto offset = temp_node_->CopyHighKeyFrom(r_node);

    // copy consolidated records to the original node
    offset = temp_node_->CopyRecordsFrom(this, 0, record_count_, offset);
    const auto rec_count = temp_node_->record_count_;
    if (!is_leaf_) {
      offset = temp_node_->CopyKeyFrom(this, high_meta_, offset);
      const auto key_len = high_meta_.key_len;
      temp_node_->meta_array_[rec_count - 1] = Metadata{offset, key_len, key_len + kPtrLen};
    }

    mutex_.UpgradeToX();

    record_count_ = rec_count;
    memcpy(&high_meta_, &(temp_node_->high_meta_), kMetaLen * (record_count_ + 1));
    memcpy(ShiftAddr(this, offset), ShiftAddr(temp_node_.get(), offset), kPageSize - offset);

    // copy right records to this nodes
    offset = CopyRecordsFrom(r_node, 0, r_node->record_count_, offset);

    // update a header
    block_size_ = kPageSize - offset;
    deleted_size_ = 0;
    next_ = r_node->next_;

    const auto new_ver = mutex_.UnlockX();
    r_node->mutex_.UpgradeToX();

    // update a header of a right node
    r_node->is_removed_ = 1;
    r_node->next_ = this;

    r_node->mutex_.UnlockX();

    // reset a temp node
    temp_node_->record_count_ = 0;

    return new_ver;
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
    size_t node_size = kHeaderLen;
    auto offset = kPageSize;
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
      const auto high_key_len = meta_array_[record_count_ - 1].key_len;
      high_meta_ = Metadata{offset, high_key_len, high_key_len};
      offset -= high_key_len;  // reserve space for a highest key
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
      typename std::vector<Node *>::const_iterator &iter,
      const typename std::vector<Node *>::const_iterator &iter_end)
  {
    size_t node_size = kHeaderLen;
    auto offset = kPageSize;
    for (; iter < iter_end; ++iter) {
      const auto *child = *iter;
      const auto meta = child->high_meta_;
      const auto high_key_len = meta.key_len;
      const auto rec_len = high_key_len + kPtrLen;

      // check whether the node has sufficient space
      node_size += rec_len + kMetaLen;
      if (node_size > kPageSize - kMinFreeSpaceSize) break;

      // insert an entry to the inner node
      offset = SetPayload(offset, &child, kPtrLen);
      offset = CopyKeyFrom(child, meta, offset);
      meta_array_[record_count_++] = Metadata{offset, high_key_len, rec_len};
    }

    // move a highest key to header
    const auto rightmost_pos = record_count_ - 1;
    const auto high_key_len = meta_array_[rightmost_pos].key_len;
    high_meta_ = Metadata{offset, high_key_len, high_key_len};
    meta_array_[rightmost_pos] = Metadata{offset + high_key_len, 0, kPtrLen};

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

  /**
   * @brief Remove this node from a tree and return a new root node.
   *
   * @return a new root node.
   */
  auto
  RemoveRoot()  //
      -> Node *
  {
    mutex_.UpgradeToX();

    // get a child node as a new root node
    auto *child = GetPayload<Node *>(0);

    // remove this node from a tree
    is_removed_ = 1;
    next_ = nullptr;

    mutex_.UnlockX();
    return child;
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
  IsRightmostOf(const std::optional<std::pair<const Key &, bool>> &end_key) const  //
      -> bool
  {
    if (!next_) return true;     // the rightmost node
    if (!end_key) return false;  // perform full scan
    const auto &high_key = GetKey(high_meta_);
    return !Comp{}(high_key, end_key->first);
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

  /**
   * @param new_rec_len the length of a new record.
   * @param ver an expected version value.
   * @retval kNeedRetry if a version check fails.
   * @retval kNeedSplit if this node requires splitting.
   * @retval kCompleted otherwise.
   */
  [[nodiscard]] auto
  NeedSplit(  //
      const size_t new_rec_len,
      uint64_t &ver)  //
      -> NodeRC
  {
    // check whether the node has space for a new record
    const auto total_size = kMetaLen * (record_count_ + 1) + block_size_ + new_rec_len;
    if (total_size <= kPageSize - kHeaderLen) {
      if (mutex_.HasSameVersion(ver)) return kCompleted;
      return kNeedRetry;
    }

    // check the node requires splitting
    if (total_size - deleted_size_ > kMaxUsedSpaceSize) {
      if (mutex_.HasSameVersion(ver)) return kNeedSplit;
      return kNeedRetry;
    }

    // this node has enough space but cleaning up is required
    if (!mutex_.TryLockX(ver)) return kNeedRetry;
    CleanUp();
    ver = mutex_.UnlockX();
    return kCompleted;
  }

  /**
   * @param ver an expected version value.
   * @retval kNeedRetry if a version check fails.
   * @retval kNeedMerge if this node requires merging.
   * @retval kCompleted otherwise.
   */
  [[nodiscard]] auto
  NeedMerge(uint64_t &ver)  //
      -> NodeRC
  {
    // check this node uses enough space
    if (GetUsedSize() < kMinUsedSpaceSize) {
      if (mutex_.HasSameVersion(ver)) return kNeedMerge;
      return kNeedRetry;
    }

    // check this node has a lot of dead space
    if (deleted_size_ <= kMaxDeletedSpaceSize) {
      if (mutex_.HasSameVersion(ver)) return kCompleted;
      return kNeedRetry;
    }

    // this node has a lot of dead space
    if (!mutex_.TryLockX(ver)) return kNeedRetry;
    CleanUp();
    ver = mutex_.UnlockX();
    return kCompleted;
  }

  /**
   * @param new_rec_len the length of a new record.
   * @param ver an expected version value.
   * @retval kNeedRetry if a version check fails.
   * @retval kNeedSplit if this node requires splitting.
   * @retval kNeedMerge if this node requires merging.
   * @retval kCompleted otherwise.
   */
  [[nodiscard]] auto
  CheckSMOs(  //
      const size_t new_rec_len,
      uint64_t &ver)  //
      -> NodeRC
  {
    // check splitting is needed
    auto rc = NeedSplit(new_rec_len, ver);
    if (rc != kCompleted) return rc;

    // check merging is needed
    rc = NeedMerge(ver);
    return rc;
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
  ::dbgroup::lock::OptimisticLock mutex_{};

  /// the pointer to the next node.
  Node *next_{nullptr};

  /// the metadata of a highest key.
  Metadata high_meta_{kPageSize, 0, 0};

  /// an actual data block (it starts with record metadata).
  Metadata meta_array_[0];

  // a temporary node for SMOs.
  static thread_local inline std::unique_ptr<Node>          //
      temp_node_{new (::operator new(kPageSize)) Node{0}};  // NOLINT
};

}  // namespace dbgroup::index::b_tree::component::oml

#endif  // B_TREE_COMPONENT_OML_NODE_VARLEN_HPP
