
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

#ifndef B_TREE_COMPONENT_VERSIONED_OSL_NODE_FIXLEN_HPP
#define B_TREE_COMPONENT_VERSIONED_OSL_NODE_FIXLEN_HPP

#include <atomic>
#include <optional>
#include <utility>
#include <vector>

// organization libraries
#include "lock/optimistic_lock.hpp"
#include "memory/epoch_based_gc.hpp"
#include "memory/epoch_manager.hpp"

// local sources
#include "b_tree/component/common.hpp"
#include "version_record.hpp"
namespace dbgroup::index::b_tree::component::osl
{

/**
 * @brief A class for representing nodes in B+trees.
 *
 * This class uses optimistic single-layer locking for concurrency controls and
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
  using ScanKey = std::optional<std::tuple<const Key &, size_t, bool>>;
  template <class Entry>
  using BulkIter = typename std::vector<Entry>::const_iterator;
  using NodeEntry = std::tuple<Key, Node *, size_t>;
  using GC_t = ::dbgroup::memory::EpochBasedGC<Node>;
  using EpochManager_t = ::dbgroup::memory::EpochManager;
  using Timestamp_t = size_t;

  /*####################################################################################
   * Public constructors and assignment operators
   *##################################################################################*/

  /**
   * @brief Construct an empty node object.
   *
   * @param is_inner a flag to indicate whether a inner node is constructed.
   */
  constexpr explicit NodeFixLen(const uint32_t is_inner)
      : is_inner_{is_inner}, is_removed_{0}, block_size_{0}, has_low_key_{0}, has_high_key_{0}
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
      : is_inner_{1},
        is_removed_{0},
        block_size_{2 * kPtrLen},
        record_count_{2},
        has_low_key_{0},
        has_high_key_{0}
  {
    keys_[1] = l_key;
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
  GetValidSplitNode(const Key &key)  //
      -> std::tuple<Node *, Key, size_t>
  {
    auto sep_key = GetHighKey();

    auto *node = this;
    if (!Comp{}(key, sep_key)) {
      node = next_;
      mutex_.UnlockSIX();
    } else {
      next_->mutex_.UnlockSIX();
    }

    return {node, std::move(sep_key), kKeyLen};
  }

  /**
   * @brief Get a mergeable-sibling node if exist.
   *
   * The returned node is locked with a SIX lock. If there is no mergeable node, this
   * node is unlocked.
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
      mutex_.UnlockSIX();
      return nullptr;
    }

    // check the right-sibling node has enough capacity for merging
    next_->mutex_.LockSIX();
    if (GetUsedSize() + next_->GetUsedSize() >= kMaxUsedSpaceSize) {
      next_->mutex_.UnlockSIX();
      mutex_.UnlockSIX();
      return nullptr;
    }

    return next_;
  }

  /**
   * @retval a lowest key in this node if exist.
   * @retval std::nullopt otherwise.
   */
  [[nodiscard]] auto
  GetLowKey() const  //
      -> std::optional<Key>
  {
    while (true) {
      const auto ver = mutex_.GetVersion();

      std::optional<Key> low_key = std::nullopt;
      if (has_low_key_) {
        low_key = keys_[record_count_ + has_high_key_];
      }

      if (mutex_.HasSameVersion(ver)) return low_key;
    }
  }

  /**
   * @retval 1st: a highest key.
   * @retval 2nd: the length of the highest key.
   */
  [[nodiscard]] auto
  GetHighKeyForSMOs() const  //
      -> std::pair<Key, size_t>
  {
    return {GetHighKey(), kKeyLen};
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
   * @brief Get a leftmost child node.
   *
   * @return the child node.
   */
  [[nodiscard]] auto
  GetLeftmostChild()  //
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
   * @param pos the position of a record.
   * @retval true if the record is deleted.
   * @retval false otherwise.
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
    return {keys_[pos], GetPayload<Payload>(pos)};
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
      const auto &index_key = keys_[pos];

      if (Comp{}(key, index_key)) {  // a target key is in a left side
        end_pos = pos - 1;
      } else if (Comp{}(index_key, key)) {  // a target key is in a right side
        begin_pos = pos + 1;
      } else {  // find an equivalent key
        return {kKeyAlreadyInserted, pos};
      }
    }

    return {kKeyNotInserted, begin_pos - is_inner_};
  }

  /**
   * @brief Get a child node of a specified key by using binary search.
   *
   * If there is no specified key in this node, this returns the minimum position that
   * is greater than the specified key.
   *
   * @param node a current node to be searched.
   * @param key a search key.
   * @return the child node that includes the given key.
   */
  [[nodiscard]] static auto
  SearchChild(  //
      Node *&node,
      const Key &key)  //
      -> Node *
  {
    Node *child{};
    while (true) {
      const auto ver = CheckKeyRange(node, key);
      if (node == nullptr) break;  // a root node was removed

      auto pos = node->SearchRecord(key).second;

      child = node->GetPayload<Node *>(pos);
      if (node->mutex_.HasSameVersion(ver)) break;
    }

    return child;
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

  /**
   * @brief Check the key range of a given node and traverse side links if needed.
   *
   * @param node a current node to be checked.
   * @param key a search key.
   * @return a node whose key range includes the search key.
   */
  static auto
  CheckKeyRange(  //
      Node *&node,
      const Key &key)  //
      -> uint64_t
  {
    uint64_t ver{};
    while (true) {
      ver = node->mutex_.GetVersion();

      // check the node is not removed
      if (node->is_removed_ == 0) {
        // check the node includes a target key
        if (node->has_high_key_ == 0) {
          if (node->mutex_.HasSameVersion(ver)) break;
          continue;
        }
        const auto &high_key = node->GetHighKey();
        if (Comp{}(key, high_key)) {
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
   * @return a node whose key range includes the search key.
   */
  static void
  CheckKeyRangeAndLockForRead(  //
      Node *&node,
      const Key &key)
  {
    while (true) {
      const auto ver = node->mutex_.GetVersion();

      // check the node is not removed
      if (node->is_removed_ == 0) {
        // check the node includes a target key
        if (node->has_high_key_ == 0) {
          if (node->mutex_.TryLockS(ver)) return;
          continue;
        }
        const auto &high_key = node->GetHighKey();
        if (Comp{}(key, high_key)) {
          if (node->mutex_.TryLockS(ver)) return;
          continue;
        }
      }

      // go to the next node
      auto *next = node->next_;
      if (!node->mutex_.HasSameVersion(ver)) continue;

      node = next;
      if (node == nullptr) return;
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
  static void
  CheckKeyRangeAndLockForWrite(  //
      Node *&node,
      const Key &key)
  {
    while (true) {
      const auto ver = node->mutex_.GetVersion();

      // check the node is not removed
      if (node->is_removed_ == 0) {
        // check the node includes a target key
        if (node->has_high_key_ == 0 || Comp{}(key, node->GetHighKey())) {
          if (node->mutex_.TryLockSIX(ver)) return;
          continue;
        }
      }

      // go to the next node
      auto *next = node->next_;
      if (!node->mutex_.HasSameVersion(ver)) continue;

      node = next;
      if (node == nullptr) return;
    }
  }

  /*####################################################################################
   * Leaf read operations
   *##################################################################################*/

  /**
   * @brief Read a payload of a specified key if it exists.
   *
   * @tparam Payload a class of payload.
   * @param node a current node to be read.
   * @param key a target key.
   * @param out_payload a reference to be stored a target payload.
   * @retval kKeyAlreadyInserted if a target record is read.
   * @retval kKeyAlreadyDeleted if a target record is deleted.
   * @retval kKeyNotInserted otherwise.
   */
  template <class Payload>
  static auto
  Read(  //
      Node *&node,
      const Key &key,
      Payload &out_payload)  //
      -> NodeRC
  {
    auto current_version = VersionRecord<Payload>{};
    while (true) {
      const auto ver = CheckKeyRange(node, key);

      auto [existence, pos] = node->SearchRecord(key);
      if (existence == kKeyAlreadyInserted) {
        memcpy(&current_version, node->GetPayloadAddr(pos), sizeof(VersionRecord<Payload>));
        if (current_version.IsDeleted()) {
          existence = kKeyAlreadyDeleted;
        } else {
          out_payload = current_version.GetPayload();
        }
      }

      if (node->mutex_.HasSameVersion(ver)) return existence;
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
   * @tparam Payload a class of payload.
   * @param key a target key to be written.
   * @param key_len the length of a target key.
   * @param payload a target payload to be written.
   * @param pay_len the length of a target payload.
   * @retval kCompleted if a record is written.
   * @retval kNeedSplit if this node should be split before inserting a record.
   */
  template <class Payload>
  auto
  Write(  //
      const Key &key,
      [[maybe_unused]] const size_t key_len,
      const Payload &payload,
      [[maybe_unused]] const size_t pay_len,
      [[maybe_unused]] GC_t &gc,
      EpochManager_t &epoch_manager)  //
      -> NodeRC
  {
    const auto rec_len = kKeyLen + pay_len_;

    // search position where this key has to be set
    const auto [existence, pos] = SearchRecord(key);
    if (existence == kKeyNotInserted) {
      if (GetUsedSize() + rec_len > kPageSize - kHeaderLen) return kNeedSplit;

      // insert a new version record
      InsertRecord(key, kKeyLen, payload, pay_len_, pos, epoch_manager);
      return kCompleted;
    }

    // there is a record with the same key, so reuse it
    AppendNewVersionRecord(pos, payload, pay_len_, mutex_, epoch_manager);

    return kCompleted;
  }

  /**
   * @brief Insert a specified kay/payload pair.
   *
   * If a specified key does not exist, this function insert a target payload into a
   * target leaf node. If a specified key exists in a target leaf node, this function
   * does nothing and returns kKeyExist as a return code.
   *
   * @param node a current node to be updated.
   * @param key a target key to be written.
   * @param key_len the length of a target key.
   * @param payload a target payload to be written.
   * @param pay_len the length of a target payload.
   * @retval kCompleted if a record is inserted.
   * @retval kKeyAlreadyInserted if there is a record with the same key.
   * @retval kNeedSplit if this node should be split before inserting a record.
   */
  template <class Payload>
  static auto
  Insert(  //
      Node *&node,
      const Key &key,
      [[maybe_unused]] const size_t key_len,
      Payload &payload,
      [[maybe_unused]] const size_t pay_len,
      [[maybe_unused]] GC_t &gc,
      EpochManager_t &epoch_manager)  //
      -> NodeRC
  {
    const auto rec_len = kKeyLen + node->pay_len_;

    while (true) {
      const auto ver = CheckKeyRange(node, key);

      // search position where this key has to be set
      const auto [existence, pos] = node->SearchRecord(key);
      if (existence == kKeyAlreadyInserted) {
        if (!node->mutex_.HasSameVersion(ver)) continue;
        auto current_ver = VersionRecord<Payload>{};
        memcpy(&current_ver, node->GetPayloadAddr(pos), sizeof(VersionRecord<Payload>));
        if (current_ver.IsDeleted()) {
          // if the key is already inserted and deleted, Insert() behaves as Update()
          node->AppendNewVersionRecord(pos, payload, node->pay_len_, node->mutex_, epoch_manager);
          return kCompleted;
        } else {  //
        return kKeyAlreadyInserted;
      }
      } else {
        // if the key has never been inserted
      // inserting a new record is required
      if (!node->mutex_.TryLockSIX(ver)) continue;
      if (node->GetUsedSize() + rec_len > kPageSize - kHeaderLen) return kNeedSplit;

        // insert a new version record
        auto current_epoch = epoch_manager.GetCurrentEpoch();
        auto new_version = VersionRecord<Payload>{current_epoch, payload};
        node->InsertRecord(key, kKeyLen, payload, node->pay_len_, pos, epoch_manager);
      return kCompleted;
      }
    }
  }

  /**
   * @brief Insert a given record into this node.
   *
   * @param key a target key to be set.
   * @param key_len the length of the key.
   * @param payload a target payload to be written.
   * @param pay_len the length of a target payload.
   * @param pos an insertion position.
   */
  template <class Payload>
  void
  InsertRecord(  //
      const Key &key,
      [[maybe_unused]] const size_t key_len,
      const Payload &payload,
      [[maybe_unused]] const size_t pay_len,
      const size_t pos,
      EpochManager_t &epoch_manager)  //
  {
    mutex_.UpgradeToX();

    // insert a new key
    const auto move_num = record_count_ - pos;
    memmove(&(keys_[pos + 1]), &(keys_[pos]), kKeyLen * (move_num + has_low_key_ + has_high_key_));
    keys_[pos] = key;

    // insert a new version record
    auto current_epoch = epoch_manager.GetCurrentEpoch();
    auto new_version = VersionRecord<Payload>{current_epoch, payload};

    const auto top_offset = kPageSize - block_size_;
    const auto move_size = pay_len_ * move_num;
    memmove(ShiftAddr(this, top_offset - pay_len_), ShiftAddr(this, top_offset), move_size);
    SetPayload(top_offset + move_size, &new_version);

    // update header information
    ++record_count_;
    block_size_ += pay_len_;

    mutex_.UnlockX();
  }

  /**
   * @brief Update a target kay with a specified payload.
   *
   * If a specified key exist, this function update a target payload. If a specified key
   * does not exist in a target leaf node, this function does nothing and returns
   * kKeyNotExist as a return code.
   *
   * @param node a current node to be updated.
   * @param key a target key to be written.
   * @param payload a target payload to be written.
   * @param pay_len the length of a target payload.
   * @retval kSuccess if a key/payload pair is written.
   * @retval kKeyNotExist otherwise.
   */
  template <class Payload>
  static auto
  Update(  //
      Node *&node,
      const Key &key,
      Payload &payload,
      [[maybe_unused]] const size_t pay_len,
      [[maybe_unused]] GC_t &gc,
      EpochManager_t &epoch_manager)  //
      -> ReturnCode
  {
    while (true) {
      const auto ver = CheckKeyRange(node, key);

      // check this node has a target record
      const auto [existence, pos] = node->SearchRecord(key);
      if (existence != kKeyAlreadyInserted) {
        if (node->mutex_.HasSameVersion(ver)) return kKeyNotExist;
        continue;
      }

      // a target record exists, so try to acquire an exclusive lock
      if (!node->mutex_.TryLockX(ver)) continue;

      // perform update operation
      node->AppendNewVersionRecord(pos, payload, node->pay_len_, node->mutex_, epoch_manager);

      return kSuccess;
    }
  }

  /**
   * @brief Delete a target key from the index.
   *
   * If a specified key exist, this function overwrites it with deleted version record. If a
   * specified key does not exist in a leaf node, this function does nothing and returns
   * kKeyNotExist as a return code.
   *
   * @param node a current node to be updated.
   * @param key a target key to be written.
   * @retval kSuccess if a key/payload pair is written by deleted version record.
   * @retval kKeyNotExist otherwise.
   */
  template <class Payload>
  static auto
  Delete(  //
      Node *&node,
      const Key &key,
      [[maybe_unused]] const size_t pay_len,
      [[maybe_unused]] GC_t &gc,
      EpochManager_t &epoch_manager)  //
      -> ReturnCode
  {
    while (true) {
      const auto ver = CheckKeyRange(node, key);

      // check this node has a target record
      const auto [existence, pos] = node->SearchRecord(key);
      if (existence != kKeyAlreadyInserted) {
        if (node->mutex_.HasSameVersion(ver)) return kKeyNotExist;
        continue;
      }

      // a target record exists, so try to acquire an exclusive lock
      if (!node->mutex_.TryLockX(ver)) continue;

      // perform delete operation
      node->AppendNewVersionRecord(pos, Payload{}, node->pay_len_, node->mutex_, epoch_manager,
                                   true);

      return kSuccess;
    }
  }

  /**
   * @brief Insert a new child node into this node.
   *
   * @param r_node a right child (i.e., new) node.
   * @param sep_key a separator key.
   * @param sep_key_len the length of the separator key.
   * @retval kCompleted if a new entry is inserted.
   * @retval kNeedRetry if previous merging has not been finished.
   * @retval kNeedSplit if this node should be split before inserting an index entry.
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

    // check an inserting position and concurrent SMOs
    const auto [existence, pos] = SearchRecord(sep_key);
    if (existence == kKeyAlreadyInserted) {
      // previous merging has not been applied, so unlock and retry
      mutex_.UnlockSIX();
      return kNeedRetry;
    }

    // recheck free space in this node
    if (GetUsedSize() + rec_len > kPageSize - kHeaderLen) return kNeedSplit;

    mutex_.UpgradeToX();

    // insert a split-right child node
    const auto move_num = record_count_ - 1 - pos;
    const auto move_size = kPtrLen * move_num;
    const auto top_offset = kPageSize - block_size_;
    memmove(ShiftAddr(this, top_offset - kPtrLen), ShiftAddr(this, top_offset), move_size);
    SetPayload(top_offset + move_size, &r_node);

    // insert a separator key
    memmove(&(keys_[pos + 2]), &(keys_[pos + 1]),
            kKeyLen * (move_num + has_low_key_ + has_high_key_));
    keys_[pos + 1] = sep_key;

    // update header information
    ++record_count_;
    block_size_ += kPtrLen;

    mutex_.UnlockX();
    return kCompleted;
  }

  /**
   * @brief Delete a child node from this node.
   *
   * @param r_child a right child (i.e., removed) node.
   * @param del_key a separater key to be deleted.
   * @retval kCompleted if the entry is deleted.
   * @retval kNeedWaitAndRetry if previous splitting has not been finished.
   * @retval kNeedMerge if this node should be merged.
   * @retval kAbortMerge if this merge operation should be aborted.
   */
  auto
  DeleteChild(const Key &del_key)  //
      -> NodeRC
  {
    // check a child node to be deleted is not rightmost
    if (record_count_ <= 1 || Comp{}(del_key, GetKey(1))) {
      mutex_.UnlockSIX();
      return kAbortMerge;
    }

    // check a position to be deleted and concurrent SMOs
    const auto [existence, r_pos] = SearchRecord(del_key);
    if (existence == kKeyNotInserted) {
      // previous splitting has not been applied, so unlock and retry
      mutex_.UnlockSIX();
      return kNeedRetry;
    }

    mutex_.UpgradeToX();

    // delete a separator key
    const auto move_num = record_count_ - 1 - r_pos;
    memmove(&(keys_[r_pos]), &(keys_[r_pos + 1]),
            kKeyLen * (move_num + has_low_key_ + has_high_key_));

    // delete a right child
    auto *top_addr = ShiftAddr(this, kPageSize - block_size_);
    memmove(ShiftAddr(top_addr, kPtrLen), top_addr, kPtrLen * move_num);

    // update this header
    --record_count_;
    block_size_ -= kPtrLen;

    if (GetUsedSize() < kMinUsedSpaceSize) {
      mutex_.DowngradeToSIX();
      return kNeedMerge;
    }

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
    const auto l_count = record_count_ / 2;
    const auto r_count = record_count_ - l_count;

    r_node->mutex_.LockX();

    // copy right half records to a right node
    r_node->pay_len_ = pay_len_;
    auto r_offset = r_node->CopyRecordsFrom(this, l_count, record_count_, kPageSize);
    r_node->keys_[r_count - 1 + has_high_key_] =
        keys_[record_count_ - 1 + has_high_key_];               // a highest key
    r_node->keys_[r_count + has_high_key_] = r_node->keys_[0];  // a lowest key

    // update a right header
    r_node->block_size_ = kPageSize - r_offset;
    r_node->next_ = next_;
    r_node->has_low_key_ = 1;
    r_node->has_high_key_ = has_high_key_;

    r_node->mutex_.DowngradeToSIX();
    mutex_.UpgradeToX();

    // update lowest/highest keys
    keys_[l_count + has_low_key_] = keys_[0];  // a lowest key
    keys_[l_count] = r_node->keys_[0];         // a highest key

    // update a header
    record_count_ = l_count;
    block_size_ -= r_node->block_size_;
    next_ = r_node;
    has_high_key_ = 1;

    mutex_.DowngradeToSIX();
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

    // copy right records to this nodes
    const auto lowest_key = keys_[record_count_ - 1 + has_high_key_ + has_low_key_];
    auto offset = CopyRecordsFrom(r_node, 0, r_node->record_count_, kPageSize - block_size_);
    if (r_node->has_high_key_) {
      keys_[record_count_ + has_low_key_] = std::move(lowest_key);
      keys_[record_count_] = r_node->keys_[r_node->record_count_];
    } else {
      keys_[record_count_] = std::move(lowest_key);
    }

    // update a header
    block_size_ = kPageSize - offset;
    next_ = r_node->next_;
    has_high_key_ = r_node->has_high_key_;

    mutex_.UnlockX();
    r_node->mutex_.UpgradeToX();

    // update a header of a right node
    r_node->is_removed_ = 1;
    r_node->next_ = this;

    r_node->mutex_.UnlockX();
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
    constexpr auto kKeyLen = sizeof(Key);
    constexpr auto kRecLen = kKeyLen + kPtrLen;

    // extract and insert entries into this node
    auto offset = kPageSize;
    auto node_size = kHeaderLen;
    for (; iter < iter_end; ++iter) {
      // check whether the node has sufficent space
      node_size += kRecLen;
      if (node_size + 2 * kKeyLen > kPageSize) break;

      // insert an entry into this node
      const auto &[key, payload, key_len] = ParseEntry(*iter);
      offset = SetPayload(offset, &payload);
      keys_[record_count_++] = key;
    }

    block_size_ = kPageSize - offset;

    // link the sibling nodes if exist
    if (prev_node != nullptr) {
      prev_node->LinkNext(this);
    }

    // set a lowest key
    has_low_key_ = 1;
    keys_[record_count_] = keys_[0];

    nodes.emplace_back(keys_[0], this, kKeyLen);
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
      node->keys_[0] = Key{};

      // go down to the lower level
      node = node->template GetPayload<Node *>(0);
    }
  }

 private:
  /*####################################################################################
   * Internal constants
   *##################################################################################*/

  /// the length of keys.
  static constexpr size_t kKeyLen = sizeof(Key);

  /// the length of node pointers.
  static constexpr size_t kPtrLen = sizeof(uintptr_t);

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
    return Comp{}(std::get<0>(*end_key), GetHighKey());
  }

  /**
   * @return Current usage of the record block.
   */
  [[nodiscard]] auto
  GetUsedSize() const  //
      -> size_t
  {
    return kKeyLen * (record_count_ + has_low_key_ + has_high_key_) + block_size_;
  }

  /**
   * @retval a highest key in this node.
   */
  [[nodiscard]] auto
  GetHighKey() const  //
      -> const Key &
  {
    return keys_[record_count_];
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
   */
  void
  LinkNext(Node *r_node)
  {
    // set a highest key in a left node
    has_high_key_ = 1;
    keys_[record_count_ + 1] = keys_[record_count_];
    keys_[record_count_] = r_node->keys_[0];

    // set a sibling link in a left node
    next_ = r_node;
  }


  /**
   * @brief Update the record, and append the new version on the chain.
   *
   * @param pos the position of the record
   * @param payload a payload of the version to write
   * @param pay_len the length of payload of the node which will be written
   * @param mutex an optimistic lock
   * @param current_epoch the timestamp of the version to write
   * @param is_delete a flag that indicates whether this function is called by Delete() or others
   *
   */
  template <class Payload>
  auto
  AppendNewVersionRecord(  //
      const size_t pos,
      const Payload &payload,
      const size_t &pay_len,
      ::dbgroup::lock::OptimisticLock &mutex,
      EpochManager_t &epoch_manager,
      bool is_delete = false)  //
      -> void
  {
    const auto &protected_epochs = epoch_manager.GetProtectedEpochs();
    auto current_epoch = protected_epochs->front();
    // TODO: GC
    auto &&new_version = VersionRecord<Payload>{current_epoch, payload, is_delete};
    auto *old_version_ptr = new VersionRecord<Payload>{};  // This is the chain's head before update,
                                                          // and will be filled by memcpy.
    new_version.SetNextPtr(old_version_ptr);              // New head's next version is old one.

    mutex.UpgradeToX();
    memcpy(old_version_ptr, GetPayloadAddr(pos), pay_len);  // Copy old one.
    memcpy(GetPayloadAddr(pos), &new_version, pay_len);
    mutex.UnlockX();
  }
  /*####################################################################################
   * Internal member variables
   *##################################################################################*/

  /// a flag for indicating this node is a leaf or internal node.
  uint32_t is_inner_ : 1;

  /// a flag for indicating this node is removed from a tree.
  uint32_t is_removed_ : 1;

  /// the total byte length of records in a node.
  uint32_t block_size_ : 30;

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
};

}  // namespace dbgroup::index::b_tree::component::osl

#endif  // B_TREE_COMPONENT_VERSIONED_OSL_NODE_FIXLEN_HPP
