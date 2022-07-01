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

#ifndef B_TREE_B_TREE_HPP
#define B_TREE_B_TREE_HPP

#include <optional>
#include <vector>

#include "component/node.hpp"

namespace dbgroup::index::b_tree
{
/**
 * @brief A class for representing Btree with variable-length keys.
 *
 * This implementation can store variable-length keys (i.e., 'text' type in PostgreSQL).
 *
 * @tparam Node a class of stored nodes.
 * @tparam Key a class of stored keys.
 * @tparam Payload a class of stored payloads (only fixed-length data for simplicity).
 * @tparam Comp a class for ordering keys.
 */
template <template <class Key, class Comp> class Node,
          class Key,
          class Payload,
          class Comp = ::std::less<Key>>
class BTreePCL
{
 public:
  /*####################################################################################
   * Type aliases
   *##################################################################################*/

  using Node_t = Node<Key, Comp>;
  using NodeRC = component::NodeRC;
  using NodeStack = std::vector<std::pair<Node_t *, size_t>>;
  using Lock = ::dbgroup::lock::PessimisticLock;

  /**
   * @brief A class for representing an iterator of scan results.
   *
   */
  class RecordIterator
  {
   public:
    /*##################################################################################
     * Public constructors and assignment operators
     *################################################################################*/

    RecordIterator(  //
        Node_t *node,
        size_t count,
        size_t pos,
        const std::optional<std::pair<const Key &, bool>> end_key,
        bool is_end)
        : node_{node},
          record_count_{count},
          current_pos_{pos},
          end_key_{std::move(end_key)},
          is_end_{is_end}
    {
    }

    RecordIterator(const RecordIterator &) = delete;
    RecordIterator(RecordIterator &&) = delete;

    auto operator=(const RecordIterator &) -> RecordIterator & = delete;
    auto operator=(RecordIterator &&) -> RecordIterator & = delete;

    /*##################################################################################
     * Public destructors
     *################################################################################*/

    ~RecordIterator() = default;

    /*##################################################################################
     * Public operators for iterators
     *################################################################################*/

    auto
    operator*() const  //
        -> std::pair<Key, Payload>
    {
      return node_->template GetRecord<Payload>(current_pos_);
    }

    constexpr void
    operator++()
    {
      ++current_pos_;
    }

    /*##################################################################################
     * Public getters
     *################################################################################*/

    /**
     * @brief Check if there are any records left.
     *
     * @retval true if there are any records or next node left.
     * @retval false otherwise.
     */
    [[nodiscard]] auto
    HasNext()  //
        -> bool
    {
      while (true) {
        if (current_pos_ < record_count_) return true;
        if (is_end_) {
          node_->ReleaseSharedLock();
          return false;
        }
        node_ = node_->GetNextNodeForRead();
        current_pos_ = 0;
        std::tie(is_end_, record_count_) = node_->SearchEndPositionFor(end_key_);
      }
    }

    /**
     * @return a Key of a current record
     */
    [[nodiscard]] auto
    GetKey() const  //
        -> Key
    {
      return node_->GetKey(current_pos_);
    }

    /**
     * @return a payload of a current record
     */
    [[nodiscard]] auto
    GetPayload() const  //
        -> Payload
    {
      return node_->template GetPayload<Payload>(current_pos_);
    }

   private:
    /// the pointer to a node that includes partial scan results.
    Node_t *node_{nullptr};

    /// the number of records in this node.
    size_t record_count_{0};

    /// the position of a current record.
    size_t current_pos_{0};

    /// the end key given from a user.
    std::optional<std::pair<const Key &, bool>> end_key_{};

    /// a flag for indicating a current node is rightmost in scan-range.
    bool is_end_{false};
  };

  /*####################################################################################
   * Public constructors and assignment operators
   *##################################################################################*/

  /**
   * @brief Construct a new BTreePCL object.
   *
   */
  explicit BTreePCL()
  {
    if constexpr (IsVariableLengthData<Key>()) {
      static_assert(sizeof(component::Metadata) + kMaxVarDataSize + sizeof(Payload)
                    < (kPageSize / 2));
    } else {
      static_assert(sizeof(component::Metadata) + sizeof(Key) + sizeof(Payload) < (kPageSize / 2));
    }
  };

  BTreePCL(const BTreePCL &) = delete;
  BTreePCL(BTreePCL &&) = delete;

  BTreePCL &operator=(const BTreePCL &) = delete;
  BTreePCL &operator=(BTreePCL &&) = delete;

  /*####################################################################################
   * Public destructors
   *##################################################################################*/

  /**
   * @brief Destroy the BTreePCL object.
   *
   */
  ~BTreePCL() = default;

  /*####################################################################################
   * Public read APIs
   *##################################################################################*/

  /**
   * @brief Read the payload corresponding to a given key if it exists.
   *
   * @param key a target key.
   * @retval the payload of a given key wrapped with std::optional if it is in this tree.
   * @retval std::nullopt otherwise.
   */
  auto
  Read(const Key &key)  //
      -> std::optional<Payload>
  {
    auto *node = SearchLeafNodeForRead(key, kClosed);

    Payload payload{};
    const auto rc = node->Read(key, payload);

    if (rc == NodeRC::kCompleted) return payload;
    return std::nullopt;
  }

  /**
   * @brief Perform a range scan with given keys.
   *
   * @param begin_key a pair of a begin key and its openness (true=closed).
   * @param end_key a pair of an end key and its openness (true=closed).
   * @return an iterator to access scanned records.
   */
  auto
  Scan(  //
      const std::optional<std::pair<const Key &, bool>> &begin_key = std::nullopt,
      const std::optional<std::pair<const Key &, bool>> &end_key = std::nullopt)  //
      -> RecordIterator
  {
    Node_t *node{};
    size_t begin_pos = 0;

    if (begin_key) {
      const auto &[key, is_closed] = begin_key.value();
      node = SearchLeafNodeForRead(key, is_closed);
      const auto [rc, pos] = node->SearchRecord(key);
      begin_pos = (rc == NodeRC::kKeyAlreadyInserted && !is_closed) ? pos + 1 : pos;
    } else {
      node = SearchLeftmostLeaf();
    }

    const auto [is_end, end_pos] = node->SearchEndPositionFor(end_key);
    return RecordIterator{node, end_pos, begin_pos, end_key, is_end};
  }

  /*####################################################################################
   * Public write APIs
   *##################################################################################*/

  /**
   * @brief Write (i.e., put) a given key/payload pair.
   *
   * If a given key does not exist in this tree, this function performs an insert
   * operation. If a given key has been already inserted, this function perfroms an
   * update operation. Thus, this function always returns kSuccess as a return code.
   *
   * @param key a target key to be written.
   * @param payload a target payload to be written.
   * @param key_len the length of a target key.
   * @return kSuccess.
   */
  auto
  Write(  //
      const Key &key,
      const Payload &payload,
      const size_t key_len = sizeof(Key))  //
      -> ReturnCode
  {
    auto &&stack = SearchLeafNodeForWrite(key, !kDelOps);
    auto *node = stack.back().first;

    const auto rc = node->template Write<Payload>(key, key_len, payload);
    if (rc != NodeRC::kCompleted) {
      if (rc == NodeRC::kNeedSplit) {
        Split<Payload>(key, stack);
        node = node->GetValidSplitNode(key);
      } else {
        stack.pop_back();
        ReleaseExclusiveLocks(stack);  // unlock ancestor nodes
        node->template Consolidate<Payload>();
      }

      // insert a record again
      node->template Write<Payload>(key, key_len, payload);
      stack.emplace_back(node, 0);  // add a current node to release its lock
    }

    ReleaseExclusiveLocks(stack);
    return kSuccess;
  }

  /**
   * @brief Insert a given key/payload pair.
   *
   * This function performs a uniqueness check in its processing. If a given key does
   * not exist in this tree, this function inserts a target payload to this tree. If
   * there is a given key in this tree, this function does nothing and returns kKeyExist
   * as a return code.
   *
   * @param key a target key to be inserted.
   * @param payload a target payload to be inserted.
   * @param key_len the length of a target key.
   * @retval kSuccess if inserted.
   * @retval kKeyExist otherwise.
   */
  auto
  Insert(  //
      const Key &key,
      const Payload &payload,
      const size_t key_len = sizeof(Key))  //
      -> ReturnCode
  {
    auto &&stack = SearchLeafNodeForWrite(key, !kDelOps);
    auto *node = stack.back().first;

    auto rc = node->template Insert<Payload>(key, key_len, payload);
    if (rc != NodeRC::kCompleted && rc != NodeRC::kKeyAlreadyInserted) {
      if (rc == NodeRC::kNeedSplit) {
        Split<Payload>(key, stack);
        node = node->GetValidSplitNode(key);
      } else {
        stack.pop_back();
        ReleaseExclusiveLocks(stack);  // unlock ancestor nodes
        node->template Consolidate<Payload>();
      }

      // insert a record again
      rc = node->template Insert<Payload>(key, key_len, payload);
      stack.emplace_back(node, 0);  // add a current node to release its lock
    }

    ReleaseExclusiveLocks(stack);
    return (rc == NodeRC::kCompleted) ? kSuccess : kKeyExist;
  }

  /**
   * @brief Update the record corresponding to a given key with a given payload.
   *
   * This function performs a uniqueness check in its processing. If there is a given
   * key in this tree, this function updates the corresponding record. If a given key
   * does not exist in this tree, this function does nothing and returns kKeyNotExist as
   * a return code.
   *
   * @param key a target key to be updated.
   * @param payload a payload for updating.
   * @param key_len the length of a target key.
   * @retval kSuccess if updated.
   * @retval kKeyNotExist otherwise.
   */
  auto
  Update(  //
      const Key &key,
      const Payload &payload)  //
      -> ReturnCode
  {
    auto &&stack = SearchLeafNodeForWrite(key, !kDelOps);
    auto *node = stack.back().first;

    const auto rc = node->template Update<Payload>(key, payload);

    ReleaseExclusiveLocks(stack);
    return (rc == NodeRC::kCompleted) ? kSuccess : kKeyNotExist;
  }

  /**
   * @brief Delete the record corresponding to a given key from this tree.
   *
   * This function performs a uniqueness check in its processing. If there is a given
   * key in this tree, this function deletes it. If a given key does not exist in this
   * tree, this function does nothing and returns kKeyNotExist as a return code.
   *
   * @param key a target key to be deleted.
   * @retval kSuccess if deleted.
   * @retval kKeyNotExist otherwise.
   */
  auto
  Delete(const Key &key)  //
      -> ReturnCode
  {
    auto &&stack = SearchLeafNodeForWrite(key, kDelOps);
    auto *node = stack.back().first;

    const auto rc = node->Delete(key);
    if (rc == NodeRC::kNeedMerge) {
      Merge<Payload>(stack);
    }

    ReleaseExclusiveLocks(stack);
    return (rc == NodeRC::kKeyNotInserted) ? kKeyNotExist : kSuccess;
  }

 private:
  /*################################################################################################
   * Internal constants
   *##############################################################################################*/

  // an expected maximum height of a tree.
  static constexpr size_t kExpectedTreeHeight = 8;

  // a flag for indicating closed-interval
  static constexpr bool kClosed = true;

  // a flag for indicating delete operations
  static constexpr bool kDelOps = true;

  /*####################################################################################
   * Internal utility functions
   *##################################################################################*/

  /**
   * @brief Get root node with shared lock
   * @return root
   */
  [[nodiscard]] auto
  GetRootForRead()  //
      -> Node_t *
  {
    mutex_.LockShared();         // tree-latch
    root_->AcquireSharedLock();  // root-latch
    mutex_.UnlockShared();       // tree-unlatch
    return root_;
  }

  /**
   * @brief Get root node with exclusive lock
   * @return root
   */
  [[nodiscard]] auto
  GetRootForWrite()  //
      -> Node_t *
  {
    mutex_.Lock();  // tree-latch
    has_tree_lock_ = true;
    root_->AcquireExclusiveLock();  // root-latch
    return root_;
  }

  /**
   * @brief Search a leaf node that may have a target key.
   *
   * @param key a search key.
   * @param ops_is_del a flag for indicating a target operation.
   * @return a stack of traversed nodes.
   */
  [[nodiscard]] auto
  SearchLeafNodeForWrite(  //
      const Key &key,
      const bool ops_is_del)  //
      -> NodeStack
  {
    NodeStack stack{};
    stack.reserve(kExpectedTreeHeight);

    auto *node = GetRootForWrite();
    stack.emplace_back(node, 0);
    while (!node->IsLeaf()) {
      const auto pos = node->SearchChild(key, kClosed);

      bool keep_lock{};
      std::tie(node, keep_lock) = node->GetChildForWrite(pos, ops_is_del);
      if (!keep_lock) {
        ReleaseExclusiveLocks(stack);
      }

      stack.emplace_back(node, pos);
    }

    return stack;
  }

  /**
   * @brief Search a leaf node that may have a target key.
   *
   * @param key a search key.
   * @param closed a flag for indicating closed/open-interval.
   * @return a stack of traversed nodes.
   */
  [[nodiscard]] auto
  SearchLeafNodeForRead(  //
      const Key &key,
      const bool range_is_closed)  //
      -> Node_t *
  {
    auto *node = GetRootForRead();
    while (!node->IsLeaf()) {
      const auto pos = node->SearchChild(key, range_is_closed);
      node = node->GetChildForRead(pos);
    }

    return node;
  }

  void
  ReleaseExclusiveLocks(NodeStack &stack)
  {
    if (has_tree_lock_) {
      mutex_.Unlock();  // unlatch tree-latch
      has_tree_lock_ = false;
    }

    for (auto [node, pos] : stack) {
      node->ReleaseExclusiveLock();
    }
    stack.clear();
  }

  /**
   * @brief Search a leftmost leaf node.
   *
   * @return a stack of traversed nodes.
   */
  [[nodiscard]] auto
  SearchLeftmostLeaf()  //
      -> Node_t *
  {
    auto *node = GetRootForRead();
    while (!node->IsLeaf()) {
      node = node->GetChildForRead(0);
    }

    return node;
  }

  /**
   * @brief Split a given node
   *
   * @tparam Value a payload type of given node
   * @param key a search key.
   * @param stack a stack of nodes in the path to given node
   */
  template <class Value>
  void
  Split(  //
      const Key &key,
      NodeStack &stack)
  {
    // perform splitting for a current node
    auto [l_node, pos] = stack.back();
    stack.pop_back();

    auto *r_node = new Node_t{l_node->IsLeaf()};
    r_node->AcquireExclusiveLock();
    l_node->template Split<Value>(r_node);

    if (stack.empty()) {
      // if node is root, create a new root
      root_ = new Node_t{l_node, r_node};
    } else {
      auto *parent = stack.back().first;
      const auto rc = parent->InsertChild(l_node, r_node, pos);
      if (rc != NodeRC::kCompleted) {
        if (rc == NodeRC::kNeedSplit) {
          Split<Node_t *>(key, stack);
          parent = parent->GetValidSplitNode(key);
        } else {  // NodeRC::kNeedConsolidation
          stack.pop_back();
          ReleaseExclusiveLocks(stack);  // unlock ancestor nodes
          parent->template Consolidate<Node_t *>();
        }

        // update the child position and insert again
        pos = parent->SearchChild(key, kClosed);
        parent->InsertChild(l_node, r_node, pos);
        stack.emplace_back(parent, 0);  // add a parent to release its lock
      }
    }
  }

  /**
   * @brief Merge given nodes
   *
   * @param stack a stack of nodes in the path to given node
   */
  template <class Value>
  void
  Merge(NodeStack &stack)
  {
    auto [node, pos] = stack.back();
    stack.pop_back();

    if (stack.empty()) {
      // a root node cannot be merged
      if (!root_->IsLeaf() && node->GetRecordCount() == 1) {
        // if a root node has only one child, shrink a tree
        root_ = node->template GetPayload<Node_t *>(0);
        delete node;
      } else {
        node->ReleaseExclusiveLock();
      }
      return;
    }

    auto *parent = stack.back().first;

    // check there is a right-sibling node
    if (pos == parent->GetRecordCount() - 1) {
      node->ReleaseExclusiveLock();
      return;
    }

    // check the right-sibling node has enough capacity for merging
    auto *right_node = parent->GetChildForWrite(pos + 1, kDelOps).first;
    if (!node->CanMerge(right_node)) return;

    // perform merging
    node->template Merge<Value>(right_node);
    const auto rc = parent->DeleteChild(node, pos + 1);
    delete right_node;

    if (rc == NodeRC::kNeedMerge) {
      // perform merging recursively
      Merge<Node_t *>(stack);
    }
  }

  /*####################################################################################
   * Internal member variables
   *##################################################################################*/

  /// root node of B+tree
  Node_t *root_ = new Node_t{true};

  /// mutex for root
  Lock mutex_{};

  /// thread local flags for managing the tree lock
  static inline thread_local bool has_tree_lock_{false};  // NOLINT
};
}  // namespace dbgroup::index::b_tree

#endif  // B_TREE_B_TREE_HPP
