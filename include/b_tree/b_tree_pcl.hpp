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
    auto *node = SearchLeafNodeForWrite(key, !kDelOps);
    node->template Write<Payload>(key, key_len, payload);
    UnlockTree();
    node->ReleaseExclusiveLock();
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
    auto *node = SearchLeafNodeForWrite(key, !kDelOps);
    auto rc = node->template Insert<Payload>(key, key_len, payload);
    UnlockTree();
    node->ReleaseExclusiveLock();
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
    auto *node = SearchLeafNodeForWrite(key, !kDelOps);
    const auto rc = node->template Update<Payload>(key, payload);
    UnlockTree();
    node->ReleaseExclusiveLock();
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
    auto *node = SearchLeafNodeForWrite(key, kDelOps);
    const auto rc = node->Delete(key);
    UnlockTree();
    node->ReleaseExclusiveLock();
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
   * @return a target leaf node.
   */
  [[nodiscard]] auto
  SearchLeafNodeForWrite(  //
      const Key &key,
      const bool ops_is_del)  //
      -> Node_t *
  {
    auto *node = GetRootForWrite();

    // check if the root has sufficient space
    if (!node->HasSufficientSpace(ops_is_del)) {
      auto child = node;
      const size_t pos = 0;
      if (ops_is_del) {
        if (!node->IsLeaf()) {
          Merge<Node_t *>(child, node, pos);
          node = root_;
        }
      } else {
        if (node->IsLeaf()) {
          Split<Payload>(child, node, pos);
        } else {
          Split<Node_t *>(child, node, pos);
        }
        child = child->GetValidSplitNode(key);
        node = child;
      }
    }

    while (!node->IsLeaf()) {
      const auto pos = node->SearchChild(key, kClosed);

      bool should_smo{};
      Node_t *child;
      std::tie(child, should_smo) = node->GetChildForWrite(pos, ops_is_del);

      // SMO eagerly
      if (should_smo) {
        if (ops_is_del) {
          if (child->IsLeaf()) {
            Merge<Payload>(child, node, pos);
          } else {
            Merge<Node_t *>(child, node, pos);
          }
        } else {
          if (child->IsLeaf()) {
            Split<Payload>(child, node, pos);
          } else {
            Split<Node_t *>(child, node, pos);
          }
          child = child->GetValidSplitNode(key);
        }
      }
      UnlockTree();
      node->ReleaseExclusiveLock();
      node = child;
    }

    return node;
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
  UnlockTree()
  {
    if (has_tree_lock_) {
      mutex_.Unlock();
      has_tree_lock_ = false;
    }
    return;
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
   * @param l_node a node to be split
   * @param parent a parent node of the l_node
   * @param pos the position of the l_node
   */
  template <class Value>
  void
  Split(  //
      Node_t *l_node,
      Node_t *parent,
      const size_t pos)
  {
    auto *r_node = new Node_t{l_node->IsLeaf()};
    r_node->AcquireExclusiveLock();
    l_node->template Split<Value>(r_node);

    if (l_node == parent) {
      // if node is root, create a new root
      root_ = new Node_t{l_node, r_node};
      UnlockTree();
    } else {
      parent->InsertChild(l_node, r_node, pos);
    }
  }

  /**
   * @brief Merge given nodes
   *
   * @param node a node to be merged
   * @param parent a parent node of the node
   * @param pos the position of the node
   */
  template <class Value>
  void
  Merge(  //
      Node_t *node,
      Node_t *parent,
      const size_t pos)
  {
    if (node == parent) {
      // a root node cannot be merged
      if (!root_->IsLeaf() && node->GetRecordCount() == 1) {
        // if a root node has only one child, shrink a tree
        auto child = node->template GetPayload<Node_t *>(0);
        child->AcquireExclusiveLock();
        root_ = child;
        delete node;
      }
      return;
    }

    // check there is a right-sibling node
    if (pos == parent->GetRecordCount() - 1) return;

    // check the right-sibling node has enough capacity for merging
    auto *right_node = parent->GetChildForWrite(pos + 1, kDelOps).first;
    if (!node->CanMerge(right_node)) return;

    // perform merging
    node->template Merge<Value>(right_node);
    parent->DeleteChild(node, pos + 1);
    delete right_node;
  }

  /*####################################################################################
   * Internal member variables
   *##################################################################################*/

  /// root node of B+tree
  Node_t *root_ = new Node_t{true};

  /// mutex for root
  ::dbgroup::lock::PessimisticLock mutex_{};

  /// thread local flags for managing the tree lock
  static inline thread_local bool has_tree_lock_{false};  // NOLINT
};
}  // namespace dbgroup::index::b_tree

#endif  // B_TREE_B_TREE_HPP
