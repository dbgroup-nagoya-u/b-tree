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

#ifndef B_TREE_B_TREE_HPP
#define B_TREE_B_TREE_HPP

#include <future>
#include <optional>
#include <thread>
#include <vector>

#include "component/node.hpp"
#include "component/record_iterator.hpp"

namespace dbgroup::index::b_tree
{
/**
 * @brief A class for representing Btree with variable-length keys.
 *
 * This implementation can store variable-length keys (i.e., 'text' type in PostgreSQL).
 *
 * @tparam Key a class of stored keys.
 * @tparam Payload a class of stored payloads (only fixed-length data for simplicity).
 * @tparam Comp a class for ordering keys.
 */
template <class Key, class Payload, class Comp = ::std::less<Key>>
class BTreePCL
{
 public:
  /*####################################################################################
   * Type aliases
   *##################################################################################*/

  using K = Key;
  using V = Payload;
  using Node_t = component::PessimisticNode<Key, Comp>;
  using BTreePCL_t = BTreePCL<Key, Payload, Comp>;
  using RecordIterator_t = component::RecordIterator<BTreePCL_t>;
  using NodeRC = component::NodeRC;

  /*####################################################################################
   * Public constructors and assignment operators
   *##################################################################################*/

  /**
   * @brief Construct a new BTreePCL object.
   *
   */
  BTreePCL() = default;

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
  ~BTreePCL()  //
  {
    DeleteChildren(root_);
  }

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

    if (rc == kSuccess) return payload;
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
      -> RecordIterator_t
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
    return RecordIterator_t{node, begin_pos, end_pos, end_key, is_end};
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
    auto *node = SearchLeafNodeForWrite(key);
    node->Write(key, key_len, payload);
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
    auto *node = SearchLeafNodeForWrite(key);
    return node->Insert(key, key_len, payload);
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
      const Payload &payload,
      [[maybe_unused]] const size_t key_len = sizeof(Key))  //
      -> ReturnCode
  {
    auto *node = SearchLeafNodeForWrite(key);
    return node->Update(key, payload);
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
  Delete(  //
      const Key &key,
      [[maybe_unused]] const size_t key_len = sizeof(Key))  //
      -> ReturnCode
  {
    auto *node = SearchLeafNodeForWrite(key);
    return node->Delete(key);
  }

  /*####################################################################################
   * Public bulkload API
   *##################################################################################*/

  template <class Entry>
  auto
  Bulkload(  //
      std::vector<Entry> &entries,
      const size_t thread_num = 1)  //
      -> ReturnCode
  {
    assert(thread_num > 0);

    if (entries.empty()) return kSuccess;

    Node_t *new_root{};
    auto &&iter = entries.cbegin();
    if (thread_num == 1) {
      // bulkloading with a single thread
      new_root = BulkloadWithSingleThread<Entry>(iter, entries.cend());
    } else {
      // bulkloading with multi-threads
      std::vector<std::future<Node_t *>> threads{};
      threads.reserve(thread_num);

      // prepare a lambda function for bulkloading
      auto loader = [&](std::promise<Node_t *> p,  //
                        size_t n,                  //
                        typename std::vector<Entry>::const_iterator iter, bool is_rightmost) {
        auto *partial_root = BulkloadWithSingleThread<Entry>(iter, iter + n, is_rightmost);
        p.set_value(partial_root);
      };

      // create threads to construct partial BTrees
      const size_t rec_num = entries.size();
      const auto rightmost_id = thread_num - 1;
      for (size_t i = 0; i < thread_num; ++i) {
        // create a partial BTree
        std::promise<Node_t *> p{};
        threads.emplace_back(p.get_future());
        const size_t n = (rec_num + i) / thread_num;
        std::thread{loader, std::move(p), n, iter, i == rightmost_id}.detach();

        // forward the iterator to the next begin position
        iter += n;
      }

      // wait for the worker threads to create BTrees
      std::vector<Node_t *> partial_trees{};
      partial_trees.reserve(thread_num);
      for (auto &&future : threads) {
        partial_trees.emplace_back(future.get());
      }
      new_root = ConstructUpperLayer(partial_trees);
    }

    // set a new root
    root_ = new_root;
    return kSuccess;
  }

 private:
  /*####################################################################################
   * Internal constants
   *##################################################################################*/

  /// an expected maximum height of a tree.
  static constexpr size_t kExpectedTreeHeight = 8;

  /// a flag for indicating closed-interval
  static constexpr bool kClosed = true;

  /// a flag for indicating delete operations
  static constexpr bool kDelOps = true;

  /// the length of payloads.
  static constexpr size_t kPayLen = sizeof(Payload);

  /// the length of child pointers.
  static constexpr size_t kPtrLen = sizeof(Node_t *);

  /// the length of metadata.
  static constexpr size_t kMetaLen = sizeof(component::Metadata);

  /// the length of a header in each node page.
  static constexpr size_t kHeaderLen = sizeof(Node_t);

  /// the expected maximum length of keys.
  static constexpr size_t kExpMaxKeyLen = (IsVarLenData<Key>()) ? kMaxVarLenDataSize : sizeof(Key);

  /// the expected maximum length of payloads (including child pointers).
  static constexpr size_t kExpMaxPayLen = (kPayLen < kPtrLen) ? kPtrLen : kPayLen;

  /// the expected maximum length of records.
  static constexpr size_t kExpMaxRecLen = kExpMaxKeyLen + kExpMaxPayLen + kMetaLen;

  /// the expected minimum block size in a certain node.
  static constexpr size_t kExpMinBlockSize = kPageSize - kHeaderLen + kExpMaxKeyLen;

  /*####################################################################################
   * Static assertions
   *##################################################################################*/

  // Each node must have space for at least two records.
  static_assert(2 * kExpMaxRecLen <= kExpMinBlockSize);

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
    mutex_.LockShared();  // tree-latch

    auto *node = root_;
    node->AcquireSharedLock();  // root-latch

    mutex_.UnlockShared();  // tree-unlatch
    return node;
  }

  /**
   * @brief Get root node with exclusive lock
   * @return root
   */
  [[nodiscard]] auto
  GetRootForWrite(const Key &key)  //
      -> Node_t *
  {
    mutex_.Lock();                  // tree-latch
    root_->AcquireExclusiveLock();  // root-latch

    // check this tree requires SMOs
    ShrinkTreeIfNeeded();
    auto *node = (root_->NeedSplit(kExpMaxRecLen)) ? RootSplit(key) : root_;

    mutex_.Unlock();
    return node;
  }

  /**
   * @brief Search a leaf node that may have a target key.
   *
   * @param key a search key.
   * @param is_closed a flag for indicating closed/open-interval.
   * @return a stack of traversed nodes.
   */
  [[nodiscard]] auto
  SearchLeafNodeForRead(  //
      const Key &key,
      const bool is_closed)  //
      -> Node_t *
  {
    auto *node = GetRootForRead();
    while (!node->IsLeaf()) {
      const auto pos = node->SearchChild(key, is_closed);
      node = node->GetChildForRead(pos);
    }

    return node;
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
   * @brief Search a leaf node that may have a target key.
   *
   * @param key a search key.
   * @return a target leaf node.
   */
  [[nodiscard]] auto
  SearchLeafNodeForWrite(const Key &key)  //
      -> Node_t *
  {
    auto *node = GetRootForWrite(key);
    while (!node->IsLeaf()) {
      // search a child node
      const auto pos = node->SearchChild(key, kClosed);
      auto *child = node->GetChildForWrite(pos);

      // perform internal SMOs eagerly
      if (child->NeedSplit(kExpMaxRecLen)) {
        Split(child, node, pos);
        child = child->GetValidSplitNode(key);
      } else if (child->NeedMerge()) {
        Merge(child, node, pos);
      } else {
        node->ReleaseExclusiveLock();
      }

      // go down to the next level
      node = child;
    }

    return node;
  }

  /**
   * @brief Delete child nodes recursively.
   *
   * Note that this function assumes that there are no other threads in operation.
   *
   * @param node a target node.
   */
  static void
  DeleteChildren(Node_t *node)
  {
    if (!node->IsLeaf()) {
      // delete children nodes recursively
      for (size_t i = 0; i < node->GetRecordCount(); ++i) {
        auto *child_node = node->template GetPayload<Node_t *>(i);
        DeleteChildren(child_node);
      }
    }

    delete node;
  }

  /*####################################################################################
   * Internal structure modification operations
   *##################################################################################*/

  /**
   * @brief Split a given node
   *
   * @param l_node a node to be split
   * @param parent a parent node of the l_node
   * @param pos the position of the l_node
   */
  void
  Split(  //
      Node_t *l_node,
      Node_t *parent,
      const size_t pos)
  {
    auto *r_node = new Node_t{l_node->IsLeaf()};
    l_node->Split(r_node);
    parent->InsertChild(l_node, r_node, pos);
  }

  [[nodiscard]] auto
  RootSplit(const Key &key)  //
      -> Node_t *
  {
    auto *l_node = root_;
    auto *r_node = new Node_t{l_node->IsLeaf()};
    l_node->Split(r_node);
    root_ = new Node_t{l_node, r_node};

    return l_node->GetValidSplitNode(key);
  }

  /**
   * @brief Merge a given node with its right-sibling node if possible.
   *
   * @param l_node a node to be merged.
   * @param parent a parent node of `l_node`.
   * @param pos the position of `l_node`.
   */
  void
  Merge(  //
      Node_t *l_node,
      Node_t *parent,
      const size_t l_pos)
  {
    // check there is a mergeable sibling node
    auto *r_node = parent->GetMergeableRightChild(l_node, l_pos);
    if (r_node == nullptr) return;

    // perform merging
    l_node->Merge(r_node);
    parent->DeleteChild(l_node, l_pos);
    delete r_node;
  }

  void
  ShrinkTreeIfNeeded()
  {
    while (root_->GetRecordCount() == 1 && !root_->IsLeaf()) {
      // if a root node has only one child, shrink a tree
      auto *child = root_->template GetPayload<Node_t *>(0);
      child->AcquireExclusiveLock();
      delete root_;
      root_ = child;
    }
  }

  /*####################################################################################
   * Internal bulkload utilities
   *##################################################################################*/

  /**
   * @brief Bulkload specified kay/payload pairs with a single thread.
   *
   * @param iter the begin position of target records.
   * @param iter_end the end position of target records.
   * @return the root node of a created BTree.
   */
  template <class Entry>
  auto
  BulkloadWithSingleThread(  //
      typename std::vector<Entry>::const_iterator &iter,
      const typename std::vector<Entry>::const_iterator &iter_end,
      const bool is_rightmost = true)  //
      -> Node_t *
  {
    std::vector<Node_t *> nodes{};
    Node_t *l_node = nullptr;
    while (iter < iter_end) {
      // load records into a leaf node
      auto *node = new Node_t{true};
      node->template Bulkload<Entry, Payload>(iter, iter_end, is_rightmost, l_node);
      nodes.emplace_back(node);
      l_node = node;
    }
    if (nodes.size() == 1) return nodes[0];
    return ConstructUpperLayer(nodes);
  }

  auto
  ConstructUpperLayer(std::vector<Node_t *> &entries)  //
      -> Node_t *
  {
    std::vector<Node_t *> nodes;
    auto &&iter = entries.cbegin();
    const auto &iter_end = entries.cend();
    Node_t *l_node = nullptr;
    while (iter < iter_end) {
      // load records into a inner node
      auto *node = new Node_t{false};
      node->Bulkload(iter, iter_end, l_node);
      nodes.emplace_back(node);
      l_node = node;
    }
    if (nodes.size() == 1) return nodes[0];
    return ConstructUpperLayer(nodes);
  }

  /*####################################################################################
   * Internal member variables
   *##################################################################################*/

  /// root node of B+tree
  Node_t *root_ = new Node_t{true};

  /// mutex for root
  ::dbgroup::lock::PessimisticLock mutex_{};
};
}  // namespace dbgroup::index::b_tree

#endif  // B_TREE_B_TREE_HPP
