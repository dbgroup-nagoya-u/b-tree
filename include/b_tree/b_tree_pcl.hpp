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

#include <future>
#include <optional>
#include <thread>
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
template <class Key, class Payload, class Comp = ::std::less<Key>>
class BTreePCL
{
 public:
  /*####################################################################################
   * Type aliases
   *##################################################################################*/

  using Node_t = component::PessimisticNode<Key, Comp>;
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
        : node_{node}, end_pos_{count}, pos_{pos}, end_key_{std::move(end_key)}, is_end_{is_end}
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
      return node_->template GetRecord<Payload>(pos_);
    }

    constexpr void
    operator++()
    {
      ++pos_;
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
        while (pos_ < end_pos_ && node_->GetMetadata(pos_).is_deleted) {
          ++pos_;
        }
        if (pos_ < end_pos_) return true;
        if (is_end_) {
          node_->ReleaseSharedLock();
          return false;
        }
        node_ = node_->GetNextNodeForRead();
        pos_ = 0;
        std::tie(is_end_, end_pos_) = node_->SearchEndPositionFor(end_key_);
      }
    }

    /**
     * @return a Key of a current record
     */
    [[nodiscard]] auto
    GetKey() const  //
        -> Key
    {
      return node_->GetKey(pos_);
    }

    /**
     * @return a payload of a current record
     */
    [[nodiscard]] auto
    GetPayload() const  //
        -> Payload
    {
      return node_->template GetPayload<Payload>(pos_);
    }

   private:
    /// the pointer to a node that includes partial scan results.
    Node_t *node_{nullptr};

    /// the number of records in this node.
    size_t end_pos_{0};

    /// the position of a current record.
    size_t pos_{0};

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
      const Payload &payload,
      [[maybe_unused]] const size_t key_len = sizeof(Key))  //
      -> ReturnCode
  {
    auto *node = SearchLeafNodeForWrite(key, !kDelOps);
    const auto rc = node->template Update<Payload>(key, payload);
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
  Delete(  //
      const Key &key,
      [[maybe_unused]] const size_t key_len = sizeof(Key))  //
      -> ReturnCode
  {
    auto *node = SearchLeafNodeForWrite(key, kDelOps);
    const auto rc = node->Delete(key);
    return (rc == NodeRC::kKeyNotInserted) ? kKeyNotExist : kSuccess;
  }

  template <class Entry>
  auto
  Bulkload(  //
      std::vector<Entry> &entries,
      const size_t thread_num = 1)  //
      -> ReturnCode
  {
    assert(thread_num > 0);

    if (entries.empty()) return kSuccess;

    Node_t *new_root;
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
  GetRootForWrite(  //
      const Key &key,
      const bool ops_is_del)  //
      -> Node_t *
  {
    mutex_.Lock();                  // tree-latch
    root_->AcquireExclusiveLock();  // root-latch

    // check if the root has sufficient space
    Node_t *node{};
    if (root_->HasSufficientSpace(ops_is_del)) {
      node = root_;
    } else if (ops_is_del) {
      node = TryShrinkTree();
    } else {
      node = RootSplit(key);
    }

    mutex_.Unlock();

    return node;
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
    auto *node = GetRootForWrite(key, ops_is_del);
    while (!node->IsLeaf()) {
      const auto pos = node->SearchChild(key, kClosed);
      auto [child, should_smo] = node->GetChildForWrite(pos, ops_is_del);

      // SMO eagerly
      if (should_smo) {
        if (ops_is_del) {
          Merge(child, node, pos);
        } else {
          Split(child, node, pos);
          child = child->GetValidSplitNode(key);
        }
      }
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

  [[nodiscard]] auto
  TryShrinkTree() -> Node_t *
  {
    // a root node cannot be merged
    while (!root_->IsLeaf() && root_->GetRecordCount() == 1) {
      // if a root node has only one child, shrink a tree
      auto *child = root_->template GetPayload<Node_t *>(0);
      child->AcquireExclusiveLock();
      delete root_;
      root_ = child;
    }
    return root_;
  }

  [[nodiscard]] auto
  RootSplit(const Key &key) -> Node_t *
  {
    auto *l_node = root_;
    auto *r_node = new Node_t{l_node->IsLeaf()};
    if (l_node->IsLeaf()) {
      l_node->template Split<Payload>(r_node);
    } else {
      l_node->template Split<Node_t *>(r_node);
    }
    root_ = new Node_t{l_node, r_node};
    auto *node = l_node->GetValidSplitNode(key);
    return node;
  }

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
  void
  Split(  //
      Node_t *l_node,
      Node_t *parent,
      const size_t pos)
  {
    auto *r_node = new Node_t{l_node->IsLeaf()};
    if (l_node->IsLeaf()) {
      l_node->template Split<Payload>(r_node);
    } else {
      l_node->template Split<Node_t *>(r_node);
    }
    parent->InsertChild(l_node, r_node, pos);
  }

  /**
   * @brief Merge given nodes
   *
   * @param node a node to be merged
   * @param parent a parent node of the node
   * @param pos the position of the node
   */
  void
  Merge(  //
      Node_t *node,
      Node_t *parent,
      const size_t pos)
  {
    // check there is a right-sibling node
    if (pos == parent->GetRecordCount() - 1) {
      parent->ReleaseExclusiveLock();
      return;
    }

    // check the right-sibling node has enough capacity for merging
    auto *r_node = parent->GetChildForWrite(pos + 1, kDelOps, false).first;
    if (!node->CanMerge(r_node)) {
      parent->ReleaseExclusiveLock();
      return;
    }

    // perform merging
    if (node->IsLeaf()) {
      node->template Merge<Payload>(r_node);
    } else {
      node->template Merge<Node_t *>(r_node);
    }
    parent->DeleteChild(node, pos + 1);
    delete r_node;
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
