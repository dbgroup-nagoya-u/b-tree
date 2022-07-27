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

// local sources
#include "component/node.hpp"
#include "component/record_iterator.hpp"

namespace dbgroup::index::b_tree
{
/**
 * @brief A class for representing B+trees with pessimistic coarse-grained locking.
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

  // aliases for bulkloading
  template <class Entry>
  using BulkIter = typename std::vector<Entry>::const_iterator;
  using BulkResult = std::pair<size_t, std::vector<Node_t *>>;
  using BulkPromise = std::promise<BulkResult>;
  using BulkFuture = std::future<BulkResult>;

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

  /**
   * @brief Bulkload specified kay/payload pairs.
   *
   * This function bulkloads given entries into this index. The entries are assumed to
   * be given as a vector of pairs of Key and Payload (or key/payload/key-length for
   * variable-length keys). Note that keys in records are assumed to be unique and
   * sorted.
   *
   * @param entries vector of entries to be bulkloaded.
   * @param thread_num the number of threads to perform bulkloading.
   * @return kSuccess.
   */
  template <class Entry>
  auto
  Bulkload(  //
      std::vector<Entry> &entries,
      const size_t thread_num = 1)  //
      -> ReturnCode
  {
    assert(thread_num > 0);
    assert(entries.size() >= thread_num);

    if (entries.empty()) return kSuccess;

    std::vector<Node_t *> nodes{};
    auto &&iter = entries.cbegin();
    if (thread_num == 1) {
      // bulkloading with a single thread
      nodes = BulkloadWithSingleThread<Entry>(iter, entries.size()).second;
    } else {
      // bulkloading with multi-threads
      std::vector<BulkFuture> futures{};
      futures.reserve(thread_num);

      // prepare a lambda function for bulkloading
      auto loader = [&](BulkPromise p, BulkIter<Entry> iter, size_t n, bool is_rightmost) {
        p.set_value(BulkloadWithSingleThread<Entry>(iter, n, is_rightmost));
      };

      // create threads to construct partial BTrees
      const size_t rec_num = entries.size();
      const auto rightmost_id = thread_num - 1;
      for (size_t i = 0; i < thread_num; ++i) {
        // create a partial BTree
        BulkPromise p{};
        futures.emplace_back(p.get_future());
        const size_t n = (rec_num + i) / thread_num;
        std::thread{loader, std::move(p), iter, n, i == rightmost_id}.detach();

        // forward the iterator to the next begin position
        iter += n;
      }

      // wait for the worker threads to create BTrees
      std::vector<BulkResult> partial_trees{};
      partial_trees.reserve(thread_num);
      size_t height = 1;
      for (auto &&future : futures) {
        partial_trees.emplace_back(future.get());
        const auto partial_height = partial_trees.back().first;
        height = (partial_height > height) ? partial_height : height;
      }

      // align the height of partial trees
      nodes.reserve(kInnerNodeCap * thread_num);
      for (auto &&[p_height, p_nodes] : partial_trees) {
        while (p_height < height) {  // NOLINT
          ConstructUpperLayer(p_nodes);
          ++p_height;
        }
        nodes.insert(nodes.end(), p_nodes.begin(), p_nodes.end());
      }
    }

    // create upper layers until a root node is created
    while (nodes.size() > 1) {
      ConstructUpperLayer(nodes);
    }
    root_ = nodes.front();

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

  /// the expected capacity of leaf nodes for bulkloading.
  static constexpr size_t kLeafNodeCap = kExpMinBlockSize / (kMetaLen + sizeof(Key) + kPayLen);

  /// the expected capacity of internal nodes for bulkloading.
  static constexpr size_t kInnerNodeCap = kExpMinBlockSize / (kMetaLen + sizeof(Key) + kPtrLen);

  /*####################################################################################
   * Static assertions
   *##################################################################################*/

  // target keys must be trivially copyable.
  static_assert(IsValidKeyType());

  // target payloads must be trivially copyable.
  static_assert(std::is_trivially_copyable_v<Payload>);

  // Each node must have space for at least two records.
  static_assert(2 * kExpMaxRecLen <= kExpMinBlockSize);

  /*####################################################################################
   * Internal utility functions
   *##################################################################################*/

  /**
   * @retval true if a target key class is trivially copyable.
   * @retval false otherwise.
   */
  [[nodiscard]] static constexpr auto
  IsValidKeyType()  //
      -> bool
  {
    if constexpr (IsVarLenData<Key>()) {
      return std::is_trivially_copyable_v<std::remove_pointer_t<Key>>;
    } else {
      return std::is_trivially_copyable_v<Key>;
    }
  }

  /**
   * @return a root node with a shared lock.
   */
  [[nodiscard]] auto
  GetRootForRead()  //
      -> Node_t *
  {
    mutex_.LockS();  // tree-latch

    auto *node = root_;
    node->LockS();  // root-latch

    mutex_.UnlockS();  // tree-unlatch
    return node;
  }

  /**
   * @brief Get a root node with a shared lock with an intent-exclusive lock.
   *
   * This function performs SMOs for a root node if required.
   *
   * @return a root node or a valid child node if a root node is split.
   */
  [[nodiscard]] auto
  GetRootForWrite(const Key &key)  //
      -> Node_t *
  {
    mutex_.LockSIX();  // tree-latch
    root_->LockSIX();  // root-latch

    // check this tree requires SMOs
    ShrinkTreeIfNeeded();
    auto *node = (root_->NeedSplit(kExpMaxRecLen)) ? RootSplit(key) : root_;

    mutex_.UnlockSIX();
    return node;
  }

  /**
   * @brief Search a leaf node that may have a target key.
   *
   * This function uses only shared locks, and a returned leaf node is also locked with
   * a shared lock.
   *
   * @param key a search key.
   * @param is_closed a flag for indicating closed/open-interval.
   * @return a leaf node that may have a target key.
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
   * This function uses only shared locks, and a returned leaf node is also locked with
   * a shared lock.
   *
   * @return a leftmost leaf node.
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
   * This function performs SMOs if internal nodes on the way should be modified. This
   * function uses SIX and X locks for modifying trees, and a returned leaf node is
   * locked by an exclusive lock.
   *
   * @param key a search key.
   * @return a leaf node that may have a target key.
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
        node->UnlockSIX();
      }

      // go down to the next level
      node = child;
    }

    node->UpgradeToX();
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
   * @brief Split a given node.
   *
   * @param l_node a node to be split.
   * @param parent a parent node of `l_node`.
   * @param pos the position of `l_node` in its parent node.
   */
  void
  Split(  //
      Node_t *l_node,
      Node_t *parent,
      const size_t pos)
  {
    parent->UpgradeToX();
    auto *r_node = new Node_t{l_node->IsLeaf()};
    l_node->Split(r_node);
    parent->InsertChild(l_node, r_node, pos);
  }

  /**
   * @brief Split a root node and install a new root node into this tree.
   *
   * @param key a search key.
   * @return one of the split child nodes that includes a given key.
   */
  [[nodiscard]] auto
  RootSplit(const Key &key)  //
      -> Node_t *
  {
    mutex_.UpgradeToX();

    auto *l_node = root_;
    auto *r_node = new Node_t{l_node->IsLeaf()};
    l_node->Split(r_node);
    root_ = new Node_t{l_node, r_node};

    mutex_.DowngradeToSIX();
    return l_node->GetValidSplitNode(key);
  }

  /**
   * @brief Merge a given node with its right-sibling node if possible.
   *
   * @param l_node a node to be merged.
   * @param parent a parent node of `l_node`.
   * @param pos the position of `l_node` in its parent node.
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
    parent->UpgradeToX();
    l_node->Merge(r_node);
    parent->DeleteChild(l_node, l_pos);
    delete r_node;
  }

  /**
   * @brief Remove a root node that has only one child node.
   *
   */
  void
  ShrinkTreeIfNeeded()
  {
    while (root_->GetRecordCount() == 1 && !root_->IsLeaf()) {
      mutex_.UpgradeToX();

      // if a root node has only one child, shrink a tree
      auto *child = root_->template GetPayload<Node_t *>(0);
      child->LockSIX();
      delete root_;
      root_ = child;

      mutex_.DowngradeToSIX();
    }
  }

  /*####################################################################################
   * Internal bulkload utilities
   *##################################################################################*/

  /**
   * @brief Bulkload specified kay/payload pairs with a single thread.
   *
   * Note that this function does not create a root node. The main process must create a
   * root node by using the nodes constructed by this function.
   *
   * @param iter the begin position of target records.
   * @param n the number of entries to be bulkloaded.
   * @param is_rightmost a flag for indicating a constructed tree is rightmost one.
   * @retval 1st: the height of a constructed tree.
   * @retval 2nd: constructed nodes in the top layer.
   */
  template <class Entry>
  auto
  BulkloadWithSingleThread(  //
      BulkIter<Entry> &iter,
      const size_t n,
      const bool is_rightmost = true)  //
      -> BulkResult
  {
    // reserve space for nodes in the leaf layer
    std::vector<Node_t *> nodes{};
    nodes.reserve((n / kLeafNodeCap) + 1);

    // load records into a leaf node
    const auto &iter_end = iter + n;
    Node_t *l_node = nullptr;
    while (iter < iter_end) {
      auto *node = new Node_t{true};
      node->template Bulkload<Entry, Payload>(iter, iter_end, is_rightmost, l_node);
      nodes.emplace_back(node);
      l_node = node;
    }

    // construct index layers
    size_t height = 1;
    if (nodes.size() > 1) {
      // continue until the number of internal nodes is sufficiently small
      do {
        ++height;
      } while (ConstructUpperLayer(nodes));
    }

    return {height, std::move(nodes)};
  }

  /**
   * @brief Construct internal nodes based on given child nodes.
   *
   * @param child_nodes child nodes in a lower layer.
   * @retval true if constructed nodes cannot be contained in a single node.
   * @retval false otherwise.
   */
  auto
  ConstructUpperLayer(std::vector<Node_t *> &child_nodes)  //
      -> bool
  {
    // reserve space for nodes in the upper layer
    std::vector<Node_t *> nodes{};
    nodes.reserve((child_nodes.size() / kInnerNodeCap) + 1);

    // load child nodes into a inner node
    auto &&iter = child_nodes.cbegin();
    const auto &iter_end = child_nodes.cend();
    while (iter < iter_end) {
      auto *node = new Node_t{false};
      node->Bulkload(iter, iter_end);
      nodes.emplace_back(node);
    }

    child_nodes = std::move(nodes);
    return child_nodes.size() > kInnerNodeCap;
  }

  /*####################################################################################
   * Internal member variables
   *##################################################################################*/

  /// a root node of this tree.
  Node_t *root_ = new Node_t{true};

  /// mutex for managing a tree lock.
  ::dbgroup::lock::PessimisticLock mutex_{};
};
}  // namespace dbgroup::index::b_tree

#endif  // B_TREE_B_TREE_HPP
