/*
 * Copyright 2023 Database Group, Nagoya University
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

#ifndef B_TREE_COMPONENT_OML_B_TREE_HPP
#define B_TREE_COMPONENT_OML_B_TREE_HPP

// C++ standard libraries
#include <future>
#include <optional>
#include <thread>
#include <vector>

// external sources
#include "memory/epoch_based_gc.hpp"

// local sources
#include "b_tree/component/oml/node_fixlen.hpp"
#include "b_tree/component/oml/node_varlen.hpp"
#include "b_tree/component/optimistic_record_iterator.hpp"

namespace dbgroup::index::b_tree::component::oml
{
/**
 * @brief A class for representing B+trees with optimistic multi-layer locking.
 *
 * This implementation can store variable-length keys (i.e., 'text' type in PostgreSQL).
 *
 * @tparam Key a class of stored keys.
 * @tparam Payload a class of stored payloads (only fixed-length data for simplicity).
 * @tparam Comp a class for ordering keys.
 * @tparam kUseVarLenLayout a flag for indicating variable-length keys.
 */
template <class Key, class Payload, class Comp, bool kUseVarLenLayout>
class BTree
{
 public:
  /*####################################################################################
   * Type aliases
   *##################################################################################*/

  using K = Key;
  using V = Payload;
  using NodeVarLen_t = NodeVarLen<Key, Comp>;
  using NodeFixLen_t = NodeFixLen<Key, Comp>;
  using Node_t = std::conditional_t<kUseVarLenLayout, NodeVarLen_t, NodeFixLen_t>;
  using BTree_t = BTree<Key, Payload, Comp, kUseVarLenLayout>;
  using RecordIterator_t = OptimisticRecordIterator<BTree_t>;
  using ScanKey = std::optional<std::tuple<const Key &, size_t, bool>>;
  using GC_t = ::dbgroup::memory::EpochBasedGC<Page>;

  // aliases for bulkloading
  template <class Entry>
  using BulkIter = typename std::vector<Entry>::const_iterator;
  using NodeEntry = std::tuple<Key, Node_t *, size_t>;
  using BulkResult = std::pair<size_t, std::vector<NodeEntry>>;
  using BulkPromise = std::promise<BulkResult>;
  using BulkFuture = std::future<BulkResult>;

  /*####################################################################################
   * Public constructors and assignment operators
   *##################################################################################*/

  /**
   * @brief Construct a new BTree object.
   *
   * @param gc_interval_micro time interval for garbage collection [us].
   * @param gc_thread_num the number of worker threads for garbage collection.
   */
  BTree(  //
      const size_t gc_interval_micro,
      const size_t gc_thread_num)
      : gc_{gc_interval_micro, gc_thread_num}
  {
    auto *root = new (GetNodePage()) Node_t{kLeafFlag};
    if constexpr (!kUseVarLenLayout) {
      root->SetPayloadLength(kPayLen);
    }
    root_.store(root, std::memory_order_release);

    gc_.StartGC();
  }

  BTree(const BTree &) = delete;
  BTree(BTree &&) = delete;

  BTree &operator=(const BTree &) = delete;
  BTree &operator=(BTree &&) = delete;

  /*####################################################################################
   * Public destructors
   *##################################################################################*/

  /**
   * @brief Destroy the BTree object.
   *
   */
  ~BTree()  //
  {
    DeleteChildren(root_.load(std::memory_order_acquire));
  }

  /*####################################################################################
   * Public read APIs
   *##################################################################################*/

  /**
   * @brief The entity of a function for reading records.
   *
   * @param key a target key.
   * @param key_len the length of the target key.
   * @retval the payload of a given key wrapped with std::optional if it is in this tree.
   * @retval std::nullopt otherwise.
   */
  auto
  Read(  //
      const Key &key,
      [[maybe_unused]] const size_t key_len)  //
      -> std::optional<Payload>
  {
    [[maybe_unused]] const auto &guard = gc_.CreateEpochGuard();

    auto *node = SearchLeafNodeForRead(key);
    Payload payload{};
    const auto rc = Node_t::Read(node, key, payload);
    if (rc == kCompleted) return payload;
    return std::nullopt;
  }

  /**
   * @brief The entity of a function for scanning records.
   *
   * @param begin_key a pair of a begin key and its openness (true=closed).
   * @param end_key a pair of an end key and its openness (true=closed).
   * @return an iterator to access scanned records.
   */
  auto
  Scan(  //
      const ScanKey &begin_key = std::nullopt,
      const ScanKey &end_key = std::nullopt)  //
      -> RecordIterator_t
  {
    auto &&guard = gc_.CreateEpochGuard();
    auto *node = (begin_key) ? SearchLeafNodeForRead(std::get<0>(*begin_key))  //
                             : SearchLeftmostLeaf();

    return RecordIterator_t{node, begin_key, end_key, std::move(guard)};
  }

  /*####################################################################################
   * Public write APIs
   *##################################################################################*/

  /**
   * @brief The entity of a function for putting records.
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
    [[maybe_unused]] const auto &guard = gc_.CreateEpochGuard();

    while (true) {
      auto *node = SearchLeafNodeForWrite(key);
      const auto rc = Node_t::CheckKeyRangeAndLockForWrite(node, key, key_len + kPayLen);
      if (rc == kNeedRetry) continue;
      node->Write(key, key_len, &payload, kPayLen);
      return kSuccess;
    }
  }

  /**
   * @brief The entity of a function for inserting records.
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
    [[maybe_unused]] const auto &guard = gc_.CreateEpochGuard();

    while (true) {
      auto *node = SearchLeafNodeForWrite(key);
      const auto rc = Node_t::Insert(node, key, key_len, &payload, kPayLen);
      if (rc == kNeedRetry) continue;
      return (rc == kCompleted) ? kSuccess : kKeyExist;
    }
  }

  /**
   * @brief The entity of a function for updating records.
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
    [[maybe_unused]] const auto &guard = gc_.CreateEpochGuard();

    auto *node = SearchLeafNodeForWrite(key);
    return Node_t::Update(node, key, &payload, kPayLen);
  }

  /**
   * @brief The entity of a function for deleting records.
   *
   * @param key a target key to be deleted.
   * @param key_len the length of a target key.
   * @retval kSuccess if deleted.
   * @retval kKeyNotExist otherwise.
   */
  auto
  Delete(  //
      const Key &key,
      [[maybe_unused]] const size_t key_len = sizeof(Key))  //
      -> ReturnCode
  {
    [[maybe_unused]] const auto &guard = gc_.CreateEpochGuard();

    auto *node = SearchLeafNodeForWrite(key);
    return Node_t::Delete(node, key);
  }

  /*####################################################################################
   * Public bulkload API
   *##################################################################################*/

  /**
   * @brief The entity of a function for bulkinserting records.
   *
   * @tparam Entry a container of a key/payload pair.
   * @param entries vector of entries to be bulkloaded.
   * @param thread_num the number of threads to perform bulkloading.
   * @return kSuccess.
   */
  template <class Entry>
  auto
  Bulkload(  //
      const std::vector<Entry> &entries,
      const size_t thread_num = 1)  //
      -> ReturnCode
  {
    if (entries.empty()) return kSuccess;

    std::vector<NodeEntry> nodes{};
    auto &&iter = entries.cbegin();
    const auto rec_num = entries.size();
    if (thread_num <= 1 || rec_num < thread_num) {
      // bulkloading with a single thread
      nodes = BulkloadWithSingleThread<Entry>(iter, rec_num).second;
    } else {
      // bulkloading with multi-threads
      std::vector<BulkFuture> futures{};
      futures.reserve(thread_num);

      // a lambda function for bulkloading with multi-threads
      auto loader = [&](BulkPromise p, BulkIter<Entry> iter, size_t n) {
        p.set_value(BulkloadWithSingleThread<Entry>(iter, n));
      };

      // create threads to construct partial BzTrees
      for (size_t i = 0; i < thread_num; ++i) {
        // create a partial B+Tree
        BulkPromise p{};
        futures.emplace_back(p.get_future());
        const size_t n = (rec_num + i) / thread_num;
        std::thread{loader, std::move(p), iter, n}.detach();

        // forward the iterator to the next begin position
        iter += n;
      }

      // wait for the worker threads to create partial trees
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
      Node_t *prev_node = nullptr;
      for (auto &&[p_height, p_nodes] : partial_trees) {
        while (p_height < height) {  // NOLINT
          p_nodes = ConstructSingleLayer<NodeEntry>(p_nodes.cbegin(), p_nodes.size());
          ++p_height;
        }
        nodes.insert(nodes.end(), p_nodes.begin(), p_nodes.end());

        // link partial trees
        if (prev_node != nullptr) {
          Node_t::LinkVerticalBorderNodes(prev_node, std::get<1>(p_nodes.front()));
        }
        prev_node = std::get<1>(p_nodes.back());
      }
    }

    // create upper layers until a root node is created
    while (nodes.size() > 1) {
      nodes = ConstructSingleLayer<NodeEntry>(nodes.cbegin(), nodes.size());
    }
    root_ = std::get<1>(nodes.front());
    Node_t::RemoveLeftmostKeys(root_);

    return kSuccess;
  }

  /**
   * @brief Collect statistical data of this tree.
   *
   * @retval 1st: the number of nodes.
   * @retval 2nd: the actual usage in bytes.
   * @retval 3rd: the virtual usage in bytes.
   */
  auto
  CollectStatisticalData()  //
      -> std::vector<std::tuple<size_t, size_t, size_t>>
  {
    std::vector<std::tuple<size_t, size_t, size_t>> stat_data{};
    auto *node = root_.load(std::memory_order_acquire);

    CollectStatisticalData(node, 0, stat_data);

    return stat_data;
  }

 private:
  /*####################################################################################
   * Internal constants
   *##################################################################################*/

  /// the length of payloads.
  static constexpr size_t kPayLen = sizeof(Payload);

  /// the length of child pointers.
  static constexpr size_t kPtrLen = sizeof(Node_t *);

  /// the length of metadata.
  static constexpr size_t kMetaLen = sizeof(Metadata);

  /// the length of a header in each node page.
  static constexpr size_t kHeaderLen = sizeof(Node_t);

  /// the maximum length of keys.
  static constexpr size_t kMaxKeyLen = (IsVarLenData<Key>()) ? kMaxVarLenDataSize : sizeof(Key);

  /// the maximum length of payloads (including child pointers).
  static constexpr size_t kMaxPayLen = (kPayLen < kPtrLen) ? kPtrLen : kPayLen;

  /// the maximum length of records.
  static constexpr size_t kMaxRecLen =
      kMaxKeyLen + kMaxPayLen + ((kUseVarLenLayout) ? kMetaLen : 0);

  /// the minimum block size in a certain node.
  static constexpr size_t kMinBlockSize = kPageSize - kHeaderLen - kMaxKeyLen;

  /// the expected length of leaf records.
  static constexpr size_t kExpLeafRecLen =
      sizeof(Key) + kPayLen + ((kUseVarLenLayout) ? kMetaLen : 0);

  /// the expected length of internal records.
  static constexpr size_t kExpInnerRecLen =
      sizeof(Key) + kPtrLen + ((kUseVarLenLayout) ? kMetaLen : 0);

  /// the expected capacity of leaf nodes for bulkloading.
  static constexpr size_t kLeafNodeCap = (kMinBlockSize - kMinFreeSpaceSize) / kExpLeafRecLen;

  /// the expected capacity of internal nodes for bulkloading.
  static constexpr size_t kInnerNodeCap = (kMinBlockSize - kMinFreeSpaceSize) / kExpInnerRecLen;

  /*####################################################################################
   * Internal utility functions
   *##################################################################################*/

  /**
   * @brief Allocate or reuse a memory region for a base node.
   *
   * @returns the reserved memory page.
   */
  [[nodiscard]] auto
  GetNodePage()  //
      -> void *
  {
    auto *page = gc_.template GetPageIfPossible<Page>();
    return (page == nullptr) ? (::dbgroup::memory::Allocate<Page>()) : page;
  }

  /**
   * @brief Get a root node with a shared lock with an intent-exclusive lock.
   *
   * This function performs SMOs for a root node if required.
   *
   * @param key a search key.
   * @retval 1st: a root node or a valid child node if a root node is split.
   * @return 2nd: version value of the node.
   */
  [[nodiscard]] auto
  GetRootForWrite(const Key &key)  //
      -> std::pair<Node_t *, uint64_t>
  {
    while (true) {
      auto *node = root_.load(std::memory_order_acquire);

      // check this tree requires SMOs
      auto [rc, ver] = node->CheckNodeStatus(kMaxKeyLen + kMaxPayLen);
      if (rc == kNeedRetry) continue;  // a root node is removed
      if (rc == kNeedSplit) {
        if (!TryRootSplit(key, node, ver)) continue;
        // a valid child node and its version value are selected in TryRootSplit
      } else if (rc == kNeedMerge) {
        if (!TryShrinkTree(node, ver)) continue;  // retry to get a new root node
      }

      return {node, ver};
    }
  }

  /**
   * @brief Search a leaf node that may have a target key.
   *
   * This function uses only shared locks, and a returned leaf node is also locked with
   * a shared lock.
   *
   * @param key a search key.
   * @return a leaf node that may have a target key.
   */
  [[nodiscard]] auto
  SearchLeafNodeForRead(const Key &key)  //
      -> Node_t *
  {
    auto *node = root_.load(std::memory_order_acquire);
    while (node->IsInner()) {
      auto *child = node->SearchChild(key);
      if (child == nullptr) {
        node = root_.load(std::memory_order_acquire);
        continue;
      }
      node = child;
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
    auto *node = root_.load(std::memory_order_acquire);
    while (node->IsInner()) {
      node = node->GetLeftmostChild();
      if (node == nullptr) {
        node = root_.load(std::memory_order_acquire);
      }
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
    auto [node, ver] = GetRootForWrite(key);
    while (node->IsInner()) {
      auto [pos, child] = node->SearchChild(key, ver, kMaxKeyLen + kMaxPayLen);
      if (child == nullptr) {
        std::tie(node, ver) = GetRootForWrite(key);
        continue;
      }

      // perform internal SMOs eagerly
      auto [rc, child_ver] = child->CheckNodeStatus(kMaxKeyLen + kMaxPayLen);
      if (rc == kNeedRetry) continue;  // the child node was removed
      if (rc == kNeedSplit) {
        if (!TrySplit(key, child, child_ver, node, ver, pos)) continue;
        // a valid child node and its version value are selected in TrySplit
      } else if (rc == kNeedMerge) {
        if (!TryMerge(child, child_ver, node, ver, pos)) continue;
        // a valid version value is selected in TryMerge
      }

      // go down to the next level
      node = child;
      ver = child_ver;
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
    if (node->IsInner()) {
      // delete children nodes recursively
      for (size_t i = 0; i < node->GetRecordCount(); ++i) {
        auto *child_node = node->template GetPayload<Node_t *>(i);
        DeleteChildren(child_node);
      }
    }

    ::dbgroup::memory::Release<Page>(node);
  }

  /**
   * @brief Collect statistical data recursively.
   *
   * @param node a target node.
   * @param level the current level in the tree.
   * @param stat_data an output statistical data.
   */
  static void
  CollectStatisticalData(  //
      Node_t *node,
      const size_t level,
      std::vector<std::tuple<size_t, size_t, size_t>> &stat_data)
  {
    node->LockS();

    // add an element for a new level
    if (stat_data.size() <= level) {
      stat_data.emplace_back(0, 0, 0);
    }

    // add statistical data of this node
    auto &[node_num, actual_usage, virtual_usage] = stat_data.at(level);
    ++node_num;
    actual_usage += node->GetNodeUsage();
    virtual_usage += kPageSize;

    // collect data recursively
    if (node->IsInner()) {
      for (size_t i = 0; i < node->GetRecordCount(); ++i) {
        auto *child = node->template GetPayload<Node_t *>(i);
        CollectStatisticalData(child, level + 1, stat_data);
      }
    }

    node->UnlockS();
  }

  /*####################################################################################
   * Internal structure modification operations
   *##################################################################################*/

  /**
   * @brief Split a given node.
   *
   * @param key a search key.
   * @param child a node to be split.
   * @param c_ver an expected version value of `child`.
   * @param parent a parent node of `child`.
   * @param p_ver an expected version value of `parent`.
   * @param pos the position of `child` in its parent node.
   * @retval true if succeeded.
   * @retval false if failed.
   */
  auto
  TrySplit(  //
      const Key &key,
      Node_t *&child,
      uint64_t &c_ver,
      Node_t *parent,
      uint64_t p_ver,
      const size_t pos)  //
      -> bool
  {
    // try to acquire locks
    if (!parent->TryLockSIX(p_ver)) return false;
    if (!child->TryLockSIX(c_ver)) {
      parent->UnlockSIX();
      return false;
    }

    // perform splitting
    auto *l_node = child;
    auto *r_node = new (GetNodePage()) Node_t{l_node->IsInner()};
    l_node->Split(r_node);
    auto &&[new_child, sep_key, sep_key_len, new_c_ver] = l_node->GetValidSplitNode(key, r_node);
    parent->InsertChild(r_node, sep_key, sep_key_len, pos);

    child = new_child;
    c_ver = new_c_ver;
    return true;
  }

  /**
   * @brief Split a root node and install a new root node into this tree.
   *
   * @param key a search key.
   * @param node a node expected to be the root node.
   * @param ver an expected version value of `node`.
   * @retval true if succeeded.
   * @retval false if failed.
   */
  [[nodiscard]] auto
  TryRootSplit(  //
      const Key &key,
      Node_t *&node,
      uint64_t &ver)  //
      -> bool
  {
    // try to acquire a lock and check the given node is a root node
    if (!node->TryLockSIX(ver)) return false;
    if (node != root_.load(std::memory_order_relaxed)) {
      // other threads modified a root node
      node->UnlockSIX();
      return false;
    }

    // perform splitting the root node
    auto *l_node = node;
    auto *r_node = new (GetNodePage()) Node_t{l_node->IsInner()};
    l_node->Split(r_node);

    // install a new root node
    auto *new_root = new (GetNodePage()) Node_t{l_node, r_node};
    root_.store(new_root, std::memory_order_release);
    auto &&[new_node, sep_key, sep_key_len, new_ver] = l_node->GetValidSplitNode(key, r_node);

    node = new_node;
    ver = new_ver;
    return true;
  }

  /**
   * @brief Merge a given node with its right-sibling node if possible.
   *
   * @param l_node a node to be merged.
   * @param c_ver an expected version value of `l_node`.
   * @param parent a parent node of `l_node`.
   * @param p_ver an expected version value of `parent`.
   * @param l_pos the position of `l_node` in its parent node.
   * @retval true if succeeded.
   * @retval false if failed.
   */
  auto
  TryMerge(  //
      Node_t *l_node,
      uint64_t &c_ver,
      Node_t *parent,
      uint64_t p_ver,
      const size_t l_pos)  //
      -> bool
  {
    // try to acquire locks
    if (!parent->TryLockSIX(p_ver)) return false;
    if (!l_node->TryLockSIX(c_ver)) {
      parent->UnlockSIX();
      return false;
    }

    // check there is a mergeable sibling node
    auto *r_node = parent->GetMergeableRightChild(l_node, l_pos);
    if (r_node == nullptr) return true;

    // perform merging
    c_ver = l_node->Merge(r_node);
    parent->DeleteChild(l_pos);
    gc_.AddGarbage<Page>(r_node);

    return true;
  }

  /**
   * @brief Remove a root node that has only one child node.
   *
   * @param node a node expected to be the root node.
   * @param ver an expected version value of `node`.
   * @retval true if succeeded.
   * @retval false if failed.
   */
  auto
  TryShrinkTree(  //
      Node_t *node,
      uint64_t ver)  //
      -> bool
  {
    if (node->GetRecordCount() > 1 || !node->IsInner()) return true;

    // if a internal-root node has only one child, try to shrink a tree
    if (!node->TryLockSIX(ver)) return false;
    if (node != root_.load(std::memory_order_relaxed)) {
      node->UnlockSIX();
      return false;
    }

    // remove the root node
    gc_.AddGarbage<Page>(node);
    node = node->RemoveRoot();
    root_.store(node, std::memory_order_release);

    return false;
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
   * @tparam Entry a container of a key/payload pair.
   * @param iter the begin position of target records.
   * @param n the number of entries to be bulkloaded.
   * @retval 1st: the height of a constructed tree.
   * @retval 2nd: constructed nodes in the top layer.
   */
  template <class Entry>
  auto
  BulkloadWithSingleThread(  //
      BulkIter<Entry> iter,
      const size_t n)  //
      -> BulkResult
  {
    // construct a data layer (leaf nodes)
    auto &&nodes = ConstructSingleLayer<Entry>(iter, n);

    // construct index layers (inner nodes)
    size_t height = 1;
    for (auto n = nodes.size(); n > kInnerNodeCap; n = nodes.size(), ++height) {
      // continue until the number of inner nodes is sufficiently small
      nodes = ConstructSingleLayer<NodeEntry>(nodes.cbegin(), n);
    }

    return {height, std::move(nodes)};
  }

  /**
   * @brief Construct internal nodes based on given child nodes.
   *
   * @tparam Entry a container of a key/payload pair.
   * @param iter the begin position of target records.
   * @param n the number of entries to be bulkloaded.
   * @return constructed nodes.
   */
  template <class Entry>
  auto
  ConstructSingleLayer(  //
      BulkIter<Entry> iter,
      const size_t n)  //
      -> std::vector<NodeEntry>
  {
    using T = std::tuple_element_t<1, Entry>;
    constexpr auto kIsInner = std::is_same_v<T, Node_t *>;

    // reserve space for nodes in the upper layer
    std::vector<NodeEntry> nodes{};
    nodes.reserve((n / (kIsInner ? kInnerNodeCap : kLeafNodeCap)) + 1);

    // load child nodes into parent nodes
    const auto &iter_end = iter + n;
    for (Node_t *prev_node = nullptr; iter < iter_end;) {
      auto *node = new (GetNodePage()) Node_t{kIsInner};
      if constexpr (!kUseVarLenLayout && !kIsInner) {
        node->SetPayloadLength(kPayLen);
      }
      node->template Bulkload<Entry>(iter, iter_end, prev_node, nodes);
      prev_node = node;
    }

    return nodes;
  }

  /*####################################################################################
   * Static assertions
   *##################################################################################*/

  // Each node must have space for at least two records.
  static_assert(2 * kMaxRecLen <= kMinBlockSize);

  /*####################################################################################
   * Internal member variables
   *##################################################################################*/

  /// a garbage collector for node pages.
  GC_t gc_{};

  /// a root node of this tree.
  std::atomic<Node_t *> root_{nullptr};
};
}  // namespace dbgroup::index::b_tree::component::oml

#endif  // B_TREE_COMPONENT_OML_B_TREE_HPP
