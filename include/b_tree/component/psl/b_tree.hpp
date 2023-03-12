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

#ifndef B_TREE_COMPONENT_PSL_B_TREE_HPP
#define B_TREE_COMPONENT_PSL_B_TREE_HPP

// C++ standard libraries
#include <future>
#include <optional>
#include <thread>
#include <vector>

// external sources
#include "memory/epoch_based_gc.hpp"

// local sources
#include "b_tree/component/psl/node_fixlen.hpp"
#include "b_tree/component/psl/node_varlen.hpp"
#include "b_tree/component/record_iterator.hpp"

namespace dbgroup::index::b_tree::component::psl
{
/**
 * @brief A class for representing B+trees with pessimistic single-layer locking.
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
  using RecordIterator_t = RecordIterator<BTree_t>;
  using ScanKey = std::optional<std::tuple<const Key &, size_t, bool>>;
  using GC_t = ::dbgroup::memory::EpochBasedGC<PageTarget>;

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
    const auto rc = node->Read(key, payload);

    if (rc == kSuccess) return payload;
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

    Node_t *node{};
    size_t begin_pos = 0;

    if (begin_key) {
      const auto &[key, key_len, is_closed] = begin_key.value();
      node = SearchLeafNodeForRead(key);
      const auto [rc, pos] = node->SearchRecord(key);
      begin_pos = (rc == NodeRC::kKeyAlreadyInserted && !is_closed) ? pos + 1 : pos;
    } else {
      node = SearchLeftmostLeaf();
    }

    const auto [is_end, end_pos] = node->SearchEndPositionFor(end_key);
    return RecordIterator_t{node, begin_pos, end_pos, end_key, is_end, std::move(guard)};
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

    auto &&stack = SearchLeafNodeForWrite(key);
    auto *node = stack.back();
    const auto rc = node->Write(key, key_len, &payload, kPayLen);

    if (rc == NodeRC::kNeedSplit) {
      // perform splitting if needed
      auto *r_node = HalfSplit(node);
      auto &&[target_node, sep_key, sep_key_len] = node->GetValidSplitNode(key);
      target_node->Write(key, key_len, &payload, kPayLen);

      // complete splitting by inserting a new entry
      CompleteSplit(stack, node, r_node, sep_key, sep_key_len);
      if constexpr (IsVarLenData<Key>()) {
        ::operator delete(sep_key);
      }
    }

    return kSuccess;
  }

  /**
   * @brief The entity of a function for inserting records.
   *
   * @param key a target key to be inserted.
   * @param payload a target payload to be inserted.
   * @param key_len the length of a target key.
   * @retval.kSuccess if inserted.
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

    auto &&stack = SearchLeafNodeForWrite(key);
    auto *node = stack.back();
    auto rc = node->Insert(key, key_len, &payload, kPayLen);
    if (rc == NodeRC::kKeyAlreadyInserted) return kKeyExist;

    if (rc == NodeRC::kNeedSplit) {
      // perform splitting if needed
      auto *r_node = HalfSplit(node);
      auto &&[target_node, sep_key, sep_key_len] = node->GetValidSplitNode(key);
      target_node->Insert(key, key_len, &payload, kPayLen);

      // complete splitting by inserting a new entry
      CompleteSplit(stack, node, r_node, sep_key, sep_key_len);
      if constexpr (IsVarLenData<Key>()) {
        ::operator delete(sep_key);
      }
    }

    return kSuccess;
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

    auto &&stack = SearchLeafNodeForWrite(key);
    auto *node = stack.back();
    return node->Update(key, &payload, kPayLen);
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

    auto &&stack = SearchLeafNodeForWrite(key);
    auto *node = stack.back();
    const auto rc = node->Delete(key);
    if (rc == NodeRC::kKeyNotInserted) return kKeyNotExist;

    if (rc == NodeRC::kNeedMerge) {
      Merge(stack, node);
    }

    return kSuccess;
  }

  /*####################################################################################
   * Public bulkload API
   *##################################################################################*/

  /**
   * @brief The entity of a function for bulkinserting records.
   *
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
    auto *page = gc_.template GetPageIfPossible<PageTarget>();
    return (page == nullptr) ? (::operator new(kPageSize, component::kCacheAlignVal)) : page;
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
    while (true) {
      // traverse side links if needed
      node = Node_t::CheckKeyRangeAndLockForRead(node, key);
      if (node == nullptr) {
        // a root node is removed
        node = root_.load(std::memory_order_acquire);
        continue;
      }
      if (!node->IsInner()) return node;  // reach a valid leaf node

      // go down to the next level
      node = node->SearchChild(key);
    }
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
    node->LockS();

    return node;
  }

  /**
   * @brief Search a leaf node that may have a target key.
   *
   * This function uses shared locks while traversing a tree, and a returned leaf node
   * is locked by a shared with intent-exclusive lock.
   *
   * @param key a search key.
   * @return a stack of nodes that may have a target key.
   */
  [[nodiscard]] auto
  SearchLeafNodeForWrite(const Key &key)  //
      -> std::vector<Node_t *>
  {
    std::vector<Node_t *> stack{};
    stack.reserve(kExpectedTreeHeight);

    auto *node = root_.load(std::memory_order_acquire);
    while (true) {
      if (!node->IsInner()) {
        // require an exclusive lock for write operations
        node = Node_t::CheckKeyRangeAndLockForWrite(node, key);
        stack.emplace_back(node);
        return stack;
      }

      // traverse side links if needed
      node = Node_t::CheckKeyRangeAndLockForRead(node, key);
      if (node == nullptr) {
        // a root node is removed
        stack.clear();
        node = root_.load(std::memory_order_acquire);
        continue;
      }
      stack.emplace_back(node);

      // go down to the next level
      node = node->SearchChild(key);
    }
  }

  /**
   * @brief Search a parent node of a given one to construct a valid node stack.
   *
   * @param stack the instance of a stack to be reused.
   * @param target_node a child node to be searched.
   */
  void
  SearchParentNode(  //
      std::vector<Node_t *> &stack,
      Node_t *target_node)
  {
    auto *node = root_.load(std::memory_order_acquire);
    const auto key = target_node->GetLowKey();
    // search a target node with its lowest key
    while (true) {
      node = Node_t::CheckKeyRangeAndLockForRead(node, *key);
      if (node == nullptr) {
        // a root node is removed
        stack.clear();
        node = root_.load(std::memory_order_acquire);
        continue;
      }

      stack.emplace_back(node);
      if (node == target_node) {
        // reach a target node
        node->UnlockS();
        stack.pop_back();
        return;
      }

      // go down to the next level
      node = node->SearchChild(*key);
    }
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

    DeleteAlignedPtr(node);
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
   * @param l_node a node to be split.
   * @return a split-right node.
   */
  [[nodiscard]] auto
  HalfSplit(Node_t *l_node)  //
      -> Node_t *
  {
    auto *r_node = new (GetNodePage()) Node_t{l_node->IsInner()};
    l_node->Split(r_node);

    return r_node;
  }

  /**
   * @brief Complete a split operation by inserting a new entry.
   *
   * @param stack nodes in the searched path.
   * @param l_child a left child node.
   * @param r_child a right child node.
   * @param l_key a separator key.
   * @param l_key_len the length of the separator key.
   */
  void
  CompleteSplit(  //
      std::vector<Node_t *> &stack,
      Node_t *l_child,
      Node_t *r_child,
      const Key &l_key,
      const size_t l_key_len)
  {
    stack.pop_back();  // remove a child node from a stack
    Node_t *node = nullptr;
    while (true) {
      if (stack.empty()) {
        // splitting a root node if possible
        if (TryRootSplit(l_child, r_child, l_key, l_key_len)) return;

        // other threads have modified this tree concurrently, so retry
        SearchParentNode(stack, r_child);
        continue;
      }

      if (node == nullptr) {
        // set a parent node if needed
        node = stack.back();
      }

      // traverse horizontally to reach a valid parent node
      node = Node_t::CheckKeyRangeAndLockForWrite(node, l_key);
      if (node == nullptr) {
        // a root node is removed
        if (TryRootSplit(l_child, r_child, l_key, l_key_len)) return;
        SearchParentNode(stack, r_child);
        continue;
      }

      // insert a new entry into a tree
      const auto rc = node->InsertChild(r_child, l_key, l_key_len);
      if (rc == NodeRC::kCompleted) return;
      if (rc == NodeRC::kNeedSplit) {
        // since a parent node is full, perform splitting
        auto *r_node = HalfSplit(node);
        auto &&[target_node, sep_key, sep_key_len] = node->GetValidSplitNode(l_key);
        target_node->InsertChild(r_child, l_key, l_key_len);

        // complete splitting by inserting a new entry
        CompleteSplit(stack, node, r_node, sep_key, sep_key_len);
        if constexpr (IsVarLenData<Key>()) {
          ::operator delete(sep_key);
        }
        return;
      }

      // previous merging has not finished, so wait and retry
      std::this_thread::sleep_for(kRetryWait);
    }
  }

  /**
   * @brief Try to split a root node.
   *
   * @param l_child a split-left (i.e., an old root) node.
   * @param r_child a split-right node.
   * @param l_key a highest key of `l_child`.
   * @param l_key_len the length of a highest key.
   * @retval true if a root node is split.
   * @retval false otherwise.
   */
  [[nodiscard]] auto
  TryRootSplit(  //
      const Node_t *l_child,
      const Node_t *r_child,
      const Key &l_key,
      const size_t l_key_len)  //
      -> bool
  {
    auto *cur_root = root_.load(std::memory_order_relaxed);
    if (cur_root != l_child) return false;

    // install a new root node
    auto *new_root = new (GetNodePage()) Node_t{l_key, l_key_len, l_child, r_child};
    root_.store(new_root, std::memory_order_release);
    return true;
  }

  void
  Merge(  //
      std::vector<Node_t *> &stack,
      Node_t *l_child)
  {
    stack.pop_back();  // remove a child node from a stack
    Node_t *node = nullptr;
    while (true) {
      // check there is a mergeable sibling node
      auto *r_child = l_child->GetMergeableSiblingNode();
      if (r_child == nullptr) return;

      if (stack.empty()) {
        // other threads have modified this tree concurrently, so retry
        SearchParentNode(stack, r_child);
        continue;
      }

      if (node == nullptr) {
        // set a parent node if needed
        node = stack.back();
        stack.pop_back();
      }

      const auto &del_key = r_child->GetLowKey();

      // traverse horizontally to reach a valid parent node
      node = Node_t::CheckKeyRangeAndLockForWrite(node, *del_key);
      if (node == nullptr) {
        // a root node is removed
        SearchParentNode(stack, r_child);
        continue;
      }

      // delete an entry from a tree
      switch (node->DeleteChild(*del_key)) {
        case NodeRC::kCompleted:
          l_child->Merge(r_child);
          gc_.AddGarbage<PageTarget>(r_child);
          return;

        case NodeRC::kAbortMerge:
          l_child->UnlockSIX();
          r_child->UnlockSIX();
          return;

        case NodeRC::kNeedRetry:
          // previous splitting has not finished, so wait and retry
          std::this_thread::sleep_for(kRetryWait);
          continue;

        case NodeRC::kNeedMerge:
        default:
          l_child->Merge(r_child);
          gc_.AddGarbage<PageTarget>(r_child);

          if (stack.empty()) {
            TryShrinkTree(node);
            return;
          }

          l_child = node;
          node = nullptr;
      }
    }
  }

  /**
   * @brief Try to remove a root node.
   *
   * @param node an expected old root node.
   */
  void
  TryShrinkTree(Node_t *node)
  {
    if (node == root_.load(std::memory_order_relaxed) && node->GetRecordCount() == 1) {
      do {
        gc_.AddGarbage<PageTarget>(node);
        node = node->RemoveRoot();
      } while (node->GetRecordCount() == 1 && node->IsInner());
      root_.store(node, std::memory_order_relaxed);
    }

    node->UnlockSIX();
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
   * @param child_nodes child nodes in a lower layer.
   * @retval true if constructed nodes cannot be contained in a single node.
   * @retval false otherwise.
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
}  // namespace dbgroup::index::b_tree::component::psl

#endif  // B_TREE_COMPONENT_PSL_B_TREE_HPP
