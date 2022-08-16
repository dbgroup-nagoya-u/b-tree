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

#ifndef B_TREE_COMPONENT_OML_B_TREE_HPP
#define B_TREE_COMPONENT_OML_B_TREE_HPP

#include <future>
#include <optional>
#include <thread>
#include <vector>

// organization libraries
#include "memory/epoch_based_gc.hpp"

// local sources
#include "b_tree/component/record_iterator.hpp"
// #include "node_fixlen.hpp"
#include "node_varlen.hpp"

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
 * @tparam kIsVarLen a flag for indicating variable-length keys.
 */
template <class Key, class Payload, class Comp, bool kIsVarLen>
class BTree
{
 public:
  /*####################################################################################
   * Type aliases
   *##################################################################################*/

  using K = Key;
  using V = Payload;
  using NodeVarLen_t = NodeVarLen<Key, Comp>;
  using NodeFixLen_t = NodeVarLen<Key, Comp>;  // NodeFixLen<Key, Comp>;
  using Node_t = std::conditional_t<kIsVarLen, NodeVarLen_t, NodeFixLen_t>;
  using BTree_t = BTree<Key, Payload, Comp, kIsVarLen>;
  using RecordIterator_t = RecordIterator<BTree_t>;
  using GC_t = ::dbgroup::memory::EpochBasedGC<Node_t>;

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
   * @brief Construct a new BTree object.
   *
   */
  BTree(  //
      const size_t gc_interval_micro,
      const size_t gc_thread_num)
      : gc_{gc_interval_micro, gc_thread_num, true}
  {
    auto *root = new Node_t{kLeafFlag};
    if constexpr (!kIsVarLen) {
      root->SetPayloadLength(kPayLen);
    }
    root_.store(root, std::memory_order_release);
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
  Read(const Key &key)  //
      -> std::optional<Payload>
  {
    [[maybe_unused]] const auto &guard = gc_.CreateEpochGuard();

    auto &&[node, ver] = SearchLeafNodeForRead(key, kClosed);

    Payload payload{};
    while (true) {
      const auto rc = node->Read(key, payload);
      const auto retry_rc = node->NeedRetryLeaf(key, ver);

      if (retry_rc == kNeedNextRetry) {
        const auto [next_node, new_ver] = node->GetNextNode();
        if (new_ver == ver) {
          node = next_node;
        }
      } else if (retry_rc == kCompleted) {
        if (rc == kSuccess) return payload;
        return std::nullopt;
      } else {
        ver = node->GetVersion();
      }
    }
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
      const std::optional<std::pair<const Key &, bool>> &begin_key = std::nullopt,
      const std::optional<std::pair<const Key &, bool>> &end_key = std::nullopt)  //
      -> RecordIterator_t
  {
    [[maybe_unused]] const auto &guard = gc_.CreateEpochGuard();

    Node_t *node{};
    size_t begin_pos = 0;

    if (begin_key) {
      const auto &[key, is_closed] = begin_key.value();
      node = SearchLeafNodeForRead(key, is_closed).first;
      node->LockS();
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

    auto *node = SearchLeafNodeForWrite(key);
    node->Write(key, key_len, &payload, kPayLen);
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

    auto *node = SearchLeafNodeForWrite(key);
    return node->Insert(key, key_len, &payload, kPayLen);
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

    auto *node = SearchLeafNodeForWrite(key);
    return node->Delete(key);
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
      Node_t *prev_node = nullptr;
      for (auto &&[p_height, p_nodes] : partial_trees) {
        while (p_height < height) {  // NOLINT
          ConstructUpperLayer(p_nodes);
          ++p_height;
        }
        nodes.insert(nodes.end(), p_nodes.begin(), p_nodes.end());

        // link partial trees
        if (prev_node != nullptr) {
          Node_t::LinkVerticalBorderNodes(prev_node, p_nodes.front());
        }
        prev_node = p_nodes.back();
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

  /// a flag for indicating leaf nodes.
  static constexpr uint32_t kLeafFlag = 1;

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

  /// the maximum length of keys.
  static constexpr size_t kMaxKeyLen = (kIsVarLen) ? kMaxVarLenDataSize : sizeof(Key);

  /// the maximum length of payloads (including child pointers).
  static constexpr size_t kMaxPayLen = (kPayLen < kPtrLen) ? kPtrLen : kPayLen;

  /// the maximum length of records.
  static constexpr size_t kMaxRecLen = kMaxKeyLen + kMaxPayLen + ((kIsVarLen) ? kMetaLen : 0);

  /// the minimum block size in a certain node.
  static constexpr size_t kMinBlockSize = kPageSize - kHeaderLen - kMaxKeyLen;

  /// the expected length of leaf records.
  static constexpr size_t kExpLeafRecLen = sizeof(Key) + kPayLen + ((kIsVarLen) ? kMetaLen : 0);

  /// the expected length of internal records.
  static constexpr size_t kExpInnerRecLen = sizeof(Key) + kPtrLen + ((kIsVarLen) ? kMetaLen : 0);

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
    auto *page = gc_.template GetPageIfPossible<Node_t>();
    return (page == nullptr) ? (::operator new(kPageSize)) : page;
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
    Node_t *node;
    while (true) {
      node = root_.load(std::memory_order_acquire);

      // check this tree requires SMOs
      auto rc = ShrinkTreeIfNeeded(node);
      if (rc == kNeedRootRetry) {
        continue;
      }

      node = root_.load(std::memory_order_acquire);
      auto [next_node, next_rc] = RootSplit(node, key);
      if (next_rc == kNeedRootRetry) {
        continue;
      }
      node = next_node;
      break;
    }

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
      -> std::pair<Node_t *, uint64_t>
  {
    while (true) {  // for root retry
      auto *node = root_.load(std::memory_order_acquire);
      while (!node->IsLeaf()) {
        const auto [pos, ver] = node->SearchChild(key, is_closed);
        auto *child = node->GetChild(pos);
        const auto rc = node->NeedRetryInner(key, ver);
        if (rc == kNeedRetry) {
          continue;
        } else if (rc == kNeedRootRetry) {
          break;
        } else {
          node = child;
        }
      }
      const auto ver = node->GetVersion();
      return {node, ver};
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
    while (!node->IsLeaf()) {
      node = node->GetChild(0);
    }
    node->LockS();
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
    Node_t *node;
    while (true) {
      node = GetRootForWrite(key);
      while (!node->IsLeaf()) {
        // search a child node
        const auto [pos, ver] = node->SearchChild(key, kClosed);
        auto *child = node->GetChild(pos);
        const auto rc = node->NeedRetryInner(key, ver);
        if (rc == kNeedRetry) {
          continue;
        } else if (rc == kNeedRootRetry) {
          break;
        }

        // perform internal SMOs eagerly
        auto retry_rc = CheckSMO(node, child, ver, true);  // check split
        if (retry_rc == kNeedRetry) {
          continue;
        } else if (retry_rc == kNeedSplit) {
          Split(child, node, pos);
          node = child->GetValidSplitNode(key);
          continue;
        }

        retry_rc = CheckSMO(node, child, ver, false);  // check merge
        if (retry_rc == kNeedRetry) {
          continue;
        } else if (retry_rc == kNeedMerge) {
          Merge(child, node, pos);
        }

        node = child;
      }
      break;
    }

    // search for a valid node in leaf nodes
    while (true) {
      const auto ver = node->GetVersion();
      node->LockX();
      const auto retry_rc = node->NeedRetryLeaf(key, ver);
      if (retry_rc == kCompleted) {
        break;
      } else if (retry_rc == kNeedRetry) {
        node->UnlockX();
        continue;
      } else {
        node = node->GetValidNextNode(key);
      }
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

    ::operator delete(node);
  }

  [[nodiscard]] auto
  CheckSMO(Node_t *parent,
           Node_t *child,
           const uint64_t parent_ver,
           const bool is_split)  //
      -> NodeRC
  {
    while (true) {
      auto [should_smo, child_ver] = (is_split) ? child->NeedSplit(kMaxRecLen) : child->NeedMerge();
      (is_split) ? child->NeedSplit(kMaxRecLen) : child->NeedMerge();
      if (!should_smo) {
        if (!child->HasSameVersion(child_ver)) {
          continue;
        } else {
          return kCompleted;
        }
      } else {
        parent->LockX();
        child->LockX();
        if (!parent->CheckMyVersion(parent_ver) || !child->CheckMyVersion(child_ver)) {
          parent->UnlockX();
          child->UnlockX();
          return kNeedRetry;
        }
        return (is_split) ? kNeedSplit : kNeedMerge;
      }
    }
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
      const size_t pos)  //
  {
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
  [[nodiscard]] auto  //
  RootSplit(Node_t *node,
            const Key &key)  //
      -> std::pair<Node_t *, NodeRC>
  {
    while (true) {
      auto [need_split, ver] = node->NeedSplit(kMaxRecLen);
      if (!need_split) {
        if (!node->HasSameVersion(ver)) {
          continue;
        } else {
          return {node, kCompleted};
        }
      } else {
        node->LockX();
        if (node != root_) {
          node->UnlockX();
          return {nullptr, kNeedRootRetry};
        } else if (!node->CheckMyVersion(ver)) {
          node->UnlockX();
          continue;
        } else {
          auto *l_node = node;
          auto *r_node = new Node_t{l_node->IsLeaf()};
          l_node->Split(r_node);
          auto *new_root = new Node_t{l_node, r_node};
          root_.store(new_root, std::memory_order_release);
          return {l_node->GetValidSplitNode(key), kCompleted};
        }
      }
    }
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
    l_node->Merge(r_node);
    parent->DeleteChild(l_node, l_pos);
    gc_.AddGarbage(r_node);
  }

  /**
   * @brief Remove a root node that has only one child node.
   *
   */
  auto
  ShrinkTreeIfNeeded(Node_t *node)  //
      -> NodeRC
  {
    while (true) {
      auto ver = node->GetVersion();
      if (node->GetRecordCount() == 1 && !node->IsLeaf()) {
        // if a root node has only one child, shrink a tree
        auto *child = node->template GetPayload<Node_t *>(0);

        node->LockX();
        if (node != root_) {  // if the root node has been moved by another thread
          node->UnlockX();
          return kNeedRootRetry;
        } else if (!node->CheckMyVersion(ver)) {  // TODO: use check ver
          node->UnlockX();
          continue;
        } else {
          node->SetRemoved();
          gc_.AddGarbage(node);
          node->UnlockX();
          node = child;
          root_.store(node, std::memory_order_relaxed);
        }
      } else {
        return kCompleted;
      }
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
  std::atomic<Node_t *> root_ = new Node_t{kLeafFlag};
};
}  // namespace dbgroup::index::b_tree::component::oml

#endif  // B_TREE_COMPONENT_OML_B_TREE_HPP
