/*
 * Copyright 2024 Database Group, Nagoya University
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

#ifndef B_TREE_DBGROUP_B_TREE_COMPONENT_B_TREE_SL_HPP_
#define B_TREE_DBGROUP_B_TREE_COMPONENT_B_TREE_SL_HPP_

// C++ standard libraries
#include <cstddef>
#include <cstdint>
#include <functional>
#include <future>
#include <memory>
#include <optional>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

// external sources
#include "dbgroup/constants.hpp"
#include "dbgroup/index/utility.hpp"
#include "dbgroup/lock/utility.hpp"
#include "dbgroup/memory/utility.hpp"

// local sources
#include "dbgroup/b_tree/component/common.hpp"
#include "dbgroup/b_tree/component/node.hpp"
#include "dbgroup/b_tree/record_iterator.hpp"
#include "dbgroup/b_tree/utility.hpp"

// macros for simplifying OCC-related constexpr.  // NOLINTBEGIN

/// @brief A macro for verifying a procedure reads stable data.
/// @note This macro calls `continue` when version verification failed.
#define DBGROUP_B_TREE_VERIFY_LEAF_VER(node, guard)         \
  if constexpr (kUseOptCCForLeaf) {                         \
    if (!(guard).VerifyVersion(kInsDelMask)) {              \
      tls_stack.emplace_back(std::bit_cast<INode *>(node)); \
      continue;                                             \
    }                                                       \
  }

/// @brief A macro for acquiring an exclusive lock.
/// @note This macro calls `continue` when version verification failed.
#define DBGROUP_B_TREE_TRY_LOCK_LEAF_X(node, guard, x_guard) \
  if constexpr (kUseOptCCForLeaf) {                          \
    (x_guard) = (guard).TryLockX(kInsDelMask);               \
    if (!(x_guard)) {                                        \
      tls_stack.emplace_back(std::bit_cast<INode *>(node));  \
      continue;                                              \
    }                                                        \
  } else {                                                   \
    (x_guard) = (guard).UpgradeToX();                        \
  }

/// @brief A macro for incrementing a version for insert/delete operations.
#define DBGROUP_B_TREE_INCREMENT_LEAF_VER(x_guard) \
  if constexpr (kUseOptCCForLeaf) {                \
    VerIncrement<kInsDelMask>(x_guard);            \
  }

// NOLINTEND

namespace dbgroup::index::b_tree::component
{
/**
 * @brief A class for representing B+trees based on single-level SMOs.
 *
 * @tparam Key A class of stored keys.
 * @tparam Payload A class of stored payloads.
 * @tparam Comp A comparator class for keys.
 * @tparam InnerLock A class for locking inner nodes.
 * @tparam LeafLock A class for locking leaf nodes.
 * @tparam kUseOptCCForLeaf A flag for using optimistic CC in leaf nodes.
 */
template <class Key,
          class Payload,
          class Comp,
          class GC,
          ::dbgroup::lock::OptimisticallyLockable InnerLock,
          ::dbgroup::lock::Lockable LeafLock,
          bool kUseOptCCForLeaf>
class BTreeSL
{
  /*##########################################################################*
   * Type aliases
   *##########################################################################*/

  using INode = Node<Key, Comp, InnerLock, kUseOptimisticCC>;
  using IOptGuard = typename INode::ReadGuard;
  using ISIXGuard = typename INode::SIXGuard;
  using IXGuard = typename INode::XGuard;

  using LNode = Node<Key, Comp, LeafLock, kUseOptCCForLeaf>;
  using LReadGuard = typename LNode::ReadGuard;
  using LCheckGuard = typename LNode::CheckGuard;
  using LOptGuard = typename LNode::OptGuard;
  using LSGuard = typename LNode::SGuard;
  using LXGuard = typename LNode::XGuard;

  using ScanKey = std::optional<std::tuple<const Key &, size_t, bool>>;
  using SIterator = RecordIterator<BTreeSL, LSGuard>;
  using OptIterator = RecordIterator<BTreeSL, LOptGuard>;
  friend SIterator;    // call sibling scan from iterators
  friend OptIterator;  // call sibling scan from iterators

  template <class Entry>
  using BulkIter = typename std::vector<Entry>::const_iterator;
  using NodeEntry = std::tuple<Key, INode *, size_t>;
  using BulkResult = std::pair<size_t, std::vector<NodeEntry>>;
  using BulkPromise = std::promise<BulkResult>;
  using BulkFuture = std::future<BulkResult>;

 public:
  /*##########################################################################*
   * Public types
   *##########################################################################*/

  using Key_t = Key;
  using Payload_t = Payload;
  using Node_t = LNode;

  /*##########################################################################*
   * Public constructors and assignment operators
   *##########################################################################*/

  /**
   * @brief Construct a new tree.
   *
   * @param gc_interval_micro GC internal [us] (default: 10ms).
   * @param gc_thread_num the number of GC threads (default: 1).
   */
  explicit BTreeSL(  //
      const size_t gc_interval_micro = kDefaultGCTime,
      const size_t gc_thread_num = kDefaultGCThreadNum)
      : gc_{std::make_shared<GC>(gc_interval_micro, gc_thread_num)}
  {
  }

  /**
   * @brief Construct a new tree.
   *
   * @param gc A shared garbage collector for reusing internal pages.
   */
  explicit BTreeSL(  //
      const std::shared_ptr<GC> &gc)
      : gc_{gc}
  {
  }

  BTreeSL(const BTreeSL &) = delete;
  BTreeSL(BTreeSL &&) = delete;

  auto operator=(const BTreeSL &) -> BTreeSL & = delete;
  auto operator=(BTreeSL &&) -> BTreeSL & = delete;

  /*##########################################################################*
   * Public destructors
   *##########################################################################*/

  /**
   * @brief Destroy the object.
   *
   */
  ~BTreeSL()
  {
    std::vector<std::pair<INode *, size_t>> stack{};
    stack.reserve(kInitialHeight);
    stack.emplace_back(root_.load(kAcquire), 0);
    while (!stack.empty()) {
      auto &[node, pos] = stack.back();
      if (node->GetLevel() > 0 && pos < node->GetRecNum()) {
        INode *child;
        node->CopyPayloadTo(pos++, child);
        stack.emplace_back(child, 0);
        continue;
      }
      ::dbgroup::memory::Release<Page>(node);
      stack.pop_back();
    }
  }

  /*##########################################################################*
   * Public read APIs
   *##########################################################################*/

  /**
   * @brief Read the payload corresponding to a given key if it exists.
   *
   * @param key A search key.
   * @param key_len The length of a target key.
   * @retval The corresponding payload of a given key if exist.
   * @retval std::nullopt otherwise.
   */
  auto
  Read(  //
      const Key &key,
      [[maybe_unused]] const size_t key_len)  //
      -> std::optional<Payload>
  {
    [[maybe_unused]] const auto &gc_guard = gc_->CreateEpochGuard();

    tls_stack.clear();
    while (true) {
      auto &&[node, guard] = SearchNode<LReadGuard>(tls_stack, key);
      const auto [found, deleted, pos] = node->CheckUniqueness(key);
      if constexpr (!kUseOptCCForLeaf) {
        if (!found || deleted) return std::nullopt;
        node->CopyPayloadTo(pos, tls_pay);
        return tls_pay;
      } else {
        if (found && !deleted) {
          node->CopyPayloadTo(pos, tls_pay);
          if (guard.VerifyVersion(kNoMask, kMaxRetry)) return tls_pay;
        } else if (guard.VerifyVersion(kNoMask, kMaxRetry)) {
          return std::nullopt;
        }
        tls_stack.emplace_back(node);
      }
    }
  }

  /**
   * @brief Perform a range scan with given keys.
   *
   * @tparam kUseOCCIfPossible A flag for using OCC when verifying records.
   * @param begin_key A pair of a begin key and its openness (true=closed).
   * @param end_key A pair of an end key and its openness (true=closed).
   * @return An iterator for accessing scanned records.
   */
  template <bool kUseOCCIfPossible = true>
  auto
  Scan(  //
      const ScanKey &begin_key = std::nullopt,
      const ScanKey &end_key = std::nullopt)
  {
    using Guard = std::conditional_t<kUseOCCIfPossible, LReadGuard, LSGuard>;

    auto &&gc_guard = gc_->CreateEpochGuard();
    LNode *node;
    Guard guard;
    size_t begin_pos;

    tls_stack.clear();
    if constexpr (kUseOptCCForLeaf) {
      thread_local Page tls_page{};  // retain the copy of a target node
      if (begin_key) {
        const auto &[key, _, closed] = *begin_key;
        while (true) {
          std::tie(node, guard) = SearchNode<Guard>(tls_stack, key);
          std::memcpy(static_cast<void *>(&tls_page), node, kPageSize);
          if (guard.VerifyVersion(kNoMask, kMaxRetry)) break;
          tls_stack.emplace_back(node);
        }
        node = std::bit_cast<LNode *>(&tls_page);
        const auto [found, deleted, pos] = node->CheckUniqueness(key);
        begin_pos = pos + static_cast<size_t>(!found || deleted || !closed);
      } else {
        std::tie(node, guard) = SearchLeftmostLeaf(tls_stack);
        do {
          std::memcpy(static_cast<void *>(&tls_page), node, kPageSize);
        } while (!guard.VerifyVersion(kNoMask, kMaxRetry));
        node = std::bit_cast<LNode *>(&tls_page);
        begin_pos = 0;
      }
    } else {
      if (begin_key) {
        const auto &[key, _, closed] = *begin_key;
        std::tie(node, guard) = SearchNode<Guard>(tls_stack, key);
        const auto [found, deleted, pos] = node->CheckUniqueness(key);
        begin_pos = pos + static_cast<size_t>(!found || deleted || !closed);
      } else {
        std::tie(node, guard) = SearchLeftmostLeaf(tls_stack);
        begin_pos = 0;
      }
    }

    return RecordIterator{this, node, std::move(guard), begin_pos, end_key, std::move(gc_guard)};
  }

  /*##########################################################################*
   * Public write APIs
   *##########################################################################*/

  /**
   * @brief Write (i.e., put) a given key/payload pair.
   *
   * @param key A target key.
   * @param payload A target payload.
   * @param key_len The length of a target key.
   * @note This function is the alias of `Upsert`.
   */
  void
  Write(  //
      const Key &key,
      const Payload &payload,
      const size_t key_len = sizeof(Key))
  {
    Upsert(key, payload, key_len);
  }

  /**
   * @brief Upsert (update or insert) a given key/payload pair.
   *
   * @param key A target key.
   * @param payload A target payload.
   * @param key_len The length of a target key.
   * @retval The previous payload if exist (i.e., the function acts as update).
   * @retval std::nullopt otherwise (i.e., the function acts as insert).
   */
  auto
  Upsert(  //
      const Key &key,
      const Payload &payload,
      const size_t key_len = sizeof(Key))  //
      -> std::optional<Payload>
  {
    [[maybe_unused]] const auto &gc_guard = gc_->CreateEpochGuard();

    tls_stack.clear();
    while (true) {
      auto &&[node, guard] = SearchNode<LCheckGuard>(tls_stack, key);
      const auto [found, deleted, pos] = node->CheckUniqueness(key);

      LXGuard x_guard;
      DBGROUP_B_TREE_TRY_LOCK_LEAF_X(node, guard, x_guard);
      if (found && !deleted) {
        node->Update(pos, payload, tls_pay, merger_);
        return tls_pay;
      }
      const auto rc = node->Insert(pos, found, key, key_len, &payload, kPayLen);
      if (rc == kCompleted) {
        DBGROUP_B_TREE_INCREMENT_LEAF_VER(x_guard);
        return std::nullopt;
      }

      Split(tls_stack, node, std::move(x_guard));
      tls_stack.emplace_back(std::bit_cast<INode *>(node));
    }
  }

  /**
   * @brief Insert a given key/payload pair.
   *
   * @param key A target key.
   * @param payload A target payload.
   * @param key_len The length of a target key.
   * @retval The current payload if exist (i.e., failed).
   * @retval std::nullopt otherwise (i.e., succeeded).
   */
  auto
  Insert(  //
      const Key &key,
      const Payload &payload,
      const size_t key_len = sizeof(Key))  //
      -> std::optional<Payload>
  {
    [[maybe_unused]] const auto &gc_guard = gc_->CreateEpochGuard();

    tls_stack.clear();
    while (true) {
      auto &&[node, guard] = SearchNode<LCheckGuard>(tls_stack, key);
      const auto [found, deleted, pos] = node->CheckUniqueness(key);
      if (found && !deleted) {
        node->CopyPayloadTo(pos, tls_pay);
        DBGROUP_B_TREE_VERIFY_LEAF_VER(node, guard);
        return tls_pay;
      }

      LXGuard x_guard;
      DBGROUP_B_TREE_TRY_LOCK_LEAF_X(node, guard, x_guard);
      const auto rc = node->Insert(pos, found, key, key_len, &payload, kPayLen);
      if (rc == kCompleted) {
        DBGROUP_B_TREE_INCREMENT_LEAF_VER(x_guard);
        return std::nullopt;
      }

      Split(tls_stack, node, std::move(x_guard));
      tls_stack.emplace_back(std::bit_cast<INode *>(node));
    }
  }

  /**
   * @brief Update a record using a given kay/payload pair.
   *
   * @param key A target key.
   * @param payload A target payload.
   * @param key_len The length of a target key.
   * @retval The previous payload if exist (i.e., succeeded).
   * @retval std::nullopt otherwise (i.e., failed).
   */
  auto
  Update(  //
      const Key &key,
      const Payload &payload,
      [[maybe_unused]] const size_t key_len = sizeof(Key))  //
      -> std::optional<Payload>
  {
    [[maybe_unused]] const auto &gc_guard = gc_->CreateEpochGuard();

    tls_stack.clear();
    while (true) {
      auto &&[node, guard] = SearchNode<LCheckGuard>(tls_stack, key);
      const auto [found, deleted, pos] = node->CheckUniqueness(key);
      if (!found || deleted) {
        DBGROUP_B_TREE_VERIFY_LEAF_VER(node, guard);
        return std::nullopt;
      }

      LXGuard x_guard;
      DBGROUP_B_TREE_TRY_LOCK_LEAF_X(node, guard, x_guard);
      node->Update(pos, payload, tls_pay, merger_);
      return tls_pay;
    }
  }

  /**
   * @brief Delete a record using a given kay.
   *
   * @param key A target key.
   * @param key_len The length of a target key.
   * @retval The deleted payload if exist (i.e., succeeded).
   * @retval std::nullopt otherwise (i.e., failed).
   */
  auto
  Delete(  //
      const Key &key,
      [[maybe_unused]] const size_t key_len = sizeof(Key))  //
      -> std::optional<Payload>
  {
    [[maybe_unused]] const auto &gc_guard = gc_->CreateEpochGuard();

    tls_stack.clear();
    while (true) {
      auto &&[node, guard] = SearchNode<LCheckGuard>(tls_stack, key);
      const auto [found, deleted, pos] = node->CheckUniqueness(key);
      if (!found || deleted) {
        DBGROUP_B_TREE_VERIFY_LEAF_VER(node, guard);
        return std::nullopt;
      }

      LXGuard x_guard;
      DBGROUP_B_TREE_TRY_LOCK_LEAF_X(node, guard, x_guard);
      const auto rc = node->Delete(pos, tls_pay);
      DBGROUP_B_TREE_INCREMENT_LEAF_VER(x_guard);
      if (rc == kNeedMerge) {
        TryMerge(tls_stack, node, x_guard.DowngradeToSIX());
      }
      return tls_pay;
    }
  }

  /*##########################################################################*
   * Public bulkload API
   *##########################################################################*/

  /**
   * @brief Bulkload specified kay/payload pairs.
   *
   * @tparam Entry A container for a key/payload pair.
   * @param entries All entries for bulkload.
   * @param thread_num The number of threads used for bulkload.
   */
  template <class Entry>
  void
  Bulkload(  //
      const std::vector<Entry> &entries,
      const size_t thread_num = 1)
  {
    if (entries.empty()) return;

    std::vector<NodeEntry> nodes{};
    auto &&iter = entries.cbegin();
    const auto rec_num = entries.size();
    size_t level = 1;
    if (thread_num <= 1 || rec_num < thread_num) {
      std::tie(level, nodes) = BulkloadWithSingleThread<Entry>(iter, rec_num);
    } else {
      // construct partial trees
      std::vector<BulkFuture> futures{};
      futures.reserve(thread_num);
      for (size_t i = 0; i < thread_num; ++i) {
        BulkPromise p{};
        futures.emplace_back(p.get_future());
        const size_t n = (rec_num + i) / thread_num;
        std::thread t{[this](BulkPromise p, BulkIter<Entry> iter, size_t n) {
                        p.set_value(BulkloadWithSingleThread<Entry>(iter, n));
                      },
                      std::move(p), iter, n};
        t.detach();
        iter += n;
      }

      // wait for the worker threads to construct partial trees
      std::vector<BulkResult> trees{};
      trees.reserve(thread_num);
      for (auto &&future : futures) {
        trees.emplace_back(future.get());
        const auto partial_height = trees.back().first;
        level = (partial_height > level) ? partial_height : level;
      }

      // align the height of partial trees and link their rightmost/leftmost nodes
      nodes.reserve(2 * kNodeCapacity * thread_num);
      INode *prev_node = nullptr;
      for (auto &&[p_level, p_nodes] : trees) {
        while (p_level < level) {
          p_nodes = ConstructSingleLevel<NodeEntry>(p_nodes.cbegin(), p_nodes.size(), p_level++);
        }
        nodes.insert(nodes.end(), p_nodes.begin(), p_nodes.end());
        INode::LinkVerticalBorderNodes(prev_node, std::get<1>(p_nodes.front()));
        prev_node = std::get<1>(p_nodes.back());
      }
    }

    // create upper layers until a root node is created
    while (nodes.size() > 1) {
      nodes = ConstructSingleLevel<NodeEntry>(nodes.cbegin(), nodes.size(), level++);
    }
    auto *new_root = std::get<1>(nodes.front());
    INode::RemoveLeftmostKeys(new_root);

    // set a new root
    auto *old_root = root_.exchange(new_root, kRelease);
    gc_->template AddGarbage<Page>(old_root);
  }

  /*##########################################################################*
   * Public APIs
   *##########################################################################*/

  /**
   * @brief Set a merger function for performing read-modify-write operations.
   *
   * @param merger A function for merging payloads.
   */
  constexpr void SetRecordMerger(  //
      Payload (*merger)(const Payload &, const Payload &))
  {
    merger_ = merger;
  }

 private:
  /*##########################################################################*
   * Internal types
   *##########################################################################*/

  /// @brief A type for allocating variable-length data.
  using KeyWOPtr = std::remove_pointer_t<Key>;

  /*##########################################################################*
   * Internal constants
   *##########################################################################*/

  /// @brief The initial capacity of a node stack.
  static constexpr size_t kInitialHeight = 8;

  /// @brief A flag for indicating inner nodes.
  static constexpr bool kInnerFlag = true;

  /// @brief The size of payloads.
  static constexpr size_t kPayLen = sizeof(Payload);

  /// @brief The expected capacity of a node when performing bulkload.
  static constexpr size_t kNodeCapacity = kPageSize / (sizeof(Key) + kPayLen + kMetaSize);

  /*##########################################################################*
   * Internal utility functions
   *##########################################################################*/

  /**
   * @return A page for nodes.
   */
  [[nodiscard]] auto
  GetNodePage()  //
      -> void *
  {
    auto *page = gc_->template GetPageIfPossible<Page>();
    return page ? page : ::dbgroup::memory::Allocate<Page>();
  }

  /**
   * @brief Search a target node horizontally.
   *
   * @tparam Guard A guard class for reading a target node.
   * @tparam NodeT A target node class.
   * @param[in,out] node A target node.
   * @param[in,out] guard The guard instance of a target node.
   * @param[in] key A search key.
   * @retval true if a container node has been found.
   * @retval false otherwise.
   */
  template <class Guard, class NodeT>
  [[nodiscard]] static auto
  SearchHorizontally(  //
      NodeT *&node,
      Guard &guard,
      const Key &key)  //
      -> bool
  {
    size_t i = 0;
    while (true) {
      auto *sib_node = node->GetSibNode();
      const auto removed = node->Removed();
      const auto included = node->Include(key);
      if constexpr (std::is_same_v<Guard, typename NodeT::OptGuard>) {
        if (!guard.VerifyVersion(kSMOMask)) continue;
      }

      if (!removed && included) return true;
      if (i++ > 0) return false;  // the traversal path may be too old

      node = sib_node;
      guard = node->template GetGuard<Guard>();
    }
  }

  /**
   * @brief Search a node that may have a target key.
   *
   * @tparam Guard A desired guard class.
   * @tparam kIsInner A flag for indicating whether a target is inner nodes.
   * @param[in,out] stack A stack of traversed nodes.
   * @param[in] key A search key.
   * @param[in] level The bottom level of search space.
   * @retval 1st: A found node.
   * @retval 2nd: The lock guard of a found node.
   */
  template <class Guard, class NodeT = LNode>
  [[nodiscard]] auto
  SearchNode(  //
      std::vector<INode *> &stack,
      const Key &key,
      const size_t level = 0) const  //
      -> std::pair<NodeT *, Guard>
  {
    INode *node;
    while (true) {
      if (stack.empty()) {
        node = root_.load(kAcquire);
      } else {
        node = stack.back();
        stack.pop_back();
      }
      if (node->GetLevel() < level) {
        stack.clear();  // the stack has expired nodes
        continue;
      }

      while (true) {
        if (node->GetLevel() == level) {
          auto *out_node = std::bit_cast<NodeT *>(node);
          auto &&guard = out_node->template GetGuard<Guard>();
          if (!SearchHorizontally(out_node, guard, key)) break;
          return {out_node, std::move(guard)};
        }

        INode *child{};
        auto &&guard = node->template GetGuard<IOptGuard>();
        while (true) {
          if (!SearchHorizontally(node, guard, key)) goto out;
          child = node->SearchChild(key);
          if (guard.VerifyVersion(kInsDelMask)) break;
        }
        stack.emplace_back(std::exchange(node, child));
      }
    out:;  // the current node has been removed or too old, so retry
    }
  }

  /**
   * @brief Search the leftmost leaf node.
   *
   * @param[in,out] stack A stack of traversed nodes.
   * @retval 1st: The leftmost node.
   * @retval 2nd: The lock guard of the leftmost node.
   */
  [[nodiscard]] auto
  SearchLeftmostLeaf(                     //
      std::vector<INode *> &stack) const  //
      -> std::pair<LNode *, LReadGuard>
  {
    auto *node = root_.load(kAcquire);
    while (true) {
      if (node->GetLevel() == 0) {
        auto *leaf = std::bit_cast<LNode *>(node);
        return {leaf, leaf->template GetGuard<LReadGuard>()};
      }

      INode *child;
      auto &&guard = node->template GetGuard<IOptGuard>();
      do {
        node->CopyPayloadTo(0, child);
      } while (!guard.VerifyVersion(kInsDelMask));
      stack.emplace_back(std::exchange(node, child));
    }
  }

  /**
   * @brief Go to the next node for scanning.
   *
   * @tparam Guard A class for representing lock guards.
   * @param[in,out] node A target node.
   * @param[in,out] guard The guard instance of a target node.
   * @return The begin position for scanning.
   */
  template <class Guard>
  [[nodiscard]] auto
  SiblingScan(  //
      LNode *&node,
      Guard &guard) const  //
      -> size_t
  {
    size_t begin_pos{};
    if constexpr (std::is_same_v<Guard, LOptGuard>) {
      auto *sib_node = node->GetSibNode();
      const auto &key = node->GetSeparatorKey().first;

      tls_stack.clear();
      do {
        tls_stack.emplace_back(sib_node);
        std::tie(sib_node, guard) = SearchNode<Guard>(tls_stack, key);
        std::memcpy(static_cast<void *>(node), sib_node, kPageSize);
      } while (!guard.VerifyVersion(kNoMask, kMaxRetry));

      auto [found, deleted, pos] = node->CheckUniqueness(key);
      begin_pos = pos + static_cast<size_t>(!found || deleted);
      if constexpr (IsVarLenData<Key>()) {
        ::dbgroup::memory::Release<KeyWOPtr>(key);
      }
    } else {
      node = node->GetSibNode();
      guard = node->template GetGuard<Guard>();
    }
    return begin_pos;
  }

  /*##########################################################################*
   * Internal structure modification operations
   *##########################################################################*/

  /**
   * @brief Split a given node and insert a new down link.
   *
   * @tparam NodeT A target node class.
   * @param stack A stack of traversed nodes.
   * @param l_child A split-left node.
   * @param l_guard The guard instance of a split-left node.
   * @param level A level where a parent node is.
   */
  template <class NodeT>
  void
  Split(  //
      std::vector<INode *> stack,
      NodeT *l_child,
      typename NodeT::XGuard l_guard,
      const size_t level = 1)
  {
    auto *r_child = new (GetNodePage()) NodeT{l_child};
    if constexpr (NodeT::kUseOptCC) {
      VerIncrement<kSMOMask>(l_guard);
    }
    l_guard = typename NodeT::XGuard{};

    const auto &[key, key_len] = l_child->GetSeparatorKey();
    while (true) {
      if (stack.empty() && TryAddRoot(level, key, key_len, l_child, r_child)) break;

      // insert a new down link
      auto &&[node, guard] = SearchNode<IOptGuard, INode>(stack, key, level);
      const auto [found, _, pos] = node->CheckUniqueness(key);  // `found` must be false
      auto &&x_guard = guard.TryLockX(kInsDelMask);
      if (!x_guard) continue;  // another thread may insert the key
      const auto rc = node->Insert(pos, found, key, key_len, &r_child, kPtrSize);
      if (rc == kCompleted) {
        VerIncrement<kInsDelMask>(x_guard);
        break;
      }

      Split(stack, node, std::move(x_guard), level + 1);
      stack.emplace_back(node);
    }
    if constexpr (IsVarLenData<Key>()) {
      ::dbgroup::memory::Release<KeyWOPtr>(key);
    }
  }

  /**
   * @brief Try to install a new root node.
   *
   * @param level The level of a new root node.
   * @param key A separator key.
   * @param key_len The length of a separator key.
   * @param l_child A left child.
   * @param r_child A right child.
   * @retval true if a new root node has been created.
   * @retval false otherwise.
   */
  [[nodiscard]] auto
  TryAddRoot(  //
      const size_t level,
      const Key &key,
      const size_t key_len,
      const void *l_child,
      const void *r_child)  //
      -> bool
  {
    auto *old_root = root_.load(kRelaxed);
    if (old_root == l_child) {
      auto *root = new (GetNodePage()) INode{level, key, key_len, l_child, r_child};
      if (root_.compare_exchange_strong(old_root, root, kRelease, kRelaxed)) return true;
      ::dbgroup::memory::Release<Page>(root);
    }
    return false;
  }

  /**
   * @brief Try to remove a down link and merge a given node.
   *
   * @tparam NodeT A target node class.
   * @tparam SIXGuard A class for representing SIX-lock guards.
   * @param stack A stack of traversed nodes.
   * @param l_child A left child (to be merged) node.
   * @param l_guard The SIX guard instance of a given child node.
   * @param level A level where a target node is.
   */
  template <class NodeT, class SIXGuard>
  void
  TryMerge(  //
      std::vector<INode *> &stack,
      NodeT *l_child,
      SIXGuard l_guard,
      const size_t level = 1)
  {
    auto &&r_guard = l_child->GetMergeableSib();
    if (!r_guard) {
      if constexpr (std::is_same_v<NodeT, INode>) {
        TryRemoveRoot(l_child, std::move(l_guard));
      }
      return;
    }

    const auto &key = l_child->GetSeparatorKey().first;
    while (true) {
      auto &&[node, guard] = SearchNode<IOptGuard, INode>(stack, key, level);
      const auto [found, _, pos] = node->CheckUniqueness(key);
      if (!found) {
        if (!guard.VerifyVersion(kInsDelMask)) continue;
        break;  // a downlink has not been inserted yet or is leftmost in a node
      }
      auto &&x_guard = guard.TryLockX(kInsDelMask);
      if (!x_guard) continue;  // another thread may modify a node

      NodeT *r_child;
      const auto rc = node->Delete(pos, r_child);
      VerIncrement<kInsDelMask>(x_guard);
      if (rc == kCompleted) {
        x_guard = IXGuard{};
        l_child->Merge(std::move(l_guard), r_child, std::move(r_guard));
      } else {
        auto &&six_guard = x_guard.DowngradeToSIX();
        l_child->Merge(std::move(l_guard), r_child, std::move(r_guard));
        TryMerge(stack, node, std::move(six_guard), level);
      }
      gc_->template AddGarbage<Page>(r_child);
      break;
    }
    if constexpr (IsVarLenData<Key>()) {
      ::dbgroup::memory::Release<KeyWOPtr>(key);
    }
  }

  /**
   * @brief Try to remove the top level of this tree.
   *
   * @param node A target inner node.
   * @param guard The SIX guard instance of a given node.
   */
  void
  TryRemoveRoot(  //
      INode *node,
      ISIXGuard guard)
  {
    auto *root = root_.load(kRelaxed);
    if (node == root && !node->GetSibNode() && node->GetRecNum() == 1) {
      node->CopyPayloadTo(0, root);
      if (root_.compare_exchange_strong(node, root, kRelease, kRelaxed)) {
        node->Remove(std::move(guard), root);
        gc_->template AddGarbage<Page>(node);
      }
    }
  }

  /*##########################################################################*
   * Internal utilities for bulkload
   *##########################################################################*/

  /**
   * @brief Bulkload specified kay/payload pairs with a single thread.
   *
   * @tparam Entry A container of a key/payload pair.
   * @param iter The begin position of target records.
   * @param n The number of entries.
   * @retval 1st: The height of a constructed tree.
   * @retval 2nd: Constructed nodes in the top level.
   */
  template <class Entry>
  auto
  BulkloadWithSingleThread(  //
      BulkIter<Entry> iter,
      size_t n)  //
      -> BulkResult
  {
    size_t level = 0;
    auto &&nodes = ConstructSingleLevel<Entry>(iter, n, level++);
    while (true) {
      n = nodes.size();
      if (n < 2 * kNodeCapacity) break;
      nodes = ConstructSingleLevel<NodeEntry>(nodes.cbegin(), n, level++);
    }
    return {level, std::move(nodes)};
  }

  /**
   * @brief Construct nodes based on given entries.
   *
   * @tparam Entry A container of a key/payload pair.
   * @param iter The begin position of target records.
   * @param n The number of entries.
   * @param level A target level to be constructed.
   * @return Constructed nodes.
   */
  template <class Entry>
  auto
  ConstructSingleLevel(  //
      BulkIter<Entry> iter,
      const size_t n,
      const size_t level)  //
      -> std::vector<NodeEntry>
  {
    std::vector<NodeEntry> nodes{};
    nodes.reserve((n / kNodeCapacity) + 1);
    const auto &iter_end = std::next(iter, n);
    for (INode *prev_node = nullptr; iter < iter_end;) {
      const auto &[key, key_len] = ParseKey(*iter);
      auto *node = new (GetNodePage()) INode{level};
      nodes.emplace_back(key, node, key_len);
      if (prev_node) {
        prev_node->LinkSiblingNode(node, key, key_len);
      }

      node->Bulkload(iter, iter_end);
      prev_node = node;
    }
    return nodes;
  }

  /*##########################################################################*
   * Static assertions
   *##########################################################################*/

  static_assert(  //
      IsTriviallyCopyable<Key>(),
      "A key type must be trivially copyable (i.e., copyable with std::memcpy).");

  static_assert(  //
      IsTriviallyCopyable<Payload>(),
      "A payload type must be trivially copyable (i.e., copyable with std::memcpy).");

  /// @brief The expected maximum size after node split.
  static constexpr size_t kMaxSplitSize = kHeaderSize  //
                                          + (kPageSize - kHeaderSize) * 3 / 4
                                          + (MaxSize<Key>() + kPayLen + kMetaSize) / 2
                                          + MaxSize<Key>();

  static_assert(  //
      kMaxSplitSize <= kPageSize,
      "The page size is too small to store given key/payload pairs.");

  static_assert(  //
      GC::template HasTarget<Page>(),
      "Garbage collector must reuse node pages.");

  /*##########################################################################*
   * Internal member variables
   *##########################################################################*/

  /// @brief A garbage collector.
  std::shared_ptr<GC> gc_{};

  /// @brief The root node of this tree.
  std::atomic<INode *> root_{new (GetNodePage()) INode{}};

  /// @brief A function for merging payloads.
  Payload (*merger_)(const Payload &, const Payload &){};

  /*##########################################################################*
   * Thread local class variables
   *##########################################################################*/

  /// @brief A stack for retaining traversed nodes.
  static thread_local inline std::vector<INode *> tls_stack{kInitialHeight};

  /// @brief A temporary payload.
  static thread_local inline Payload tls_pay{};
};

}  // namespace dbgroup::index::b_tree::component

#endif  // B_TREE_DBGROUP_B_TREE_COMPONENT_B_TREE_SL_HPP_
