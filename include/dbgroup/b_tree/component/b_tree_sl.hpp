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
#include "dbgroup/memory/epoch_based_gc.hpp"
#include "dbgroup/memory/utility.hpp"

// local sources
#include "dbgroup/b_tree/component/common.hpp"
#include "dbgroup/b_tree/component/node.hpp"
#include "dbgroup/b_tree/record_iterator.hpp"

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
  using LXGuard = typename LNode::XGuard;

  using ScanKey = std::optional<std::tuple<const Key &, size_t, bool>>;
  using Iterator = RecordIterator<BTreeSL>;
  friend Iterator;  // call sibling scan from iterators

  template <class Entry>
  using BulkIter = typename std::vector<Entry>::const_iterator;
  using NodeEntry = std::tuple<Key, INode *, size_t>;
  using BulkResult = std::pair<size_t, std::vector<NodeEntry>>;
  using BulkPromise = std::promise<BulkResult>;
  using BulkFuture = std::future<BulkResult>;

  using GC = ::dbgroup::memory::EpochBasedGC<Page>;

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
   */
  BTreeSL() = default;

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
        stack.emplace_back(node->GetChild(pos++), 0);
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
   * @note If you store variable-length payloads, thread local storage
   * temporarily retains the returned payload. You need to create a deep copy to
   * refer to the returned payload safely whenever you want.
   */
  auto
  Read(  //
      const Key &key,
      [[maybe_unused]] const size_t key_len)  //
      -> std::optional<Payload>
  {
    [[maybe_unused]] const auto &gc_guard = gc_.CreateEpochGuard();

    std::vector<INode *> stack{};
    stack.reserve(kInitialHeight);
    while (true) {
      auto &&[node, guard] = SearchNode<LReadGuard>(stack, key);
      const auto &payload = node->template Read<Payload>(key);
      if constexpr (kUseOptCCForLeaf) {
        if (!guard.VerifyVersion()) {
          stack.emplace_back(node);
          continue;
        }
      }
      return payload;
    }
  }

  /**
   * @brief Perform a range scan with given keys.
   *
   * @param begin_key A pair of a begin key and its openness (true=closed).
   * @param end_key A pair of an end key and its openness (true=closed).
   * @return An iterator for accessing scanned records.
   */
  auto
  Scan(  //
      const ScanKey &begin_key = std::nullopt,
      const ScanKey &end_key = std::nullopt)  //
      -> Iterator
  {
    auto &&gc_guard = gc_.CreateEpochGuard();

    LNode *node;
    LReadGuard guard;
    size_t begin_pos;
    std::vector<INode *> stack{};
    stack.reserve(kInitialHeight);
    if constexpr (kUseOptCCForLeaf) {
      thread_local Page tls_page{};  // retain the copy of a target node
      if (begin_key) {
        const auto &[key, _, closed] = *begin_key;
        while (true) {
          std::tie(node, guard) = SearchNode<LReadGuard>(stack, key);
          std::memcpy(static_cast<void *>(&tls_page), node, kPageSize);
          if (guard.VerifyVersion()) break;
          stack.emplace_back(node);
        }
        node = std::bit_cast<LNode *>(&tls_page);
        const auto [found, deleted, pos] = node->CheckUniqueness(key);
        begin_pos = pos + static_cast<size_t>(!found || deleted || !closed);
      } else {
        std::tie(node, guard) = SearchLeftmostLeaf(stack);
        do {
          std::memcpy(static_cast<void *>(&tls_page), node, kPageSize);
        } while (!guard.VerifyVersion());
        node = std::bit_cast<LNode *>(&tls_page);
        begin_pos = 0;
      }
    } else {
      if (begin_key) {
        const auto &[key, _, closed] = *begin_key;
        std::tie(node, guard) = SearchNode<LReadGuard>(stack, key);
        const auto [found, deleted, pos] = node->CheckUniqueness(key);
        begin_pos = pos + static_cast<size_t>(!found || deleted || !closed);
      } else {
        std::tie(node, guard) = SearchLeftmostLeaf(stack);
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
   * @param pay_len The length of a target payload.
   * @return kSuccess.
   * @note This function always overwrites a payload and can be optimized for
   * that purpose; the procedure can omit the key uniqueness check.
   */
  auto
  Write(  //
      const Key &key,
      const Payload &payload,
      const size_t key_len = sizeof(Key),
      const size_t pay_len = sizeof(Payload))  //
      -> ReturnCode
  {
    [[maybe_unused]] const auto &gc_guard = gc_.CreateEpochGuard();

    const auto *pay_addr = GetSrcAddr(payload);
    std::vector<INode *> stack{};
    stack.reserve(kInitialHeight);
    while (true) {
      auto &&[node, x_guard] = SearchNode<LXGuard>(stack, key);
      const auto rc = node->Write(x_guard, key, key_len, pay_addr, pay_len);
      if (rc == kCompleted) return kSuccess;

      const auto &[sep_key, sep_key_len, r_node] = Split(node, std::move(x_guard));
      AddDownLink(stack, sep_key, sep_key_len, node, r_node);
      stack.emplace_back(std::bit_cast<INode *>(node));
    }
  }

  /**
   * @brief Insert a given key/payload pair.
   *
   * @param key A target key.
   * @param payload A target payload.
   * @param key_len The length of a target key.
   * @param pay_len The length of a target payload.
   * @retval kSuccess if a given record is inserted.
   * @retval kKeyExist if a target key has been already inserted.
   */
  auto
  Insert(  //
      const Key &key,
      const Payload &payload,
      const size_t key_len = sizeof(Key),
      const size_t pay_len = sizeof(Payload))  //
      -> ReturnCode
  {
    [[maybe_unused]] const auto &gc_guard = gc_.CreateEpochGuard();

    const auto *pay_addr = GetSrcAddr(payload);
    std::vector<INode *> stack{};
    stack.reserve(kInitialHeight);
    while (true) {
      auto &&[node, guard] = SearchNode<LCheckGuard>(stack, key);
      const auto [found, deleted, pos] = node->CheckUniqueness(key);
      if (found && !deleted) {
        if constexpr (kUseOptCCForLeaf) {
          if (!guard.VerifyVersion(kInsDelMask)) continue;
        }
        return kKeyExist;
      }

      LXGuard x_guard;
      if constexpr (kUseOptCCForLeaf) {
        x_guard = guard.TryLockX(kInsDelMask);
        if (!x_guard) continue;  // another thread may insert the key
      } else {
        x_guard = guard.UpgradeToX();
      }
      const auto rc = node->Insert(x_guard, pos, found, key, key_len, pay_addr, pay_len);
      if (rc == kCompleted) return kSuccess;

      const auto &[sep_key, sep_key_len, r_node] = Split(node, std::move(x_guard));
      AddDownLink(stack, sep_key, sep_key_len, node, r_node);
      stack.emplace_back(std::bit_cast<INode *>(node));
    }
  }

  /**
   * @brief Update a record using a given kay/payload pair.
   *
   * @param key A target key.
   * @param payload A target payload.
   * @param key_len The length of a target key.
   * @param pay_len The length of a target payload.
   * @retval kSuccess if a given record is updated.
   * @retval kKeyNotExist if a target key is not in this tree.
   */
  auto
  Update(  //
      const Key &key,
      const Payload &payload,
      const size_t key_len = sizeof(Key),
      const size_t pay_len = sizeof(Payload))  //
      -> ReturnCode
  {
    [[maybe_unused]] const auto &gc_guard = gc_.CreateEpochGuard();

    const auto *pay_addr = GetSrcAddr(payload);
    std::vector<INode *> stack{};
    stack.reserve(kInitialHeight);
    while (true) {
      auto &&[node, guard] = SearchNode<LCheckGuard>(stack, key);
      const auto [found, deleted, pos] = node->CheckUniqueness(key);
      if (!found || deleted) {
        if constexpr (kUseOptCCForLeaf) {
          if (!guard.VerifyVersion(kInsDelMask)) continue;
        }
        return kKeyNotExist;
      }

      LXGuard x_guard;
      if constexpr (kUseOptCCForLeaf) {
        x_guard = guard.TryLockX(kInsDelMask);
        if (!x_guard) continue;  // another thread may insert the key
      } else {
        x_guard = guard.UpgradeToX();
      }
      const auto rc = node->Update(x_guard, pos, key, key_len, pay_addr, pay_len);
      if (rc == kCompleted) return kSuccess;

      const auto &[sep_key, sep_key_len, r_node] = Split(node, std::move(x_guard));
      AddDownLink(stack, sep_key, sep_key_len, node, r_node);
      stack.emplace_back(std::bit_cast<INode *>(node));
    }
  }

  /**
   * @brief Delete a record using a given kay.
   *
   * @param key A target key.
   * @param key_len The length of a target key.
   * @retval kSuccess if a given record is deleted.
   * @retval kKeyNotExist if a target key is not in this tree.
   */
  auto
  Delete(  //
      const Key &key,
      [[maybe_unused]] const size_t key_len = sizeof(Key))  //
      -> ReturnCode
  {
    [[maybe_unused]] const auto &gc_guard = gc_.CreateEpochGuard();

    std::vector<INode *> stack{};
    stack.reserve(kInitialHeight);
    while (true) {
      auto &&[node, guard] = SearchNode<LCheckGuard>(stack, key);
      const auto [found, deleted, pos] = node->CheckUniqueness(key);
      if (!found || deleted) {
        if constexpr (kUseOptCCForLeaf) {
          if (!guard.VerifyVersion(kInsDelMask)) continue;
        }
        return kKeyNotExist;
      }

      LXGuard x_guard;
      if constexpr (kUseOptCCForLeaf) {
        x_guard = guard.TryLockX(kInsDelMask);
        if (!x_guard) continue;  // another thread may insert the key
      } else {
        x_guard = guard.UpgradeToX();
      }
      const auto rc = node->Delete(x_guard, pos);
      if (rc == kNeedMerge) {
        TryMerge(stack, node, x_guard.DowngradeToSIX());
      }
      return kSuccess;
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
    gc_.AddGarbage<Page>(old_root);

    return ReturnCode::kSuccess;
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

  /// @brief The expected capacity of a node when performing bulkload.
  static constexpr size_t kNodeCapacity = kPageSize / (sizeof(Key) + sizeof(Payload) + kMetaSize);

  /*##########################################################################*
   * Internal utility functions
   *##########################################################################*/

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

      if constexpr (std::is_same_v<Guard, typename NodeT::SIXGuard>) {
        return !removed && included;  // do not use side links with the SIX lock
      } else {
        if (!removed && included) return true;
        if (i++ > 0) return false;  // the traversal path may be too old
      }

      node = sib_node;
      guard = node->template GetGuard<Guard>();
    }
  }

  /**
   * @brief Acquire a guard instance for a target node.
   *
   * @tparam Guard A desired guard class.
   * @tparam NodeT A target node class.
   * @param[in,out] node A target node.
   * @param[in] key A search key.
   * @return A guard instance for reading/modifying a target node.
   */
  template <class Guard, class NodeT>
  [[nodiscard]] static auto
  AcquireGuard(  //
      NodeT *&node,
      const Key &key)  //
      -> Guard
  {
    constexpr auto kNeedX = std::is_same_v<Guard, typename NodeT::XGuard>;
    using TmpGuard = std::conditional_t<kNeedX, typename NodeT::CheckGuard, Guard>;

    auto &&guard = node->template GetGuard<TmpGuard>();
    while (SearchHorizontally(node, guard, key)) {
      if constexpr (!std::is_same_v<Guard, typename NodeT::XGuard>) {
        return guard;
      } else if constexpr (NodeT::kUseOptCC) {
        auto &&x_guard = guard.TryLockX(kSMOMask);
        if (x_guard) return x_guard;
      } else {
        return guard.UpgradeToX();
      }
    }
    return Guard{};
  }

  /**
   * @return A page for nodes.
   */
  [[nodiscard]] auto
  GetNodePage()  //
      -> void *
  {
    auto *page = gc_.GetPageIfPossible<Page>();
    return page ? page : ::dbgroup::memory::Allocate<Page>();
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
  template <class Guard, bool kIsInner = false>
  [[nodiscard]] auto
  SearchNode(  //
      std::vector<INode *> &stack,
      const Key &key,
      const size_t level = 0) const  //
      -> std::pair<std::conditional_t<kIsInner, INode *, LNode *>, Guard>
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
          auto *out_node = std::bit_cast<std::conditional_t<kIsInner, INode *, LNode *>>(node);
          auto &&guard = AcquireGuard<Guard>(out_node, key);
          if (guard) return {out_node, std::move(guard)};
          break;
        }

        INode *child{};
        auto &&guard = node->template GetGuard<IOptGuard>();
        while (true) {
          if (!SearchHorizontally(node, guard, key)) goto out;
          child = node->SearchChild(key);
          if (guard.VerifyVersion(kInsDelMask)) break;
        }
        stack.emplace_back(node);
        node = child;
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

      INode *child{};
      auto &&guard = node->template GetGuard<IOptGuard>();
      do {
        child = node->GetChild(0);
      } while (!guard.VerifyVersion(kInsDelMask));
      stack.emplace_back(node);
      node = child;
    }
  }

  /**
   * @brief Go to the next node for scanning.
   *
   * @param[in,out] node A target node.
   * @param[in,out] guard The guard instance of a target node.
   * @return The begin position for scanning.
   */
  [[nodiscard]] auto
  SiblingScan(  //
      LNode *&node,
      LReadGuard &guard) const  //
      -> size_t
  {
    size_t begin_pos{};
    if constexpr (kUseOptCCForLeaf) {
      auto *sib_node = node->GetSibNode();
      const auto &key = node->GetSeparatorKey().first;

      std::vector<INode *> stack{};
      stack.reserve(kInitialHeight);
      do {
        stack.emplace_back(sib_node);
        std::tie(sib_node, guard) = SearchNode<LReadGuard>(stack, key);
        std::memcpy(static_cast<void *>(node), sib_node, kPageSize);
      } while (!guard.VerifyVersion());

      auto [found, deleted, pos] = node->CheckUniqueness(key);
      begin_pos = pos + static_cast<size_t>(!found || deleted);
      if constexpr (IsVarLenData<Key>()) {
        ::dbgroup::memory::Release<KeyWOPtr>(key);
      }
    } else {
      node = node->GetSibNode();
      guard = node->template GetGuard<LReadGuard>();
    }
    return begin_pos;
  }

  /*##########################################################################*
   * Internal structure modification operations
   *##########################################################################*/

  /**
   * @brief Split a given node.
   *
   * @tparam NodeT A target node class.
   * @param stack A stack of traversed nodes.
   * @param l_node A split-left node.
   * @param l_guard The guard instance of a split-left node.
   * @retval 1st: A separator key.
   * @retval 2nd: The size of a separator key.
   * @retval 3rd: A split-right node.
   */
  template <class NodeT>
  [[nodiscard]] auto
  Split(  //
      NodeT *l_node,
      typename NodeT::XGuard l_guard)  //
      -> std::tuple<Key, size_t, NodeT *>
  {
    auto *r_node = new (GetNodePage()) NodeT{l_node};
    if constexpr (NodeT::kUseOptCC) {
      VerIncrement<kSMOMask>(l_guard);
    }
    auto &&[sep_key, sep_key_len] = l_node->GetSeparatorKey();
    return {sep_key, sep_key_len, r_node};
  }

  /**
   * @brief Add a down link to complete split operations.
   *
   * @param key A new key.
   * @param key_len The length of a new key.
   * @param l_child A source (left-split) node.
   * @param r_child A new (right-split) node to be inserted.
   * @param level A level where a target node is.
   */
  void
  AddDownLink(  //
      std::vector<INode *> stack,
      const Key &key,
      const size_t key_len,
      const void *l_child,
      const void *r_child,
      const size_t level = 1)
  {
    while (true) {
      if (stack.empty() && TryAddNewRoot(level, key, key_len, l_child, r_child)) break;

      // insert a new down link
      auto &&[node, guard] = SearchNode<IOptGuard, kInnerFlag>(stack, key, level);
      const auto [found, _, pos] = node->CheckUniqueness(key);  // `found` must be false
      auto &&x_guard = guard.TryLockX(kInsDelMask);
      if (!x_guard) continue;  // another thread may insert the key
      const auto rc = node->Insert(x_guard, pos, found, key, key_len, &r_child, kPtrSize);
      if (rc == kCompleted) break;

      const auto &[sep_key, sep_key_len, r_node] = Split(node, std::move(x_guard));
      AddDownLink(stack, sep_key, sep_key_len, node, r_node, level + 1);
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
  TryAddNewRoot(  //
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
    auto &&[r_child, r_guard] = l_child->GetMergeableSib();
    if (!r_child) {
      if constexpr (std::is_same_v<NodeT, INode>) {
        TryRemoveRoot(l_child, std::move(l_guard));
      }
      return;
    }

    const auto &key = l_child->GetSeparatorKey().first;
    while (true) {
      auto &&[node, guard] = SearchNode<IOptGuard, kInnerFlag>(stack, key, level);
      const auto [found, _, pos] = node->CheckUniqueness(key);
      if (!found) {
        if (!guard.VerifyVersion(kInsDelMask)) continue;
        break;  // a downlink has not been inserted yet or is leftmost in a node
      }
      auto &&x_guard = guard.TryLockX(kInsDelMask);
      if (!x_guard) continue;  // another thread may modify a node

      const auto rc = node->Delete(x_guard, pos);
      if (rc == kCompleted) {
        static_cast<void>(std::move(x_guard));
        l_child->Merge(std::move(l_guard), r_child, std::move(r_guard));
      } else {
        auto &&six_guard = x_guard.DowngradeToSIX();
        l_child->Merge(std::move(l_guard), r_child, std::move(r_guard));
        TryMerge(stack, node, std::move(six_guard), level);
      }
      gc_.AddGarbage<Page>(r_child);
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
      auto *new_root = node->GetChild(0);
      if (root_.compare_exchange_strong(root, new_root, kRelease, kRelaxed)) {
        node->Remove(std::move(guard), new_root);
        gc_.AddGarbage<Page>(root);
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
   * Internal member variables
   *##########################################################################*/

  /// @brief A garbage collector.
  GC gc_{};

  /// @brief The root node of this tree.
  std::atomic<INode *> root_{new (GetNodePage()) INode{}};
};

}  // namespace dbgroup::index::b_tree::component

#endif  // B_TREE_DBGROUP_B_TREE_COMPONENT_B_TREE_SL_HPP_
