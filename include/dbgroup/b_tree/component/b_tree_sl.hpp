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
  using IXGuard = typename INode::XGuard;

  using LNode = Node<Key, Comp, LeafLock, kUseOptCCForLeaf>;
  using LReadGuard = typename LNode::ReadGuard;
  using LCheckGuard = typename LNode::CheckGuard;
  using LXGuard = typename LNode::XGuard;

 public:
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
    const auto *pay_addr = GetSrcAddr(payload);
    std::vector<INode *> stack{};
    stack.reserve(kInitialHeight);
    while (true) {
      auto &&[node, guard] = SearchNode<LCheckGuard>(stack, key);
      const auto [found, pos] = node->SearchRecord(key);
      if (found) {
        if constexpr (kUseOptCCForLeaf) {
          if (!guard.VerifyVersion()) continue;
        }
        return kKeyExist;
      }

      LXGuard x_guard;
      if constexpr (kUseOptCCForLeaf) {
        x_guard = guard.TryLockX();
        if (!x_guard) continue;  // another thread may insert the key
      } else {
        x_guard = guard.UpgradeToX();
      }
      const auto rc = node->Insert(x_guard, pos, key, key_len, pay_addr, pay_len);
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
   * @retval kSuccess if a given record is inserted.
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
    const auto *pay_addr = GetSrcAddr(payload);
    std::vector<INode *> stack{};
    stack.reserve(kInitialHeight);
    while (true) {
      auto &&[node, guard] = SearchNode<LCheckGuard>(stack, key);
      const auto [found, pos] = node->SearchRecord(key);
      if (!found) {
        if constexpr (kUseOptCCForLeaf) {
          if (!guard.VerifyVersion()) continue;
        }
        return kKeyNotExist;
      }

      LXGuard x_guard;
      if constexpr (kUseOptCCForLeaf) {
        x_guard = guard.TryLockX();
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
  static constexpr bool kIsInner = true;

  /*##########################################################################*
   * Internal utility functions
   *##########################################################################*/

  /**
   * @brief Search a target node horizontally.
   *
   * @tparam kIsInner A flag for indicating whether a target is inner nodes.
   * @tparam Guard A guard class for reading a target node.
   * @tparam NodeT A target node class.
   * @param[in,out] node A target node.
   * @param[in,out] guard The guard instance of a target node.
   * @param[in] key A search key.
   * @retval true if a container node has been found.
   * @retval false otherwise.
   */
  template <bool kIsInner, class Guard, class NodeT>
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
      if constexpr (kIsInner || kUseOptCCForLeaf) {
        if (!guard.VerifyVersion()) continue;
      }

      if (!removed && included) return true;
      if ((removed && !sib_node) || i++ > 1) return false;

      node = sib_node;
      guard = node->template GetGuard<Guard>();
    }
  }

  /**
   * @brief Acquire a guard instance for a target node.
   *
   * @tparam kIsInner A flag for indicating whether a target is inner nodes.
   * @tparam Guard A desired guard class.
   * @tparam NodeT A target node class.
   * @param[in,out] node A target node.
   * @param[in] key A search key.
   * @return A guard instance for reading/modifying a target node.
   */
  template <bool kIsInner, class Guard, class NodeT>
  [[nodiscard]] static auto
  AcquireGuard(  //
      NodeT *&node,
      const Key &key)  //
      -> Guard
  {
    constexpr auto kRequireRead = std::is_same_v<Guard, LReadGuard>;
    using CheckGuardForL = std::conditional_t<kRequireRead, LReadGuard, LCheckGuard>;
    using CheckGuard = std::conditional_t<kIsInner, IOptGuard, CheckGuardForL>;

    auto &&guard = node->template GetGuard<CheckGuard>();
    while (SearchHorizontally<kIsInner>(node, guard, key)) {
      if constexpr (kIsInner || (kUseOptCCForLeaf && std::is_same_v<Guard, LXGuard>)) {
        auto &&x_guard = guard.TryLockX();
        if (x_guard) return x_guard;
      } else if constexpr (std::is_same_v<Guard, LXGuard>) {
        return guard.UpgradeToX();
      } else {
        return guard;
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
    return ::dbgroup::memory::Allocate<Page>();
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

      while (true) {
        if (node->GetLevel() == level) {  // a level region is immutable
          auto *out_node = std::bit_cast<std::conditional_t<kIsInner, INode *, LNode *>>(node);
          auto &&guard = AcquireGuard<kIsInner, Guard>(out_node, key);
          if (guard) return {out_node, std::move(guard)};
          break;
        }

        INode *child{};
        auto &&guard = node->template GetGuard<IOptGuard>();
        while (true) {
          if (!SearchHorizontally<kIsInner>(node, guard, key)) goto out;
          child = node->SearchChild(key);
          if (guard.VerifyVersion()) break;
        }
        stack.emplace_back(node);
        node = child;
      }
    out:;  // the current node has been removed or too old, so retry
    }
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
    if constexpr (!std::is_same_v<typename NodeT::OptGuard, void>) {
      l_guard.SetVersion((l_guard.GetVersion() & kSMOMask) + kSMOVerUnit);
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
      auto &&[node, x_guard] = SearchNode<IXGuard, kIsInner>(stack, key, level);
      const auto rc = node->Write(x_guard, key, key_len, &r_child, kPtrSize);
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

  /*##########################################################################*
   * Internal member variables
   *##########################################################################*/

  /// @brief The root node of this tree.
  std::atomic<INode *> root_{new (::dbgroup::memory::Allocate<Page>()) INode{}};
};

}  // namespace dbgroup::index::b_tree::component

#endif  // B_TREE_DBGROUP_B_TREE_COMPONENT_B_TREE_SL_HPP_
