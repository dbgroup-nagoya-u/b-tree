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

#ifndef B_TREE_DBGROUP_B_TREE_B_TREE_HPP_
#define B_TREE_DBGROUP_B_TREE_B_TREE_HPP_

// C++ standard libraries
#include <cstddef>
#include <cstdint>
#include <functional>
#include <future>
#include <memory>
#include <optional>
#include <tuple>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <vector>

// external libraries
#include <dbgroup/constants.hpp>
#include <dbgroup/index/utility.hpp>
#include <dbgroup/lock/mcs_lock.hpp>
#include <dbgroup/lock/optimistic_lock.hpp>
#include <dbgroup/lock/pessimistic_lock.hpp>
#include <dbgroup/lock/utility.hpp>
#include <dbgroup/memory/epoch_based_gc.hpp>
#include <dbgroup/memory/utility.hpp>
#include <dbgroup/thread/epoch_guard.hpp>

// local sources
#include "dbgroup/b_tree/component/common.hpp"
#include "dbgroup/b_tree/component/node.hpp"
#include "dbgroup/b_tree/utility.hpp"

namespace dbgroup::index::b_tree
{
/**
 * @brief A class for representing B+trees based on single-level SMOs.
 *
 * @tparam Key A class of stored keys.
 * @tparam Payload A class of stored payloads.
 * @tparam Comp A comparator class for keys.
 * @tparam GC A garbage collector for reusing node pages.
 * @tparam Lock A class for locking each node.
 * @tparam Delta A class for representing delta updates.
 */
template <class Key,
          class Payload,
          class Comp = std::less<Key>,
          class GC = ::dbgroup::memory::EpochBasedGC<>,
          ::dbgroup::lock::OptimisticallyLockable Lock = ::dbgroup::lock::OptimisticLock,
          class Delta = Payload>
class BTree
{
  /*##########################################################################*
   * Internal types
   *##########################################################################*/

  using KeyWOPtr = std::remove_pointer_t<Key>;
  using ScanKey = std::optional<std::tuple<Key, size_t, bool>>;

  using GCGuard = dbgroup::thread::EpochGuard;
  using OptGuard = typename Lock::OptGuard;
  using SIXGuard = typename Lock::SIXGuard;
  using XGuard = typename Lock::XGuard;

  using RC = component::NodeRC;
  using Node = component::Node<Key, Comp, Lock>;

  template <class Entry>
  using BulkIter = typename std::vector<Entry>::const_iterator;
  using NodeEntry = std::tuple<Key, Node *, size_t>;
  using BulkResult = std::pair<size_t, std::vector<NodeEntry>>;
  using BulkPromise = std::promise<BulkResult>;
  using BulkFuture = std::future<BulkResult>;

 public:
  /*##########################################################################*
   * Public types
   *##########################################################################*/

  /**
   * @brief A class for iterating through scan results.
   *
   * @tparam kForward A flag for indicating forward or backward order.
   * @note The APIs of this class are not thread-safe.
   */
  template <bool kForward = true>
  class Iterator
  {
   public:
    /*########################################################################*
     * Public constructors and assignment operators
     *########################################################################*/

    constexpr Iterator() = default;

    /**
     * @brief Construct a new iterator object.
     *
     * @param index A source index structure.
     * @param node The current scanning node.
     * @param begin_pos The begin position for scanning.
     * @param end_key  The end key given from a user.
     * @param gc_grd  The guard instance for preventing GC from reclaiming nodes.
     * @param stack A stack of traversed nodes.
     */
    Iterator(  //
        BTree *index,
        Node *node,
        const int32_t begin_pos,
        const ScanKey &end_key,
        GCGuard gc_grd,
        std::vector<std::tuple<Node *, int32_t, OptGuard>> stack = {})
        : node_{node},
          pos_{begin_pos},
          index_{index},
          end_key_{end_key},
          gc_grd_{std::move(gc_grd)},
          stack_{std::move(stack)}
    {
      node_->template SearchEndPosition<kForward>(end_key_, is_end_, end_pos_);
      FetchRecord();
    }

    Iterator(  //
        Iterator &&obj) noexcept
        : node_{std::exchange(obj.node_, nullptr)},
          pos_{obj.pos_},
          end_pos_{obj.end_pos_},
          is_end_{obj.is_end_},
          payload_{obj.payload_},
          index_{obj.index_},
          end_key_{obj.end_key_},
          gc_grd_{std::move(obj.gc_grd_)},
          stack_{std::move(obj.stack_)}
    {
      std::memcpy(key_, obj.key_, kMaxKeySize);
    }

    auto
    operator=(                    //
        Iterator &&rhs) noexcept  //
        -> Iterator &
    {
      node_ = std::exchange(rhs.node_, nullptr);
      pos_ = rhs.pos_;
      end_pos_ = rhs.end_pos_;
      is_end_ = rhs.is_end_;
      payload_ = rhs.payload_;
      index_ = rhs.index_;
      end_key_ = rhs.end_key_;
      gc_grd_ = std::move(rhs.gc_grd_);
      stack_ = std::move(rhs.stack_);
      std::memcpy(key_, rhs.key_, kMaxKeySize);
      return *this;
    }

    // forbit copying
    Iterator(const Iterator &) = delete;
    auto operator=(const Iterator &) -> Iterator & = delete;

    /*########################################################################*
     * Public destructors
     *########################################################################*/

    ~Iterator()
    {
      if (index_ && node_) {
        index_->FreePage(node_);
      }
    }

    /*########################################################################*
     * Public operators for iterators
     *########################################################################*/

    /**
     * @retval true if this iterator indicates a live record.
     * @retval false otherwise.
     */
    [[nodiscard]]
    constexpr explicit
    operator bool() const noexcept
    {
      return node_;
    }

    /**
     * @retval 1st: A key indicated by the iterator.
     * @retval 2nd: A payload indicated by the iterator.
     */
    [[nodiscard]]
    constexpr auto
    operator*() const noexcept  //
        -> std::pair<Key, Payload>
    {
      return {GetKey(), payload_};
    }

    /**
     * @brief Forward the iterator.
     *
     */
    void
    operator++()
    {
      if constexpr (kForward) {
        ++pos_;
      } else {
        --pos_;
      }
      FetchRecord();
    }

    /*########################################################################*
     * Public APIs
     *########################################################################*/

    /**
     * @return A key indicated by the iterator.
     */
    [[nodiscard]]
    constexpr auto
    GetKey() const noexcept  //
        -> Key
    {
      auto *addr = std::bit_cast<void *>(&key_[0]);
      if constexpr (IsVarLenData<Key>()) {
        return std::bit_cast<Key>(addr);
      } else {
        return *std::bit_cast<Key *>(addr);
      }
    }

    /**
     * @return A payload indicated by the iterator.
     */
    [[nodiscard]]
    constexpr auto
    GetPayload() const noexcept  //
        -> Payload
    {
      return payload_;
    }

   private:
    /*########################################################################*
     * Internal constants
     *########################################################################*/

    /// @brief The maximum size of keys.
    static constexpr size_t kMaxKeySize = component::MaxSize<Key>();

    /*########################################################################*
     * Internal utilities
     *########################################################################*/

    /**
     * @brief Fetch the record of the current position from a node.
     *
     */
    void
    FetchRecord()
    {
      while (true) {
        if constexpr (kForward) {
          if (pos_ < end_pos_) {
            if (node_->CopyRecordTo(pos_, key_, &payload_)) return;
            ++pos_;
            continue;
          }
        } else {
          if (pos_ >= end_pos_) {
            if (node_->CopyRecordTo(pos_, key_, &payload_)) return;
            --pos_;
            continue;
          }
        }

        if (is_end_) {
          index_->FreePage(node_);
          node_ = nullptr;
          return;
        }

        if constexpr (kForward) {
          pos_ = index_->SiblingScan(node_);
        } else {
          pos_ = index_->SiblingScanReverse(node_, stack_);
        }
        node_->template SearchEndPosition<kForward>(end_key_, is_end_, end_pos_);
      }
    }

    /*########################################################################*
     * Internal member variables
     *########################################################################*/

    /// @brief The current scanning node.
    Node *node_{};

    /// @brief The position of the current record.
    int32_t pos_{};

    /// @brief The end position of records in the node.
    int32_t end_pos_{};

    /// @brief A flag for indicating a current node is rightmost in scan-range.
    bool is_end_{};

    /// @brief The key of the current record.
    alignas(alignof(KeyWOPtr)) std::byte key_[kMaxKeySize]{};

    /// @brief The payload of the current record.
    Payload payload_{};

    /// @brief A source index structure.
    const BTree *index_{};

    /// @brief The end key given from a user.
    ScanKey end_key_{};

    /// @brief The guard instance for preventing GC from reclaiming nodes.
    GCGuard gc_grd_{};

    /// @brief A stack of traversed nodes.
    /// @note This filed is only used for backward scan.
    std::vector<std::tuple<Node *, int32_t, OptGuard>> stack_{};
  };

  /*##########################################################################*
   * Public constructors and assignment operators
   *##########################################################################*/

  /**
   * @brief Construct a new tree.
   *
   * @param page_size The page size for each node.
   * @param gc_interval_micro GC internal [us] (default: 10ms).
   * @param gc_thread_num the number of GC threads (default: 1).
   */
  explicit BTree(  //
      const size_t page_size = k1Ki,
      const size_t gc_interval_micro = kDefaultGCTime,
      const size_t gc_thread_num = kDefaultGCThreadNum)
      : page_size_{static_cast<uint16_t>(page_size)},
        gc_{std::make_shared<GC>(gc_interval_micro, gc_thread_num)}
  {
    VerifyPageSize();
  }

  /**
   * @brief Construct a new tree.
   *
   * @param gc A shared garbage collector for reusing internal pages.
   * @param page_size The page size for each node.
   */
  explicit BTree(  //
      const std::shared_ptr<GC> &gc,
      const size_t page_size = k1Ki)
      : page_size_{static_cast<uint16_t>(page_size)}, gc_{gc}
  {
    VerifyPageSize();
  }

  BTree(const BTree &) = delete;
  BTree(BTree &&) = delete;

  auto operator=(const BTree &) -> BTree & = delete;
  auto operator=(BTree &&) -> BTree & = delete;

  /*##########################################################################*
   * Public destructors
   *##########################################################################*/

  /**
   * @brief Destroy the object.
   *
   */
  ~BTree()
  {
    std::vector<std::pair<Node *, size_t>> stack{};
    stack.reserve(kInitialHeight);
    stack.emplace_back(root_.load(kAcquire), 0);
    while (!stack.empty()) {
      auto &[node, pos] = stack.back();
      if (node->GetLevel() > 0 && pos < node->GetRecNum()) {
        Node *child{};
        node->CopyPayloadTo(pos++, &child);
        stack.emplace_back(child, 0);
        continue;
      }
      FreePage(node);
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
      [[maybe_unused]] const size_t key_len = sizeof(Key))  //
      -> std::optional<Payload>
  {
    [[maybe_unused]] const auto &gc_grd = gc_->CreateEpochGuard();

    std::optional<Payload> out_pay{};
    std::vector<Node *> stack{};
    stack.reserve(kInitialHeight);
    while (true) {
      auto &&[node, grd] = SearchNode(stack, key);
      const auto [found, deleted, pos] = node->CheckUniqueness(key);
      if (found && !deleted) {
        out_pay = node->template GetPayload<Payload>(pos);
      } else {
        out_pay = std::nullopt;
      }
      if (grd.VerifyVersion(kNoMask)) return out_pay;
      stack.emplace_back(node);
    }
  }

  /**
   * @brief Perform a range scan in forward order with given keys.
   *
   * @param begin_key A pair of a begin key and its openness (true=closed).
   * @param end_key A pair of an end key and its openness (true=closed).
   * @return An iterator for accessing the records in forward order.
   */
  auto
  Scan(  //
      const ScanKey &begin_key = std::nullopt,
      const ScanKey &end_key = std::nullopt)
  {
    Node *node;
    int32_t begin_pos;
    auto &&gc_grd = gc_->CreateEpochGuard();

    std::vector<Node *> stack{};
    stack.reserve(kInitialHeight);
    auto *tmp_page = AllocPage();
    if (begin_key) {
      const auto &[key, _, closed] = *begin_key;
      while (true) {
        auto &&[tmp_node, grd] = SearchNode(stack, key);
        std::memcpy(tmp_page, tmp_node, page_size_);
        if (grd.VerifyVersion(kNoMask)) break;
        stack.emplace_back(tmp_node);
      }
      node = std::bit_cast<Node *>(tmp_page);
      const auto [found, deleted, pos] = node->CheckUniqueness(key);
      begin_pos = pos + static_cast<int32_t>(!found || deleted || !closed);
    } else {
      auto &&[tmp_node, grd] = SearchLeftmostLeaf(stack);
      do {
        std::memcpy(tmp_page, tmp_node, page_size_);
      } while (!grd.VerifyVersion(kNoMask));
      node = std::bit_cast<Node *>(tmp_page);
      begin_pos = 0;
    }

    return Iterator{this, node, begin_pos, end_key, std::move(gc_grd)};
  }

  /**
   * @brief Perform a range scan in backward order with given keys.
   *
   * @param begin_key A pair of a begin key and its openness (true=closed).
   * @param end_key A pair of an end key and its openness (true=closed).
   * @return An iterator for accessing the records in backward order.
   */
  auto
  ScanBackward(  //
      const ScanKey &begin_key = std::nullopt,
      const ScanKey &end_key = std::nullopt)
  {
    Node *node;
    int32_t begin_pos;
    auto &&gc_grd = gc_->CreateEpochGuard();

    std::vector<std::tuple<Node *, int32_t, OptGuard>> stack{};
    stack.reserve(kInitialHeight);
    auto *tmp_page = AllocPage();
    if (end_key) {
      const auto &[key, _, closed] = *end_key;
      while (true) {
        SearchNodeBackward(stack, key);
        auto &[tmp_node, _, grd] = stack.back();
        std::memcpy(tmp_page, tmp_node, page_size_);
        if (grd.VerifyVersion(kNoMask)) break;
      }
      node = std::bit_cast<Node *>(tmp_page);
      const auto [found, deleted, pos] = node->CheckUniqueness(key);
      begin_pos = pos - static_cast<int32_t>(!found || deleted || !closed);
    } else {
      while (true) {
        SearchRightmostLeaf(stack);
        auto &[tmp_node, _, grd] = stack.back();
        std::memcpy(tmp_page, tmp_node, page_size_);
        if (grd.VerifyVersion(kNoMask)) break;
      }
      node = std::bit_cast<Node *>(tmp_page);
      begin_pos = node->GetRecNum() - 1;
    }

    constexpr bool kBackward = false;
    return Iterator<kBackward>{
        this, node, begin_pos, begin_key, std::move(gc_grd), std::move(stack)};
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
   * @param merger A function for merging payloads.
   * @note This function is the alias of `Upsert`.
   */
  void
  Write(  //
      const Key &key,
      const Payload &payload,
      const size_t key_len = sizeof(Key),
      Payload (*merger)(const Payload &, const Payload &) = nullptr)
  {
    Upsert(key, payload, key_len, merger);
  }

  /**
   * @brief Upsert (update or insert) a given key/payload pair.
   *
   * @param key A target key.
   * @param payload A target payload.
   * @param key_len The length of a target key.
   * @param merger A function for merging payloads.
   * @retval The previous payload if exist (i.e., the function acts as update).
   * @retval std::nullopt otherwise.
   * @note If a given verifier finds inconsistent record modifications, this
   * function does nothing and returns std::nullopt.
   */
  auto
  Upsert(  //
      const Key &key,
      const Payload &payload,
      const size_t key_len = sizeof(Key),
      Payload (*merger)(const Payload &, const Payload &) = nullptr)  //
      -> std::optional<Payload>
  {
    [[maybe_unused]] const auto &gc_grd = gc_->CreateEpochGuard();

    std::optional<Payload> out_pay{};
    std::vector<Node *> stack{};
    stack.reserve(kInitialHeight);
    while (true) {
      auto &&[node, grd] = SearchNode(stack, key);
      const auto [found, deleted, pos] = node->CheckUniqueness(key);
      auto &&x_grd = grd.TryLockX(kInsDelMask);
      if (x_grd) {
        if (found && !deleted) {
          out_pay = node->Update(pos, payload, merger);
          break;
        }
        const auto rc = node->Insert(pos, found, key, key_len, &payload, kPayLen);
        if (rc == RC::kCompleted) {
          VerIncrement<kInsDelMask>(x_grd);
          break;
        }
        Split(stack, node, std::move(x_grd));
      }
      stack.emplace_back(std::bit_cast<Node *>(node));
    }
    return out_pay;
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
    [[maybe_unused]] const auto &gc_grd = gc_->CreateEpochGuard();

    std::optional<Payload> out_pay{};
    std::vector<Node *> stack{};
    stack.reserve(kInitialHeight);
    while (true) {
      auto &&[node, grd] = SearchNode(stack, key);
      const auto [found, deleted, pos] = node->CheckUniqueness(key);
      if (found && !deleted) {
        out_pay = node->template GetPayload<Payload>(pos);
        if (grd.VerifyVersion(kInsDelMask)) break;
      } else {
        auto &&x_grd = grd.TryLockX(kInsDelMask);
        if (x_grd) {
          const auto rc = node->Insert(pos, found, key, key_len, &payload, kPayLen);
          if (rc == RC::kCompleted) {
            VerIncrement<kInsDelMask>(x_grd);
            out_pay = std::nullopt;
            break;
          }
          Split(stack, node, std::move(x_grd));
        }
      }
      stack.emplace_back(std::bit_cast<Node *>(node));
    }
    return out_pay;
  }

  /**
   * @brief Update a record using a given kay/payload pair.
   *
   * @tparam T A delta-payload class (`Payload` or `Delta`).
   * @param key A target key.
   * @param payload A target payload.
   * @param key_len The length of a target key.
   * @param merger A function for merging payloads.
   * @param pay_len The length of a target payload.
   * @retval The previous payload if exist (i.e., succeeded).
   * @retval std::nullopt otherwise (i.e., failed).
   */
  template <class T = Payload>
  auto
  Update(  //
      const Key &key,
      const std::type_identity_t<T> &payload,
      [[maybe_unused]] const size_t key_len = sizeof(Key),
      Payload (*merger)(const Payload &, const T &) = nullptr,
      [[maybe_unused]] const size_t pay_len = sizeof(T))  //
      -> std::optional<Payload>
  {
    static_assert(                                               //
        std::is_same_v<T, Payload> || std::is_same_v<T, Delta>,  //
        "The update API requires the pre-defined Payload or Delta types.");

    [[maybe_unused]] const auto &gc_grd = gc_->CreateEpochGuard();

    std::optional<Payload> out_pay{};
    std::vector<Node *> stack{};
    stack.reserve(kInitialHeight);
    while (true) {
      auto &&[node, grd] = SearchNode(stack, key);
      const auto [found, deleted, pos] = node->CheckUniqueness(key);
      if (!found || deleted) {
        if (grd.VerifyVersion(kInsDelMask)) break;
      } else {
        auto &&x_grd = grd.TryLockX(kInsDelMask);
        if (x_grd) {
          out_pay = node->Update(pos, payload, merger);
          break;
        }
      }
      stack.emplace_back(std::bit_cast<Node *>(node));
    }
    return out_pay;
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
    [[maybe_unused]] const auto &gc_grd = gc_->CreateEpochGuard();

    std::optional<Payload> out_pay{};
    std::vector<Node *> stack{};
    stack.reserve(kInitialHeight);
    while (true) {
      auto &&[node, grd] = SearchNode(stack, key);
      const auto [found, deleted, pos] = node->CheckUniqueness(key);
      if (!found || deleted) {
        if (grd.VerifyVersion(kInsDelMask)) break;
      } else {
        auto &&x_grd = grd.TryLockX(kInsDelMask);
        if (x_grd) {
          Payload tmp_pay{};
          const auto rc = node->Delete(pos, &tmp_pay);
          out_pay = tmp_pay;
          VerIncrement<kInsDelMask>(x_grd);
          if (rc == RC::kNeedMerge) {
            TryMerge(stack, node, x_grd.DowngradeToSIX());
          }
          break;
        }
      }
      stack.emplace_back(std::bit_cast<Node *>(node));
    }
    return out_pay;
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
      nodes.reserve(2 * (page_size_ / kRecLen) * thread_num);
      Node *prev_node = nullptr;
      for (auto &&[p_level, p_nodes] : trees) {
        while (p_level < level) {
          p_nodes = ConstructSingleLevel<NodeEntry>(p_nodes.cbegin(), p_nodes.size(), p_level++);
        }
        nodes.insert(nodes.end(), p_nodes.begin(), p_nodes.end());
        Node::LinkVerticalBorderNodes(prev_node, std::get<1>(p_nodes.front()));
        prev_node = std::get<1>(p_nodes.back());
      }
    }

    // create upper layers until a root node is created
    while (nodes.size() > 1) {
      nodes = ConstructSingleLevel<NodeEntry>(nodes.cbegin(), nodes.size(), level++);
    }
    auto *new_root = std::get<1>(nodes.front());
    Node::RemoveLeftmostKeys(new_root);

    // set a new root
    auto *old_root = root_.exchange(new_root, kRelease);
    gc_->AddGarbage(old_root, page_size_);
  }

 private:
  /*##########################################################################*
   * Internal constants
   *##########################################################################*/

  /// @brief The initial capacity of a node stack.
  static constexpr size_t kInitialHeight = 8;

  /// @brief The size of payloads.
  static constexpr size_t kPayLen = sizeof(Payload);

  /// @brief The expected size of records.
  static constexpr size_t kRecLen = (sizeof(Key) + kPayLen + kWordSize);

  /*##########################################################################*
   * Internal utility functions
   *##########################################################################*/

  /**
   * @return A page for a new node.
   */
  [[nodiscard]]
  auto
  AllocPage()  //
      -> void *
  {
    auto *page = gc_->GetPageIfPossible(page_size_);
    if (!page) {
      page = ::operator new(page_size_, static_cast<std::align_val_t>(align_val_));
    }
    return page;
  }

  /**
   * @param page A page for a node.
   */
  void
  FreePage(  //
      void *page) const
  {
    ::operator delete(page, static_cast<std::align_val_t>(align_val_));
  }

  /**
   * @brief Search a target node horizontally.
   *
   * @param[in,out] node A target node.
   * @param[in,out] grd The guard instance of a target node.
   * @param[in] key A search key.
   * @param[in] reverse A flag for searching in reverse order.
   * @retval true if a container node has been found.
   * @retval false otherwise.
   */
  [[nodiscard]]
  static auto
  SearchHorizontally(  //
      Node *&node,
      OptGuard &grd,
      const Key &key,
      const bool reverse = false)  //
      -> bool
  {
    size_t i = 0;
    while (true) {
      auto *sib_node = node->GetSibNode();
      const auto removed = node->Removed();
      const auto included = node->Include(key, reverse);
      if (!grd.VerifyVersion(kSMOMask)) continue;

      if (!removed && included) return true;
      if (!sib_node || i++ > 0) return false;  // the traversed path may be too old

      node = sib_node;
      grd = node->GetGuard();
    }
  }

  /**
   * @brief Search a node that may have a target key.
   *
   * @param[in,out] stack A stack of traversed nodes.
   * @param[in] key A search key.
   * @param[in] level The bottom level of search space.
   * @retval 1st: A found node.
   * @retval 2nd: The lock guard of a found node.
   */
  [[nodiscard]]
  auto
  SearchNode(  //
      std::vector<Node *> &stack,
      const Key &key,
      const size_t level = 0) const  //
      -> std::pair<Node *, OptGuard>
  {
    Node *node;
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
        Node *child{};
        auto &&grd = node->GetGuard();
        while (true) {
          if (!SearchHorizontally(node, grd, key)) goto out;
          if (node->GetLevel() == level) return {node, std::move(grd)};

          child = node->SearchChild(key);
          if (grd.VerifyVersion(kInsDelMask)) break;
        }
        stack.emplace_back(std::exchange(node, child));
      }
    out:;  // the current node has been removed or too old, so retry
    }
  }

  /**
   * @brief Search a node that may have a target key.
   *
   * @param[in,out] stack A stack of traversed nodes.
   * @param[in] key A search key.
   * @note A found node is retained in the stack.
   */
  void
  SearchNodeBackward(  //
      std::vector<std::tuple<Node *, int32_t, OptGuard>> &stack,
      const Key &key) const
  {
    constexpr bool kBackward = true;
    Node *child{};
    while (true) {
      if (stack.empty()) {
        auto *root = root_.load(kAcquire);
        stack.emplace_back(root, 0, root->GetGuard());
      }

      while (true) {
        auto &[node, pos, grd] = stack.back();
        while (true) {
          if (!SearchHorizontally(node, grd, key, kBackward)) goto out;
          if (node->GetLevel() == 0) return;

          pos = node->SearchRecord(key).second - 1;
          node->CopyPayloadTo(pos--, &child);
          if (grd.VerifyVersion(kInsDelMask)) break;
        }
        stack.emplace_back(child, 0, child->GetGuard());
      }
    out:  // the current node has been removed or too old, so retry
      stack.pop_back();
    }
  }

  /**
   * @brief Search the leftmost leaf node.
   *
   * @param[in,out] stack A stack of traversed nodes.
   * @retval 1st: The leftmost node.
   * @retval 2nd: The lock guard of the leftmost node.
   */
  [[nodiscard]]
  auto
  SearchLeftmostLeaf(                    //
      std::vector<Node *> &stack) const  //
      -> std::pair<Node *, OptGuard>
  {
    auto *node = root_.load(kAcquire);
    while (true) {
      if (node->GetLevel() == 0) return {node, node->GetGuard()};

      Node *child{};
      auto &&grd = node->GetGuard();
      do {
        node->CopyPayloadTo(0, &child);
      } while (!grd.VerifyVersion(kInsDelMask));
      stack.emplace_back(std::exchange(node, child));
    }
  }

  /**
   * @brief Search the rightmost leaf node.
   *
   * @param[in,out] stack A stack of traversed nodes.
   */
  void
  SearchRightmostLeaf(  //
      std::vector<std::tuple<Node *, int32_t, OptGuard>> &stack) const
  {
    while (true) {
      if (stack.empty()) {
        auto *root = root_.load(kAcquire);
        stack.emplace_back(root, 0, root->GetGuard());
      }

      while (true) {
        Node *child{};
        auto cnt = 0;
        auto &[node, pos, grd] = stack.back();
        while (true) {
          const auto removed = node->Removed();
          auto *sib = node->GetSibNode();
          if (!grd.VerifyVersion(kSMOMask)) continue;

          if (!removed && !sib) {
            if (node->GetLevel() == 0) return;
            pos = node->GetRecNum() - 1;
            node->CopyPayloadTo(pos--, &child);
            if (grd.VerifyVersion(kInsDelMask)) break;
            continue;
          }

          if ((removed && !sib) || (sib && cnt++ > 0)) goto out;
          node = sib;
          grd = node->GetGuard();
        }
        stack.emplace_back(child, 0, child->GetGuard());
      }
    out:  // the current node has been removed or too old, so retry
      stack.pop_back();
    }
  }

  /**
   * @brief Go to the next node for scanning.
   *
   * @param[in,out] node A target node.
   * @return The begin position for scanning.
   */
  [[nodiscard]]
  auto
  SiblingScan(            //
      Node *&node) const  //
      -> size_t
  {
    OptGuard grd{};
    auto *sib_node = node->GetSibNode();
    const auto &key = node->GetHighKey().first;

    std::vector<Node *> stack{};
    stack.reserve(kInitialHeight);
    do {
      stack.emplace_back(sib_node);
      std::tie(sib_node, grd) = SearchNode(stack, key);
      std::memcpy(static_cast<void *>(node), sib_node, page_size_);
    } while (!grd.VerifyVersion(kNoMask));

    auto [found, deleted, pos] = node->CheckUniqueness(key);
    if constexpr (IsVarLenData<Key>()) {
      ::dbgroup::memory::Release<KeyWOPtr>(key);
    }
    return pos + static_cast<size_t>(!found || deleted);
  }

  /**
   * @brief Go to the previous node for scanning.
   *
   * @param[in,out] node A target node.
   * @param[in,out] stack A stack of traversed nodes.
   * @return The begin position for scanning.
   */
  [[nodiscard]]
  auto
  SiblingScanReverse(  //
      Node *&node,
      std::vector<std::tuple<Node *, int32_t, OptGuard>> &stack) const  //
      -> size_t
  {
    stack.pop_back();
    while (!stack.empty()) {
      auto &[node, pos, grd] = stack.back();
      if (pos >= 0) {
        Node *child{};
        node->CopyPayloadTo(pos--, &child);
        if (grd.VerifyVersion(kInsDelMask)) {
          stack.emplace_back(child, 0, child->GetGuard());
        }
        break;
      }
      stack.pop_back();
    }

    const auto &key = node->GetLowKey();
    while (true) {
      SearchNodeBackward(stack, key);
      auto &[tmp_node, _, grd] = stack.back();
      std::memcpy(static_cast<void *>(node), tmp_node, page_size_);
      if (grd.VerifyVersion(kNoMask)) break;
    }

    const auto pos = node->SearchRecord(key).second - 1;
    if constexpr (IsVarLenData<Key>()) {
      ::dbgroup::memory::Release<KeyWOPtr>(key);
    }
    return pos;
  }

  /*##########################################################################*
   * Internal structure modification operations
   *##########################################################################*/

  /**
   * @brief Split a given node and insert a new down link.
   *
   * @param stack A stack of traversed nodes.
   * @param l_child A split-left node.
   * @param l_grd The guard instance of a split-left node.
   */
  void
  Split(  //
      std::vector<Node *> stack,
      Node *l_child,
      XGuard l_grd)
  {
    auto *r_child = new (AllocPage()) Node{l_child};
    VerIncrement<kSMOMask>(l_grd);
    const auto &[key, key_len] = l_child->GetHighKey();
    l_grd = XGuard{};

    const auto level = l_child->GetLevel() + 1;
    while (true) {
      if (stack.empty()) {  // create a new root
        auto *old_root = root_.load(kRelaxed);
        if (old_root == std::bit_cast<Node *>(l_child)) {
          auto *root = new (AllocPage()) Node{level, key, key_len, l_child, r_child};
          if (root_.compare_exchange_strong(old_root, root, kRelease, kRelaxed)) break;
          FreePage(root);
        }
      }

      // insert a new down link
      auto &&[node, grd] = SearchNode(stack, key, level);
      const auto [found, pos] = node->SearchRecord(key);
      assert(!found);
      auto &&x_grd = grd.TryLockX(kInsDelMask);
      if (!x_grd) continue;  // another thread may insert the key
      const auto rc = node->Insert(pos, found, key, key_len, &r_child, kWordSize);
      if (rc == RC::kCompleted) {
        VerIncrement<kInsDelMask>(x_grd);
        break;
      }

      Split(stack, node, std::move(x_grd));
      stack.emplace_back(node);
    }
    if constexpr (IsVarLenData<Key>()) {
      ::dbgroup::memory::Release<KeyWOPtr>(key);
    }
  }

  /**
   * @brief Try to remove a down link and merge a given node.
   *
   * @param stack A stack of traversed nodes.
   * @param l_child A left child (to be merged) node.
   * @param l_grd The SIX guard instance of a given child node.
   */
  void
  TryMerge(  //
      std::vector<Node *> &stack,
      Node *l_child,
      SIXGuard l_grd)
  {
    const auto level = l_child->GetLevel() + 1;
    auto &&r_grd = l_child->GetMergeableSib();
    if (!r_grd) {  // remove a root node if possible
      auto *root = root_.load(kRelaxed);
      if (level > 1 && l_child == root && !l_child->GetSibNode() && l_child->GetRecNum() == 1) {
        l_child->CopyPayloadTo(0, &root);  // `root` has become a new root
        if (root_.compare_exchange_strong(l_child, root, kRelease, kRelaxed)) {
          l_child->Remove(std::move(l_grd));
          gc_->AddGarbage(l_child, page_size_);
        }
      }
      return;
    }

    const auto &key = l_child->GetHighKey().first;
    while (true) {
      auto &&[node, grd] = SearchNode(stack, key, level);
      const auto [found, pos] = node->SearchRecord(key);
      if (!found) {
        if (!grd.VerifyVersion(kInsDelMask)) continue;
        break;  // a downlink has not been inserted yet or is leftmost in a node
      }
      auto &&x_grd = grd.TryLockX(kInsDelMask);
      if (!x_grd) continue;  // another thread may modify a node

      Node *r_child{};
      const auto rc = node->Delete(pos, &r_child);
      VerIncrement<kInsDelMask>(x_grd);
      if (rc == RC::kCompleted) {
        x_grd = XGuard{};
        l_child->Merge(std::move(l_grd), r_child, std::move(r_grd));
      } else {
        auto &&six_grd = x_grd.DowngradeToSIX();
        l_child->Merge(std::move(l_grd), r_child, std::move(r_grd));
        TryMerge(stack, node, std::move(six_grd));
      }
      gc_->AddGarbage(r_child, page_size_);
      break;
    }
    if constexpr (IsVarLenData<Key>()) {
      ::dbgroup::memory::Release<KeyWOPtr>(key);
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
      if (n < 2 * (page_size_ / kRecLen)) break;
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
    nodes.reserve((n / (page_size_ / kRecLen)) + 1);
    const auto &iter_end = std::next(iter, n);
    for (Node *prev_node = nullptr; iter < iter_end;) {
      const auto &[key, key_len] = ParseKey(*iter);
      auto *node = new (AllocPage()) Node{page_size_, align_val_, level};
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

  void
  VerifyPageSize()
  {
    constexpr size_t kMaxKeySize = component::MaxSize<Key>();
    const size_t max_split_size = component::kHeaderSize  // NOLINTBEGIN
                                  + (page_size_ - component::kHeaderSize) * 0.75
                                  + (kMaxKeySize + kPayLen + kWordSize) * 0.5
                                  + kMaxKeySize;  // NOLINTEND
    if (max_split_size > page_size_) {
      throw std::runtime_error{"The page size is too small to store records."};
    }
    if (page_size_ < k512 || k64Ki < page_size_ || std::popcount(page_size_) > 1) {
      throw std::runtime_error{"The page size must be 512, 1024, 2048, ..., or 65536."};
    }
  }

  /*##########################################################################*
   * Internal member variables
   *##########################################################################*/

  /// @brief The page size for each node.
  const uint16_t page_size_{};

  /// @brief The alignment size for each node.
  const uint16_t align_val_{static_cast<uint16_t>(GetAlignValOnVirtualPages(page_size_))};

  /// @brief A garbage collector.
  std::shared_ptr<GC> gc_{};

  /// @brief The root node of this tree.
  std::atomic<Node *> root_{new (AllocPage()) Node{page_size_, align_val_}};
};

}  // namespace dbgroup::index::b_tree

#endif  // B_TREE_DBGROUP_B_TREE_B_TREE_HPP_
