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

#ifndef B_TREE_DBGROUP_B_TREE_RECORD_ITERATOR_HPP_
#define B_TREE_DBGROUP_B_TREE_RECORD_ITERATOR_HPP_

// C++ standard libraries
#include <deque>
#include <optional>
#include <utility>

// external sources
#include "dbgroup/constants.hpp"
#include "dbgroup/index/utility.hpp"
#include "dbgroup/lock/utility.hpp"
#include "dbgroup/memory/utility.hpp"

// local sources
#include "dbgroup/b_tree/component/common.hpp"
#include "dbgroup/b_tree/component/node.hpp"

namespace dbgroup::index::b_tree
{
/**
 * @brief A class for representing an iterator of scan results.
 *
 * @tparam Index A source index class.
 */
template <class Index, class Guard>
class RecordIterator
{
  /*##########################################################################*
   * Type aliases
   *##########################################################################*/

  using Key = typename Index::Key_t;
  using Payload = typename Index::Payload_t;
  using Node = typename Index::Node_t;
  using VerGuard = typename Node::VerGuard;
  using ScanKey = std::optional<std::tuple<Key, size_t, bool>>;
  using EpochGuard = ::dbgroup::thread::EpochGuard;

 public:
  /*##########################################################################*
   * Public constructors and assignment operators
   *##########################################################################*/

  constexpr RecordIterator() = default;

  /**
   * @brief Construct a new iterator object.
   *
   * @param index A source index structure.
   * @param node The current scanning node.
   * @param guard A guard instance for locking the current node.
   * @param begin_pos The begin position for scanning.
   * @param end_key  The end key given from a user.
   * @param gc_guard  The guard instance for preventing GC from reclaiming nodes.
   */
  RecordIterator(  //
      Index *index,
      Node *node,
      Guard guard,
      const size_t begin_pos,
      const ScanKey &end_key,
      EpochGuard gc_guard)
      : node_{node},
        guard_{std::move(guard)},
        pos_{begin_pos},
        index_{index},
        end_key_{end_key},
        gc_guard_{std::move(gc_guard)}
  {
    std::tie(is_end_, end_pos_) = node_->SearchEndPosition(end_key_);
    FetchRecord();
  }

  RecordIterator(  //
      RecordIterator &&obj) noexcept
      : node_{obj.node_},
        guard_{std::move(obj.guard_)},
        pos_{obj.pos_},
        end_pos_{obj.end_pos_},
        is_end_{obj.is_end_},
        payload_{obj.payload_},
        verifier_{std::move(obj.verifier_)},
        index_{obj.index_},
        end_key_{obj.end_key_},
        gc_guard_{std::move(obj.gc_guard_)}
  {
    std::memcpy(key_, obj.key_, kMaxKeySize);
  }

  auto
  operator=(                          //
      RecordIterator &&rhs) noexcept  //
      -> RecordIterator &
  {
    node_ = rhs.node_;
    guard_ = std::move(rhs.guard_);
    pos_ = rhs.pos_;
    end_pos_ = rhs.end_pos_;
    is_end_ = rhs.is_end_;
    payload_ = rhs.payload_;
    verifier_ = std::move(rhs.verifier_);
    index_ = rhs.index_;
    end_key_ = rhs.end_key_;
    gc_guard_ = std::move(rhs.gc_guard_);
    std::memcpy(key_, rhs.key_, kMaxKeySize);
    return *this;
  }

  // forbit copying
  RecordIterator(const RecordIterator &) = delete;
  auto operator=(const RecordIterator &) -> RecordIterator & = delete;

  /*##########################################################################*
   * Public destructors
   *##########################################################################*/

  ~RecordIterator() = default;

  /*##########################################################################*
   * Public operators for iterators
   *##########################################################################*/

  /**
   * @retval true if this iterator indicates a live record.
   * @retval false otherwise.
   */
  [[nodiscard]] constexpr explicit
  operator bool() const
  {
    return node_;
  }

  /**
   * @retval 1st: A key indicated by the iterator.
   * @retval 2nd: A payload indicated by the iterator.
   */
  [[nodiscard]] constexpr auto
  operator*() const  //
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
    ++pos_;
    FetchRecord();
  }

  /*##########################################################################*
   * Public APIs
   *##########################################################################*/

  /**
   * @return A key indicated by the iterator.
   */
  [[nodiscard]] constexpr auto
  GetKey() const  //
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
  [[nodiscard]] constexpr auto
  GetPayload() const  //
      -> Payload
  {
    return payload_;
  }

  /**
   * @brief Construct a verifier for snapshot reads or phantom avoidance.
   *
   */
  constexpr void
  PrepareVerifier()
  {
    verifier_ = std::make_optional<ScanVerifier>();
  }

  /**
   * @retval true if scanned records are a snapshot at some point.
   * @retval false otherwise.
   * @note If you have not prepared a verifier, this function always returns
   * true.
   */
  [[nodiscard]] auto
  VerifySnapshot()  //
      -> bool
  {
    if (!verifier_) return true;
    return (*verifier_).Verify(kNoMask);
  }

  /**
   * @retval true if scanned records do not have phantoms (insertion/deletion).
   * @retval false otherwise.
   * @note If you have not prepared a verifier, this function always returns
   * true.
   */
  [[nodiscard]] auto
  VerifyNoPhantom()  //
      -> bool
  {
    if (!verifier_) return true;
    return (*verifier_).Verify(kInsDelMask);
  }

 private:
  /*##########################################################################*
   * Internal types
   *##########################################################################*/

  /// @brief A type for allocating variable-length keys.
  using KeyWOPtr = std::remove_pointer_t<Key>;

  /**
   * @brief A class for verifying scan results.
   *
   */
  class ScanVerifier
  {
   public:
    /*########################################################################*
     * Public constructors and assignment operators
     *########################################################################*/

    constexpr ScanVerifier() = default;

    ScanVerifier(  //
        ScanVerifier &&obj) noexcept
        : guards_{std::move(obj.guards_)}
    {
    }

    auto
    operator=(                        //
        ScanVerifier &&rhs) noexcept  //
        -> ScanVerifier &
    {
      guards_ = std::move(rhs.guards_);
      return *this;
    }

    // forbit copying
    ScanVerifier(const ScanVerifier &) = delete;
    auto operator=(const ScanVerifier &) -> ScanVerifier & = delete;

    /*########################################################################*
     * Public destructors
     *########################################################################*/

    ~ScanVerifier() = default;

    /*########################################################################*
     * Public APIs
     *########################################################################*/

    /**
     * @brief Add a guard instance to this verifier.
     *
     * @param guard A guard instance to be added.
     */
    void
    Add(  //
        Guard guard)
    {
      guards_.emplace_back(VerGuard{std::move(guard)});
    }

    /**
     * @brief Perform verification based on registered guards.
     *
     * @param mask A bitmask for indicating target bits.
     * @retval true if this does not find changes using a given bitmask.
     * @retval false otherwise.
     */
    [[nodiscard]] auto
    Verify(                   //
        const uint32_t mask)  //
        -> bool
    {
      auto verified = true;
      for (auto &guard : guards_) {
        verified = guard.ImmediateVerify(mask);
        if (!verified) break;
      }
      return verified;
    }

   private:
    /*########################################################################*
     * Internal member variables
     *########################################################################*/

    /// @brief Guard instances for verification.
    std::deque<VerGuard> guards_{};
  };

  /*##########################################################################*
   * Internal constants
   *##########################################################################*/

  /// @brief The maximum size of keys.
  static constexpr size_t kMaxKeySize = component::MaxSize<Key>();

  /*##########################################################################*
   * Internal utilities
   *##########################################################################*/

  /**
   * @brief Fetch the record of the current position from a node.
   *
   */
  void
  FetchRecord()
  {
    while (true) {
      if (pos_ < end_pos_) {
        if (node_->CopyRecordTo(pos_, key_, &payload_)) return;
        ++pos_;
        continue;
      }

      if (is_end_) {
        node_ = nullptr;
        return;
      }

      if (verifier_) {
        (*verifier_).Add(std::move(guard_));
      }
      pos_ = index_->SiblingScan(node_, guard_);
      std::tie(is_end_, end_pos_) = node_->SearchEndPosition(end_key_);
    }
  }

  /*##########################################################################*
   * Internal member variables
   *##########################################################################*/

  /// @brief The current scanning node.
  Node *node_{};

  /// @brief A guard instance for locking the current node.
  Guard guard_{};

  /// @brief The position of the current record.
  size_t pos_{};

  /// @brief The end position of records in the node.
  size_t end_pos_{};

  /// @brief A flag for indicating a current node is rightmost in scan-range.
  bool is_end_{};

  /// @brief The key of the current record.
  alignas(alignof(KeyWOPtr)) std::byte key_[kMaxKeySize]{};

  /// @brief The payload of the current record.
  Payload payload_{};

  /// @brief A verifier for snapshot read or phantom avoidance.
  std::optional<ScanVerifier> verifier_{};

  /// @brief A source index structure.
  const Index *index_{};

  /// @brief The end key given from a user.
  ScanKey end_key_{};

  /// @brief The guard instance for preventing GC from reclaiming nodes.
  EpochGuard gc_guard_{};
};

}  // namespace dbgroup::index::b_tree

#endif  // B_TREE_DBGROUP_B_TREE_RECORD_ITERATOR_HPP_
