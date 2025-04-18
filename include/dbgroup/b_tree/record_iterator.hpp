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
template <class Index>
class RecordIterator
{
 public:
  /*##########################################################################*
   * Type aliases
   *##########################################################################*/

  using Key = typename Index::Key_t;
  using Payload = typename Index::Payload_t;
  using Node = typename Index::Node_t;
  using Guard = typename Node::ReadGuard;
  using ScanKey = std::optional<std::tuple<const Key &, size_t, bool>>;
  using EpochGuard = ::dbgroup::thread::EpochGuard;

  /*##########################################################################*
   * Public constructors and assignment operators
   *##########################################################################*/

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
      ScanKey end_key,
      EpochGuard gc_guard)
      : node_{node},
        guard_{std::move(guard)},
        pos_{begin_pos},
        index_{index},
        end_key_{std::move(end_key)},
        gc_guard_{std::move(gc_guard)}
  {
    std::tie(is_end_, end_pos_) = node_->SearchEndPosition(end_key_);
    FetchRecord();
  }

  RecordIterator(const RecordIterator &) = delete;
  RecordIterator(RecordIterator &&) = delete;

  auto operator=(const RecordIterator &) -> RecordIterator & = delete;
  auto operator=(RecordIterator &&) -> RecordIterator & = delete;

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
    return {GetKey(), GetPayload()};
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
    auto *addr = std::bit_cast<void *>(&payload_[0]);
    if constexpr (IsVarLenData<Payload>()) {
      return std::bit_cast<Payload>(addr);
    } else {
      return *std::bit_cast<Payload *>(addr);
    }
  }

 private:
  /*##########################################################################*
   * Internal types
   *##########################################################################*/

  /// @brief A type for allocating variable-length keys.
  using KeyWOPtr = std::remove_pointer_t<Key>;

  /// @brief A type for allocating variable-length payloads.
  using PayWOPtr = std::remove_pointer_t<Payload>;

  /*##########################################################################*
   * Internal constants
   *##########################################################################*/

  /// @brief A flag for using optimistic CC.
  static constexpr bool kUseOptCC = Node::kUseOptCC;

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
        if (node_->CopyRecordTo(pos_, key_, payload_)) return;
        ++pos_;
        continue;
      }

      if (is_end_) {
        node_ = nullptr;
        return;
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
  alignas(alignof(KeyWOPtr)) std::byte key_[component::MaxSize<Key>()]{};

  /// @brief The payload of the current record.
  alignas(alignof(PayWOPtr)) std::byte payload_[component::MaxSize<Payload>()]{};

  /// @brief A source index structure.
  const Index *index_{};

  /// @brief The end key given from a user.
  const ScanKey end_key_{};

  /// @brief The guard instance for preventing GC from reclaiming nodes.
  const EpochGuard gc_guard_{};
};

}  // namespace dbgroup::index::b_tree

#endif  // B_TREE_DBGROUP_B_TREE_RECORD_ITERATOR_HPP_
