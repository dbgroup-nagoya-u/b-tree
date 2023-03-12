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

#ifndef B_TREE_COMPONENT_OPTIMISTIC_RECORD_ITERATOR_HPP
#define B_TREE_COMPONENT_OPTIMISTIC_RECORD_ITERATOR_HPP

// C++ standard libraries
#include <optional>
#include <utility>

// external sources
#include "memory/epoch_based_gc.hpp"

// local sources
#include "b_tree/component/common.hpp"

namespace dbgroup::index::b_tree::component
{
/**
 * @brief A class for representing an iterator of scan results.
 *
 */
template <class Index>
class OptimisticRecordIterator
{
 public:
  /*####################################################################################
   * Type aliases
   *##################################################################################*/

  using Key = typename Index::K;
  using Payload = typename Index::V;
  using Node = typename Index::Node_t;
  using ScanKey = std::optional<std::tuple<const Key &, size_t, bool>>;
  using EpochGuard = ::dbgroup::memory::component::EpochGuard;

  /*####################################################################################
   * Public constructors and assignment operators
   *##################################################################################*/

  /**
   * @brief Construct a new iterator object.
   *
   * @param node a node to be scanned.
   * @param begin_key the begin key of this scan operation.
   * @param end_key a copied end key of this scan operation.
   * @param guard a guard instance to protect nodes from GC.
   */
  OptimisticRecordIterator(  //
      Node *node,
      const ScanKey &begin_key,
      ScanKey end_key,
      EpochGuard &&guard)
      : node_{node},
        end_key_{std::move(end_key)},
        has_key_{begin_key.has_value()},
        guard_{std::move(guard)}
  {
    if (has_key_) {
      // copy the given key
      const auto &[begine_k, len, closed] = begin_key.value();
      if constexpr (IsVarLenData<Key>()) {
        memcpy(key_.get(), begine_k, len);
      } else {
        key_ = begine_k;
      }
      is_closed_ = closed;
    }

    ResetIterator();
  }

  OptimisticRecordIterator(const OptimisticRecordIterator &) = delete;
  OptimisticRecordIterator(OptimisticRecordIterator &&) = delete;

  auto operator=(const OptimisticRecordIterator &) -> OptimisticRecordIterator & = delete;
  auto operator=(OptimisticRecordIterator &&) -> OptimisticRecordIterator & = delete;

  /*####################################################################################
   * Public destructors
   *##################################################################################*/

  /**
   * @brief Destroy the iterator object.
   *
   */
  ~OptimisticRecordIterator() = default;

  /*####################################################################################
   * Public operators for iterators
   *##################################################################################*/

  /**
   * @retval true if this iterator indicates a live record.
   * @retval false otherwise.
   */
  explicit operator bool() { return HasRecord(); }

  /**
   * @retval 1st: a key indicated by the iterator.
   * @retval 2nd: a payload indicated by the iterator.
   */
  auto
  operator*() const  //
      -> std::pair<Key, Payload>
  {
    return {GetKey(), GetPayload()};
  }

  /**
   * @brief Forward the iterator.
   *
   */
  constexpr void
  operator++()
  {
    ++pos_;
  }

  /*####################################################################################
   * Public getters
   *##################################################################################*/

  /**
   * @brief Check if there are any records left.
   *
   * @retval true if there are still records to be scanned.
   * @retval false otherwise.
   */
  [[nodiscard]] auto
  HasRecord()  //
      -> bool
  {
    thread_local auto tmp_key = ReserveDataRegion<Key>();
    [[maybe_unused]] size_t len{};

    while (node_ != nullptr) {
      // check this node has unread records
      if (pos_ < end_pos_) {
        // copy the current key to a temporary region
        if constexpr (IsVarLenData<Key>()) {
          len = node_->GetKeySize(pos_);
          memcpy(tmp_key.get(), node_->GetKey(pos_), len);
        } else {
          tmp_key = node_->GetKey(pos_);
        }

        // check the current record is active
        if (node_->RecordIsDeleted(pos_)) {
          if (!node_->HasSameVersion(ver_)) {
            ResetIterator();  // version check failed, so retry
            continue;
          }

          // the current record was deleted, so skip
          ++pos_;
          if constexpr (IsVarLenData<Key>()) {
            memcpy(key_.get(), tmp_key.get(), len);
          } else {
            key_ = std::move(tmp_key);
          }
          has_key_ = true;
          is_closed_ = false;
          continue;
        }

        // copy the current payload and check the version value
        payload_ = node_->template GetPayload<Payload>(pos_);
        if (!node_->HasSameVersion(ver_)) {
          ResetIterator();  // version check failed, so retry
          continue;
        }

        // copy the current key from the temporary region
        if constexpr (IsVarLenData<Key>()) {
          memcpy(key_.get(), tmp_key.get(), len);
        } else {
          key_ = std::move(tmp_key);
        }
        has_key_ = true;
        is_closed_ = false;

        return true;
      }

      // check this node is rightmost for a given end key
      if (is_end_) {
        node_ = nullptr;
        break;
      }

      // go to the next node
      auto *next = node_->GetNextNode();
      if (!node_->HasSameVersion(ver_)) {
        ResetIterator();  // version check failed, so retry
        continue;
      }

      node_ = next;
      ResetIterator();
    }

    return false;
  }

  /**
   * @return a key indicated by the iterator.
   */
  [[nodiscard]] auto
  GetKey() const  //
      -> Key
  {
    if constexpr (IsVarLenData<Key>()) {
      return key_.get();
    } else {
      return key_;
    }
  }

  /**
   * @return a payload indicated by the iterator.
   */
  [[nodiscard]] auto
  GetPayload() const  //
      -> Payload
  {
    return payload_;
  }

 private:
  /*####################################################################################
   * Internal utilities
   *##################################################################################*/

  /**
   * @tparam T a target class.
   * @return The allocated region for target data.
   */
  template <class T>
  static auto
  ReserveDataRegion()
  {
    if constexpr (IsVarLenData<T>()) {
      using WOPtr = std::remove_pointer_t<T>;
      auto *page = new (::operator new(kMaxVarLenDataSize)) WOPtr{};
      return std::unique_ptr<WOPtr>{page};
    } else {
      return T{};
    }
  }

  /**
   * @brief Initialize the iterator by using the current search key.
   *
   */
  void
  ResetIterator()
  {
    if (has_key_) {
      const auto &cur_key = GetKey();
      ver_ = Node::CheckKeyRange(node_, cur_key);
      const auto [rc, pos] = node_->SearchRecord(cur_key);
      pos_ = (rc == NodeRC::kKeyAlreadyInserted && !is_closed_) ? pos + 1 : pos;
    } else {
      ver_ = node_->GetVersion();
      pos_ = 0;
    }

    std::tie(is_end_, end_pos_) = node_->SearchEndPositionFor(end_key_);
  }

  /*####################################################################################
   * Internal member variables
   *##################################################################################*/

  /// a node that includes partial scan results.
  Node *node_{nullptr};

  /// a current version value for lock-free scanning.
  uint64_t ver_{};

  /// the position of a current record.
  size_t pos_{0};

  /// the end position of records in the node.
  size_t end_pos_{0};

  /// the end key given from a user.
  ScanKey end_key_{};

  /// a flag for indicating a current node is rightmost in scan-range.
  bool is_end_{false};

  /// a flag for indicating this iterator has a current search key.
  bool has_key_{true};

  /// a flag for indicating the current search key is closed-interval.
  bool is_closed_{false};

  /// the guard instance for protecting nodes from GC.
  EpochGuard guard_{};

  /// the payload of the current record.
  Payload payload_{};

  /// the key of the current record.
  static thread_local inline auto key_ = ReserveDataRegion<Key>();  // NOLINT
};

}  // namespace dbgroup::index::b_tree::component

#endif  // B_TREE_COMPONENT_OPTIMISTIC_RECORD_ITERATOR_HPP
