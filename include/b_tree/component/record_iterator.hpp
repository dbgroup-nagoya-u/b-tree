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

#ifndef B_TREE_COMPONENT_RECORD_ITERATOR_HPP
#define B_TREE_COMPONENT_RECORD_ITERATOR_HPP

#include <optional>
#include <utility>

// local sources
#include "common.hpp"

namespace dbgroup::index::b_tree::component
{
/**
 * @brief A class for representing an iterator of scan results.
 *
 */
template <class Index>
class RecordIterator
{
 public:
  /*####################################################################################
   * Type aliases
   *##################################################################################*/

  using Key = typename Index::K;
  using Payload = typename Index::V;
  using Node = typename Index::Node_t;

  /*####################################################################################
   * Public constructors and assignment operators
   *##################################################################################*/

  /**
   * @brief Construct a new iterator object.
   *
   * @param node a node to be scanned.
   * @param begin_pos the begin position for scanning in the node.
   * @param end_pos the end position for scanning in the node.
   * @param end_key a copied end key of this scan operation.
   * @param is_end a flag for indicating the node is rightmost in this scan operation.
   */
  RecordIterator(  //
      Node *node,
      size_t begin_pos,
      size_t end_pos,
      std::optional<std::pair<const Key &, bool>> end_key,
      bool is_end)
      : node_{node},
        pos_{begin_pos},
        end_pos_{end_pos},
        end_key_{std::move(end_key)},
        is_end_{is_end}
  {
  }

  RecordIterator(const RecordIterator &) = delete;
  RecordIterator(RecordIterator &&) = delete;

  auto operator=(const RecordIterator &) -> RecordIterator & = delete;
  auto operator=(RecordIterator &&) -> RecordIterator & = delete;

  /*####################################################################################
   * Public destructors
   *##################################################################################*/

  /**
   * @brief Destroy the iterator object.
   *
   */
  ~RecordIterator() = default;

  /*####################################################################################
   * Public operators for iterators
   *##################################################################################*/

  /**
   * @retval 1st: a key indicated by the iterator.
   * @retval 2nd: a payload indicated by the iterator.
   */
  auto
  operator*() const  //
      -> std::pair<Key, Payload>
  {
    return node_->template GetRecord<Payload>(pos_);
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
  HasNext()  //
      -> bool
  {
    while (true) {
      // check records remain in this node
      while (pos_ < end_pos_ && node_->RecordIsDeleted(pos_)) {
        ++pos_;  // skip deleted records
      }
      if (pos_ < end_pos_) return true;

      // check this node is rightmost for a given end key
      if (is_end_) {
        node_->UnlockS();
        return false;
      }

      // go to the next node
      node_ = node_->GetNextNodeForRead();
      pos_ = 0;
      std::tie(is_end_, end_pos_) = node_->SearchEndPositionFor(end_key_);
    }
  }

  /**
   * @return a key indicated by the iterator.
   */
  [[nodiscard]] auto
  GetKey() const  //
      -> Key
  {
    return node_->GetKey(pos_);
  }

  /**
   * @return a payload indicated by the iterator.
   */
  [[nodiscard]] auto
  GetPayload() const  //
      -> Payload
  {
    return node_->template GetPayload<Payload>(pos_);
  }

 private:
  /*####################################################################################
   * Internal member variables
   *##################################################################################*/

  /// a node that includes partial scan results.
  Node *node_{nullptr};

  /// the position of a current record.
  size_t pos_{0};

  /// the end position of records in the node.
  size_t end_pos_{0};

  /// the end key given from a user.
  std::optional<std::pair<const Key &, bool>> end_key_{};

  /// a flag for indicating a current node is rightmost in scan-range.
  bool is_end_{false};
};

}  // namespace dbgroup::index::b_tree::component

#endif  // B_TREE_COMPONENT_RECORD_ITERATOR_HPP
