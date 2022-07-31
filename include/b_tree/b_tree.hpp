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

#ifndef B_TREE_B_TREE_HPP
#define B_TREE_B_TREE_HPP

// local sources
#include "component/pcl/b_tree.hpp"

namespace dbgroup::index::b_tree
{
/**
 * @brief A class for representing B+trees.
 *
 * @tparam Key a class of stored keys.
 * @tparam Payload a class of stored payloads (only fixed-length data for simplicity).
 * @tparam Comp a class for ordering keys.
 * @tparam kUseOptimisticLock a flag for using optimistic locking.
 * @tparam kUseFineGrainedSMOs a flag for using fine-grained SMOs.
 * @tparam kIsVarLen a flag for indicating variable-length keys.
 */
template <class Key,
          class Payload,
          class Comp = std::less<Key>,
          bool kUseOptimisticLock = true,
          bool kUseFineGrainedSMOs = true,
          bool kIsVarLen = IsVarLenData<Key>()>
class BTree
{
 public:
  /*####################################################################################
   * Type aliases
   *##################################################################################*/

  using BTree_t = component::pcl::BTree<Key, Payload, Comp, kIsVarLen>;
  using RecordIterator_t = component::RecordIterator<BTree_t>;

  /*####################################################################################
   * Public constructors and assignment operators
   *##################################################################################*/

  /**
   * @brief Construct a new BTree object.
   *
   */
  BTree() = default;

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
  ~BTree() = default;

  /*####################################################################################
   * Public read APIs
   *##################################################################################*/

  /**
   * @brief Read the payload corresponding to a given key if it exists.
   *
   * @param key a target key.
   * @retval the payload of a given key wrapped with std::optional if it is in this tree.
   * @retval std::nullopt otherwise.
   */
  auto
  Read(const Key &key)  //
      -> std::optional<Payload>
  {
    return b_tree_.Read(key);
  }

  /**
   * @brief Perform a range scan with given keys.
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
    return b_tree_.Scan(begin_key, end_key);
  }

  /*####################################################################################
   * Public write APIs
   *##################################################################################*/

  /**
   * @brief Write (i.e., put) a given key/payload pair.
   *
   * If a given key does not exist in this tree, this function performs an insert
   * operation. If a given key has been already inserted, this function perfroms an
   * update operation. Thus, this function always returns kSuccess as a return code.
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
    return b_tree_.Write(key, payload, key_len);
  }

  /**
   * @brief Insert a given key/payload pair.
   *
   * This function performs a uniqueness check in its processing. If a given key does
   * not exist in this tree, this function inserts a target payload to this tree. If
   * there is a given key in this tree, this function does nothing and returns kKeyExist
   * as a return code.
   *
   * @param key a target key to be inserted.
   * @param payload a target payload to be inserted.
   * @param key_len the length of a target key.
   * @retval kSuccess if inserted.
   * @retval kKeyExist otherwise.
   */
  auto
  Insert(  //
      const Key &key,
      const Payload &payload,
      const size_t key_len = sizeof(Key))  //
      -> ReturnCode
  {
    return b_tree_.Insert(key, payload, key_len);
  }

  /**
   * @brief Update the record corresponding to a given key with a given payload.
   *
   * This function performs a uniqueness check in its processing. If there is a given
   * key in this tree, this function updates the corresponding record. If a given key
   * does not exist in this tree, this function does nothing and returns kKeyNotExist as
   * a return code.
   *
   * @param key a target key to be updated.
   * @param payload a payload for updating.
   * @param key_len the length of a target key (maybe unused in B+trees).
   * @retval kSuccess if updated.
   * @retval kKeyNotExist otherwise.
   */
  auto
  Update(  //
      const Key &key,
      const Payload &payload,
      const size_t key_len = sizeof(Key))  //
      -> ReturnCode
  {
    return b_tree_.Update(key, payload, key_len);
  }

  /**
   * @brief Delete the record corresponding to a given key from this tree.
   *
   * This function performs a uniqueness check in its processing. If there is a given
   * key in this tree, this function deletes it. If a given key does not exist in this
   * tree, this function does nothing and returns kKeyNotExist as a return code.
   *
   * @param key a target key to be deleted.
   * @param key_len the length of a target key (maybe unused in B+trees).
   * @retval kSuccess if deleted.
   * @retval kKeyNotExist otherwise.
   */
  auto
  Delete(  //
      const Key &key,
      const size_t key_len = sizeof(Key))  //
      -> ReturnCode
  {
    return b_tree_.Delete(key, key_len);
  }

  /*####################################################################################
   * Public bulkload API
   *##################################################################################*/

  /**
   * @brief Bulkload specified kay/payload pairs.
   *
   * This function bulkloads given entries into this index. The entries are assumed to
   * be given as a vector of pairs of Key and Payload (or key/payload/key-length for
   * variable-length keys). Note that keys in records are assumed to be unique and
   * sorted.
   *
   * @param entries vector of entries to be bulkloaded.
   * @param thread_num the number of threads to perform bulkloading.
   * @return kSuccess.
   */
  template <class Entry>
  auto
  Bulkload(  //
      std::vector<Entry> &entries,
      const size_t thread_num = 1)  //
      -> ReturnCode
  {
    assert(thread_num > 0);
    assert(entries.size() >= thread_num);

    return b_tree_.Bulkload(entries, thread_num);
  }

 private:
  /*####################################################################################
   * Internal utility functions
   *##################################################################################*/

  /**
   * @retval true if a target key class is trivially copyable.
   * @retval false otherwise.
   */
  [[nodiscard]] static constexpr auto
  IsValidKeyType()  //
      -> bool
  {
    if constexpr (IsVarLenData<Key>()) {
      if constexpr (!kIsVarLen) {
        // tried to construct B+trees for fixed-length data, but variable-length keys were input
        return false;
      } else {
        return std::is_trivially_copyable_v<std::remove_pointer_t<Key>>;
      }
    } else {
      return std::is_trivially_copyable_v<Key>;
    }
  }

  /*####################################################################################
   * Static assertions
   *##################################################################################*/

  // target keys must be trivially copyable.
  static_assert(IsValidKeyType());

  // target payloads must be trivially copyable.
  static_assert(std::is_trivially_copyable_v<Payload>);

  /*####################################################################################
   * Internal member variables
   *##################################################################################*/

  /// an actual B+tree instance.
  BTree_t b_tree_{};
};

/*######################################################################################
 * Aliases for convenience
 *####################################################################################*/

/// an alias for using a B+tree with pessimistic coarse-grained locking.
template <class Key, class Payload, class Comp = std::less<Key>>
using BTreePCL = BTree<Key, Payload, Comp, false, false>;

/// an alias for using a B+tree with pessimistic coarse-grained locking and variable-length keys.
template <class Key, class Payload, class Comp = std::less<Key>>
using BTreePCLVarLen = BTree<Key, Payload, Comp, false, false, true>;

/// an alias for using a B+tree with pessimistic coarse-grained locking and fixed-length keys.
template <class Key, class Payload, class Comp = std::less<Key>>
using BTreePCLFixLen = BTree<Key, Payload, Comp, false, false, false>;

}  // namespace dbgroup::index::b_tree

#endif  // B_TREE_B_TREE_HPP
