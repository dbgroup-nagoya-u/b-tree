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

#ifndef B_TREE_B_TREE_HPP
#define B_TREE_B_TREE_HPP

// local sources
#include "utility.hpp"

// actual B+tree implementations
#include "b_tree/component/oml/b_tree.hpp"
#include "b_tree/component/osl/b_tree.hpp"
#include "b_tree/component/pml/b_tree.hpp"
#include "b_tree/component/psl/b_tree.hpp"

namespace dbgroup::index::b_tree
{
/**
 * @brief A class for representing B+trees.
 *
 * @tparam Key a class of stored keys.
 * @tparam Payload a class of stored payloads (only fixed-length data for simplicity).
 * @tparam Comp a class for ordering keys.
 * @tparam kUseOptimisticLock a flag for using optimistic locking.
 * @tparam kUseSingleLayerLock a flag for using fine-grained SMOs.
 * @tparam kUseVarLenLayout a flag for indicating variable-length keys.
 */
template <class Key,
          class Payload,
          class Comp = std::less<Key>,
          bool kUseOptimisticLock = true,
          bool kUseSingleLayerLock = true,
          bool kUseVarLenLayout = IsVarLenData<Key>()>
class BTree
{
 public:
  /*####################################################################################
   * Type aliases
   *##################################################################################*/

  template <class K, class V, class C, bool VAR>
  using BTreePML = component::pml::BTree<K, V, C, VAR>;

  template <class K, class V, class C, bool VAR>
  using BTreePSL = component::psl::BTree<K, V, C, VAR>;

  template <class K, class V, class C, bool VAR>
  using BTreeOML = component::oml::BTree<K, V, C, VAR>;

  template <class K, class V, class C, bool VAR>
  using BTreeOSL = component::osl::BTree<K, V, C, VAR>;

  using BTree_t =
      std::conditional_t<kUseOptimisticLock,
                         std::conditional_t<kUseSingleLayerLock,
                                            BTreeOSL<Key, Payload, Comp, kUseVarLenLayout>,
                                            BTreeOML<Key, Payload, Comp, kUseVarLenLayout>>,
                         std::conditional_t<kUseSingleLayerLock,
                                            BTreePSL<Key, Payload, Comp, kUseVarLenLayout>,
                                            BTreePML<Key, Payload, Comp, kUseVarLenLayout>>>;

  using ScanKey = std::optional<std::tuple<const Key &, size_t, bool>>;

  /*####################################################################################
   * Public constructors and assignment operators
   *##################################################################################*/

  /**
   * @brief Construct a new BTree object.
   *
   */
  explicit BTree(  //
      const size_t gc_interval_micro = kDefaultGCTime,
      const size_t gc_thread_num = kDefaultGCThreadNum)
      : b_tree_{gc_interval_micro, gc_thread_num}
  {
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
  ~BTree() = default;

  /*####################################################################################
   * Public read APIs
   *##################################################################################*/

  /**
   * @brief Read the payload corresponding to a given key if it exists.
   *
   * @param key a target key.
   * @param key_len the length of the target key.
   * @retval the payload of a given key wrapped with std::optional if it is in this tree.
   * @retval std::nullopt otherwise.
   */
  auto
  Read(  //
      const Key &key,
      const size_t key_len = sizeof(Key))  //
      -> std::optional<Payload>
  {
    return b_tree_.Read(key, key_len);
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
      const ScanKey &begin_key = std::nullopt,
      const ScanKey &end_key = std::nullopt)
  {
    return b_tree_.Scan(begin_key, end_key);
  }

  /*####################################################################################
   * Public write APIs
   *##################################################################################*/

  /**
   * @brief Write (i.e., put) a given key/payload pair.
   *
   * This function always overwrites a payload and can be optimized for that purpose;
   * the procedure may omit the key uniqueness check.
   *
   * @param key a target key.
   * @param payload a target payload.
   * @param key_len the length of the target key.
   * @param pay_len the length of the target payload.
   * @return kSuccess.
   */
  auto
  Write(  //
      const Key &key,
      const Payload &payload,
      const size_t key_len = sizeof(Key),
      [[maybe_unused]] const size_t pay_len = sizeof(Payload))  //
      -> ReturnCode
  {
    return b_tree_.Write(key, payload, key_len);
  }

  /**
   * @brief Insert a given key/payload pair.
   *
   * This function performs a uniqueness check on its processing. If the given key does
   * not exist in this tree, this function inserts a target payload into this tree. If
   * the given key exists in this tree, this function does nothing and returns kKeyExist.
   *
   * @param key a target key.
   * @param payload a target payload.
   * @param key_len the length of the target key.
   * @param pay_len the length of the target payload.
   * @retval kSuccess if inserted.
   * @retval kKeyExist otherwise.
   */
  auto
  Insert(  //
      const Key &key,
      const Payload &payload,
      const size_t key_len = sizeof(Key),
      [[maybe_unused]] const size_t pay_len = sizeof(Payload))  //
      -> ReturnCode
  {
    return b_tree_.Insert(key, payload, key_len);
  }

  /**
   * @brief Update the record corresponding to a given key with a given payload.
   *
   * This function performs a uniqueness check on its processing. If the given key
   * exists in this tree, this function updates the corresponding payload. If the given
   * key does not exist in this tree, this function does nothing and returns
   * kKeyNotExist.
   *
   * @param key a target key.
   * @param payload a target payload.
   * @param key_len the length of the target key.
   * @param pay_len the length of the target payload.
   * @retval kSuccess if updated.
   * @retval kKeyNotExist otherwise.
   */
  auto
  Update(  //
      const Key &key,
      const Payload &payload,
      const size_t key_len = sizeof(Key),
      [[maybe_unused]] const size_t pay_len = sizeof(Payload))  //
      -> ReturnCode
  {
    return b_tree_.Update(key, payload, key_len);
  }

  /**
   * @brief Delete the record corresponding to a given key from this tree.
   *
   * This function performs a uniqueness check on its processing. If the given key
   * exists in this tree, this function deletes it. If the given key does not exist in
   * this tree, this function does nothing and returns kKeyNotExist.
   *
   * @param key a target key.
   * @param key_len the length of the target key.
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
   * This function loads the given entries into this index, assuming that the entries
   * are given as a vector of key/payload pairs (or the tuples key/payload/key-length
   * for variable-length keys). Note that keys in records are assumed to be unique and
   * sorted.
   *
   * @tparam Entry a container of a key/payload pair.
   * @param entries the vector of entries to be bulkloaded.
   * @param thread_num the number of threads used for bulk loading.
   * @return kSuccess.
   */
  template <class Entry>
  auto
  Bulkload(  //
      const std::vector<Entry> &entries,
      const size_t thread_num = 1)  //
      -> ReturnCode
  {
    assert(thread_num > 0);
    assert(entries.size() >= thread_num);

    return b_tree_.Bulkload(entries, thread_num);
  }

  /*####################################################################################
   * Public utilities
   *##################################################################################*/

  /**
   * @brief Collect statistical data of this tree.
   *
   * @retval 1st: the number of nodes.
   * @retval 2nd: the actual usage in bytes.
   * @retval 3rd: the virtual usage (i.e., reserved memory) in bytes.
   */
  auto
  CollectStatisticalData()  //
      -> std::vector<std::tuple<size_t, size_t, size_t>>
  {
    return b_tree_.CollectStatisticalData();
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
    if constexpr (kUseVarLenLayout && IsVarLenData<Key>()) {
      // check a base type is trivially copyable
      return std::is_trivially_copyable_v<std::remove_pointer_t<Key>>;
    } else if constexpr (IsVarLenData<Key>()) {
      // cannot use optimized page layouts with variable-length data
      return false;
    } else {
      // check a given key type is trivially copyable
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

/// a B+tree based on pessimistic multi-layer locking (PML).
template <class Key, class Payload, class Comp = std::less<Key>>
using BTreePML = BTree<Key, Payload, Comp, kPessimisticLock, kMultiLayerLock>;

/// a B+tree based on PML with generic page layouts.
template <class Key, class Payload, class Comp = std::less<Key>>
using BTreePMLVarLen =
    BTree<Key, Payload, Comp, kPessimisticLock, kMultiLayerLock, !kOptimizeForFixLenData>;

/// a B+tree based on PML with optimized page layouts for fixed-length data.
template <class Key, class Payload, class Comp = std::less<Key>>
using BTreePMLFixLen =
    BTree<Key, Payload, Comp, kPessimisticLock, kMultiLayerLock, kOptimizeForFixLenData>;

/// a B+tree with pessimistic single-layer locking (PSL).
template <class Key, class Payload, class Comp = std::less<Key>>
using BTreePSL = BTree<Key, Payload, Comp, kPessimisticLock, kSingleLayerLock>;

/// a B+tree based on PSL with generic page layouts.
template <class Key, class Payload, class Comp = std::less<Key>>
using BTreePSLVarLen =
    BTree<Key, Payload, Comp, kPessimisticLock, kSingleLayerLock, !kOptimizeForFixLenData>;

/// a B+tree based on PSL with generic page layouts.
template <class Key, class Payload, class Comp = std::less<Key>>
using BTreePSLFixLen =
    BTree<Key, Payload, Comp, kPessimisticLock, kSingleLayerLock, kOptimizeForFixLenData>;

/// a B+tree based on optimistic multi-layer locking (OML).
template <class Key, class Payload, class Comp = std::less<Key>>
using BTreeOML = BTree<Key, Payload, Comp, kOptimisticLock, kMultiLayerLock>;

/// a B+tree based on OML with generic page layouts.
template <class Key, class Payload, class Comp = std::less<Key>>
using BTreeOMLVarLen =
    BTree<Key, Payload, Comp, kOptimisticLock, kMultiLayerLock, !kOptimizeForFixLenData>;

/// a B+tree based on OML with optimized page layouts for fixed-length data.
template <class Key, class Payload, class Comp = std::less<Key>>
using BTreeOMLFixLen =
    BTree<Key, Payload, Comp, kOptimisticLock, kMultiLayerLock, kOptimizeForFixLenData>;

/// a B+tree with optimistic single-layer locking (OSL).
template <class Key, class Payload, class Comp = std::less<Key>>
using BTreeOSL = BTree<Key, Payload, Comp, kOptimisticLock, kSingleLayerLock>;

/// a B+tree based on OSL with generic page layouts.
template <class Key, class Payload, class Comp = std::less<Key>>
using BTreeOSLVarLen =
    BTree<Key, Payload, Comp, kOptimisticLock, kSingleLayerLock, !kOptimizeForFixLenData>;

/// a B+tree based on OSL with generic page layouts.
template <class Key, class Payload, class Comp = std::less<Key>>
using BTreeOSLFixLen =
    BTree<Key, Payload, Comp, kOptimisticLock, kSingleLayerLock, kOptimizeForFixLenData>;

}  // namespace dbgroup::index::b_tree

#endif  // B_TREE_B_TREE_HPP
