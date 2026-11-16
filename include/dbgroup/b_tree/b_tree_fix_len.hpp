/*
 * Copyright 2026 Database Group, Nagoya University
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

#ifndef B_TREE_DBGROUP_B_TREE_B_TREE_FIX_LEN_HPP_
#define B_TREE_DBGROUP_B_TREE_B_TREE_FIX_LEN_HPP_

// C++ standard libraries
#include <functional>

// our libraries
#include <dbgroup/index/fix_len_wrapper.hpp>
#include <dbgroup/lock/omcs_lock.hpp>
#include <dbgroup/lock/utility.hpp>
#include <dbgroup/memory/epoch_based_gc.hpp>

// local sources
#include "dbgroup/b_tree/b_tree.hpp"

namespace dbgroup::index::b_tree
{
/**
 * @brief B+trees for fixed-length payloads with APIs returning `std::optional`.
 *
 * @tparam Key A class of stored keys.
 * @tparam Payload A class of stored payloads.
 * @tparam Comp A comparator class for keys.
 * @tparam GC A garbage collector for reusing node pages.
 * @tparam Lock A class for locking each node.
 * @tparam Delta A class for representing delta updates.
 * @note This class is defined by inheritance instead of an alias template so
 * that it can receive pack expansions as its template arguments.
 */
template <class Key,
          class Payload,
          class Comp = std::less<Key>,
          class GC = ::dbgroup::memory::EpochBasedGC<>,
          ::dbgroup::lock::OptimisticallyLockable Lock = ::dbgroup::lock::OMCSLock,
          class Delta = Payload>
class BTreeFixLen : public FixLenWrapper<BTree<Key, Payload, Comp, GC, Lock, Delta>>
{
  /*##########################################################################*
   * Type aliases
   *##########################################################################*/

  using Wrapper = FixLenWrapper<BTree<Key, Payload, Comp, GC, Lock, Delta>>;

 public:
  /*##########################################################################*
   * Public constructors and assignment operators
   *##########################################################################*/

  using Wrapper::Wrapper;
};

}  // namespace dbgroup::index::b_tree

#endif  // B_TREE_DBGROUP_B_TREE_B_TREE_FIX_LEN_HPP_
