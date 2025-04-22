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
#include <functional>

// external libraries
#include "dbgroup/lock/mcs_lock.hpp"
#include "dbgroup/lock/optimistic_lock.hpp"
#include "dbgroup/lock/pessimistic_lock.hpp"
#include "dbgroup/memory/epoch_based_gc.hpp"

// local sources
#include "dbgroup/b_tree/component/b_tree_sl.hpp"
#include "dbgroup/b_tree/utility.hpp"

namespace dbgroup::index::b_tree
{
/**
 * @brief A B<sup>link</sup>-tree with optimistic concurrency controls.
 *
 * @tparam Key A key type.
 * @tparam Payload A payload type.
 * @tparam Comp A comparator class for keys.
 * @tparam GC A garbage collector for internal pages.
 */
template <class Key,
          class Payload,
          class Comp = std::less<Key>,
          class GC = ::dbgroup::memory::EpochBasedGC<Page>>
using BTreeOSL = component::BTreeSL<Key,
                                    Payload,
                                    Comp,
                                    GC,
                                    ::dbgroup::lock::OptimisticLock,
                                    ::dbgroup::lock::OptimisticLock,
                                    kUseOptimisticCC>;

/**
 * @brief A B<sup>link</sup>-tree with OCC (inner nodes) and pessimistic locks
 * (leaf nodes).
 *
 * @tparam Key A key type.
 * @tparam Payload A payload type.
 * @tparam Comp A comparator class for keys.
 * @tparam GC A garbage collector for internal pages.
 */
template <class Key,
          class Payload,
          class Comp = std::less<Key>,
          class GC = ::dbgroup::memory::EpochBasedGC<Page>>
using BTreePSL = component::BTreeSL<Key,
                                    Payload,
                                    Comp,
                                    GC,
                                    ::dbgroup::lock::OptimisticLock,
                                    ::dbgroup::lock::PessimisticLock,
                                    !kUseOptimisticCC>;

}  // namespace dbgroup::index::b_tree

#endif  // B_TREE_DBGROUP_B_TREE_B_TREE_HPP_
