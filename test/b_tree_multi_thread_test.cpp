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

// the corresponding header
#include "dbgroup/b_tree/b_tree.hpp"

// external sources
#include <dbgroup/index_fixtures/index_fixture_multi_thread.hpp>

namespace dbgroup::index::test
{
/*############################################################################*
 * Preparation for typed testing
 *############################################################################*/

template <class K, class V, class C>
using Index = ::dbgroup::index::b_tree::BTree<K, V, C>;

using TestTargets = ::testing::Types<  //
    IndexInfo<Index, UInt8, UInt8>,    // fixed-length keys
    IndexInfo<Index, UInt4, UInt8>,    // small keys
    IndexInfo<Index, Var, UInt8>,      // varlen keys
    IndexInfo<Index, UInt8, UInt4>,    // fixed-length keys/small payloads
    IndexInfo<Index, UInt4, UInt4>,    // small keys/small payloads
    IndexInfo<Index, Var, UInt4>       // varlen keys/small payloads
    >;
TYPED_TEST_SUITE(IndexMultiThreadFixture, TestTargets);

/*############################################################################*
 * Unit test definitions
 *############################################################################*/

#include <dbgroup/index_fixtures/index_fixture_multi_thread_test_definitions.hpp>

}  // namespace dbgroup::index::test
