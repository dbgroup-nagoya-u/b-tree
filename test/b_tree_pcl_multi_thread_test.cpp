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

#include "b_tree/b_tree_pcl.hpp"

// organization libraries
#include "external/index-fixtures/index_fixture_multi_thread.hpp"

namespace dbgroup::index::b_tree
{
/**
 * @brief Use CString as variable-length data in tests.
 *
 */
template <>
constexpr auto
IsVarLenData<char *>()  //
    -> bool
{
  return true;
}

}  // namespace dbgroup::index::b_tree

namespace dbgroup::index::test
{
/*######################################################################################
 * Preparation for typed testing
 *####################################################################################*/

template <class K, class V, class C>
using BTreePCL = ::dbgroup::index::b_tree::BTreePCL<K, V, C>;

using TestTargets = ::testing::Types<        //
    IndexInfo<BTreePCL, UInt8, UInt8>,       // fixed-length keys
    IndexInfo<BTreePCL, UInt4, UInt8>,       // small keys
    IndexInfo<BTreePCL, UInt8, UInt4>,       // small payloads
    IndexInfo<BTreePCL, UInt4, UInt4>,       // small keys/payloads
    IndexInfo<BTreePCL, Var, UInt8>,         // variable-length keys
    IndexInfo<BTreePCL, Ptr, Ptr>,           // pointer keys/payloads
    IndexInfo<BTreePCL, Original, Original>  // original class keys/payloads
    >;
TYPED_TEST_SUITE(IndexMultiThreadFixture, TestTargets);

/*######################################################################################
 * Unit test definitions
 *####################################################################################*/

#include "external/index-fixtures/index_fixture_multi_thread_test_definitions.hpp"

}  // namespace dbgroup::index::test
