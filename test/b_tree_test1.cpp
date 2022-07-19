#include "b_tree/b_tree_pcl.hpp"
#include "external/index-fixtures/index_fixture.hpp"

namespace dbgroup::index::b_tree
{
/**
 * @brief Use CString as variable-length data in tests.
 *
 */
template <>
constexpr auto
IsVariableLengthData<char *>()  //
    -> bool
{
  return true;
}

}  // namespace dbgroup::index::b_tree

namespace dbgroup::index::test
{
/*######################################################################################
 * API implementation status
 *####################################################################################*/

/*######################################################################################
 * Preparation for typed testing
 *####################################################################################*/
template <class K, class V, class C>
using BTreePCL = ::dbgroup::index::b_tree::BTreePCL<K, V, C>;

using TestTargets = ::testing::Types<   //
    IndexInfo<BTreePCL, UInt8, UInt8>,  // fixed-length keys
    IndexInfo<BTreePCL, UInt8, Int8>,   // fixed-length keys with append-mode
    IndexInfo<BTreePCL, UInt4, UInt8>,  // small keys
    IndexInfo<BTreePCL, UInt4, Int8>,   // small keys with append-mode
    IndexInfo<BTreePCL, UInt8, UInt4>   // small payloads with append-mode
    >;
TYPED_TEST_SUITE(IndexFixture, TestTargets);

/*######################################################################################
 * Unit test definitions
 *####################################################################################*/

#include "external/index-fixtures/index_fixture_test_definitions.hpp"
}  // namespace dbgroup::index::test
