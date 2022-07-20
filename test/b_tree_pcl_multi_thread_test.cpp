#include "b_tree/b_tree_pcl.hpp"
#include "external/index-fixtures/index_fixture_multi_thread.hpp"

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
 * Preparation for typed testing
 *####################################################################################*/

template <class K, class V, class C>
using BTreePCL = ::dbgroup::index::b_tree::BTreePCL<K, V, C>;

using TestTargets = ::testing::Types<  //
    IndexInfo<BTreePCL, UInt8, UInt8>  // fixed-length keys
    >;
TYPED_TEST_SUITE(IndexMultiThreadFixture, TestTargets);

/*######################################################################################
 * Unit test definitions
 *####################################################################################*/

#include "external/index-fixtures/index_fixture_multi_thread_test_definitions.hpp"

}  // namespace dbgroup::index::test
