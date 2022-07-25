/*
 * Copyright 2021 Database Group, Nagoya University
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

#include "b_tree/component/node.hpp"

#include "common.hpp"
#include "gtest/gtest.h"

namespace dbgroup::index::b_tree::component::test
{
/*######################################################################################
 * Global constants
 *####################################################################################*/

constexpr size_t kMaxRecSize = 24 + sizeof(Metadata);
constexpr size_t kKeyNumForTest = 1e5;
constexpr bool kLeafFlag = true;
constexpr bool kExpectSuccess = true;
constexpr bool kExpectFailed = false;
constexpr bool kExpectKeyExist = true;
constexpr bool kExpectKeyNotExist = false;

template <class KeyType, class PayloadType>
struct KeyPayload {
  using Key = KeyType;
  using Payload = PayloadType;
};

template <class KeyPayload>
class NodeFixture : public testing::Test
{
  // extract key-payload types
  using Key = typename KeyPayload::Key::Data;
  using Payload = typename KeyPayload::Payload::Data;
  using KeyComp = typename KeyPayload::Key::Comp;
  using PayloadComp = typename KeyPayload::Payload::Comp;

  // define type aliases for simplicity
  using Node_t = PessimisticNode<Key, KeyComp>;

 protected:
  /*################################################################################################
   * Internal constants
   *##############################################################################################*/

  static constexpr size_t kKeyLen = GetDataLength<Key>();
  static constexpr size_t kPayLen = GetDataLength<Payload>();
  static constexpr size_t kRecLen = kKeyLen + kPayLen;
  static constexpr size_t kRecNumInNode =
      (kPageSize - kHeaderLength) / (kRecLen + sizeof(Metadata));

  /*################################################################################################
   * Setup/Teardown
   *##############################################################################################*/

  void
  SetUp() override
  {
    static_assert(kPageSize > kMaxRecSize + kHeaderLength,
                  "The page size is too small to perform unit tests.");

    node_.reset(new Node_t{kLeafFlag});

    PrepareTestData(keys_, kKeyNumForTest);      // NOLINT
    PrepareTestData(payloads_, kKeyNumForTest);  // NOLINT
  }

  void
  TearDown() override
  {
    ReleaseTestData(keys_, kKeyNumForTest);      // NOLINT
    ReleaseTestData(payloads_, kKeyNumForTest);  // NOLINT
  }

  /*####################################################################################
   * Operation wrappers
   *##################################################################################*/

  auto
  Write(  //
      const size_t key_id,
      const size_t payload_id)
  {
    node_->AcquireExclusiveLock();
    const auto rc = node_->Write(keys_[key_id], kKeyLen, payloads_[payload_id]);
    node_->ReleaseExclusiveLock();
    return rc;
  }

  auto
  Insert(  //
      const size_t key_id,
      const size_t payload_id)
  {
    node_->AcquireExclusiveLock();
    const auto rc = node_->Insert(keys_[key_id], kKeyLen, payloads_[payload_id]);
    node_->ReleaseExclusiveLock();
    return rc;
  }

  auto
  Update(  //
      const size_t key_id,
      const size_t payload_id)
  {
    node_->AcquireExclusiveLock();
    const auto rc = node_->Update(keys_[key_id], payloads_[payload_id]);
    node_->ReleaseExclusiveLock();
    return rc;
  }

  auto
  Delete(const size_t key_id)
  {
    node_->AcquireExclusiveLock();
    const auto rc = node_->Delete(keys_[key_id]);
    node_->ReleaseExclusiveLock();
    return rc;
  }

  void
  Consolidate()
  {
    node_->AcquireExclusiveLock();
    node_->template Consolidate<Payload>();
    node_->ReleaseExclusiveLock();
  }

  auto
  Split(Node_t *right_node)
  {
    node_->AcquireExclusiveLock();
    node_->template Split<Payload>(right_node);

    node_->ReleaseExclusiveLock();
    right_node->ReleaseExclusiveLock();
  }

  auto
  Merge(Node_t *right_node)
  {
    node_->AcquireExclusiveLock();
    right_node->AcquireExclusiveLock();

    node_->template Merge<Payload>(right_node);

    node_->ReleaseExclusiveLock();
    right_node->ReleaseExclusiveLock();
  }

  /*####################################################################################
   * Utility functions
   *##################################################################################*/

  void
  PrepareNode()
  {
    const size_t max_num = kRecNumInNode;

    for (size_t i = 0; i < max_num; ++i) {
      Write(i, i);
    }
  }

  /*####################################################################################
   * Functions for verification
   *##################################################################################*/

  void
  VerifyRead(  //
      const size_t key_id,
      const size_t expected_id,
      const bool expect_success)
  {
    const NodeRC expected_rc = (expect_success) ? kCompleted : kKeyNotInserted;

    Payload payload{};
    node_->AcquireSharedLock();
    const auto rc = node_->Read(keys_[key_id], payload);

    EXPECT_EQ(expected_rc, rc);
    if (expect_success) {
      EXPECT_TRUE(IsEqual<PayloadComp>(payloads_[expected_id], payload));
      if constexpr (IsVarLenData<Payload>()) {
        ::operator delete(payload);
      }
    }
  }

  void
  VerifyWrite(  //
      const size_t key_id,
      const size_t payload_id,
      const bool expect_success)
  {
    const NodeRC expected_rc = (expect_success) ? kCompleted : kNeedConsolidation;
    auto rc = Write(key_id, payload_id);

    EXPECT_EQ(expected_rc, rc);
  }

  void
  VerifyInsert(  //
      const size_t key_id,
      const size_t payload_id,
      const bool expect_success,
      const bool expect_key_exist = false)
  {
    NodeRC expected_rc = kCompleted;
    if (!expect_success) {
      expected_rc = (expect_key_exist) ? kKeyAlreadyInserted : kNeedSplit;
    }
    auto rc = Insert(key_id, payload_id);

    EXPECT_EQ(expected_rc, rc);
  }

  void
  VerifyUpdate(  //
      const size_t key_id,
      const size_t payload_id,
      const bool expect_success,
      const bool expect_key_exist = false)
  {
    NodeRC expected_rc = kCompleted;
    if (!expect_success) {
      expected_rc = (expect_key_exist) ? kNeedConsolidation : kKeyNotInserted;
    }
    auto rc = Update(key_id, payload_id);

    EXPECT_EQ(expected_rc, rc);
  }

  void
  VerifyDelete(  //
      const size_t key_id,
      const bool expect_success,
      const bool expect_invoking_merge = false)
  {
    NodeRC expected_rc = kCompleted;
    if (!expect_success) {
      expected_rc = kKeyNotInserted;
    }
    if (expect_invoking_merge) {
      expected_rc = kNeedMerge;
    }

    auto rc = Delete(key_id);

    EXPECT_EQ(expected_rc, rc);
  }

  void
  VerifyConsolidation()
  {
    PrepareNode();

    // node requires consolidation
    VerifyInsert(kRecNumInNode, kRecNumInNode, kExpectFailed);

    for (size_t i = 0; i < kRecNumInNode / 2; ++i) {
      Delete(i);
    }

    // consolidate node
    Consolidate();

    VerifyInsert(kRecNumInNode, kRecNumInNode, kExpectSuccess);

    for (size_t i = 0; i < kRecNumInNode / 2; ++i) {
      VerifyRead(i, i, kExpectFailed);
    }

    for (size_t i = kRecNumInNode / 2; i < kRecNumInNode + 1; ++i) {
      VerifyRead(i, i, kExpectSuccess);
    }
  }

  void
  VerifySplit()
  {
    PrepareNode();

    // node requires consolidation
    VerifyInsert(kRecNumInNode, kRecNumInNode, kExpectFailed, kExpectKeyNotExist);

    auto *right_node = new Node_t{kLeafFlag};
    Split(right_node);

    const auto l_count = (kRecNumInNode - 1) / 2 + 1;  // ceiling

    for (size_t i = 0; i < l_count; ++i) {
      VerifyRead(i, i, kExpectSuccess);
    }

    node_.reset(right_node);
    for (size_t i = l_count; i < kRecNumInNode; ++i) {
      VerifyRead(i, i, kExpectSuccess);
    }
  }

  void
  VerifyMerge()
  {
    PrepareNode();

    auto *right_node = new Node_t{kLeafFlag};
    Split(right_node);

    // to avoid overflow of merged node
    Delete(0);
    right_node->Delete(keys_[kRecNumInNode - 1]);

    Merge(right_node);

    VerifyRead(0, 0, kExpectFailed);
    VerifyRead(kRecNumInNode - 1, kRecNumInNode - 1, kExpectFailed);

    for (size_t i = 1; i < kRecNumInNode - 1; ++i) {
      VerifyRead(i, i, kExpectSuccess);
    }

    delete right_node;
  }

  /*################################################################################################
   * Internal member variables
   *##############################################################################################*/

  // actual keys and payloads
  Key keys_[kKeyNumForTest]{};
  Payload payloads_[kKeyNumForTest]{};

  std::unique_ptr<Node_t> node_{nullptr};
};

/*######################################################################################
 * Preparation for typed testing
 *####################################################################################*/

using KeyPayloadPairs = ::testing::Types<  //
    KeyPayload<UInt8, UInt8>,              // fixed keys and in-place payloads
    KeyPayload<Var, UInt8>,                // variable keys and in-place payloads
    KeyPayload<Ptr, Ptr>,                  // pointer keys/payloads
    KeyPayload<UInt8, Original>,           // original class payloads
    KeyPayload<UInt8, Int8>,               // fixed keys and appended payloads
    KeyPayload<Var, Int8>                  // variable keys and appended payloads
    >;
TYPED_TEST_SUITE(NodeFixture, KeyPayloadPairs);

/*######################################################################################
 * Write operation
 *####################################################################################*/

TYPED_TEST(NodeFixture, WriteWithUniqueKeysReadWrittenValues)
{
  const size_t max_num = TestFixture::kRecNumInNode;

  // write records
  for (size_t i = 0; i < max_num; ++i) {
    TestFixture::VerifyWrite(i, i, kExpectSuccess);
  }

  // the node has the written records
  for (size_t i = 0; i < max_num; ++i) {
    TestFixture::VerifyRead(i, i, kExpectSuccess);
  }
}

TYPED_TEST(NodeFixture, WriteWithDuplicateKeysReadLatestValues)
{
  const size_t max_num = TestFixture::kRecNumInNode / 2;

  // write base records
  for (size_t i = 0; i < max_num; ++i) {
    TestFixture::VerifyWrite(i, i, kExpectSuccess);
  }

  // overwrite the records with different values
  for (size_t i = 0; i < max_num; ++i) {
    TestFixture::VerifyWrite(i, i + 1, kExpectSuccess);
  }

  // the node has the written records
  for (size_t i = 0; i < max_num; ++i) {
    TestFixture::VerifyRead(i, i + 1, kExpectSuccess);
  }
}

/*--------------------------------------------------------------------------------------
 * Insert operation
 *------------------------------------------------------------------------------------*/

TYPED_TEST(NodeFixture, InsertWithUniqueKeysReadWrittenValues)
{
  const size_t max_num = TestFixture::kRecNumInNode;

  // insert records
  for (size_t i = 0; i < max_num; ++i) {
    TestFixture::VerifyInsert(i, i, kExpectSuccess);
  }

  // the node has the inserted records
  for (size_t i = 0; i < max_num; ++i) {
    TestFixture::VerifyRead(i, i, kExpectSuccess);
  }
}

TYPED_TEST(NodeFixture, InsertWithDuplicateKeysFail)
{
  const size_t max_num = TestFixture::kRecNumInNode / 2;

  // insert base records
  for (size_t i = 0; i < max_num; ++i) {
    TestFixture::VerifyInsert(i, i, kExpectSuccess);
  }

  // insert operations will fail with inserted-keys
  for (size_t i = 0; i < max_num; ++i) {
    TestFixture::VerifyInsert(i, i + 1, kExpectFailed, kExpectKeyExist);
  }

  // the node has the inserted records
  for (size_t i = 0; i < max_num; ++i) {
    TestFixture::VerifyRead(i, i, kExpectSuccess);
  }
}

/*--------------------------------------------------------------------------------------
 * Update operation
 *------------------------------------------------------------------------------------*/

TYPED_TEST(NodeFixture, UpdateWithUniqueKeysFail)
{
  const size_t max_num = TestFixture::kRecNumInNode / 2;

  // write base records
  for (size_t i = 0; i < max_num; ++i) {
    TestFixture::Write(i, i);
  }

  // update operations will fail with not-inserted keys
  for (size_t i = max_num; i < 2 * max_num; ++i) {
    TestFixture::VerifyUpdate(i, i + 1, kExpectFailed, kExpectKeyNotExist);
  }

  // the node has the written records
  for (size_t i = 0; i < max_num; ++i) {
    TestFixture::VerifyRead(i, i, kExpectSuccess);
  }

  // the failed update operations do not modify the node
  for (size_t i = max_num; i < 2 * max_num; ++i) {
    TestFixture::VerifyRead(i, i + 1, kExpectFailed);
  }
}

TYPED_TEST(NodeFixture, UpdateWithDuplicateKeysReadLatestValues)
{
  const size_t max_num = TestFixture::kRecNumInNode / 2;

  // write base records
  for (size_t i = 0; i < max_num; ++i) {
    TestFixture::Write(i, i);
  }

  // overwrite the records with different values
  for (size_t i = 0; i < max_num; ++i) {
    TestFixture::VerifyUpdate(i, i + 1, kExpectSuccess);
  }

  // the node has the updated records
  for (size_t i = 0; i < max_num; ++i) {
    TestFixture::VerifyRead(i, i + 1, kExpectSuccess);
  }
}
/*--------------------------------------------------------------------------------------
 * Delete operation
 *------------------------------------------------------------------------------------*/

TYPED_TEST(NodeFixture, DeleteWithUniqueKeysFail)
{
  const size_t max_num = TestFixture::kRecNumInNode;

  // write base records
  for (size_t i = 0; i < max_num; ++i) {
    TestFixture::VerifyWrite(i, i, kExpectSuccess);
  }

  // delete operations will fail with not-inserted keys
  for (size_t i = max_num; i < 2 * max_num; ++i) {
    TestFixture::VerifyDelete(i, kExpectFailed, kExpectKeyNotExist);
  }

  // the node has the written records
  for (size_t i = 0; i < max_num; ++i) {
    TestFixture::VerifyRead(i, i, kExpectSuccess);
  }
}

TYPED_TEST(NodeFixture, DeleteWithDuplicateKeysSuccess)
{
  const size_t max_num = TestFixture::kRecNumInNode;
  const size_t delete_num =
      ((TestFixture::kRecLen + sizeof(Metadata)) * max_num - kMinUsedSpaceSize)
      / (TestFixture::kRecLen + sizeof(Metadata));

  // write base records
  for (size_t i = 0; i < max_num; ++i) {
    TestFixture::VerifyWrite(i, i, kExpectSuccess);
  }

  // delete inserted keys
  for (size_t i = 0; i < delete_num; ++i) {
    TestFixture::VerifyDelete(i, kExpectSuccess);
  }

  // delete and get kNeedMerge
  TestFixture::VerifyDelete(delete_num, kExpectFailed, true);

  // the node does not have deleted records
  for (size_t i = 0; i < delete_num + 1; ++i) {
    TestFixture::VerifyRead(i, i, kExpectFailed);
  }

  // the node does not have deleted records
  for (size_t i = delete_num + 1; i < max_num; ++i) {
    TestFixture::VerifyRead(i, i, kExpectSuccess);
  }
}

/*--------------------------------------------------------------------------------------
 * SMO operation
 *------------------------------------------------------------------------------------*/

TYPED_TEST(NodeFixture, SplitDivideWrittenRecordsIntoTwoNodes)
{  //
  TestFixture::VerifySplit();
}

TYPED_TEST(NodeFixture, MergeCopyWrittenRecordsIntoSingleNode)
{  //
  TestFixture::VerifyMerge();
}

TYPED_TEST(NodeFixture, ConsolidateNode)
{  //
  TestFixture::VerifyConsolidation();
}

}  // namespace dbgroup::index::b_tree::component::test
