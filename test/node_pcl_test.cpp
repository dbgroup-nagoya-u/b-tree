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

#include "b_tree/component/pcl/node_varlen.hpp"

// external libraries
#include "gtest/gtest.h"

namespace dbgroup::index::b_tree::component::pcl::test
{
/*######################################################################################
 * Global constants
 *####################################################################################*/

constexpr bool kLeafFlag = true;
constexpr bool kExpectSuccess = true;
constexpr bool kExpectFailed = false;

class NodeFixture : public testing::Test
{
 protected:
  /*####################################################################################
   * Type aliases
   *##################################################################################*/

  // extract key-payload types
  using Key = uint64_t;
  using Payload = uint64_t;
  using KeyComp = std::less<Key>;
  using PayloadComp = std::less<Payload>;

  // define type aliases for simplicity
  using Node_t = NodeVarLen<Key, KeyComp>;

  /*####################################################################################
   * Internal constants
   *##################################################################################*/

  static constexpr size_t kHeaderLen = sizeof(Node_t);
  static constexpr size_t kMetaLen = sizeof(Metadata);
  static constexpr size_t kKeyLen = sizeof(Key);
  static constexpr size_t kPayLen = sizeof(Payload);
  static constexpr size_t kRecLen = kKeyLen + kPayLen + kMetaLen;
  static constexpr size_t kBlockSize = kPageSize - kHeaderLen - kKeyLen;
  static constexpr size_t kRecNumInNode = kBlockSize / kRecLen - ((kBlockSize / kRecLen) % 2);

  /*####################################################################################
   * Setup/Teardown
   *##################################################################################*/

  void
  SetUp() override
  {
    node_ = std::make_unique<Node_t>(kLeafFlag);
  }

  void
  TearDown() override
  {
    node_.reset(nullptr);
  }

  /*####################################################################################
   * Operation wrappers
   *##################################################################################*/

  void
  LockX()
  {
    node_->LockSIX();
    node_->UpgradeToX();
  }

  auto
  Write(  //
      const Key key,
      const Payload payload)
  {
    LockX();
    node_->Write(key, kKeyLen, payload);
  }

  auto
  Insert(  //
      const size_t key,
      const size_t payload)
  {
    LockX();
    return node_->Insert(key, kKeyLen, payload);
  }

  auto
  Update(  //
      const size_t key,
      const size_t payload)
  {
    LockX();
    return node_->Update(key, payload);
  }

  auto
  Delete(const size_t key)
  {
    LockX();
    return node_->Delete(key);
  }

  /*####################################################################################
   * Functions for verification
   *##################################################################################*/

  void
  VerifyRead(  //
      const Key key,
      const Payload expected_val,
      const bool expect_success)
  {
    const auto expected_rc = (expect_success) ? kSuccess : kKeyNotExist;

    Payload payload{};
    node_->LockS();
    const auto rc = node_->Read(key, payload);

    EXPECT_EQ(expected_rc, rc);
    if (expect_success) {
      EXPECT_EQ(expected_val, payload);
    }
  }

  void
  VerifyInsert(  //
      const Key key,
      const Payload payload,
      const bool expect_success)
  {
    const auto expected_rc = (expect_success) ? kSuccess : kKeyExist;
    auto rc = Insert(key, payload);

    EXPECT_EQ(expected_rc, rc);
  }

  void
  VerifyUpdate(  //
      const Key key,
      const Payload payload,
      const bool expect_success)
  {
    const auto expected_rc = (expect_success) ? kSuccess : kKeyNotExist;
    auto rc = Update(key, payload);

    EXPECT_EQ(expected_rc, rc);
  }

  void
  VerifyDelete(  //
      const Key key,
      const bool expect_success)
  {
    const auto expected_rc = (expect_success) ? kSuccess : kKeyNotExist;
    auto rc = Delete(key);

    EXPECT_EQ(expected_rc, rc);
  }

  /*####################################################################################
   * Internal member variables
   *##################################################################################*/

  std::unique_ptr<Node_t> node_{nullptr};
};

/*######################################################################################
 * Test definitions
 *####################################################################################*/

/*--------------------------------------------------------------------------------------
 * Leaf read/write APIs
 *------------------------------------------------------------------------------------*/

TEST_F(NodeFixture, WritePerformUpsertOperation)
{
  for (size_t i = 0; i < kRecNumInNode; ++i) {
    Write(i, i);
  }
  for (size_t i = 0; i < kRecNumInNode; ++i) {
    VerifyRead(i, i, kExpectSuccess);
  }
  for (size_t i = 0; i < kRecNumInNode; ++i) {
    Write(i, i + 1);
  }
  for (size_t i = 0; i < kRecNumInNode; ++i) {
    VerifyRead(i, i + 1, kExpectSuccess);
  }
}

TEST_F(NodeFixture, InsertFailedWithDuplicateKeys)
{
  for (size_t i = 0; i < kRecNumInNode; ++i) {
    VerifyInsert(i, i, kExpectSuccess);
  }
  for (size_t i = 0; i < kRecNumInNode; ++i) {
    VerifyRead(i, i, kExpectSuccess);
  }
  for (size_t i = 0; i < kRecNumInNode; ++i) {
    VerifyInsert(i, i + 1, kExpectFailed);
  }
  for (size_t i = 0; i < kRecNumInNode; ++i) {
    VerifyRead(i, i, kExpectSuccess);
  }
}

TEST_F(NodeFixture, UpdateFailedWithNotInsertedKeys)
{
  for (size_t i = 0; i < kRecNumInNode; ++i) {
    VerifyUpdate(i, i + 1, kExpectFailed);
  }
  for (size_t i = 0; i < kRecNumInNode; ++i) {
    VerifyRead(i, i, kExpectFailed);
  }
  for (size_t i = 0; i < kRecNumInNode; ++i) {
    Write(i, i);
  }
  for (size_t i = 0; i < kRecNumInNode; ++i) {
    VerifyUpdate(i, i + 1, kExpectSuccess);
  }
  for (size_t i = 0; i < kRecNumInNode; ++i) {
    VerifyRead(i, i + 1, kExpectSuccess);
  }
}

TEST_F(NodeFixture, DeleteFailedWithNotInsertedKeys)
{
  for (size_t i = 0; i < kRecNumInNode; ++i) {
    VerifyDelete(i, kExpectFailed);
  }
  for (size_t i = 0; i < kRecNumInNode; ++i) {
    VerifyRead(i, i, kExpectFailed);
  }
  for (size_t i = 0; i < kRecNumInNode; ++i) {
    Write(i, i);
  }
  for (size_t i = 0; i < kRecNumInNode; ++i) {
    VerifyDelete(i, kExpectSuccess);
  }
  for (size_t i = 0; i < kRecNumInNode; ++i) {
    VerifyRead(i, i, kExpectFailed);
  }
}

/*--------------------------------------------------------------------------------------
 * Leaf SMO APIs
 *------------------------------------------------------------------------------------*/

TEST_F(NodeFixture, SplitDivideWrittenRecordsIntoTwoNodes)
{
  // fill the node
  for (size_t i = 0; i < kRecNumInNode; ++i) {
    Write(i, i);
  }

  // perform splitting
  auto *r_node = new Node_t{kLeafFlag};
  node_->LockSIX();
  node_->Split(r_node);
  node_->UnlockSIX();

  // check the split nodes have the same number of records
  const auto l_count = node_->GetRecordCount();
  const auto r_count = r_node->GetRecordCount();
  EXPECT_EQ(l_count, r_count);

  // check the split nodes have the written records
  for (size_t i = 0; i < kRecNumInNode; ++i) {
    if (i == l_count) {
      node_.reset(r_node);
    }
    VerifyRead(i, i, kExpectSuccess);
  }
}

TEST_F(NodeFixture, MergeTwoNodesIntoSingleNode)
{
  constexpr auto kHalfNum = kRecNumInNode / 2;

  // fill a right node
  for (size_t i = kHalfNum; i < kRecNumInNode; ++i) {
    Write(i, i);
  }
  auto *r_node = node_.release();

  // fill a left node
  node_ = std::make_unique<Node_t>(kLeafFlag);
  for (size_t i = 0; i < kHalfNum; ++i) {
    Write(i, i);
  }

  // merge the two nodes
  node_->LockSIX();
  node_->Merge(r_node);
  node_->UnlockSIX();
  delete r_node;

  // check the merged node has the written records
  for (size_t i = 0; i < kRecNumInNode; ++i) {
    VerifyRead(i, i, kExpectSuccess);
  }
}

}  // namespace dbgroup::index::b_tree::component::pcl::test
