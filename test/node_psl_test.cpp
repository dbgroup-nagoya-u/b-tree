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

#include "b_tree/component/psl/node_fixlen.hpp"
#include "b_tree/component/psl/node_varlen.hpp"

// external sources
#include "gtest/gtest.h"

namespace dbgroup::index::b_tree::component::psl::test
{
/*######################################################################################
 * Global constants
 *####################################################################################*/

constexpr uint32_t kLeafFlag = 0;
constexpr bool kExpectSuccess = true;
constexpr bool kExpectFailed = false;

/*######################################################################################
 * Type aliases
 *####################################################################################*/

// extract key-payload types
using Key = uint64_t;
using Payload = uint64_t;
using KeyComp = std::less<Key>;
using PayloadComp = std::less<Payload>;
using NodeVarLen_t = NodeVarLen<Key, KeyComp>;
using NodeFixLen_t = NodeFixLen<Key, KeyComp>;

/*######################################################################################
 * Fixture definitions
 *####################################################################################*/

template <class Node>
class NodeFixture : public testing::Test
{
 protected:
  /*####################################################################################
   * Internal constants
   *##################################################################################*/

  static constexpr size_t kHeaderLen = sizeof(Node);
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
    node_ = CreateNode(kLeafFlag);
  }

  void
  TearDown() override
  {
    ::dbgroup::memory::Release<Page>(node_);
  }

  /*####################################################################################
   * Utility functions
   *##################################################################################*/

  auto
  CreateNode(const bool is_leaf)  //
      -> Node *
  {
    auto *node = new (::dbgroup::memory::Allocate<Page>()) Node{is_leaf};
    if constexpr (std::is_same_v<Node, NodeFixLen_t>) {
      node->SetPayloadLength(kPayLen);
    }
    return node;
  }

  /*####################################################################################
   * Operation wrappers
   *##################################################################################*/

  auto
  Write(  //
      const Key key,
      const Payload payload)
  {
    node_->LockSIX();
    node_->Write(key, kKeyLen, &payload, kPayLen);
  }

  auto
  Insert(  //
      const size_t key,
      const size_t payload)
  {
    node_->LockSIX();
    return node_->Insert(key, kKeyLen, &payload, kPayLen);
  }

  auto
  Update(  //
      const size_t key,
      const size_t payload)
  {
    node_->LockSIX();
    return node_->Update(key, &payload, kPayLen);
  }

  auto
  Delete(const size_t key)
  {
    node_->LockSIX();
    auto rc = node_->Delete(key);
    if (rc == kNeedMerge) {
      node_->UnlockSIX();
    }
    return rc;
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
    auto rc = Insert(key, payload);
    if (expect_success) {
      EXPECT_NE(rc, kKeyAlreadyInserted);
    } else {
      EXPECT_EQ(rc, kKeyAlreadyInserted);
    }
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
    auto rc = Delete(key);
    if (expect_success) {
      EXPECT_NE(rc, kKeyNotInserted);
    } else {
      EXPECT_EQ(rc, kKeyNotInserted);
    }
  }

  /*####################################################################################
   * Internal test definitions
   *##################################################################################*/

  /*------------------------------------------------------------------------------------
   * Leaf read/write APIs
   *----------------------------------------------------------------------------------*/

  void
  TestWrite()
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

  void
  TestInsert()
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

  void
  TestUpdate()
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

  void
  TestDelete()
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

  /*------------------------------------------------------------------------------------
   * Leaf SMO APIs
   *----------------------------------------------------------------------------------*/

  void
  TestSplit()
  {
    // fill the node
    for (size_t i = 0; i < kRecNumInNode; ++i) {
      Write(i, i);
    }

    // perform splitting
    auto *r_node = CreateNode(kLeafFlag);
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
        ::dbgroup::memory::Release<Page>(node_);
        node_ = r_node;
      }
      VerifyRead(i, i, kExpectSuccess);
    }
  }

  void
  TestMerge()
  {
    constexpr auto kHalfNum = kRecNumInNode / 2;

    // fill a right node
    for (size_t i = kHalfNum; i < kRecNumInNode; ++i) {
      Write(i, i);
    }
    auto *r_node = node_;

    // fill a left node
    node_ = CreateNode(kLeafFlag);
    for (size_t i = 0; i < kHalfNum; ++i) {
      Write(i, i);
    }

    // merge the two nodes
    node_->LockSIX();
    r_node->LockSIX();
    node_->Merge(r_node);
    ::dbgroup::memory::Release<Page>(r_node);

    // check the merged node has the written records
    for (size_t i = 0; i < kRecNumInNode; ++i) {
      VerifyRead(i, i, kExpectSuccess);
    }
  }

  /*####################################################################################
   * Internal member variables
   *##################################################################################*/

  Node *node_{nullptr};
};

/*######################################################################################
 * Preparation for typed testing
 *####################################################################################*/

using TestTargets = ::testing::Types<NodeVarLen_t, NodeFixLen_t>;
TYPED_TEST_SUITE(NodeFixture, TestTargets);

/*######################################################################################
 * Test definitions
 *####################################################################################*/

TYPED_TEST(NodeFixture, WritePerformUpsertOperation) { TestFixture::TestWrite(); }

TYPED_TEST(NodeFixture, InsertFailedWithDuplicateKeys) { TestFixture::TestInsert(); }

TYPED_TEST(NodeFixture, UpdateFailedWithNotInsertedKeys) { TestFixture::TestUpdate(); }

TYPED_TEST(NodeFixture, DeleteFailedWithNotInsertedKeys) { TestFixture::TestDelete(); }

TYPED_TEST(NodeFixture, SplitDivideWrittenRecordsIntoTwoNodes) { TestFixture::TestSplit(); }

TYPED_TEST(NodeFixture, MergeTwoNodesIntoSingleNode) { TestFixture::TestMerge(); }

}  // namespace dbgroup::index::b_tree::component::psl::test
