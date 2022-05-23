#include <random>

#include "b_tree/b_tree.hpp"
#include "common.hpp"
#include "gtest/gtest.h"

namespace dbgroup::index::b_tree::component
{
/*######################################################################################
 * Global constants
 *####################################################################################*/

constexpr bool kLeafFlag = true;
constexpr bool kExpectSuccess = true;
constexpr bool kExpectFailed = false;
constexpr bool kExpectKeyExist = true;
constexpr bool kExpectKeyNotExist = false;
constexpr bool kWriteTwice = true;
constexpr bool kWithWrite = true;
constexpr bool kWithDelete = true;
constexpr bool kShuffled = true;
constexpr bool kRangeClosed = true;

template <class KeyType, class PayloadType>
struct KeyPayload {
  using Key = KeyType;
  using Payload = PayloadType;
};

template <class KeyPayload>
class BTreeFixture : public testing::Test  // NOLINT
{
  // extract key-payload types
  using Key = typename KeyPayload::Key::Data;
  using Payload = typename KeyPayload::Payload::Data;
  using KeyComp = typename KeyPayload::Key::Comp;
  using PayloadComp = typename KeyPayload::Payload::Comp;

  // define type aliases for simplicity
  using Metadata = component::Metadata;
  using Node_t = component::Node<Key, KeyComp>;
  using NodeRC = component::NodeRC;
  using NodeStack = std::vector<std::pair<Node_t *, size_t>>;
  using BTree_t = BTree<Key, Payload, KeyComp>;

 protected:
  /*################################################################################################
   * Internal constants
   *##############################################################################################*/

  static constexpr size_t kKeyLen = GetDataLength<Key>();
  static constexpr size_t kPayLen = GetDataLength<Payload>();
  static constexpr size_t kRecLen = kKeyLen + kPayLen;
  static constexpr size_t kRecNumInNode =
      (kPageSize - kHeaderLength) / (kRecLen + sizeof(Metadata));
  static constexpr size_t kMaxRecNumForTest = 2 * kRecNumInNode * kRecNumInNode;
  static constexpr size_t kKeyNumForTest = 2 * kRecNumInNode * kRecNumInNode + 2;

  /*################################################################################################
   * Setup/Teardown
   *##############################################################################################*/

  void
  SetUp() override
  {
    PrepareTestData(keys_, kKeyNumForTest);      // NOLINT
    PrepareTestData(payloads_, kKeyNumForTest);  // NOLINT

    b_tree_ = std::make_unique<BTree_t>();
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
  Insert(  //
      const size_t key_id,
      const size_t payload_id)  //
      -> ReturnCode
  {
    return b_tree_->Insert(keys_[key_id], payloads_[payload_id], kKeyLen);
  }

  /*####################################################################################
   * Utility functions
   *##################################################################################*/

  [[nodiscard]] auto
  CreateTargetIDs(  //
      const size_t rec_num,
      const bool is_shuffled) const  //
      -> std::vector<size_t>
  {
    std::mt19937_64 rand_engine{kRandomSeed};

    std::vector<size_t> target_ids{};
    target_ids.reserve(rec_num);
    for (size_t i = 0; i < rec_num; ++i) {
      target_ids.emplace_back(i);
    }

    if (is_shuffled) {
      std::shuffle(target_ids.begin(), target_ids.end(), rand_engine);
    }

    return target_ids;
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
    const auto read_val = b_tree_->Read(keys_[key_id]);
    if (expect_success) {
      EXPECT_TRUE(read_val);

      const auto expected_val = payloads_[expected_id];
      const auto actual_val = read_val.value();
      EXPECT_TRUE(component::IsEqual<PayloadComp>(expected_val, actual_val));
    } else {
      EXPECT_FALSE(read_val);
    }
  }

  void
  VerifyScan(  //
      const std::optional<std::pair<size_t, bool>> begin_ref,
      const std::optional<std::pair<size_t, bool>> end_ref)
  {
    std::optional<std::pair<const Key &, bool>> begin_key = std::nullopt;
    size_t begin_pos = 0;
    if (begin_ref) {
      auto &&[begin_id, begin_closed] = *begin_ref;
      begin_key.emplace(keys_[begin_id], begin_closed);
      begin_pos = (begin_closed) ? begin_id : begin_id + 1;
    }

    std::optional<std::pair<const Key &, bool>> end_key = std::nullopt;
    size_t end_pos = 0;
    if (end_ref) {
      auto &&[end_id, end_closed] = *end_ref;
      end_key.emplace(keys_[end_id], end_closed);
      end_pos = (end_closed) ? end_id + 1 : end_id;
    }

    auto iter = b_tree_->Scan(begin_key, end_key);

    for (; iter.HasNext(); ++iter, ++begin_pos) {
      auto [key, payload] = *iter;
      EXPECT_TRUE(IsEqual<KeyComp>(keys_[begin_pos], key));
      EXPECT_TRUE(IsEqual<PayloadComp>(payloads_[begin_pos], payload));
    }
    EXPECT_FALSE(iter.HasNext());

    if (end_ref) {
      EXPECT_EQ(begin_pos, end_pos);
    }
  }

  void
  VerifyWrite(  //
      const size_t key_id,
      const size_t pay_id)
  {
    auto rc = b_tree_->Write(keys_[key_id], payloads_[pay_id], kKeyLen);

    EXPECT_EQ(ReturnCode::kSuccess, rc);
  }

  void
  VerifyInsert(  //
      const size_t key_id,
      const size_t payload_id,
      const bool expect_success)
  {
    ReturnCode expected_rc = b_tree::kSuccess;
    expected_rc = (expect_success) ? expected_rc : b_tree::kKeyExist;
    auto rc = Insert(key_id, payload_id);
    EXPECT_EQ(expected_rc, rc);
  }

  void
  VerifyUpdate(  //
      const size_t key_id,
      const size_t payload_id,
      const bool expect_success)
  {
    ReturnCode expected_rc = (expect_success) ? kSuccess : kKeyNotExist;

    auto rc = b_tree_->Update(keys_[key_id], payloads_[payload_id]);
    EXPECT_EQ(expected_rc, rc);
  }

  void
  VerifyDelete(  //
      const size_t key_id,
      const bool expect_success)
  {
    ReturnCode expected_rc = (expect_success) ? kSuccess : kKeyNotExist;

    auto rc = b_tree_->Delete(keys_[key_id]);
    EXPECT_EQ(expected_rc, rc);
  }

  void
  VerifyWritesWith(  //
      const bool write_twice,
      const bool is_shuffled,
      const size_t ops_num = kMaxRecNumForTest)
  {
    const auto &target_ids = CreateTargetIDs(ops_num, is_shuffled);

    for (size_t i = 0; i < ops_num; ++i) {
      const auto id = target_ids.at(i);
      VerifyWrite(id, id);
    }
    if (write_twice) {
      for (size_t i = 0; i < ops_num; ++i) {
        const auto id = target_ids.at(i);
        VerifyWrite(id, id + 1);
      }
    }
    for (size_t i = 0; i < ops_num; ++i) {
      const auto key_id = target_ids.at(i);
      const auto val_id = (write_twice) ? key_id + 1 : key_id;
      VerifyRead(key_id, val_id, kExpectSuccess);
    }
  }

  void
  VerifyInsertsWith(  //
      const bool write_twice,
      const bool with_delete,
      const bool is_shuffled)
  {
    const auto &target_ids = CreateTargetIDs(kMaxRecNumForTest, is_shuffled);

    for (size_t i = 0; i < kMaxRecNumForTest; ++i) {
      const auto id = target_ids.at(i);
      VerifyInsert(id, id, kExpectSuccess);
    }
    if (with_delete) {
      for (size_t i = 0; i < kMaxRecNumForTest; ++i) {
        const auto id = target_ids.at(i);
        VerifyDelete(id, kExpectSuccess);
      }
    }
    if (write_twice) {
      for (size_t i = 0; i < kMaxRecNumForTest; ++i) {
        const auto key_id = target_ids.at(i);
        VerifyInsert(key_id, key_id + 1, with_delete);
      }
    }
    for (size_t i = 0; i < kMaxRecNumForTest; ++i) {
      const auto key_id = target_ids.at(i);
      const auto val_id = (write_twice && with_delete) ? key_id + 1 : key_id;
      VerifyRead(key_id, val_id, kExpectSuccess);
    }
  }

  void
  VerifyUpdatesWith(  //
      const bool with_write,
      const bool with_delete,
      const bool is_shuffled)
  {
    const auto &target_ids = CreateTargetIDs(kMaxRecNumForTest, is_shuffled);
    const auto expect_update = with_write && !with_delete;

    if (with_write) {
      for (size_t i = 0; i < kMaxRecNumForTest; ++i) {
        const auto id = target_ids.at(i);
        VerifyWrite(id, id);
      }
    }
    if (with_delete) {
      for (size_t i = 0; i < kMaxRecNumForTest; ++i) {
        const auto id = target_ids.at(i);
        VerifyDelete(id, kExpectSuccess);
      }
    }
    for (size_t i = 0; i < kMaxRecNumForTest; ++i) {
      const auto key_id = target_ids.at(i);
      VerifyUpdate(key_id, key_id + 1, expect_update);
    }
    for (size_t i = 0; i < kMaxRecNumForTest; ++i) {
      const auto key_id = target_ids.at(i);
      const auto val_id = (expect_update) ? key_id + 1 : key_id;
      VerifyRead(key_id, val_id, expect_update);
    }
  }

  void
  VerifyDeletesWith(  //
      const bool with_write,
      const bool with_delete,
      const bool is_shuffled)
  {
    const auto &target_ids = CreateTargetIDs(kMaxRecNumForTest, is_shuffled);
    const auto expect_delete = with_write && !with_delete;

    if (with_write) {
      for (size_t i = 0; i < kMaxRecNumForTest; ++i) {
        const auto id = target_ids.at(i);
        VerifyWrite(id, id);
      }
    }
    if (with_delete) {
      for (size_t i = 0; i < kMaxRecNumForTest; ++i) {
        const auto id = target_ids.at(i);
        VerifyDelete(id, kExpectSuccess);
      }
    }
    for (size_t i = 0; i < kMaxRecNumForTest; ++i) {
      const auto key_id = target_ids.at(i);
      VerifyDelete(key_id, expect_delete);
    }
    for (size_t i = 0; i < kMaxRecNumForTest; ++i) {
      const auto key_id = target_ids.at(i);
      VerifyRead(key_id, key_id, kExpectFailed);
    }
  }

  /*################################################################################################
   * Internal member variables
   *##############################################################################################*/

  // actual keys and payloads
  Key keys_[kKeyNumForTest]{};
  Payload payloads_[kKeyNumForTest]{};

  std::unique_ptr<BTree_t> b_tree_{nullptr};
};

/*######################################################################################
 * Preparation for typed testing
 *####################################################################################*/

using TestTargets = ::testing::Types<  //
    KeyPayload<UInt8, UInt8>,          // fixed-length keys
    KeyPayload<UInt4, UInt8>,          // small keys
    KeyPayload<UInt8, UInt4>,          // small payloads
    KeyPayload<UInt4, UInt4>,          // small keys/payloads
    KeyPayload<Var, UInt8>,            // variable-length keys
    KeyPayload<Ptr, Ptr>,              // pointer keys/payloads
    KeyPayload<Original, Original>     // original class payloads
    >;
TYPED_TEST_SUITE(BTreeFixture, TestTargets);

/*######################################################################################
 * Unit test definitions
 *####################################################################################*/

/*--------------------------------------------------------------------------------------
 * Structure modification operations
 *------------------------------------------------------------------------------------*/

TYPED_TEST(BTreeFixture, WriteWithConsolidationReadWrittenValues)
{
  const size_t rec_num = TestFixture::kRecNumInNode;
  TestFixture::VerifyWritesWith(!kWriteTwice, !kShuffled, rec_num);
}

TYPED_TEST(BTreeFixture, WriteWithRootLeafSplitReadWrittenValues)
{
  const size_t rec_num = TestFixture::kRecNumInNode * 1.5;
  TestFixture::VerifyWritesWith(!kWriteTwice, !kShuffled, rec_num);
}

TYPED_TEST(BTreeFixture, WriteWithRootInternalSplitReadWrittenValues)
{
  const size_t rec_num = TestFixture::kRecNumInNode * TestFixture::kRecNumInNode;
  TestFixture::VerifyWritesWith(!kWriteTwice, !kShuffled, rec_num);
}

TYPED_TEST(BTreeFixture, WriteWithInternalSplitReadWrittenValues)
{
  TestFixture::VerifyWritesWith(!kWriteTwice, !kShuffled);
}

/*--------------------------------------------------------------------------------------
 * Read operation tests
 *------------------------------------------------------------------------------------*/

TYPED_TEST(BTreeFixture, ReadWithEmptyIndexFail)
{  //
  TestFixture::VerifyRead(0, 0, kExpectFailed);
}

/*--------------------------------------------------------------------------------------
 * Scan operation
 *------------------------------------------------------------------------------------*/

TYPED_TEST(BTreeFixture, ScanWithoutKeysPerformFullScan)
{
  const size_t rec_num = TestFixture::kMaxRecNumForTest;

  for (size_t i = 0; i < rec_num; ++i) {
    TestFixture::VerifyInsert(i, i, kExpectSuccess);
  }

  TestFixture::VerifyScan(std::nullopt, std::nullopt);
}

TYPED_TEST(BTreeFixture, ScanWithClosedRangeIncludeLeftRightEnd)
{
  const size_t rec_num = TestFixture::kMaxRecNumForTest;

  for (size_t i = 0; i < rec_num; ++i) {
    TestFixture::VerifyInsert(i, i, kExpectSuccess);
  }

  TestFixture::VerifyScan(std::make_pair(0, kRangeClosed),
                          std::make_pair(rec_num - 1, kRangeClosed));
}

TYPED_TEST(BTreeFixture, ScanWithOpenedRangeExcludeLeftRightEnd)
{
  const size_t rec_num = TestFixture::kMaxRecNumForTest;

  for (size_t i = 0; i < rec_num; ++i) {
    TestFixture::VerifyInsert(i, i, kExpectSuccess);
  }

  TestFixture::VerifyScan(std::make_pair(0, !kRangeClosed),
                          std::make_pair(rec_num - 1, !kRangeClosed));
}

/*--------------------------------------------------------------------------------------
 * Write operation
 *------------------------------------------------------------------------------------*/

TYPED_TEST(BTreeFixture, WriteWithDuplicateKeysSucceed)
{
  TestFixture::VerifyWritesWith(kWriteTwice, !kShuffled);
}

TYPED_TEST(BTreeFixture, RandomWriteWithDuplicateKeysSucceed)
{
  TestFixture::VerifyWritesWith(kWriteTwice, kShuffled);
}

/*--------------------------------------------------------------------------------------
 * Insert operation
 *------------------------------------------------------------------------------------*/

TYPED_TEST(BTreeFixture, InsertWithUniqueKeysReadInsertedValues)
{
  TestFixture::VerifyInsertsWith(!kWriteTwice, !kWithDelete, !kShuffled);
}

TYPED_TEST(BTreeFixture, InsertWithDuplicateKeysFail)
{
  TestFixture::VerifyInsertsWith(kWriteTwice, !kWithDelete, !kShuffled);
}

TYPED_TEST(BTreeFixture, InsertWithDeletedKeysSucceed)
{
  TestFixture::VerifyInsertsWith(kWriteTwice, kWithDelete, !kShuffled);
}

TYPED_TEST(BTreeFixture, RandomInsertWithUniqueKeysReadInsertedValues)
{
  TestFixture::VerifyInsertsWith(!kWriteTwice, !kWithDelete, kShuffled);
}

TYPED_TEST(BTreeFixture, RandomInsertWithDuplicateKeysFail)
{
  TestFixture::VerifyInsertsWith(kWriteTwice, !kWithDelete, kShuffled);
}

TYPED_TEST(BTreeFixture, RandomInsertWithDeletedKeysSucceed)
{
  TestFixture::VerifyInsertsWith(kWriteTwice, kWithDelete, kShuffled);
}

/*--------------------------------------------------------------------------------------
 * Update operation
 *------------------------------------------------------------------------------------*/

TYPED_TEST(BTreeFixture, UpdateWithDuplicateKeysSucceed)
{
  TestFixture::VerifyUpdatesWith(kWithWrite, !kWithDelete, !kShuffled);
}

TYPED_TEST(BTreeFixture, UpdateWithNotInsertedKeysFail)
{
  TestFixture::VerifyUpdatesWith(!kWithWrite, !kWithDelete, !kShuffled);
}

TYPED_TEST(BTreeFixture, UpdateWithDeletedKeysFail)
{
  TestFixture::VerifyUpdatesWith(kWithWrite, kWithDelete, !kShuffled);
}

TYPED_TEST(BTreeFixture, RandomUpdateWithDuplicateKeysSucceed)
{
  TestFixture::VerifyUpdatesWith(kWithWrite, !kWithDelete, kShuffled);
}

TYPED_TEST(BTreeFixture, RandomUpdateWithNotInsertedKeysFail)
{
  TestFixture::VerifyUpdatesWith(!kWithWrite, !kWithDelete, kShuffled);
}

TYPED_TEST(BTreeFixture, RandomUpdateWithDeletedKeysFail)
{
  TestFixture::VerifyUpdatesWith(kWithWrite, kWithDelete, kShuffled);
}

/*--------------------------------------------------------------------------------------
 * Delete operation
 *------------------------------------------------------------------------------------*/

TYPED_TEST(BTreeFixture, DeleteWithDuplicateKeysSucceed)
{
  TestFixture::VerifyDeletesWith(kWithWrite, !kWithDelete, !kShuffled);
}

TYPED_TEST(BTreeFixture, DeleteWithNotInsertedKeysFail)
{
  TestFixture::VerifyDeletesWith(!kWithWrite, !kWithDelete, !kShuffled);
}

TYPED_TEST(BTreeFixture, DeleteWithDeletedKeysFail)
{
  TestFixture::VerifyDeletesWith(kWithWrite, kWithDelete, !kShuffled);
}

TYPED_TEST(BTreeFixture, RandomDeleteWithDuplicateKeysSucceed)
{
  TestFixture::VerifyDeletesWith(kWithWrite, !kWithDelete, kShuffled);
}

TYPED_TEST(BTreeFixture, RandomDeleteWithNotInsertedKeysFail)
{
  TestFixture::VerifyDeletesWith(!kWithWrite, !kWithDelete, kShuffled);
}

TYPED_TEST(BTreeFixture, RandomDeleteWithDeletedKeysFail)
{
  TestFixture::VerifyDeletesWith(kWithWrite, kWithDelete, kShuffled);
}

}  // namespace dbgroup::index::b_tree::component
