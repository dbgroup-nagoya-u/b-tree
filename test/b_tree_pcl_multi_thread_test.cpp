#include <random>
#include <thread>
#include <vector>

#include "b_tree/b_tree_pcl.hpp"
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
class BTreePCLFixture : public testing::Test  // NOLINT
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
  using BTreePCL_t = BTreePCL<Key, Payload, KeyComp>;

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

    b_tree_pcl_ = std::make_unique<BTreePCL_t>();
  }

  void
  TearDown() override
  {
    ReleaseTestData(keys_, kKeyNumForTest);      // NOLINT
    ReleaseTestData(payloads_, kKeyNumForTest);  // NOLINT
    b_tree_pcl_.reset(nullptr);
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
    return b_tree_pcl_->Insert(keys_[key_id], payloads_[payload_id], kKeyLen);
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
    const auto read_val = b_tree_pcl_->Read(keys_[key_id]);
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
  VerifyWrite(  //
      const size_t key_id,
      const size_t pay_id)
  {
    auto rc = b_tree_pcl_->Write(keys_[key_id], payloads_[pay_id], kKeyLen);

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
  VerifyReadMultiThreads(  //
      const bool is_shuffled)
  {
    const auto &target_ids = CreateTargetIDs(kMaxRecNumForTest, is_shuffled);
    // Insert key with single thread
    for (size_t i = 0; i < kMaxRecNumForTest; i++) {
      const auto id = target_ids.at(i);
      VerifyInsert(id, id, kExpectSuccess);
    }

    // Read payload with multi-threads
    std::vector<std::thread> threads{};
    threads.reserve(kThreadNum);

    auto f = [&]() {
      for (size_t i = 0; i < kMaxRecNumForTest; i++) {
        const auto id = target_ids.at(i);
        VerifyRead(id, id, kExpectSuccess);
      }
    };

    for (size_t i = 0; i < kThreadNum; ++i) {
      threads.emplace_back(f);
    }

    for (auto &&t : threads) {
      t.join();
    }
  }
  /*################################################################################################
   * Internal member variables
   *##############################################################################################*/

  // actual keys and payloads
  Key keys_[kKeyNumForTest]{};
  Payload payloads_[kKeyNumForTest]{};

  std::unique_ptr<BTreePCL_t> b_tree_pcl_{nullptr};
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
TYPED_TEST_SUITE(BTreePCLFixture, TestTargets);

/*######################################################################################
 * Unit test definitions
 *####################################################################################*/

/*--------------------------------------------------------------------------------------
 * Read operation tests
 *------------------------------------------------------------------------------------*/

TYPED_TEST(BTreePCLFixture, ReadWithEmptyIndexFail)
{  //
  TestFixture::VerifyRead(0, 0, kExpectFailed);
}

TYPED_TEST(BTreePCLFixture, ReadWithMultiThreadSuccess)
{  //
  TestFixture::VerifyReadMultiThreads(true);
}

}  // namespace dbgroup::index::b_tree::component
