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

#include "b_tree/component/metadata.hpp"

#include "common.hpp"
#include "gtest/gtest.h"

namespace dbgroup::index::b_tree::component::test
{
class MetadataFixture : public testing::Test
{
 protected:
  /*################################################################################################
   * Internal constants
   *##############################################################################################*/

  static constexpr size_t kExpectedOffset = 256;
  static constexpr size_t kExpectedKeyLength = 8;
  static constexpr size_t kExpectedTotalLength = 16;

  /*################################################################################################
   * Setup/Teardown
   *##############################################################################################*/

  void
  SetUp() override
  {
    meta_ = Metadata{kExpectedOffset, kExpectedKeyLength, kExpectedTotalLength};
  }

  void
  TearDown() override
  {
  }

  /*################################################################################################
   * Internal member variables
   *##############################################################################################*/

  Metadata meta_{};
};

/*--------------------------------------------------------------------------------------------------
 * Constructor tests
 *------------------------------------------------------------------------------------------------*/

TEST_F(MetadataFixture, ConstructDefaultMetadataCorrectlyInitialized)
{
  EXPECT_EQ(kExpectedOffset, meta_.offset);
  EXPECT_EQ(kExpectedKeyLength, meta_.key_length);
  EXPECT_EQ(kExpectedTotalLength, meta_.total_length);
}

/*--------------------------------------------------------------------------------------------------
 * Getter/setter tests
 *------------------------------------------------------------------------------------------------*/

TEST_F(MetadataFixture, GetPayloadLengthDefaultMetadataReturnCorrectPayloadLength)
{
  EXPECT_EQ(kExpectedTotalLength - kExpectedKeyLength, meta_.total_length - meta_.key_length);
}

TEST_F(MetadataFixture, SetOffsetDefaultMetadataGetUpdatedOffset)
{
  const size_t updated_offset = kExpectedOffset / 2;

  meta_.offset = updated_offset;
  EXPECT_EQ(updated_offset, meta_.offset);
}

TEST_F(MetadataFixture, SetTotalLengthDefaultMetadataGetUpdatedTotalLength)
{
  const size_t updated_total_length = kExpectedTotalLength / 2;

  meta_.total_length = updated_total_length;
  EXPECT_EQ(updated_total_length, meta_.total_length);
}

/*--------------------------------------------------------------------------------------------------
 * Operator tests
 *------------------------------------------------------------------------------------------------*/

TEST_F(MetadataFixture, EQWithSameMetadataReturnTrue)
{
  const Metadata meta_a{kExpectedOffset, kExpectedKeyLength, kExpectedTotalLength};
  const Metadata meta_b{kExpectedOffset, kExpectedKeyLength, kExpectedTotalLength};

  EXPECT_TRUE(meta_a == meta_b);
}

TEST_F(MetadataFixture, EQWithDifferentMetadataReturnFalse)
{
  const Metadata meta_a{kExpectedOffset, kExpectedKeyLength, kExpectedTotalLength};
  const Metadata meta_b{};

  EXPECT_FALSE(meta_a == meta_b);
}

TEST_F(MetadataFixture, NEQWithSameMetadataReturnFalse)
{
  const Metadata meta_a{kExpectedOffset, kExpectedKeyLength, kExpectedTotalLength};
  const Metadata meta_b{kExpectedOffset, kExpectedKeyLength, kExpectedTotalLength};

  EXPECT_FALSE(meta_a != meta_b);
}

TEST_F(MetadataFixture, NEQWithDifferentMetadataReturnTrue)
{
  const Metadata meta_a{kExpectedOffset, kExpectedKeyLength, kExpectedTotalLength};
  const Metadata meta_b{};

  EXPECT_TRUE(meta_a != meta_b);
}

}  // namespace dbgroup::index::b_tree::component::test
