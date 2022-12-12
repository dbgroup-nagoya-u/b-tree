
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

#ifndef B_TREE_COMPONENT_VERSIONED_OSL_VERSION_NODE
#define B_TREE_COMPONENT_VERSIONED_OSL_VERSION_NODE

namespace dbgroup::index::b_tree::component::osl
{

/**
 * @brief A class for representing version records
 *
 * @tparam Payload a target Payload class
 * @tparam Timestamp_t a timestamp class
 */
template <class Payload, class Timestamp_t>
class VersionRecord
{
 public:
  /*##################################################################################
   * Public constructors
   *################################################################################*/

  /**
   * @brief Construct a new instance.
   *
   * @param timestamp a timestamp which represents when the version was written
   * @param payload a payload of tuple
   * @param next a pointer to a next node.
   */
  VersionRecord(  //
      Timestamp_t timestamp,
      Payload payload,
      VersionRecord<Payload, Timestamp_t> *next = nullptr)
      : timestamp_{timestamp}, payload_{payload}, next_{next}
  {
  }
  /*##################################################################################
   * Public destructor
   *################################################################################*/

  /**
   * @brief Destroy the instance.
   *
   */
  ~VersionRecord() = default;

  /*##################################################################################
   * Public Getters
   *################################################################################*/

  /**
   * @return the timestamp(i.e., version) of this node
   */
  [[nodiscard]] const auto
  GetTimestamp() const  //
      -> Timestamp_t
  {
    return timestamp_;
  }

  /**
   * @return the payload of this version
   */
  [[nodiscard]] const auto
  GetPayload() const  //
      -> Payload
  {
    return payload_;
  }

  /**
   * @return a pointer to the next node on chain
   * @return nullptr when this node is the latest version
   */
  [[nodiscard]] const auto
  GetNextPtr() const  //
      -> VersionRecord<Payload, Timestamp_t> *
  {
    return next_;
  }

  /*##################################################################################
   * Public Setter
   *################################################################################*/

  /**
   * @brief Set a next node pointer directly
   *
   * @param node_ptr a pointer to the node which will be the next of this on the chain
   */
  auto
  SetNextPtr(VersionRecord *node_ptr)  //
      -> void
  {
    next_ = node_ptr;
    return;
  }

  /*##################################################################################
   * Internal member
   *################################################################################*/

  const Timestamp_t timestamp_;
  const Payload payload_;
  VersionRecord *next_;
};

}  // namespace dbgroup::index::b_tree::component::osl
#endif  // B_TREE_COMPONENT_VERSIONED_OSL_VERSION_NODE
