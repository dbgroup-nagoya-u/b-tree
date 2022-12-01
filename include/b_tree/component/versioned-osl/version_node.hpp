
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

namespace dbgroup::index::b_tree::component::versioned_osl
{

/**
 * @brief A class for representing version nodes
 *
 * @tparam Payload a target Payload class
 * @tparam Timestamp_t a timestamp class
 */
template <class Payload, class Timestamp_t>
class VersionNode
{
 public:
  /*####################################################################################
   * Type aliases
   *##################################################################################*/

  // using Timestamp_t = size_t;

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
  VersionNode(  //
      Timestamp_t timestamp,
      Payload payload,
      VersionNode<Payload> *next = nullptr)
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
  ~VersionNode() = default;

  /*##################################################################################
   * Public Getters
   *################################################################################*/

  [[nodiscard]] const auto
  GetTimestamp() const  //
      -> Timestamp_t
  {
    // TODO: memcpyが必要？要確認
    // そもそもpublicメンバにすればいいかも？
    return timestamp_;
  }

  [[nodiscard]] const auto
  GetPayload() const  //
      -> Payload
  {
    // TODO: memcpyが必要？要確認
    return payload_;
  }

  [[nodiscard]] const auto
  GetNextPtr() const  //
      -> VersionNode<Payload> *
  {
    return next_;
  }
  /*##################################################################################
   * Public Setter
   *################################################################################*/
  auto
  SetNextPtr(VersionNode *node_ptr)  //
      -> void
  {
    // TODO: memcpyが必要？要確認
    next_ = node_ptr;
    return;
  }

  /*##################################################################################
   * Public utility function
   *################################################################################*/

  /*##################################################################################
   * Internal member
   *################################################################################*/
  const Timestamp_t timestamp_;
  const Payload payload_;
  VersionNode *next_;
};

}  // namespace dbgroup::index::b_tree::component::versioned_osl
