/*
 * Copyright 2024 Database Group, Nagoya University
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

#ifndef TEST_PAGE_SIZE_WRAPPER_HPP_
#define TEST_PAGE_SIZE_WRAPPER_HPP_

// C++ standard libraries
#include <cstddef>
#include <memory>
#include <optional>
#include <tuple>
#include <type_traits>
#include <vector>

namespace dbgroup::index::test
{
/*############################################################################*
 * Global constants
 *############################################################################*/

constexpr size_t kPageSize = (DBGROUP_TEST_PAGE_SIZE);
constexpr size_t kGCIntervalMS = 10;

/*############################################################################*
 * Fixture definitions
 *############################################################################*/

template <template <class Key, class Payload, class Comp, class... Others> class IndexBase,
          class Key,
          class Payload,
          class Comp>
class PageSizeWrapper
{
  /*##########################################################################*
   * Type aliases
   *##########################################################################*/

  using ScanKey = std::optional<std::tuple<Key, size_t, bool>>;
  using Index = IndexBase<Key, Payload, Comp>;

 public:
  /*##########################################################################*
   * Public constructors and assignment operators
   *##########################################################################*/

  PageSizeWrapper() = default;

  PageSizeWrapper(const PageSizeWrapper &) = delete;
  PageSizeWrapper(PageSizeWrapper &&) = delete;

  auto operator=(const PageSizeWrapper &) -> PageSizeWrapper & = delete;
  auto operator=(PageSizeWrapper &&) -> PageSizeWrapper & = delete;

  /*##########################################################################*
   * Public destructors
   *##########################################################################*/

  ~PageSizeWrapper() = default;

  /*##########################################################################*
   * Public read/write APIs
   *##########################################################################*/

  auto
  Read(  //
      const Key &key,
      const size_t key_len = sizeof(Key))
  {
    return index_->Read(key, key_len);
  }

  auto
  Scan(  //
      const ScanKey &begin_key = std::nullopt,
      const ScanKey &end_key = std::nullopt)
  {
    return index_->Scan(begin_key, end_key);
  }

  auto
  ScanBackward(  //
      const ScanKey &begin_key = std::nullopt,
      const ScanKey &end_key = std::nullopt)
  {
    return index_->ScanBackward(begin_key, end_key);
  }

  void
  Write(  //
      const Key &key,
      const Payload &payload,
      const size_t key_len = sizeof(Key),
      Payload (*merger)(const Payload &, const Payload &) = nullptr)
  {
    index_->Write(key, payload, key_len, merger);
  }

  auto
  Upsert(  //
      const Key &key,
      const Payload &payload,
      const size_t key_len = sizeof(Key),
      Payload (*merger)(const Payload &, const Payload &) = nullptr)
  {
    return index_->Upsert(key, payload, key_len, merger);
  }

  auto
  Insert(  //
      const Key &key,
      const Payload &payload,
      const size_t key_len = sizeof(Key))
  {
    return index_->Insert(key, payload, key_len);
  }

  template <class T = Payload>
  auto
  Update(  //
      const Key &key,
      const std::type_identity_t<T> &payload,
      const size_t key_len = sizeof(Key),
      Payload (*merger)(const Payload &, const T &) = nullptr,
      const size_t pay_len = sizeof(T))
  {
    return index_->Update(key, payload, key_len, merger, pay_len);
  }

  auto
  Delete(  //
      const Key &key,
      const size_t key_len = sizeof(Key))
  {
    return index_->Delete(key, key_len);
  }

  /*##########################################################################*
   * Public utility APIs
   *##########################################################################*/

  template <class Entry>
  auto
  Bulkload(  //
      const std::vector<Entry> &entries,
      const size_t thread_num = 1)
  {
    return index_->Bulkload(entries, thread_num);
  }

 private:
  /*##########################################################################*
   * Internal member variables
   *##########################################################################*/

  std::unique_ptr<Index> index_{std::make_unique<Index>(kPageSize, kGCIntervalMS)};
};

}  // namespace dbgroup::index::test

#endif  // TEST_PAGE_SIZE_WRAPPER_HPP_
