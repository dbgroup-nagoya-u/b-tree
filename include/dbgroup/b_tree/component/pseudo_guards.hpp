/*
 * Copyright 2026 Database Group, Nagoya University
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

#ifndef B_TREE_DBGROUP_B_TREE_COMPONENT_PSEUDO_GUARDS_HPP_
#define B_TREE_DBGROUP_B_TREE_COMPONENT_PSEUDO_GUARDS_HPP_

// C++ standard libraries
#include <cstdint>
#include <utility>

namespace dbgroup::index::b_tree::component
{
// forward declarations
template <class Lock>
class PseudoXGuard;

/**
 * @brief A guard that holds an X lock in place of an SIX lock.
 *
 * This class is used for comparing lock protocols. It provides the same
 * interfaces as SIX-lock guards, but it acquires an X lock from the beginning,
 * so that optimistic readers cannot read a node until a writer releases it.
 *
 * @tparam Lock A lock class for concurrency controls.
 */
template <class Lock>
class PseudoSIXGuard
{
 public:
  /*##########################################################################*
   * Public constructors and assignment operators
   *##########################################################################*/

  constexpr PseudoSIXGuard() noexcept = default;

  /**
   * @param x_grd An X-lock guard to be held.
   */
  explicit PseudoSIXGuard(  //
      typename Lock::XGuard&& x_grd) noexcept
      : x_grd_{std::move(x_grd)}
      , ver_{x_grd_.GetVersion()}
  {
  }

  /**
   * @param x_grd An X-lock guard to be held.
   * @param ver A version that this guard regards as the current one.
   */
  PseudoSIXGuard(  //
      typename Lock::XGuard&& x_grd,
      const uint32_t ver) noexcept
      : x_grd_{std::move(x_grd)}
      , ver_{ver}
  {
  }

  PseudoSIXGuard(const PseudoSIXGuard&) = delete;
  PseudoSIXGuard(PseudoSIXGuard&&) noexcept = default;

  auto operator=(const PseudoSIXGuard&) -> PseudoSIXGuard& = delete;
  auto operator=(PseudoSIXGuard&&) noexcept -> PseudoSIXGuard& = default;

  /*##########################################################################*
   * Public destructors
   *##########################################################################*/

  /**
   * @brief Destroy this instance and release a lock if holding.
   *
   */
  ~PseudoSIXGuard() = default;

  /*##########################################################################*
   * Public APIs
   *##########################################################################*/

  /**
   * @retval true if this instance has the lock ownership.
   * @retval false otherwise.
   */
  constexpr explicit
  operator bool() const noexcept
  {
    return static_cast<bool>(x_grd_);
  }

  /**
   * @return The version when this guard was created.
   */
  [[nodiscard]]
  constexpr auto
  GetVersion() const noexcept  //
      -> uint32_t
  {
    return ver_;
  }

  /**
   * @brief Convert this guard to an X-lock guard without releasing the lock.
   *
   * @return The lock guard for an X lock.
   * @note After calling the function, this lock guard abandons the lock's
   * ownership.
   */
  [[nodiscard]]
  auto
  UpgradeToX() noexcept  //
      -> PseudoXGuard<Lock>
  {
    return PseudoXGuard<Lock>{std::move(x_grd_), ver_};
  }

 private:
  /*##########################################################################*
   * Internal member variables
   *##########################################################################*/

  /// @brief An actual X-lock guard.
  typename Lock::XGuard x_grd_{};

  /// @brief A version that this guard regards as the current one.
  uint32_t ver_{};
};

/**
 * @brief A guard that holds an X lock and is convertible to `PseudoSIXGuard`.
 *
 * @tparam Lock A lock class for concurrency controls.
 */
template <class Lock>
class PseudoXGuard
{
 public:
  /*##########################################################################*
   * Public constructors and assignment operators
   *##########################################################################*/

  constexpr PseudoXGuard() noexcept = default;

  /**
   * @param x_grd An X-lock guard to be held.
   * @param ver A version that this guard regards as the current one.
   * @note As with actual X-lock guards, the version is incremented by one when
   * this guard releases the lock unless a desired version is set.
   */
  PseudoXGuard(  //
      typename Lock::XGuard&& x_grd,
      const uint32_t ver) noexcept
      : x_grd_{std::move(x_grd)}
      , ver_{ver}
      , new_ver_{ver + 1U}
  {
    if (x_grd_) {
      x_grd_.SetVersion(new_ver_);
    }
  }

  PseudoXGuard(const PseudoXGuard&) = delete;
  PseudoXGuard(PseudoXGuard&&) noexcept = default;

  auto operator=(const PseudoXGuard&) -> PseudoXGuard& = delete;
  auto operator=(PseudoXGuard&&) noexcept -> PseudoXGuard& = default;

  /*##########################################################################*
   * Public destructors
   *##########################################################################*/

  /**
   * @brief Destroy this instance and release a lock if holding.
   *
   */
  ~PseudoXGuard() = default;

  /*##########################################################################*
   * Public APIs
   *##########################################################################*/

  /**
   * @retval true if this instance has the lock ownership.
   * @retval false otherwise.
   */
  constexpr explicit
  operator bool() const noexcept
  {
    return static_cast<bool>(x_grd_);
  }

  /**
   * @return The version when this guard was created.
   */
  [[nodiscard]]
  constexpr auto
  GetVersion() const noexcept  //
      -> uint32_t
  {
    return ver_;
  }

  /**
   * @brief Set a desired version after unlocking.
   *
   * @param ver A desired version after unlocking.
   */
  void
  SetVersion(  //
      const uint32_t ver) noexcept
  {
    new_ver_ = ver;
    x_grd_.SetVersion(ver);
  }

  /**
   * @brief Convert this guard to a `PseudoSIXGuard` without releasing the lock.
   *
   * @return The lock guard that regards a desired version as the current one.
   * @note After calling the function, this lock guard abandons the lock's
   * ownership.
   */
  [[nodiscard]]
  auto
  DowngradeToSIX() noexcept  //
      -> PseudoSIXGuard<Lock>
  {
    return PseudoSIXGuard<Lock>{std::move(x_grd_), new_ver_};
  }

 private:
  /*##########################################################################*
   * Internal member variables
   *##########################################################################*/

  /// @brief An actual X-lock guard.
  typename Lock::XGuard x_grd_{};

  /// @brief A version that this guard regards as the current one.
  uint32_t ver_{};

  /// @brief A desired version after unlocking.
  uint32_t new_ver_{};
};

}  // namespace dbgroup::index::b_tree::component

#endif  // B_TREE_DBGROUP_B_TREE_COMPONENT_PSEUDO_GUARDS_HPP_
