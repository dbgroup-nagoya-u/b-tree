cmake_minimum_required(VERSION "3.8")

# decline in-source build
if("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}")
  message(FATAL_ERROR "In-source builds are not allowed.")
endif()

# project information
project(
  b_tree
  VERSION "0.1.0"
  DESCRIPTION "An open-source B+tree library."
  LANGUAGES "CXX"
)

# define source directories
set(B_TREE_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}")

#--------------------------------------------------------------------------------------#
# B-tree's configurations
#--------------------------------------------------------------------------------------#

set(B_TREE_PAGE_SIZE "1024" CACHE STRING "Page size in bytes.")
set(B_TREE_MAX_DELETED_SPACE_SIZE "(kPageSize / 4)" CACHE STRING "Invoking a clean-up-operation if the deleted space size of a node exceeds this threshold.")
set(B_TREE_MIN_FREE_SPACE_SIZE "(kPageSize / 8)" CACHE STRING "Invoking a split operation if the free space size of a node after cleaning up becomes lower than this threshold.")
set(B_TREE_MIN_USED_SPACE_SIZE "(kPageSize / 16)" CACHE STRING "Invoking a merge operation if the used space size of a node becomes lower than this threshold.")
set(B_TREE_MAX_VARLEN_DATA_SIZE "128" CACHE STRING "The expected maximum size of a variable-length data.")

#--------------------------------------------------------------------------------------#
# Configure cpp-utility
#--------------------------------------------------------------------------------------#
add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/external/cpp-utility")

#--------------------------------------------------------------------------------------#
# Configure garbage collector
#--------------------------------------------------------------------------------------#
add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/external/memory-manager")

#--------------------------------------------------------------------------------------#
# Build targets
#--------------------------------------------------------------------------------------#

if(NOT TARGET b_tree::b_tree)
  add_library(b_tree INTERFACE)
  add_library(b_tree::b_tree ALIAS b_tree)
  target_compile_features(b_tree INTERFACE
    "cxx_std_17"
  )
  target_compile_definitions(b_tree INTERFACE
    B_TREE_PAGE_SIZE=${B_TREE_PAGE_SIZE}
    B_TREE_MAX_DELETED_SPACE_SIZE=${B_TREE_MAX_DELETED_SPACE_SIZE}
    B_TREE_MIN_FREE_SPACE_SIZE=${B_TREE_MIN_FREE_SPACE_SIZE}
    B_TREE_MIN_USED_SPACE_SIZE=${B_TREE_MIN_USED_SPACE_SIZE}
    B_TREE_MAX_VARLEN_DATA_SIZE=${B_TREE_MAX_VARLEN_DATA_SIZE}
  )
  target_include_directories(b_tree INTERFACE
    "${B_TREE_SOURCE_DIR}/include"
  )
  target_link_libraries(b_tree INTERFACE
    cpp_utility::cpp_utility
    memory_manager::memory_manager
  )
endif()

#--------------------------------------------------------------------------------------#
# Build unit tests if required
#--------------------------------------------------------------------------------------#

option(B_TREE_BUILD_TESTS, "Building unit tests for this library if ON." OFF)
if(${B_TREE_BUILD_TESTS})
  # enable CTest
  enable_testing()
  # include unit tests
  add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/test")
endif()
