cmake_minimum_required(VERSION "3.8")

# decline in-source build
if("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}")
  message(FATAL_ERROR "In-source builds are not allowed.")
endif()

# project information
project(
  hello_world
  VERSION "0.1.0"
  DESCRIPTION "A template project for C++ language."
  LANGUAGES "CXX"
)

# define source directories
set(BTREE_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}")

#--------------------------------------------------------------------------------------#
# Bw-tree's configurations
#--------------------------------------------------------------------------------------#

set(BTREE_PAGE_SIZE "1024" CACHE STRING "Page size in bytes.")
set(BTREE_MAX_DELETED_SPACE_SIZE "(kPageSize * 15 / 16)" CACHE STRING "Invoking a split-operation if the deleted space size of a node exceeds this threshold.")
set(BTREE_MIN_FREE_SPACE_SIZE "(kPageSize / 16)" CACHE STRING "Invoking a split-operation if the free space size of a node becomes lower than this threshold.")
set(BTREE_MIN_USED_SPACE_SIZE "(kPageSize / 16)" CACHE STRING "Invoking merge-operation if if the used space size of a node becomes lower than this threshold.")
set(BTREE_MAX_VARIABLE_DATA_SIZE "128" CACHE STRING "The expected maximum size of a variable-length data.")

#--------------------------------------------------------------------------------------#
# Build targets
#--------------------------------------------------------------------------------------#

add_library(btree INTERFACE)
target_compile_features(btree INTERFACE
  "cxx_std_17"
)
target_compile_definitions(btree INTERFACE
  BTREE_PAGE_SIZE=${BTREE_PAGE_SIZE}
  BTREE_MAX_DELETED_SPACE_SIZE=${BTREE_MAX_DELETED_SPACE_SIZE}
  BTREE_MIN_FREE_SPACE_SIZE=${BTREE_MIN_FREE_SPACE_SIZE}
  BTREE_MIN_USED_SPACE_SIZE=${BTREE_MIN_USED_SPACE_SIZE}
  BTREE_MAX_VARIABLE_DATA_SIZE=${BTREE_MAX_VARIABLE_DATA_SIZE}
)
target_compile_options(
  btree INTERFACE
  -Wall
  -Wextra
  $<$<STREQUAL:${CMAKE_BUILD_TYPE},"Release">:"-O2 -march=native">
  $<$<STREQUAL:${CMAKE_BUILD_TYPE},"RelWithDebInfo">:"-g3 -Og -pg">
  $<$<STREQUAL:${CMAKE_BUILD_TYPE},"Debug">:"-g3 -O0 -pg">
)
target_include_directories(
  btree INTERFACE
  "${BTREE_SOURCE_DIR}/include"
)

#--------------------------------------------------------------------------------------#
# Build unit tests if required
#--------------------------------------------------------------------------------------#

option(BUILD_TESTS, "Build unit tests" OFF)
if(${BUILD_TESTS})
  # enable CTest
  enable_testing()
  # include unit tests
  add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/test")
endif()
