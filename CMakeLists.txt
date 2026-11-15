cmake_minimum_required(VERSION "3.11")

# decline in-source build
if("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}")
  message(FATAL_ERROR "In-source builds are not allowed.")
endif()

# project information
project(
  b_tree
  VERSION "0.2.0"
  DESCRIPTION "An open-source B+tree library."
  LANGUAGES "CXX"
)

if(NOT TARGET dbgroup::${PROJECT_NAME})
  # enable features
  include(FetchContent)

  #----------------------------------------------------------------------------#
  # Build parameters
  #----------------------------------------------------------------------------#

  set(
    B_TREE_DEFAULT_PAGE_SIZE
    "1024" CACHE STRING
    "A default page size in bytes."
  )

  set(
    B_TREE_DEFAULT_GC_INTERVAL
    "100" CACHE STRING
    "A default interval for GC in milliseconds."
  )

  set(
    B_TREE_DEFAULT_GC_THREADS
    "1" CACHE STRING
    "The default number of GC threads."
  )

  set(
    B_TREE_MAX_VARLEN_KEY_SIZE
    "32" CACHE STRING
    "The expected maximum size of a variable-length key."
  )

  set(
    B_TREE_MAX_VARLEN_PAYLOAD_SIZE
    "32" CACHE STRING
    "The expected maximum size of a variable-length payload."
  )

  #----------------------------------------------------------------------------#
  # Configurations
  #----------------------------------------------------------------------------#

  # use multi-threads
  find_package(Threads)

  # use our C++ utilities
  FetchContent_Declare(
    cpp-utility
    GIT_REPOSITORY "https://github.com/dbgroup-nagoya-u/cpp-utility"
    GIT_TAG "c1b810e65630b1f9964802069c6a0249bddedfa0"
  )
  FetchContent_MakeAvailable(cpp-utility)

   # use our local GC
  FetchContent_Declare(
    memory-manager
    GIT_REPOSITORY "https://github.com/dbgroup-nagoya-u/memory-manager"
    GIT_TAG "ef90a32f78ae1bb915ce59aa1a2d0eaab1a838de" # v0.7.0
  )
  FetchContent_MakeAvailable(memory-manager)

  #----------------------------------------------------------------------------#
  # Build targets
  #----------------------------------------------------------------------------#

  add_library(${PROJECT_NAME} INTERFACE)
  add_library(dbgroup::${PROJECT_NAME} ALIAS ${PROJECT_NAME})
  target_compile_features(${PROJECT_NAME} INTERFACE
    "cxx_std_20"
  )
  target_compile_definitions(${PROJECT_NAME} INTERFACE
    B_TREE_DEFAULT_PAGE_SIZE=${B_TREE_DEFAULT_PAGE_SIZE}
    B_TREE_DEFAULT_GC_INTERVAL=${B_TREE_DEFAULT_GC_INTERVAL}
    B_TREE_DEFAULT_GC_THREADS=${B_TREE_DEFAULT_GC_THREADS}
    B_TREE_MAX_VARLEN_KEY_SIZE=${B_TREE_MAX_VARLEN_KEY_SIZE}
    B_TREE_MAX_VARLEN_PAYLOAD_SIZE=${B_TREE_MAX_VARLEN_PAYLOAD_SIZE}
  )
  target_include_directories(${PROJECT_NAME} INTERFACE
    "${CMAKE_CURRENT_SOURCE_DIR}/include"
  )
  target_link_libraries(${PROJECT_NAME} INTERFACE
    dbgroup::cpp_utility
    dbgroup::memory_manager
  )

  #----------------------------------------------------------------------------#
  # Build unit tests if required
  #----------------------------------------------------------------------------#

  option(B_TREE_BUILD_TESTS, "Building unit tests for this library if ON." OFF)
  if(${B_TREE_BUILD_TESTS})
    # enable CTest
    enable_testing()
    # include unit tests
    add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/test")
  endif()
endif()
