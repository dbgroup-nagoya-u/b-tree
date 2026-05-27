cmake_minimum_required(VERSION "3.11")

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
    B_TREE_MAX_VARLEN_DATA_SIZE
    "32" CACHE STRING
    "The expected maximum size of a variable-length data."
  )

  #----------------------------------------------------------------------------#
  # Configurations
  #----------------------------------------------------------------------------#

  # use multi-threads
  find_package(Threads)

  # use our C++ utilities
  FetchContent_Declare(
    cpp-utility
    GIT_REPOSITORY "https://github.com/dbgroup-nagoya-u/cpp-utility.git"
    GIT_TAG "8da07347db75c1498a849a8817534590ceb67738"
  )
  FetchContent_MakeAvailable(cpp-utility)

   # use our GC library
  FetchContent_Declare(
    memory-manager
    GIT_REPOSITORY "https://github.com/dbgroup-nagoya-u/memory-manager.git"
    GIT_TAG "5982e394ff15d14db39c8ecd0c4516d65faf4074"
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
    B_TREE_MAX_VARLEN_DATA_SIZE=${B_TREE_MAX_VARLEN_DATA_SIZE}
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
