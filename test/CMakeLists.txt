#--------------------------------------------------------------------------------------#
# Configure GoogleTest
#--------------------------------------------------------------------------------------#

include(FetchContent)
FetchContent_Declare(
  googletest
  GIT_REPOSITORY "https://github.com/google/googletest.git"
  GIT_TAG "58d77fa8070e8cec2dc1ed015d66b454c8d78850" # v1.12.1
)
FetchContent_MakeAvailable(googletest)

#--------------------------------------------------------------------------------------#
# Configurations for testing
#--------------------------------------------------------------------------------------#

set(DBGROUP_TEST_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}")
set(DBGROUP_TEST_THREAD_NUM "8" CACHE STRING "The number of worker threads for testing.")
set(DBGROUP_TEST_RANDOM_SEED "10" CACHE STRING "A fixed seed value to reproduce unit tests.")
set(DBGROUP_TEST_EXEC_NUM "1E5" CACHE STRING "The number of executions per a thread.")
option(DBGROUP_TEST_BUILD_MULTI_THREAD_TESTS "Building multi-threading tests if ON." ON)
option(DBGROUP_TEST_OVERRIDE_MIMALLOC "Override entire memory allocation with mimalloc" OFF)

#--------------------------------------------------------------------------------------#
# Build Unit Tests
#--------------------------------------------------------------------------------------#

if(${DBGROUP_TEST_OVERRIDE_MIMALLOC})
  find_package(mimalloc 1.7 REQUIRED)
endif()

# define function to add unit tests in the same format
function(ADD_DBGROUP_TEST DBGROUP_TEST_TARGET)
  add_executable(${DBGROUP_TEST_TARGET}
    "${CMAKE_CURRENT_SOURCE_DIR}/${DBGROUP_TEST_TARGET}.cpp"
  )
  target_compile_features(${DBGROUP_TEST_TARGET} PRIVATE
    "cxx_std_17"
  )
  target_compile_definitions(${DBGROUP_TEST_TARGET} PRIVATE
    DBGROUP_TEST_THREAD_NUM=${DBGROUP_TEST_THREAD_NUM}
    DBGROUP_TEST_RANDOM_SEED=${DBGROUP_TEST_RANDOM_SEED}
    DBGROUP_TEST_EXEC_NUM=${DBGROUP_TEST_EXEC_NUM}
    INDEX_FIXTURE_THREAD_NUM=${DBGROUP_TEST_THREAD_NUM}
    INDEX_FIXTURE_RANDOM_SEED=${DBGROUP_TEST_RANDOM_SEED}
    INDEX_FIXTURE_EXEC_NUM_PER_THREAD=${DBGROUP_TEST_EXEC_NUM}
  )
  target_compile_options(${DBGROUP_TEST_TARGET} PRIVATE
    -Wall
    -Wextra
    $<$<STREQUAL:${CMAKE_BUILD_TYPE},"Release">:"-O2 -march=native">
    $<$<STREQUAL:${CMAKE_BUILD_TYPE},"RelWithDebInfo">:"-g3 -Og -pg">
    $<$<STREQUAL:${CMAKE_BUILD_TYPE},"Debug">:"-g3 -O0 -pg">
  )
  target_link_libraries(${DBGROUP_TEST_TARGET} PRIVATE
    dbgroup::b_tree
    GTest::gtest_main
    $<$<BOOL:${DBGROUP_TEST_OVERRIDE_MIMALLOC}>:mimalloc-static>
  )
  add_test(NAME ${DBGROUP_TEST_TARGET} COMMAND $<TARGET_FILE:${DBGROUP_TEST_TARGET}>)
endfunction()

# add unit tests to build targets
ADD_DBGROUP_TEST("metadata_test")
ADD_DBGROUP_TEST("node_pml_test")
ADD_DBGROUP_TEST("node_psl_test")
ADD_DBGROUP_TEST("node_oml_test")
ADD_DBGROUP_TEST("node_osl_test")
ADD_DBGROUP_TEST("b_tree_pml_varlen_test")
ADD_DBGROUP_TEST("b_tree_pml_fixlen_test")
ADD_DBGROUP_TEST("b_tree_psl_varlen_test")
ADD_DBGROUP_TEST("b_tree_psl_fixlen_test")
ADD_DBGROUP_TEST("b_tree_oml_varlen_test")
ADD_DBGROUP_TEST("b_tree_oml_fixlen_test")
ADD_DBGROUP_TEST("b_tree_osl_varlen_test")
ADD_DBGROUP_TEST("b_tree_osl_fixlen_test")
if(${DBGROUP_TEST_BUILD_MULTI_THREAD_TESTS})
  ADD_DBGROUP_TEST("b_tree_pml_varlen_multi_thread_test")
  ADD_DBGROUP_TEST("b_tree_pml_fixlen_multi_thread_test")
  ADD_DBGROUP_TEST("b_tree_psl_varlen_multi_thread_test")
  ADD_DBGROUP_TEST("b_tree_psl_fixlen_multi_thread_test")
  ADD_DBGROUP_TEST("b_tree_oml_varlen_multi_thread_test")
  ADD_DBGROUP_TEST("b_tree_oml_fixlen_multi_thread_test")
  ADD_DBGROUP_TEST("b_tree_osl_varlen_multi_thread_test")
  ADD_DBGROUP_TEST("b_tree_osl_fixlen_multi_thread_test")
endif()
