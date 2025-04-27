#------------------------------------------------------------------------------#
# Configuration
#------------------------------------------------------------------------------#

include(FetchContent)

FetchContent_Declare(
  index_fixtures
  GIT_REPOSITORY "https://github.com/dbgroup-nagoya-u/index-fixtures.git"
  GIT_TAG "8efe82507c8189d1d801bd00c027c33447126ec5"
)
FetchContent_MakeAvailable(index_fixtures)

#------------------------------------------------------------------------------#
# Build Unit Tests
#------------------------------------------------------------------------------#

# define function to add unit tests in the same format
function(ADD_DBGROUP_TEST DBGROUP_TEST_TARGET)
  add_executable(${DBGROUP_TEST_TARGET}
    "${CMAKE_CURRENT_SOURCE_DIR}/${DBGROUP_TEST_TARGET}.cpp"
  )
  target_compile_features(${DBGROUP_TEST_TARGET} PRIVATE
    "cxx_std_20"
  )
  target_compile_definitions(${DBGROUP_TEST_TARGET} PRIVATE
    DBGROUP_TEST_THREAD_NUM=${DBGROUP_TEST_THREAD_NUM}
    DBGROUP_TEST_RANDOM_SEED=${DBGROUP_TEST_RANDOM_SEED}
    DBGROUP_TEST_EXEC_NUM=${DBGROUP_TEST_EXEC_NUM}
  )
  target_compile_options(${DBGROUP_TEST_TARGET} PRIVATE
    $<$<STREQUAL:"${CMAKE_BUILD_TYPE}","Release">:-march=native>
    $<$<STREQUAL:"${CMAKE_BUILD_TYPE}","Debug">:-g3>
    -Wall
    -Wextra
  )
  target_link_libraries(${DBGROUP_TEST_TARGET} PRIVATE
    dbgroup::b_tree
    dbgroup::index_fixtures
  )
  add_test(NAME ${DBGROUP_TEST_TARGET} COMMAND $<TARGET_FILE:${DBGROUP_TEST_TARGET}>)
endfunction()

# add unit tests to build targets
ADD_DBGROUP_TEST("b_tree_psl_test")
ADD_DBGROUP_TEST("b_tree_osl_test")
ADD_DBGROUP_TEST("b_tree_psl_multi_thread_test")
ADD_DBGROUP_TEST("b_tree_osl_multi_thread_test")
