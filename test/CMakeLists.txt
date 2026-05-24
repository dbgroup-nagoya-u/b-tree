#------------------------------------------------------------------------------#
# Configuration
#------------------------------------------------------------------------------#

include(FetchContent)

# disable tests for not implemented APIs
set(DBGROUP_TEST_DISABLE_SCAN_VERIFIER_TEST ON CACHE BOOL "" FORCE)

FetchContent_Declare(
  index_fixtures
  GIT_REPOSITORY "https://github.com/dbgroup-nagoya-u/index-fixtures.git"
  GIT_TAG "c090e574caa034c65e8b4e53374501b61a7c0066"
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
  target_compile_options(${DBGROUP_TEST_TARGET} PRIVATE
    $<$<STREQUAL:"${CMAKE_BUILD_TYPE}","Release">:-march=native>
    $<$<STREQUAL:"${CMAKE_BUILD_TYPE}","Debug">:-g3>
    -Wall
    -Wextra
  )
  target_link_libraries(${DBGROUP_TEST_TARGET} PRIVATE
    dbgroup::${PROJECT_NAME}
    dbgroup::index_fixtures
  )
  add_test(NAME ${DBGROUP_TEST_TARGET} COMMAND $<TARGET_FILE:${DBGROUP_TEST_TARGET}>)
endfunction()

# add unit tests to build targets
ADD_DBGROUP_TEST("single_thread_test")
ADD_DBGROUP_TEST("multi_thread_test")
