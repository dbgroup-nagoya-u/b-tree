#------------------------------------------------------------------------------#
# Configuration
#------------------------------------------------------------------------------#

include(FetchContent)

set(
  B_TREE_TEST_SANITIZER
  "none" CACHE STRING
  "Sanitizers for unit tests and their dependencies (none, thread, or address)."
)
set_property(CACHE B_TREE_TEST_SANITIZER PROPERTY STRINGS none thread address)

if(B_TREE_TEST_SANITIZER STREQUAL "thread")
  set(DBGROUP_SANITIZER_FLAGS -fsanitize=thread)
  if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    message(WARNING
      "[b_tree] GCC's ThreadSanitizer may abort with \"unexpected memory mapping\" "
      "under address space layout randomization. Use Clang 18 or later, or run "
      "tests via \"setarch $(uname -m) -R\".")
  endif()
elseif(B_TREE_TEST_SANITIZER STREQUAL "address")
  # make UBSan errors fail unit tests
  set(DBGROUP_SANITIZER_FLAGS
    -fsanitize=address,undefined
    -fno-sanitize-recover=undefined
    -fno-omit-frame-pointer
  )
elseif(NOT B_TREE_TEST_SANITIZER STREQUAL "none")
  message(FATAL_ERROR
    "[b_tree] B_TREE_TEST_SANITIZER must be none, thread, or address.")
endif()

if(DBGROUP_SANITIZER_FLAGS)
  if(DBGROUP_TEST_OVERRIDE_MIMALLOC)
    message(FATAL_ERROR
      "[b_tree] Sanitizers cannot be used with DBGROUP_TEST_OVERRIDE_MIMALLOC.")
  endif()

  # instrument test fixtures and libraries fetched below
  add_compile_options(${DBGROUP_SANITIZER_FLAGS})
  string(REPLACE ";" " " DBGROUP_SANITIZER_LINK_FLAGS "${DBGROUP_SANITIZER_FLAGS}")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${DBGROUP_SANITIZER_LINK_FLAGS}")

  # instrument our dependencies built in the parent directory
  foreach(DBGROUP_DEP_TARGET cpp_utility memory_manager)
    if(TARGET ${DBGROUP_DEP_TARGET})
      target_compile_options(${DBGROUP_DEP_TARGET} PRIVATE ${DBGROUP_SANITIZER_FLAGS})
    endif()
  endforeach()
endif()

# disable tests for not implemented APIs
set(DBGROUP_TEST_DISABLE_SCAN_VERIFIER_TEST ON CACHE BOOL "" FORCE)

FetchContent_Declare(
  index_fixtures
  GIT_REPOSITORY "https://github.com/dbgroup-nagoya-u/index-fixtures.git"
  GIT_TAG "36f51571d4177508e74d5b842d1365304a356bcb"
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
  if(B_TREE_TEST_SANITIZER STREQUAL "thread")
    # embed the suppression file into default ThreadSanitizer options
    target_sources(${DBGROUP_TEST_TARGET} PRIVATE
      "${CMAKE_CURRENT_SOURCE_DIR}/tsan_default_options.cpp"
    )
    target_compile_definitions(${DBGROUP_TEST_TARGET} PRIVATE
      B_TREE_TSAN_SUPPRESSIONS="${CMAKE_CURRENT_SOURCE_DIR}/tsan.supp"
    )
  endif()
  add_test(NAME ${DBGROUP_TEST_TARGET} COMMAND $<TARGET_FILE:${DBGROUP_TEST_TARGET}>)
endfunction()

# add unit tests to build targets
ADD_DBGROUP_TEST("single_thread_test")
ADD_DBGROUP_TEST("multi_thread_test")
