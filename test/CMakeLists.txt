#--------------------------------------------------------------------------------------#
# Configure GoogleTest
#--------------------------------------------------------------------------------------#

configure_file(CMakeLists-googletest.txt.in googletest-download/CMakeLists.txt)
execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
  RESULT_VARIABLE result
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download )
if(result)
  message(FATAL_ERROR "CMake step for googletest failed: ${result}")
endif()
execute_process(COMMAND ${CMAKE_COMMAND} --build .
  RESULT_VARIABLE result
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download )
if(result)
  message(FATAL_ERROR "Build step for googletest failed: ${result}")
endif()

# Prevent overriding the parent project's compiler/linker
# settings on Windows
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

# Add googletest directly to our build. This defines
# the gtest and gtest_main targets.
add_subdirectory(${CMAKE_CURRENT_BINARY_DIR}/googletest-src
                 ${CMAKE_CURRENT_BINARY_DIR}/googletest-build
                 EXCLUDE_FROM_ALL)

#--------------------------------------------------------------------------------------#
# Build Unit Tests
#--------------------------------------------------------------------------------------#

# define function to add unit tests in the same format

function(ADD_B_TREE_TEST B_TREE_TEST_TARGET)
  add_executable(${B_TREE_TEST_TARGET}
    "${CMAKE_CURRENT_SOURCE_DIR}/${B_TREE_TEST_TARGET}.cpp"
  )
  target_compile_features(${B_TREE_TEST_TARGET} PRIVATE
    "cxx_std_17"
  )
  target_compile_options(${B_TREE_TEST_TARGET} PRIVATE
    -Wall
    -Wextra
    $<$<STREQUAL:${CMAKE_BUILD_TYPE},"Release">:"-O2 -march=native">
    $<$<STREQUAL:${CMAKE_BUILD_TYPE},"RelWithDebInfo">:"-g3 -Og -pg">
    $<$<STREQUAL:${CMAKE_BUILD_TYPE},"Debug">:"-g3 -O0 -pg">
  )
  target_link_libraries(${B_TREE_TEST_TARGET} PRIVATE
    b_tree
    gtest_main
  )
  add_test(NAME ${B_TREE_TEST_TARGET} COMMAND $<TARGET_FILE:${B_TREE_TEST_TARGET}>)
endfunction()

# add unit tests to build targets
ADD_B_TREE_TEST("metadata_test")
ADD_B_TREE_TEST("node_test")
ADD_B_TREE_TEST("b_tree_test")
ADD_B_TREE_TEST("b_tree_pcl_multi_thread_test")
