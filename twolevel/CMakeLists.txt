cmake_minimum_required(VERSION 3.10)
project(twolevel CXX C)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
set(CMAKE_C_FLAGS_RELEASE "-O3 -DNDEBUG")
add_compile_options(-march=native)

# The reportable-shortcuts source, compiled here (in this build directory) with
# the same flags; nothing in its own folder is changed.
set(RS_DIR "$ENV{HOME}/shortcuts/reportable-shortcuts" CACHE PATH "reportable-shortcuts source")
add_subdirectory(${RS_DIR} rs_build EXCLUDE_FROM_ALL)
# (its build_index is used by the self-test of run_twolevel.sh)

add_executable(build_twolevel build_twolevel.cpp)
target_link_libraries(build_twolevel PRIVATE shortcuts)
add_executable(query_twolevel query_twolevel.cpp)
target_link_libraries(query_twolevel PRIVATE shortcuts)
