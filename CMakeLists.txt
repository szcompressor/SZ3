cmake_minimum_required (VERSION 3.0)
project(SZ)
set (SZ_VERSION_MAJOR 3)
set (SZ_VERSION_MINOR 0)
set (CMAKE_CXX_STANDARD 11)
set (CMAKE_CXX_STANDARD_REQUIRED ON)
set (CMAKE_CXX_FLAGS "-g -O3")
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)

add_subdirectory (src)
add_subdirectory (test)

