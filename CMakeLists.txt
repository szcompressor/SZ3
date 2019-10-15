cmake_minimum_required (VERSION 3.0)
project(SZ)
project(First VERSION 3.0.0.0)
set (CMAKE_CXX_STANDARD 11)
set (CMAKE_CXX_STANDARD_REQUIRED ON)
set (CMAKE_CXX_FLAGS "-g -O3")

add_subdirectory (src)
add_subdirectory (test)

