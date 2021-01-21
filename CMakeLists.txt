cmake_minimum_required(VERSION 3.0)
project(SZ)
#project(First VERSION 3.0.0.0)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "-O3")

find_library(ZSTD_LIBS zstd HINTS "$ENV{HOME}/tools/zstd-1.4.4/install/lib")
set(ZSTD_INCLUDES "$ENV{HOME}/tools/zstd-1.4.4/install/include")

add_subdirectory(src)
add_subdirectory(test)

