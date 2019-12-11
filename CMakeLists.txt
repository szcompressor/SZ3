cmake_minimum_required(VERSION 3.0)
project(SZ)
project(First VERSION 3.0.0.0)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_library(ZSTD_LIBS zstd HINTS "/Users/kzhao/Documents/Workspace/zstd-1.4.2/install/lib")
set(ZSTD_INCLUDES "/Users/kzhao/Documents/Workspace/zstd-1.4.2/install/include")

add_subdirectory(src)
add_subdirectory(test)

