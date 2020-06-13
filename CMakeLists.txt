cmake_minimum_required(VERSION 3.0)
project(SZ)
#project(First VERSION 3.0.0.0)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "-O3")
add_definitions(-DPOLY_REGRESSION_COEF_FOLDER="${CMAKE_CURRENT_SOURCE_DIR}/data/")

find_library(ZSTD_LIBS zstd HINTS "$ENV{HOME}/tools/zstd-1.4.4/install/lib")
set(ZSTD_INCLUDES "$ENV{HOME}/tools/zstd-1.4.4/install/include")

#download modified version of lzbench by kai zhao from https://github.com/ayzk/lzbench
find_library(lzbench_LIBS lzbench HINTS "$ENV{HOME}/code/lzbench/")
set(lzbench_INCLUDES "$ENV{HOME}/code/lzbench/")
if (lzbench_LIBS)
    add_definitions(-DENABLE_LZBENCH="TRUE")
endif ()

add_subdirectory(src)
add_subdirectory(test)

