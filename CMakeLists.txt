cmake_minimum_required(VERSION 3.0)
project(SZ3)
#project(First VERSION 3.0.0.0)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "-O3")
include(GNUInstallDirs)

include(CheckCXXCompilerFlag)
include(CTest)

find_package(PkgConfig)
pkg_search_module(ZSTD IMPORTED_TARGET libzstd)

add_library(${PROJECT_NAME} INTERFACE)
target_include_directories(
        ${PROJECT_NAME} INTERFACE
        $<BUILD_INTERFACE:${${PROJECT_NAME}_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_FULL_INCLUDEDIR}>
        PkgConfig::ZSTD)


install(DIRECTORY ${${PROJECT_NAME}_SOURCE_DIR}/include/ DESTINATION ${CMAKE_INSTALL_FULL_INCLUDEDIR})

add_subdirectory(test)
add_subdirectory(exaalt)

