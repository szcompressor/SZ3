cmake_minimum_required (VERSION 3.13)
project(SZ3)
#project(First VERSION 3.0.0.0)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "-O3")
include(GNUInstallDirs)

include(CheckCXXCompilerFlag)
include(CTest)
#by default pass no 3rd party exports
set(thirdparty_export "")

find_package(PkgConfig)


pkg_search_module(ZSTD IMPORTED_TARGET libzstd)
if(ZSTD_FOUND)
    set(ZSTD_dep PkgConfig::ZSTD)
else()
    add_subdirectory(zstd)
    set(ZSTD_dep zstd)
    list(APPEND thirdparty_export "zstd")
endif()

add_library(${PROJECT_NAME} INTERFACE)
target_include_directories(
        ${PROJECT_NAME} INTERFACE
        $<BUILD_INTERFACE:${${PROJECT_NAME}_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_FULL_INCLUDEDIR}>
        PkgConfig::ZSTD)


install(DIRECTORY ${${PROJECT_NAME}_SOURCE_DIR}/include/ DESTINATION ${CMAKE_INSTALL_FULL_INCLUDEDIR})

if(BUILD_TESTING)
  add_subdirectory(test)
endif()
#add_subdirectory(exaalt_cache)

