cmake_minimum_required(VERSION 3.0)
project(SZ3)
#project(First VERSION 3.0.0.0)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "-O3")
include(GNUInstallDirs)

find_library(ZSTD_LIBS zstd HINTS "$ENV{HOME}/tools/zstd-1.4.4/install/lib")
set(ZSTD_INCLUDES "$ENV{HOME}/tools/zstd-1.4.4/install/include")

add_library(${PROJECT_NAME} INTERFACE)
target_include_directories(
        ${PROJECT_NAME} INTERFACE
        $<BUILD_INTERFACE:${${PROJECT_NAME}_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_FULL_INCLUDEDIR}>
        ${ZSTD_INCLUDES})


install(DIRECTORY ${${PROJECT_NAME}_SOURCE_DIR}/include/  DESTINATION ${CMAKE_INSTALL_FULL_INCLUDEDIR})

add_subdirectory(test)
add_subdirectory(exaalt)

