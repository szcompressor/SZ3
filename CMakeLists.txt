cmake_minimum_required(VERSION 3.18)
project(SZ3 VERSION 3.3.1)

#data version defines the version of the compressed data format
#it is not always equal to the program version (e.g., SZ3 v3.1.0 and SZ3 v.3.1.1 may use the same data version of v.3.1.0)
#only update data version if the new version of the program changes compressed data format
set(SZ3_DATA_VERSION 3.3.0)

include(GNUInstallDirs)


option(BUILD_SHARED_LIBS "build shared libraries by default" ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
if (MSVC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /bigobj")
else () #-Wunused-parameter
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wsign-compare -Wsuggest-override -Wold-style-cast -Wcast-qual -Wzero-as-null-pointer-constant -Wextra-semi -Wnon-virtual-dtor -Wcast-align -Wunused-variable")
endif ()
set(default_build_type "Release")
if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS "No build type was set. Setting build type to ${default_build_type}.")
    set(CMAKE_BUILD_TYPE ${default_build_type} CACHE
            STRING "Choose the type to build" FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release"
            "MinSizeRel" "RelWithDebInfo")
endif ()

find_package(PkgConfig)

configure_file(
        ${CMAKE_CURRENT_SOURCE_DIR}/include/SZ3/version.hpp.in
        ${CMAKE_CURRENT_BINARY_DIR}/include/SZ3/version.hpp)

add_library(${PROJECT_NAME} INTERFACE)
target_include_directories(
        ${PROJECT_NAME} INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
target_compile_features(${PROJECT_NAME}
        INTERFACE cxx_std_17
)

find_package(OpenMP)
if (OpenMP_FOUND)
    target_link_libraries(${PROJECT_NAME} INTERFACE OpenMP::OpenMP_CXX)
endif ()

if (MSVC)
    option(SZ3_USE_BUNDLED_ZSTD "prefer the bundled version of Zstd" ON)
else ()
    option(SZ3_USE_BUNDLED_ZSTD "prefer the bundled version of Zstd" OFF)
endif ()

option(SZ3_DEBUG_TIMINGS "print debug timing information" OFF)
if (SZ3_DEBUG_TIMINGS)
    target_compile_definitions(${PROJECT_NAME} INTERFACE SZ3_DEBUG_TIMINGS=1)
else ()
    target_compile_definitions(${PROJECT_NAME} INTERFACE SZ3_DEBUG_TIMINGS=0)
endif ()

option(SZ3_USE_SKA_HASH "use ska hashmap instead of std hashmap, default is ska for better speed" ON)
if (SZ3_USE_SKA_HASH)
    target_compile_definitions(${PROJECT_NAME} INTERFACE SZ3_USE_SKA_HASH=1)
else ()
    target_compile_definitions(${PROJECT_NAME} INTERFACE SZ3_USE_SKA_HASH=0)
endif ()

pkg_search_module(ZSTD IMPORTED_TARGET libzstd)
if (ZSTD_FOUND AND NOT SZ3_USE_BUNDLED_ZSTD)
    target_link_libraries(${PROJECT_NAME} INTERFACE PkgConfig::ZSTD)
else ()
    add_subdirectory(tools/zstd)
    target_link_libraries(${PROJECT_NAME} INTERFACE zstd)
endif ()

option(BUILD_TESTING "Build tests" OFF)
if (BUILD_TESTING)
    enable_testing()
    add_subdirectory(tools/test)
endif ()

find_package(GSL)
if (GSL_FOUND)
    target_compile_definitions(${PROJECT_NAME} INTERFACE -DSZ3_ENABLE_GSL)
    target_link_libraries(${PROJECT_NAME} INTERFACE GSL::gsl)
endif ()

option(BUILD_SZ3_BINARY "build SZ3 tools" ON)
if (BUILD_SZ3_BINARY)
    add_subdirectory(tools/sz3)
    add_subdirectory(tools/sz3c)
endif ()

option(BUILD_H5Z_FILTER "build the H5Z-SZ3 filter" OFF)
if (BUILD_H5Z_FILTER)
    add_subdirectory(tools/H5Z-SZ3)
endif ()

option(BUILD_MDZ "build MDZ (MMD-SZ)" OFF)
if (BUILD_MDZ)
    add_subdirectory(tools/mdz)
endif ()


install(TARGETS ${PROJECT_NAME}
        EXPORT SZ3Targets
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(FILES
        "${CMAKE_CURRENT_BINARY_DIR}/include/SZ3/version.hpp"
        DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/SZ3/"
)
install(EXPORT SZ3Targets NAMESPACE SZ3:: DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/SZ3)
include(CMakePackageConfigHelpers)
configure_package_config_file(${CMAKE_CURRENT_SOURCE_DIR}/SZ3Config.cmake.in
        "${CMAKE_CURRENT_BINARY_DIR}/SZ3Config.cmake"
        INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/SZ3
)
write_basic_package_version_file(
        "${CMAKE_CURRENT_BINARY_DIR}/SZ3ConfigVersion.cmake"
        VERSION ${PROJECT_VERSION}
        COMPATIBILITY AnyNewerVersion
)
install(FILES
        ${CMAKE_CURRENT_BINARY_DIR}/SZ3Config.cmake
        ${CMAKE_CURRENT_BINARY_DIR}/SZ3ConfigVersion.cmake
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/SZ3
)
#export sz3 target for external use
export(TARGETS SZ3 FILE SZ3.cmake)
