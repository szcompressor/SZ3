#Either use:
#find_package(HDF5 1.12.1 REQUIRED)
#or manually set the hdf5 paths:
#set(HDF5_LIBRARIES "/path/to/hdf5-1.12.1/hdf5/lib/libhdf5.so")
#set(HDF5_INCLUDE_DIRS "/path/to/hdf5-1.12.1/hdf5/include")

find_package(HDF5 COMPONENTS C REQUIRED)

#message("HDF5LIB PATH:" ${HDF5_LIBRARIES})
#message("HDF5INCLUDE PATH:" ${HDF5_INCLUDE_DIRS})

add_library(hdf5sz3
        src/H5Z_SZ3.cpp
)
add_library(SZ3::hdf5sz3 ALIAS hdf5sz3)
# Export only what HDF5SZ3_EXPORT marks: an exported SZ3 internal can be taken over by another SZ3
# copy loaded in the same process.
set_target_properties(hdf5sz3 PROPERTIES CXX_VISIBILITY_PRESET hidden VISIBILITY_INLINES_HIDDEN ON POSITION_INDEPENDENT_CODE ON)
if (NOT BUILD_SHARED_LIBS)
    # Otherwise H5Z_SZ3.hpp marks the functions of a static hdf5sz3 dllimport on Windows.
    target_compile_definitions(hdf5sz3 PUBLIC HDF5SZ3_STATIC)
endif ()

# SZ3core only: SZ3::SZ3 has OpenMP only where the consumer finds it, and the filter must link everywhere.
target_link_libraries(hdf5sz3
        PUBLIC SZ3core HDF5::HDF5
)
# Do not add an MPI search: the filter calls no MPI function, and a parallel HDF5 brings its own.

# Reads the data SZ3 3.1.x's filter wrote (hdf5plugin 7.1.0 and earlier ship it) with the code of SZ3 3.1.7, which is
# fetched here; H5Z_SZ3_v3_1.cpp is the only file that sees it.
option(H5Z_SZ3_READ_V3_1 "Read HDF5 data written by the SZ3 3.1.x filter (fetches SZ3 3.1.7)" OFF)
if (H5Z_SZ3_READ_V3_1)
    include(FetchContent)
    FetchContent_Declare(sz3_v3_1
            URL https://github.com/szcompressor/SZ3/archive/refs/tags/v3.1.7.tar.gz
            URL_HASH SHA256=027dea310793ff1e80e595331fa56c8fd1ac9ecfe40f45e4be7d3a3e7ba7d390
            SOURCE_SUBDIR no-cmake-project)  # headers only: do not build or install SZ3 3.1.7
    FetchContent_MakeAvailable(sz3_v3_1)
    function(configure_sz3_v3_1_version)  # 3.1.7's version.hpp, in a function so its PROJECT_* stay local
        set(PROJECT_NAME SZ3)
        set(PROJECT_VERSION 3.1.7)
        set(PROJECT_VERSION_MAJOR 3)
        set(PROJECT_VERSION_MINOR 1)
        set(PROJECT_VERSION_PATCH 7)
        set(PROJECT_VERSION_TWEAK 0)
        configure_file(${sz3_v3_1_SOURCE_DIR}/include/SZ3/version.hpp.in ${CMAKE_CURRENT_BINARY_DIR}/sz3_v3_1/SZ3/version.hpp)
    endfunction()
    configure_sz3_v3_1_version()
    add_library(hdf5sz3_v3_1 OBJECT src/H5Z_SZ3_v3_1.cpp)
    set_target_properties(hdf5sz3_v3_1 PROPERTIES CXX_VISIBILITY_PRESET hidden VISIBILITY_INLINES_HIDDEN ON
            POSITION_INDEPENDENT_CODE ON)
    target_compile_features(hdf5sz3_v3_1 PRIVATE cxx_std_17)
    # SYSTEM: 3.1.7 is fetched as it was and does not meet the warnings SZ3 builds with
    target_include_directories(hdf5sz3_v3_1 SYSTEM PRIVATE ${sz3_v3_1_SOURCE_DIR}/include ${CMAKE_CURRENT_BINARY_DIR}/sz3_v3_1)
    # Zstd's headers only, not SZ3core's: they would shadow 3.1.7's
    if (TARGET sz3_zstd)
        target_include_directories(hdf5sz3_v3_1 PRIVATE $<TARGET_PROPERTY:sz3_zstd,INTERFACE_INCLUDE_DIRECTORIES>)
    else ()
        target_include_directories(hdf5sz3_v3_1 PRIVATE ${ZSTD_INCLUDE_DIRS})
    endif ()
    target_sources(hdf5sz3 PRIVATE $<TARGET_OBJECTS:hdf5sz3_v3_1>)
    target_compile_definitions(hdf5sz3 PRIVATE H5Z_SZ3_READ_V3_1)
endif ()

# Keep this: without a RUNPATH the installed library cannot find a conda, Spack or module HDF5.
set_target_properties(hdf5sz3 PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE)
# SZ3 changes these functions' signatures and SZ3::Config's layout only in a new major.minor.
set_target_properties(hdf5sz3 PROPERTIES VERSION ${SZ3_VERSION} SOVERSION ${SZ3_VERSION_MAJOR}.${SZ3_VERSION_MINOR})

target_include_directories(hdf5sz3
        PUBLIC
        # BUILD_INTERFACE only: exported, this machine's HDF5 headers would outrank a consumer's.
        # Quoted, or CMake reads the leading $ as a relative path and prepends the source directory.
        "$<BUILD_INTERFACE:${HDF5_INCLUDE_DIRS}>"
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
        "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/hdf5_sz3>"
)

if (BUILD_SZ3_BINARY)
    add_subdirectory(tools)
endif ()
# Only at top level, as SZ3's own tests: a parent's BUILD_TESTING is for its own tests.
if (BUILD_TESTING AND SZ3_IS_TOP_LEVEL)
    add_subdirectory(test)
endif ()

if (SZ3_INSTALL)
    install(TARGETS hdf5sz3 EXPORT HDF5SZ3
            ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
            LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
            RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    )

    # A copy for HDF5 to find on its own. Its own directory, because HDF5 dlopens every lib*.so on
    # HDF5_PLUGIN_PATH.
    set(H5Z_SZ3_PLUGIN_INSTALL_DIR "${CMAKE_INSTALL_LIBDIR}/plugin" CACHE STRING
            "Directory to install the HDF5 filter plugin copy into; empty installs no copy")
    if (BUILD_SHARED_LIBS AND H5Z_SZ3_PLUGIN_INSTALL_DIR)
        install(TARGETS hdf5sz3
                LIBRARY DESTINATION ${H5Z_SZ3_PLUGIN_INSTALL_DIR}
                RUNTIME DESTINATION ${H5Z_SZ3_PLUGIN_INSTALL_DIR}
                # The Windows import library belongs beside the real one, not where HDF5 scans.
                ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        )
    endif ()

    install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/hdf5_sz3)
    install(EXPORT HDF5SZ3 NAMESPACE SZ3:: DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/SZ3)
endif ()
export(TARGETS hdf5sz3 FILE HDF5SZ3.cmake NAMESPACE SZ3::)
