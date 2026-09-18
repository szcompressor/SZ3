cmake_minimum_required(VERSION 3.18)
project(SZ3 VERSION 3.3.3)

#data version defines the version of the compressed data format
#it is not always equal to the program version (e.g., SZ3 v3.1.0 and SZ3 v.3.1.1 may use the same data version of v.3.1.0)
#only update data version if the new version of the program changes compressed data format
set(SZ3_DATA_VERSION 3.3.2)

include(GNUInstallDirs)


option(BUILD_SHARED_LIBS "build shared libraries by default" ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
if (MSVC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /bigobj")
else () #-Wunused-parameter
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wreorder -Wsign-compare -Wsuggest-override -Wold-style-cast -Wcast-qual -Wextra-semi -Wnon-virtual-dtor -Wcast-align -Wunused-variable")
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

# ----------------------------------------------------------------------------------------------
# Zstd.
#
# SZ3 needs the Zstd *library* and, since Lossless_zstd.hpp declares the four Simple-API functions
# it calls, no Zstd *header* at all. Only the link has to be arranged, here and again in
# SZ3Config.cmake on the consumer's machine.
#
# Three ways of finding a system Zstd, tried in order, because no one of them is everywhere:
#   1. find_package(zstd CONFIG) -- the official CMake package. Absent entirely on RHEL/Rocky/Alma
#      8-10, Fedora <=43, Ubuntu 20.04 and 22.04, Debian bullseye, and default Spack.
#   2. pkg-config -- needs libzstd.pc, which on Debian and Ubuntu ships in libzstd-dev (Priority
#      optional, not installed by default) rather than in the libzstd1 every system already has.
#   3. find_library -- last resort, and the only one that needs neither a .cmake nor a .pc file.
# Then the vendored copy. SZ3_USE_BUNDLED_ZSTD=ON skips straight to it.
# ----------------------------------------------------------------------------------------------
if (MSVC)
    option(SZ3_USE_BUNDLED_ZSTD "prefer the bundled version of Zstd" ON)
else ()
    option(SZ3_USE_BUNDLED_ZSTD "prefer the bundled version of Zstd" OFF)
endif ()
set(SZ3_ZSTD_LIBRARY "" CACHE STRING
        "Zstd library to link: an imported target name or a full path. Empty searches for one.")

set(_sz3_zstd_link "")        # what SZ3 links in THIS build; never goes into the install export
set(SZ3_ZSTD_HEADER_DIRS "")    # where a real zstd.h lives, for the declaration-agreement test only
set(_sz3_zstd_libdirs "")     # recorded in SZ3Config.cmake as a search HINT, never as a link path

if (NOT SZ3_USE_BUNDLED_ZSTD)
    if (SZ3_ZSTD_LIBRARY)
        set(_sz3_zstd_link "${SZ3_ZSTD_LIBRARY}")
    endif ()
    if (NOT _sz3_zstd_link)
        find_package(zstd CONFIG QUIET)
        # Shared first: where both exist zstd::libzstd follows BUILD_SHARED_LIBS, so it would
        # quietly pull a static Zstd into a shared SZ3. zstd::libzstd itself only arrived in Zstd
        # 1.5.6, and a single-linkage build (Arch, conda-forge, one vcpkg triplet) has exactly one
        # of _shared/_static, so all three names have to be tried.
        foreach (_t zstd::libzstd_shared zstd::libzstd zstd::libzstd_static)
            if (TARGET ${_t} AND NOT _sz3_zstd_link)
                set(_sz3_zstd_link ${_t})
                get_target_property(SZ3_ZSTD_HEADER_DIRS ${_t} INTERFACE_INCLUDE_DIRECTORIES)
                if (NOT SZ3_ZSTD_HEADER_DIRS)
                    set(SZ3_ZSTD_HEADER_DIRS "")
                endif ()
            endif ()
        endforeach ()
    endif ()
    if (NOT _sz3_zstd_link AND PKG_CONFIG_FOUND)
        pkg_search_module(ZSTD QUIET IMPORTED_TARGET libzstd)
        if (ZSTD_FOUND)
            # The resolved libraries, not PkgConfig::ZSTD: an imported target exists only inside
            # the build that created it, so naming one in the install export made the consumer's
            # find_package(SZ3) die at generate time with "the link interface of target SZ3::SZ3
            # contains: PkgConfig::ZSTD".
            set(_sz3_zstd_link "${ZSTD_LINK_LIBRARIES}")
            set(SZ3_ZSTD_HEADER_DIRS "${ZSTD_INCLUDE_DIRS}")
            set(_sz3_zstd_libdirs "${ZSTD_LIBRARY_DIRS}")
        endif ()
    endif ()
    if (NOT _sz3_zstd_link)
        find_library(SZ3_ZSTD_LIBRARY_PATH NAMES zstd libzstd zstd_static)
        if (SZ3_ZSTD_LIBRARY_PATH)
            set(_sz3_zstd_link "${SZ3_ZSTD_LIBRARY_PATH}")
            get_filename_component(_sz3_zstd_libdirs "${SZ3_ZSTD_LIBRARY_PATH}" DIRECTORY)
        endif ()
    endif ()
endif ()

if (_sz3_zstd_link)
    # Where that Zstd actually sits, for SZ3Config.cmake to offer find_library as a HINT. The
    # CONFIG route resolves to an imported target rather than a path, so the path has to be read
    # back off the target; without this the hint came out empty exactly when find_package(zstd)
    # had succeeded, which is the common case on a conda or Homebrew prefix.
    if (NOT _sz3_zstd_libdirs AND TARGET ${_sz3_zstd_link})
        # The configuration names to try are the target's own: Debian and Ubuntu build their zstd
        # package with IMPORTED_CONFIGURATIONS=NONE, so the location is IMPORTED_LOCATION_NONE and
        # guessing RELEASE finds nothing. Ask, then fall back to the usual spellings.
        string(TOUPPER "${CMAKE_BUILD_TYPE}" _sz3_cfg)
        get_target_property(_sz3_cfgs ${_sz3_zstd_link} IMPORTED_CONFIGURATIONS)
        if (NOT _sz3_cfgs)
            set(_sz3_cfgs "")
        endif ()
        foreach (_c IN LISTS _sz3_cfgs ITEMS ${_sz3_cfg} RELEASE NOCONFIG NONE "")
            if (_c)
                set(_prop IMPORTED_LOCATION_${_c})
            else ()
                set(_prop IMPORTED_LOCATION)
            endif ()
            get_target_property(_sz3_loc ${_sz3_zstd_link} ${_prop})
            if (_sz3_loc AND EXISTS "${_sz3_loc}" AND NOT _sz3_zstd_libdirs)
                get_filename_component(_sz3_zstd_libdirs "${_sz3_loc}" DIRECTORY)
            endif ()
        endforeach ()
    endif ()
    set(SZ3_ZSTD_IS_EXTERNAL ON)
    set(SZ3_ZSTD_PROVIDER "system")
    # BUILD_INTERFACE for the link as well as for any include directory. Which Zstd this is and
    # where it lives are facts about the machine that built SZ3, and neither survives being written
    # into the install export: the resolved path left a consumer linking a libzstd that was no
    # longer there ("No rule to make target .../libzstd.so"), and an imported target name is not
    # even defined in their project. SZ3Config.cmake resolves Zstd again on their machine instead.
    target_link_libraries(${PROJECT_NAME} INTERFACE $<BUILD_INTERFACE:${_sz3_zstd_link}>)
    # No include directory is exported even under BUILD_INTERFACE beyond what the test needs:
    # nothing SZ3 installs includes zstd.h. Under conda, Spack or Homebrew the directory Zstd was
    # found in is a whole prefix holding other packages' headers, so exporting it put this
    # machine's hdf5.h ahead of the one the consumer found.
else ()
    set(SZ3_ZSTD_IS_EXTERNAL OFF)
    set(SZ3_ZSTD_PROVIDER "vendored")
    if (NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/tools/zstd/CMakeLists.txt")
        message(FATAL_ERROR
                "No system Zstd was found and tools/zstd/ is not present, so there is nothing to "
                "fall back to. Install Zstd (libzstd-dev on Debian/Ubuntu, libzstd-devel on "
                "Fedora/RHEL, zstd on Homebrew), or point SZ3_ZSTD_LIBRARY at the library.")
    endif ()
    add_subdirectory(tools/zstd)
    # Exported as SZ3::zstd with an IMPORTED_LOCATION under ${_IMPORT_PREFIX}, so unlike a system
    # Zstd this one relocates with the install prefix and needs no resolver on the consumer's side.
    target_link_libraries(${PROJECT_NAME} INTERFACE zstd)
    set(SZ3_ZSTD_HEADER_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/tools/zstd/lib")
endif ()
# The declaration-agreement test needs a real zstd.h. A system Zstd found by find_library reports
# no include directory at all -- that is the point of that fallback -- so fall back to the vendored
# header, which checks the same claim: that these four signatures do not vary.
if (NOT SZ3_ZSTD_HEADER_DIRS AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/tools/zstd/lib/zstd.h")
    set(SZ3_ZSTD_HEADER_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/tools/zstd/lib")
endif ()

message(STATUS "SZ3: Zstd provider = ${SZ3_ZSTD_PROVIDER}")
if (SZ3_ZSTD_IS_EXTERNAL)
    message(STATUS "SZ3: Zstd link = ${_sz3_zstd_link}")
endif ()

option(BUILD_TESTING "Build tests" OFF)
if (BUILD_TESTING)
    enable_testing()
    add_subdirectory(tools/test)
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

option(BUILD_PARAVIEW_PLUGIN "build the ParaView plugin" OFF)
if (BUILD_PARAVIEW_PLUGIN)
    add_subdirectory(tools/paraview)
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
