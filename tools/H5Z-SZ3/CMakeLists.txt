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
set_target_properties(hdf5sz3 PROPERTIES CXX_VISIBILITY_PRESET hidden VISIBILITY_INLINES_HIDDEN ON)

# The filter is built with SZ3's OpenMP, and its consumers are not: they get the headers without it.
# A static filter's need for the OpenMP runtime reaches them as link-only, through SZ3 below.
target_link_libraries(hdf5sz3
        PUBLIC SZ3_no_openmp HDF5::HDF5
        PRIVATE SZ3
)
# Do not add an MPI search: the filter calls no MPI function, and a parallel HDF5 brings its own.

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

add_subdirectory(tools)
if (BUILD_TESTING)
    add_subdirectory(test)
endif ()

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
export(TARGETS hdf5sz3 FILE HDF5SZ3.cmake NAMESPACE SZ3::)
