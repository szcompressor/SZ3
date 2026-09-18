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

target_link_libraries(hdf5sz3
        PUBLIC SZ3 HDF5::HDF5
)
# No MPI search here. The filter calls no MPI function, and a parallel HDF5 brings the MPI it was
# built against with it. Searching separately put MPI::MPI_CXX in the installed export, which is
# what made SZ3Config resolve MPI with a REQUIRED find_package and end its consumer's configure.
# HDF5_IS_PARALLEL gated nothing else, so it no longer goes to the parent scope either.

# Without this the installed library records no RUNPATH and ldd reports libhdf5.so.* as not found
# wherever HDF5 lives off the loader's default path -- every conda, Spack or module-file HDF5.
set_target_properties(hdf5sz3 PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE)

target_include_directories(hdf5sz3
        PUBLIC
        # BUILD_INTERFACE only: unguarded, the install export baked in this machine's HDF5 include
        # path, and a consumer who found a different HDF5 compiled against those headers. An
        # installed consumer gets them from HDF5::HDF5 instead.
        $<BUILD_INTERFACE:${HDF5_INCLUDE_DIRS}>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/hdf5_sz3>
)

add_subdirectory(test)

install(TARGETS hdf5sz3 EXPORT HDF5SZ3
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

# A second copy of the same library, in a directory that holds nothing else, for HDF5 to find by
# itself. Without one the only way to reach the filter is HDF5_PLUGIN_PATH=<prefix>/lib, and that
# is both destructive and expensive: the variable REPLACES HDF5's compiled-in default rather than
# extending it (H5PLpath.c, H5PL__create_path_table), and HDF5 dlopens every lib*.so in each
# directory on it until one answers to the filter id (H5PLpath.c, H5PL__find_plugin_in_path).
# lib/plugin is what the HDF Group's own filters use; the cache variable is H5Z-ZFP's precedent,
# for a packager who installs the copy somewhere else or not at all. HDF5's own pkg-config
# PluginDir is deliberately not consulted: it names a path outside this prefix.
set(H5Z_SZ3_PLUGIN_INSTALL_DIR "${CMAKE_INSTALL_LIBDIR}/plugin" CACHE STRING
        "Directory to install the HDF5 filter plugin copy into; empty installs no copy")
if (BUILD_SHARED_LIBS AND H5Z_SZ3_PLUGIN_INSTALL_DIR)
    install(TARGETS hdf5sz3
            LIBRARY DESTINATION ${H5Z_SZ3_PLUGIN_INSTALL_DIR}
            RUNTIME DESTINATION ${H5Z_SZ3_PLUGIN_INSTALL_DIR}
            # Windows puts the DLL under RUNTIME and an import library under ARCHIVE. The import
            # library belongs beside the other one, not in a directory HDF5 scans.
            ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    )
endif ()

install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/hdf5_sz3)
install(EXPORT HDF5SZ3 NAMESPACE SZ3:: DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/SZ3)
export(TARGETS hdf5sz3 FILE HDF5SZ3.cmake NAMESPACE SZ3::)
