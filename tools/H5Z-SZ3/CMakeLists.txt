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
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/hdf5_sz3)
install(EXPORT HDF5SZ3 NAMESPACE SZ3:: DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/SZ3)
export(TARGETS hdf5sz3 FILE HDF5SZ3.cmake NAMESPACE SZ3::)
