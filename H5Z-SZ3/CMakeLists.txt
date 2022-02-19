find_package(HDF5 REQUIRED)
add_library(
        hdf5sz3
        src/H5Z_SZ3.cpp
)
target_link_libraries(
        hdf5sz3
        PUBLIC SZ3 ${HDF5_LIBRARIES}
)
target_include_directories(
        hdf5sz3
        PUBLIC
        ${HDF5_INCLUDE_DIRS}
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:hdf5_sz3>
)

add_subdirectory(test)

install(TARGETS hdf5sz3 EXPORT HDF5SZ3
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        )
install(DIRECTORY include DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/hdf5_sz3)
export(TARGETS hdf5sz3 FILE HDF5SZ3.cmake)
