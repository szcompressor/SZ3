function(build_hdf5_test)
    get_filename_component(test_name ${ARGV0} NAME_WE)
    add_executable(${test_name} ${ARGV})
    target_link_libraries(${test_name} PUBLIC SZ3 hdf5sz3)
endfunction(build_hdf5_test)

build_hdf5_test(sz3ToHDF5.cpp)
build_hdf5_test(dsz3FromHDF5.cpp)
build_hdf5_test(convertBinToHDF5.cpp)
build_hdf5_test(cdvalueHelper.cpp)
# The access-mode suite runs against an install tree, so it installs one of its own first. Only
# registered when the HDF5 command-line tools it drives are present -- CI installs them and runs
# the script directly, so coverage there does not depend on this.
find_program(H5REPACK_EXECUTABLE h5repack)
find_program(H5DUMP_EXECUTABLE h5dump)
find_program(H5LS_EXECUTABLE h5ls)
if (BUILD_TESTING AND H5REPACK_EXECUTABLE AND H5DUMP_EXECUTABLE AND H5LS_EXECUTABLE)
    get_filename_component(H5Z_SZ3_TOOL_DIR "${H5DUMP_EXECUTABLE}" DIRECTORY)
    add_test(NAME SZ3_FilterAccessModes
            COMMAND ${CMAKE_COMMAND}
            -D SZ3_BUILD_DIR=${CMAKE_BINARY_DIR}
            -D SZ3_PREFIX=${CMAKE_CURRENT_BINARY_DIR}/_accessmodes
            -D SZ3_SCRIPT=${CMAKE_CURRENT_SOURCE_DIR}/filterAccessModes.sh
            -D SZ3_H5BIN=${H5Z_SZ3_TOOL_DIR}
            -P ${CMAKE_CURRENT_SOURCE_DIR}/runAccessModes.cmake)
else ()
    message(STATUS "SZ3_FilterAccessModes not registered: h5repack, h5dump and h5ls are not all on PATH")
endif ()
