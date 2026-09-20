add_executable(filterLargeDims filterLargeDims.cpp)
target_link_libraries(filterLargeDims PRIVATE hdf5sz3)
add_test(NAME SZ3_FilterLargeDims
        COMMAND filterLargeDims ${CMAKE_CURRENT_BINARY_DIR}/filterLargeDims.h5)

# CI installs the HDF5 tools and runs filterAccessModes.sh directly, so coverage does not depend
# on this registration.
find_program(H5REPACK_EXECUTABLE h5repack)
find_program(H5DUMP_EXECUTABLE h5dump)
find_program(H5LS_EXECUTABLE h5ls)
if (H5REPACK_EXECUTABLE AND H5DUMP_EXECUTABLE AND H5LS_EXECUTABLE)
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
