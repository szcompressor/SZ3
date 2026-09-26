enable_testing()
# Try to find GTest locally
find_package(GTest QUIET)
if(NOT GTest_FOUND)
    message(STATUS "GTest not found locally; fetching from GitHub...")
    include(FetchContent)
    FetchContent_Declare(
            googletest
            GIT_REPOSITORY https://github.com/google/googletest.git
            GIT_TAG        v1.16.0
    )
    if (MSVC)
        set(BUILD_SHARED_LIBS OFF)
    endif ()
    FetchContent_MakeAvailable(googletest)
else()
    message(STATUS "Using locally installed GTest")
endif()

include(CTest)
include(GoogleTest)

# Collect all test source files
file(GLOB TEST_SOURCES "${CMAKE_SOURCE_DIR}/tools/test/modules/*.cpp")

foreach (TEST_SRC IN LISTS TEST_SOURCES)
    get_filename_component(TEST_NAME ${TEST_SRC} NAME_WE)
    message(STATUS "Test target = ${TEST_NAME}")
    add_executable(${TEST_NAME} ${TEST_SRC})
    target_link_libraries(${TEST_NAME} PRIVATE GTest::gtest_main SZ3)
    # Automatically discover tests for this executable
    # Use PRE_TEST mode to avoid running executable at build time (important for Windows DLL path issues)
    gtest_discover_tests(${TEST_NAME} DISCOVERY_MODE PRE_TEST)
endforeach()

# A program built without OpenMP reads back what one built with it wrote with conf.openmp = true.
if (OpenMP_FOUND)
    add_executable(openmp_stream_with openmp_stream.cpp)
    target_link_libraries(openmp_stream_with PRIVATE SZ3)
    add_executable(openmp_stream_without openmp_stream.cpp)
    target_link_libraries(openmp_stream_without PRIVATE SZ3core)

    add_test(NAME SZ3_OpenMPStream_WriteWithOpenMP
            COMMAND openmp_stream_with with write ${CMAKE_CURRENT_BINARY_DIR}/stream_with)
    # More than one thread, or the stream holds one chunk and the read below proves nothing.
    set_tests_properties(SZ3_OpenMPStream_WriteWithOpenMP PROPERTIES
            ENVIRONMENT OMP_NUM_THREADS=4 FIXTURES_SETUP stream_with)
    add_test(NAME SZ3_OpenMPStream_ReadWithoutOpenMP
            COMMAND openmp_stream_without without read ${CMAKE_CURRENT_BINARY_DIR}/stream_with)
    set_tests_properties(SZ3_OpenMPStream_ReadWithoutOpenMP PROPERTIES FIXTURES_REQUIRED stream_with)
endif ()
