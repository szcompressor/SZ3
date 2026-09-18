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

# test_zstd_decl is the one test that needs a real zstd.h, because proving SZ3's declarations agree
# with the header is the whole of what it does. SZ3_ZSTD_HEADER_DIRS can be empty only if a packager
# deleted tools/zstd/ AND pkg-config named no include directory, in which case the test is dropped
# -- with a warning, because a test that quietly stops testing is worse than no test at all.
if (NOT SZ3_ZSTD_HEADER_DIRS)
    message(WARNING
            "No zstd.h was located, so test_zstd_decl (which checks that Lossless_zstd.hpp's "
            "declarations still match Zstd's own header) will not be built. Everything else is "
            "unaffected: SZ3 itself needs no zstd.h.")
    list(REMOVE_ITEM TEST_SOURCES "${CMAKE_SOURCE_DIR}/tools/test/modules/test_zstd_decl.cpp")
endif ()

foreach (TEST_SRC IN LISTS TEST_SOURCES)
    get_filename_component(TEST_NAME ${TEST_SRC} NAME_WE)
    message(STATUS "Test target = ${TEST_NAME}")
    add_executable(${TEST_NAME} ${TEST_SRC})
    target_link_libraries(${TEST_NAME} PRIVATE GTest::gtest_main SZ3)
    if (TEST_NAME STREQUAL "test_zstd_decl")
        # The second translation unit, which includes the same two headers in the opposite order.
        target_sources(${TEST_NAME} PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/zstd_decl_header_first.cpp")
        # PRIVATE, and only for this target: putting a zstd.h on SZ3's own interface is exactly
        # what this change exists to stop.
        target_include_directories(${TEST_NAME} PRIVATE ${SZ3_ZSTD_HEADER_DIRS})
    endif ()
    # Automatically discover tests for this executable
    # Use PRE_TEST mode to avoid running executable at build time (important for Windows DLL path issues)
    gtest_discover_tests(${TEST_NAME} DISCOVERY_MODE PRE_TEST)
endforeach()