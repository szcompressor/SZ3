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
    # Otherwise `make install` puts googletest into SZ3's prefix.
    set(INSTALL_GTEST OFF)
    FetchContent_MakeAvailable(googletest)
else()
    message(STATUS "Using locally installed GTest")
endif()

include(CTest)
include(GoogleTest)

# Collect all test source files
file(GLOB TEST_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/modules/*.cpp")

foreach (TEST_SRC IN LISTS TEST_SOURCES)
    get_filename_component(TEST_NAME ${TEST_SRC} NAME_WE)
    message(STATUS "Test target = ${TEST_NAME}")
    add_executable(${TEST_NAME} ${TEST_SRC})
    target_link_libraries(${TEST_NAME} PRIVATE GTest::gtest_main SZ3)
    # Automatically discover tests for this executable
    # Use PRE_TEST mode to avoid running executable at build time (important for Windows DLL path issues)
    gtest_discover_tests(${TEST_NAME} DISCOVERY_MODE PRE_TEST)
endforeach()
target_compile_definitions(test_previous_data_version PRIVATE SZ3_TEST_STREAMS="${CMAKE_CURRENT_SOURCE_DIR}/streams")

# GCC refuses -ffast-math but not these options, which together on an FMA target let it drop floating-point NaN tests.
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64|AMD64")
    add_executable(biomd_finite_math biomd_finite_math.cpp)
    target_link_libraries(biomd_finite_math PRIVATE SZ3core)
    target_compile_options(biomd_finite_math PRIVATE -funsafe-math-optimizations -ffinite-math-only -march=x86-64-v3)
    # it runs only where the machine has the instructions -march=x86-64-v3 uses
    include(CheckCXXSourceRuns)
    set(CMAKE_REQUIRED_FLAGS -march=x86-64-v3)
    check_cxx_source_runs("int main() { return __builtin_cpu_supports(\"avx2\") && __builtin_cpu_supports(\"fma\") ? 0 : 1; }"
            SZ3_HOST_HAS_X86_64_V3)
    unset(CMAKE_REQUIRED_FLAGS)
    if (SZ3_HOST_HAS_X86_64_V3)
        add_test(NAME SZ3_BioMD_NaNUnderFiniteMath COMMAND biomd_finite_math)
    endif ()
endif ()

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
    # ALGO_BIOMD coordinates, split along the atoms
    add_test(NAME SZ3_OpenMPStream_BioMDWriteWithOpenMP
            COMMAND openmp_stream_with with write ${CMAKE_CURRENT_BINARY_DIR}/stream_biomd biomd)
    set_tests_properties(SZ3_OpenMPStream_BioMDWriteWithOpenMP PROPERTIES
            ENVIRONMENT OMP_NUM_THREADS=4 FIXTURES_SETUP stream_biomd)
    add_test(NAME SZ3_OpenMPStream_BioMDReadWithoutOpenMP
            COMMAND openmp_stream_without without read ${CMAKE_CURRENT_BINARY_DIR}/stream_biomd)
    set_tests_properties(SZ3_OpenMPStream_BioMDReadWithoutOpenMP PROPERTIES FIXTURES_REQUIRED stream_biomd)

    add_executable(openmp_regions openmp_regions.cpp)
    target_link_libraries(openmp_regions PRIVATE SZ3)
    add_test(NAME SZ3_OpenMP_ChunkThrows COMMAND openmp_regions throw)
    set_tests_properties(SZ3_OpenMP_ChunkThrows PROPERTIES TIMEOUT 60)
    # MinGW's libgomp crashes in nested parallel regions (GCC 15.2 and 16.2).
    if (NOT MINGW)
        add_test(NAME SZ3_OpenMP_FewerThreadsThanAsked COMMAND openmp_regions nested)
        set_tests_properties(SZ3_OpenMP_FewerThreadsThanAsked PROPERTIES
                ENVIRONMENT "OMP_NUM_THREADS=4,4;OMP_THREAD_LIMIT=8;OMP_MAX_ACTIVE_LEVELS=2" TIMEOUT 300)
    endif ()
endif ()
