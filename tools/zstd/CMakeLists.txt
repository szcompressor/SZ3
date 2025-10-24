include(FetchContent)

if(POLICY CMP0135)
        cmake_policy(SET CMP0135 NEW)
endif()
FetchContent_Declare(
        zstdfetched
        URL "https://github.com/facebook/zstd/releases/download/v1.4.5/zstd-1.4.5.tar.gz"
        SOURCE_SUBDIR ""  # Disables automatic add_subdirectory()
)
FetchContent_MakeAvailable(zstdfetched)

if (MSVC)
    set(ZSTD_BUILD_TYPE STATIC)
else()
    set(ZSTD_BUILD_TYPE SHARED)
endif()
set(zstd_SOURCE_DIR "${zstdfetched_SOURCE_DIR}/lib")

add_library(zstd ${ZSTD_BUILD_TYPE}
        "${zstd_SOURCE_DIR}/common/entropy_common.c"
        "${zstd_SOURCE_DIR}/common/pool.c"
        "${zstd_SOURCE_DIR}/common/threading.c"
        "${zstd_SOURCE_DIR}/common/debug.c"
        "${zstd_SOURCE_DIR}/common/xxhash.c"
        "${zstd_SOURCE_DIR}/common/fse_decompress.c"
        "${zstd_SOURCE_DIR}/common/zstd_common.c"
        "${zstd_SOURCE_DIR}/common/error_private.c"

        "${zstd_SOURCE_DIR}/compress/zstd_ldm.c"
        "${zstd_SOURCE_DIR}/compress/zstd_lazy.c"
        "${zstd_SOURCE_DIR}/compress/huf_compress.c"
        "${zstd_SOURCE_DIR}/compress/zstd_opt.c"
        "${zstd_SOURCE_DIR}/compress/zstd_double_fast.c"
        "${zstd_SOURCE_DIR}/compress/zstd_compress.c"
        "${zstd_SOURCE_DIR}/compress/zstd_compress_superblock.c"
        "${zstd_SOURCE_DIR}/compress/zstd_compress_sequences.c"
        "${zstd_SOURCE_DIR}/compress/zstd_compress_literals.c"
        "${zstd_SOURCE_DIR}/compress/zstd_fast.c"
        "${zstd_SOURCE_DIR}/compress/fse_compress.c"
        "${zstd_SOURCE_DIR}/compress/zstdmt_compress.c"
        "${zstd_SOURCE_DIR}/compress/hist.c"

        "${zstd_SOURCE_DIR}/decompress/zstd_decompress.c"
        "${zstd_SOURCE_DIR}/decompress/huf_decompress.c"
        "${zstd_SOURCE_DIR}/decompress/zstd_ddict.c"
        "${zstd_SOURCE_DIR}/decompress/zstd_decompress.c"
        "${zstd_SOURCE_DIR}/decompress/zstd_decompress_block.c"

        "${zstd_SOURCE_DIR}/deprecated/zbuff_common.c"
        "${zstd_SOURCE_DIR}/deprecated/zbuff_compress.c"
        "${zstd_SOURCE_DIR}/deprecated/zbuff_decompress.c"

        "${zstd_SOURCE_DIR}/legacy/zstd_v05.c"
        "${zstd_SOURCE_DIR}/legacy/zstd_v04.c"
        "${zstd_SOURCE_DIR}/legacy/zstd_v06.c"
        "${zstd_SOURCE_DIR}/legacy/zstd_v07.c"
        "${zstd_SOURCE_DIR}/legacy/zstd_v03.c"
        "${zstd_SOURCE_DIR}/legacy/zstd_v02.c"
        "${zstd_SOURCE_DIR}/legacy/zstd_v01.c"

        "${zstd_SOURCE_DIR}/dictBuilder/cover.c"
        "${zstd_SOURCE_DIR}/dictBuilder/divsufsort.c"
        "${zstd_SOURCE_DIR}/dictBuilder/zdict.c"
        "${zstd_SOURCE_DIR}/dictBuilder/fastcover.c"
)

target_include_directories(zstd
        PUBLIC
        $<BUILD_INTERFACE:${zstd_SOURCE_DIR}/>
        PRIVATE
        ${zstd_SOURCE_DIR}/common
        ${zstd_SOURCE_DIR}/compress
        ${zstd_SOURCE_DIR}/decompress
        ${zstd_SOURCE_DIR}/deprecated
        ${zstd_SOURCE_DIR}/dictBuilder
        ${zstd_SOURCE_DIR}/dll
        ${zstd_SOURCE_DIR}/legacy
)

install(TARGETS zstd EXPORT SZ3Targets
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

install(FILES
        "${zstd_SOURCE_DIR}/deprecated/zbuff.h"
        "${zstd_SOURCE_DIR}/dictBuilder/zdict.h"
        "${zstd_SOURCE_DIR}/zstd.h"
        "${zstd_SOURCE_DIR}/common/zstd_errors.h"
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/
)

export(TARGETS zstd FILE zstd.cmake)