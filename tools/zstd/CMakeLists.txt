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

# Static everywhere, as it already was on MSVC. Built shared, this produced a libzstd.so whose
# SONAME was the unversioned development name, so everything linking it recorded
# `NEEDED libzstd.so` -- a name a real system libzstd answers to just as readily.
set(ZSTD_BUILD_TYPE STATIC)
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

# libsz3_zstd, not libzstd: under its own name this overwrote the libzstd already in the install
# prefix. Hidden visibility keeps the archive's ZSTD_* symbols inside whatever links it instead of
# re-exporting them from libhdf5sz3.so, where they interposed on the consumer's own Zstd; the two
# macros undo the visibility("default") Zstd's headers put on its public API.
set_target_properties(zstd PROPERTIES
        OUTPUT_NAME sz3_zstd
        POSITION_INDEPENDENT_CODE ON
        C_VISIBILITY_PRESET hidden)
target_compile_definitions(zstd
        PUBLIC SZ3_BUNDLED_ZSTD=1
        PRIVATE ZSTDLIB_VISIBILITY= ZSTDERRORLIB_VISIBILITY=)

# Staged where only the full SZ3/bundled_zstd/ path reaches it, so no directory SZ3 puts on a
# consumer's include path holds a zstd.h their own #include <zstd.h> would open. SZ3 already exposes
# both of these include roots, so this target needs no public include directory of its own.
configure_file("${zstd_SOURCE_DIR}/zstd.h"
        "${SZ3_BINARY_DIR}/include/SZ3/bundled_zstd/zstd.h" COPYONLY)

target_include_directories(zstd
        PRIVATE
        ${zstd_SOURCE_DIR}
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

# Only zstd.h: it is all Lossless_zstd.hpp includes, and it includes nothing else of Zstd's.
install(FILES
        "${SZ3_BINARY_DIR}/include/SZ3/bundled_zstd/zstd.h"
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/SZ3/bundled_zstd/
)

export(TARGETS zstd FILE zstd.cmake)