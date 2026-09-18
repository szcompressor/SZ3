include(FetchContent)

if(POLICY CMP0135)
        cmake_policy(SET CMP0135 NEW)
endif()
# 1.5.6 and not 1.5.7, on compression THROUGHPUT rather than output size; both are measured in
# the commit message, with their units, on Zstd alone and again end to end. 1.5.7 rewrote what
# level 3 does -- its changelog lists "compression ratio improvement for dfast, aka levels 3 and
# 4" and "better block boundaries" -- and level 3 is what Lossless_zstd uses. Measured on Zstd
# alone, 1.5.7 compresses 26-30% fewer MB/s than the 1.4.5 it would replace while 1.5.6 compresses
# 15-38% more; output size moves by under 0.6% for either. End to end that is worth ~3% of
# SZ_compress on ALGO_INTERP_LORENZO, where Zstd is only 3-5% of the work, and 9-29% on
# ALGO_LOSSLESS, where Zstd is the whole algorithm. 1.5.6 is the fastest of the three in every
# case measured, which is why it is here; 1.5.7 is not disqualified, just slower.
FetchContent_Declare(
        zstdfetched
        URL "https://github.com/facebook/zstd/releases/download/v1.5.6/zstd-1.5.6.tar.gz"
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
# visibility macros undo the visibility("default") Zstd's headers put on its public API.
# ZSTD_DISABLE_ASM, because the list above is C only: on x86-64 huf_decompress.c otherwise calls
# HUF_decompress4X{1,2}_usingDTable_internal_fast_asm_loop, which only decompress/huf_decompress_amd64.S
# defines, and the link fails. Assembling that .S would need enable_language(ASM), for a language
# MSVC does not assemble at all. What that costs was measured rather than assumed: built both ways
# on an x86-64 with BMI2, the assembly loop decompresses 0.1-2.0% faster on three of four payloads
# and 10.6% faster on the fourth. Even giving that up, 1.5.6 decompresses every one of them 2-22%
# faster than the 1.4.5 it replaces, so nothing regresses. Adding the .S is a live option if some
# day the MSVC cost is acceptable; it builds and runs correctly under GCC.
set_target_properties(zstd PROPERTIES
        OUTPUT_NAME sz3_zstd
        POSITION_INDEPENDENT_CODE ON
        C_VISIBILITY_PRESET hidden)
target_compile_definitions(zstd
        PUBLIC SZ3_BUNDLED_ZSTD=1
        PRIVATE ZSTDLIB_VISIBILITY= ZSTDERRORLIB_VISIBILITY= ZSTD_DISABLE_ASM=1)

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