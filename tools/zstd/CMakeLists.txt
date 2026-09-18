# The vendored Zstd, built when the root CMakeLists found no system Zstd or was asked for this one.
#
#   Upstream: https://github.com/facebook/zstd/releases/download/v1.5.6/zstd-1.5.6.tar.gz
#   SHA256:   8c29e06cf42aacc1eafc4077ae2ec6c6fcb96a626157e0593d5e82a34fd403c1
#   License:  LICENSE beside this file (BSD-3-Clause; upstream dual-licenses GPLv2)
#
# lib/common, lib/compress, lib/decompress and the two headers are byte-for-byte from that tarball,
# so `diff -r` shows whether this tree was modified and a re-vendor replaces them wholesale.
# lib/deprecated, lib/legacy and lib/dictBuilder are omitted: SZ3 calls none of them.
#
# Not 1.5.7: it rewrote compression level 3, the level Lossless_zstd uses, for 26-30% less
# throughput and under 0.6% of output size. Its correctness fix (facebook/zstd#4129) caps
# ZSTD_CURRENT_MAX so 32-bit arithmetic in ZSTD_window_update() cannot emit an out-of-window match
# offset, which needs a ZSTD_CCtx reused across compressions. SZ3 holds none -- one-shot
# ZSTD_compress builds and frees a CCtx per call. That stops being true if SZ3 ever keeps one:
#   grep -rE 'ZSTD_(createCCtx|CCtx|createCStream|CStream|compressStream|compress2)' include tools

set(_zstd_lib "${CMAKE_CURRENT_SOURCE_DIR}/lib")

# Listed rather than globbed: a glob would silently start compiling whatever a re-vendor drops in.
add_library(zstd STATIC
        "${_zstd_lib}/common/debug.c"
        "${_zstd_lib}/common/entropy_common.c"
        "${_zstd_lib}/common/error_private.c"
        "${_zstd_lib}/common/fse_decompress.c"
        "${_zstd_lib}/common/pool.c"
        "${_zstd_lib}/common/threading.c"
        "${_zstd_lib}/common/xxhash.c"
        "${_zstd_lib}/common/zstd_common.c"

        "${_zstd_lib}/compress/fse_compress.c"
        "${_zstd_lib}/compress/hist.c"
        "${_zstd_lib}/compress/huf_compress.c"
        "${_zstd_lib}/compress/zstd_compress.c"
        "${_zstd_lib}/compress/zstd_compress_literals.c"
        "${_zstd_lib}/compress/zstd_compress_sequences.c"
        "${_zstd_lib}/compress/zstd_compress_superblock.c"
        "${_zstd_lib}/compress/zstd_double_fast.c"
        "${_zstd_lib}/compress/zstd_fast.c"
        "${_zstd_lib}/compress/zstd_lazy.c"
        "${_zstd_lib}/compress/zstd_ldm.c"
        "${_zstd_lib}/compress/zstd_opt.c"
        "${_zstd_lib}/compress/zstdmt_compress.c"

        "${_zstd_lib}/decompress/huf_decompress.c"
        "${_zstd_lib}/decompress/zstd_ddict.c"
        "${_zstd_lib}/decompress/zstd_decompress.c"
        "${_zstd_lib}/decompress/zstd_decompress_block.c"
)

# libsz3_zstd, not libzstd, or this overwrites the libzstd already in the install prefix. Hidden
# visibility keeps the archive's ZSTD_* symbols out of libhdf5sz3.so, where they interposed on the
# consumer's own Zstd.
set_target_properties(zstd PROPERTIES
        OUTPUT_NAME sz3_zstd
        POSITION_INDEPENDENT_CODE ON
        C_VISIBILITY_PRESET hidden)

# The visibility macros undo the visibility("default") Zstd's headers put on its public API, which
# C_VISIBILITY_PRESET alone cannot. ZSTD_DISABLE_ASM because the sources above are C only: without
# it huf_decompress.c calls loops that only huf_decompress_amd64.S defines, and assembling that
# needs a language MSVC does not have. ZSTD_LEGACY_SUPPORT=0 because lib/legacy is not vendored and
# zstd_decompress.c includes it when that is >= 1.
target_compile_definitions(zstd
        PRIVATE ZSTDLIB_VISIBILITY= ZSTDERRORLIB_VISIBILITY= ZSTD_DISABLE_ASM=1 ZSTD_LEGACY_SUPPORT=0)

# PRIVATE, and no public include directory: offering one would put a zstd.h on a consumer's include
# path ahead of their own. Same reason no header is installed from here.
target_include_directories(zstd
        PRIVATE
        "${_zstd_lib}"
        "${_zstd_lib}/common"
        "${_zstd_lib}/compress"
        "${_zstd_lib}/decompress"
)

install(TARGETS zstd EXPORT SZ3Targets
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

export(TARGETS zstd FILE zstd.cmake)
