# The vendored Zstd. Built only when no system Zstd was resolved, or when SZ3_USE_BUNDLED_ZSTD is
# ON; the root CMakeLists decides, this file only builds what it is told to.
#
# ---------------------------------------------------------------------------------------------
# Provenance
#   Upstream:  https://github.com/facebook/zstd/releases/download/v1.5.6/zstd-1.5.6.tar.gz
#   Version:   1.5.6
#   SHA256:    8c29e06cf42aacc1eafc4077ae2ec6c6fcb96a626157e0593d5e82a34fd403c1
#   License:   see LICENSE beside this file (BSD-3-Clause; upstream dual-licenses GPLv2)
#   Contents:  lib/common, lib/compress, lib/decompress, lib/zstd.h and lib/zstd_errors.h,
#              copied byte for byte out of that tarball. lib/deprecated (the ZBUFF API),
#              lib/legacy (decoders for pre-1.0 frames) and lib/dictBuilder (ZDICT) are not
#              here: SZ3 calls none of them, they were 14 of the 39 files this used to compile,
#              and a frame SZ3 wrote is never a pre-1.0 frame.
#
#   To verify this copy was not modified locally, which is the point of copying whole
#   directories rather than picking files out of them:
#     curl -L -o /tmp/z.tgz https://github.com/facebook/zstd/releases/download/v1.5.6/zstd-1.5.6.tar.gz
#     tar xzf /tmp/z.tgz -C /tmp
#     diff -r /tmp/zstd-1.5.6/lib/common    tools/zstd/lib/common
#     diff -r /tmp/zstd-1.5.6/lib/compress  tools/zstd/lib/compress
#     diff -r /tmp/zstd-1.5.6/lib/decompress tools/zstd/lib/decompress
#
#   To re-vendor a different version, replace those three directories and the two headers from
#   the new tarball, update the version and hash above, and re-read the note below first.
#
# ---------------------------------------------------------------------------------------------
# Why 1.5.6 and not 1.5.7, and what would change that
#
#   On throughput. 1.5.7 rewrote what compression level 3 does, and level 3 is the level
#   Lossless_zstd uses. Measured on Zstd alone, 1.5.7 compresses 26-30% fewer MB/s than the
#   1.4.5 this replaced while 1.5.6 compresses 15-38% more; output size moves under 0.6% either
#   way. End to end that is ~3% of SZ_compress on ALGO_INTERP_LORENZO and 9-29% on ALGO_LOSSLESS.
#   The numbers, with their units and baselines, are in the commit that took this to 1.5.6.
#
#   1.5.7's one listed correctness fix does not apply to SZ3, and this is the note to re-read
#   before changing the version, because after vendoring there is no URL left to carry it.
#   That fix (facebook/zstd#4129, commit 09cb37cb) caps ZSTD_CURRENT_MAX below 2 GB in 32-bit
#   mode: through 1.5.6 the cap is (3U << 29) + (1U << ZSTD_WINDOWLOG_MAX), which on a 32-bit
#   target is 2.5 GB, and once the window index `curr` passes 2 GB the ptrdiff_t arithmetic in
#   ZSTD_window_update() goes negative and the encoder can emit a match offset outside the
#   window. `curr` only reaches 2 GB inside a ZSTD_CCtx that is reused across compressions.
#   SZ3 holds no context: it calls the one-shot ZSTD_compress/ZSTD_decompress, and ZSTD_compress
#   builds a fresh ZSTD_CCtx per call and frees it on return, so `curr` restarts at
#   ZSTD_WINDOW_START_INDEX (2) every time and ends at srcSize + 2. Measured, not assumed: 5 GB
#   pushed through one-shot calls reaches max curr = srcSize + 2 exactly, at 64 MiB, 256 MiB and
#   1 GiB per call, while the same 5 GB through one reused ZSTD_CCtx reaches ZSTD_CURRENT_MAX.
#   The same defect is in the 1.4.5 this replaced, so 1.5.6 is not new exposure.
#
#   THIS ASSESSMENT EXPIRES if SZ3 ever acquires a lossless path that holds a ZSTD_CCtx across
#   calls -- a streaming compressor, or any reuse of ZSTD_compressCCtx, which accumulates `curr`
#   just as streaming does. To re-check:
#     grep -rE 'ZSTD_(createCCtx|CCtx|createCStream|CStream|compressStream|compress2)' include tools
#   That must stay empty. If it does not, 1.5.7 or newer becomes required for 32-bit builds.
# ---------------------------------------------------------------------------------------------

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

# libsz3_zstd, not libzstd: under its own name this overwrote the libzstd already in the install
# prefix. Hidden visibility keeps the archive's ZSTD_* symbols inside whatever links it instead of
# re-exporting them from libhdf5sz3.so, where they interposed on the consumer's own Zstd; the two
# visibility macros undo the visibility("default") Zstd's headers put on its public API.
set_target_properties(zstd PROPERTIES
        OUTPUT_NAME sz3_zstd
        POSITION_INDEPENDENT_CODE ON
        C_VISIBILITY_PRESET hidden)

# ZSTD_DISABLE_ASM, because the list above is C only: on x86-64 huf_decompress.c otherwise calls
# HUF_decompress4X{1,2}_usingDTable_internal_fast_asm_loop, which only decompress/huf_decompress_amd64.S
# defines, and the link fails. Assembling that .S would need enable_language(ASM), for a language
# MSVC does not assemble at all. What that costs was measured rather than assumed: built both ways
# on an x86-64 with BMI2, the assembly loop decompresses 0.1-2.0% faster on three of four payloads
# and 10.6% faster on the fourth. Even giving that up, 1.5.6 decompresses every one of them 2-22%
# faster than the 1.4.5 it replaces, so nothing regresses. The .S is still in lib/decompress, both
# because these directories are kept byte-identical to upstream and because adding it stays a live
# option: it builds and runs correctly under GCC.
# ZSTD_LEGACY_SUPPORT=0 is stated rather than left to the default, because lib/legacy is not
# vendored: zstd_decompress.c and zstd_ddict.c include ../legacy/zstd_legacy.h when it is >= 1.
target_compile_definitions(zstd
        PRIVATE ZSTDLIB_VISIBILITY= ZSTDERRORLIB_VISIBILITY= ZSTD_DISABLE_ASM=1 ZSTD_LEGACY_SUPPORT=0)

# PRIVATE, and no public include directory at all. Nothing SZ3 installs includes zstd.h any more
# -- Lossless_zstd.hpp declares the four Simple-API functions it calls -- so this target has no
# header to offer, and offering one would put a zstd.h on a consumer's include path ahead of their
# own. That is also why there is no install(FILES) for a header here.
target_include_directories(zstd
        PRIVATE
        "${_zstd_lib}"
        "${_zstd_lib}/common"
        "${_zstd_lib}/compress"
        "${_zstd_lib}/decompress"
)

# In the export this becomes SZ3::zstd with an IMPORTED_LOCATION of
# ${_IMPORT_PREFIX}/lib/libsz3_zstd.a, so it relocates with the prefix and a consumer has nothing
# to find. That is the whole reason the vendored path needs no resolver in SZ3Config.cmake.
install(TARGETS zstd EXPORT SZ3Targets
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

export(TARGETS zstd FILE zstd.cmake)
