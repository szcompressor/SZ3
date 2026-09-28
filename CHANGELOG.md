# SZ3 version history

Oldest version first. Each entry says whether the version has a git tag and a GitHub release, and which data formats it can decompress when that changed. A version that was current for less than two weeks has no entry of its own; its changes are listed under the next version that has one.

## 3.0.0 (2021-01-22)

No tag and no GitHub release. This is commit `fb76989`, "first release of sz3".

- First release of SZ3, the C++ version of SZ. A compressor is put together from separate modules: predictor, quantizer, encoder and lossless stage.
- Header-only C++17 library. `sz_demo` compresses float data with Lorenzo and regression prediction, Huffman coding and Zstd.
- Zstd had to be pointed to by editing `ZSTD_LIBS` and `ZSTD_INCLUDES` in `CMakeLists.txt`.

## 3.0.1 (2021-04-15)

Tag `v3.0.1`. The GitHub release has no notes.

**New**
- Zstd is found through pkg-config. If it is not found, a bundled copy of Zstd 1.3.5 is built.
- Arithmetic and bypass encoders, a pre-processor stage (transpose, filter), and a truncation compressor with its example `sz_truncate`.
- `sz_fast`, a compressor built on `SZMetaFrontend`, which follows the SZauto design.

**Changes**
- CMake 3.13 or newer is required.

**Fixes**
- Faster Huffman and arithmetic decoding.

## 3.0.2 (2021-06-01)

Tag `v3.0.2`. The GitHub release has no notes.

**New**
- Point-wise relative error bound, through the separate example `sz_pw`.
- The Zstd level can be set when constructing `Lossless_zstd`.
- The truncation compressor can keep a chosen number of bytes per value.

## 3.1.3 (2022-02-04)

Tag `v3.1.3`. The GitHub release is titled "V3.1.3". This entry also covers 3.1.0 (tag `v3.1.0`, 2022-01-21), 3.1.1 (tag `v3.1.1`, 2022-01-26) and 3.1.2 (tag `v3.1.2`, 2022-02-02). Each was current for less than two weeks, and each has a GitHub release.

**New**
- Interpolation-based compression (the ICDE'21 dynamic spline interpolation). The default algorithm, `ALGO_INTERP_LORENZO`, chooses between interpolation and Lorenzo prediction and tunes its settings automatically. `ALGO_INTERP` and `ALGO_LORENZO_REG` pick one predictor family.
- A public API in `include/SZ3/api/sz.hpp`: `SZ_compress(conf, data, outSize)` and `SZ_decompress(conf, cmpData, cmpSize)`, configured with `SZ::Config`. `SZ_decompress(conf, cmpData, cmpSize, decData)` decompresses into a buffer the caller allocated.
- A command-line tool, `sz`: `-i` input, `-z` compressed file, `-o` decompressed file, `-f`/`-d` data type, `-1`…`-4` dimensions, `-M <mode> <bound>`, and `-a` for error statistics. Most SZ2 command-line arguments are still accepted. `-h2` lists them.
- INI configuration files, passed with `-c`. The example is `test/sz.config`, now `tools/sz3/sz3.config`.
- OpenMP compression and decompression for all algorithms. Set `conf.openmp = true`, or `OpenMP = YES` in a configuration file. The data is split along its slowest dimension, and only the ABS and REL modes are supported with OpenMP.
- Error-bound modes ABS, REL, PSNR, NORM (L2 norm), ABS_AND_REL and ABS_OR_REL. On the command line the bounds can also be given with `-A`, `-R`, `-S` and `-N`.
- Integer data. The command line takes 32-bit and 64-bit integers with `-I 32` or `-I 64`.
- Text output of decompressed data, with `-t`.
- Data with more than four dimensions (folded into four).
- `make install` installs a CMake package. `find_package(SZ3)` provides the `SZ3::SZ3` target.
- CMake options `SZ3_USE_BUNDLED_ZSTD` and `SZ3_DEBUG_TIMINGS` (debug timings on by default).
- A polynomial regression predictor and a run-length encoder.
- The README lists the SZ2 command-line options SZ3 does not support (`-c` in SZ2 format, `-p`, `-T`, `-P`).

**Changes**
- Headers moved from `include/` to `include/SZ3/`. Include them as `SZ3/...`.
- The compression configuration is stored at the end of the compressed data, so decompression needs only the compressed bytes.
- `sz` replaces `sz_demo`. `sz_pw` and `sz_truncate` became beta examples.
- The version macros are in the generated header `SZ3/version.hpp`.
- 3.1.1 changed the stored configuration (narrower field types and an `openmp` flag). Data compressed by 3.1.0 cannot be decompressed by 3.1.1 or later.
- 3.1.2 replaced the `#define` constants of 3.1.0 and 3.1.1 with enums in `namespace SZ`. `METHOD_*` became `ALGO_*` (the field `cmprMethod` became `cmprAlgo`), and `ABS`/`REL` became `EB_ABS`/`EB_REL`.
- 3.1.2 also renamed several `Config` fields, for example `enable_lorenzo` to `lorenzo`, `quant_state_num` to `quantbinCnt`, `block_size` to `blockSize`, and `interp_op` to `interpAlgo`.

## 3.1.3.1 (2022-03-16)

Tag `v3.1.3.1`, no GitHub release.

**Fixes**
- `SZ_compress` no longer modifies the input data.
- Fixes to interpolation tuning and data sampling.
- The build works on 32-bit Windows, and `long` was replaced by `int64_t` for portability.

## 3.1.4 (2022-04-09)

No tag and no GitHub release. This is commit `75b2c7e`.

**New**
- Builds and runs on Windows with MSYS2/MinGW (instructions: https://github.com/szcompressor/SZ3/issues/5#issuecomment-1094039224). Native Visual Studio builds came in 3.1.8.

**Fixes**
- Fixes for OpenMP, the command line and the CMake build.

## 3.1.5.1 (2022-04-29)

Tag `v3.1.5.1`, no GitHub release. This entry also covers 3.1.5 (commit `be95b02`, 2022-04-28, no tag and no GitHub release), which was current for one day.

**New**
- H5Z-SZ3, an HDF5 filter (ID 32024), built with `-DBUILD_H5Z_FILTER=ON`. `sz3ToHDF5` and `dsz3FromHDF5` are test tools. [#11](https://github.com/szcompressor/SZ3/pull/11)

**Changes**
- The bundled Zstd moved to `tools/zstd` and was updated to 1.4.5.
- The examples moved from `test/` to `examples/`, and `sz.config` was renamed `sz3.config`.
- The generated `SZ3/version.hpp` is written to the build directory and installed.

**Fixes**
- libpressio integration fixes. [#6](https://github.com/szcompressor/SZ3/pull/6)
- Explicit type in `std::multiplies`, for compilers that required it. [#10](https://github.com/szcompressor/SZ3/pull/10)
- A memory-allocation fix.

## 3.1.7 (2022-11-12)

Tag `v3.1.7`, no GitHub release. This entry also covers 3.1.5.4 (tag `v3.1.5.4`, 2022-10-24) and 3.1.6 (tag `v3.1.6`, 2022-11-01). Each was current for less than two weeks, and neither has a GitHub release. Versions 3.1.5.2 and 3.1.5.3 exist only as commits.

**New**
- A C API in `tools/sz3c` (`SZ_compress_args`, `SZ_decompress`) with SZ2's signatures and constants.
- A Python API, `tools/pysz/pysz.py`, that loads the SZ3 or SZ2 shared library through ctypes.
- MDZ, a compressor for molecular-dynamics trajectories, in `tools/mdz`, built with `-DBUILD_MDZ=ON`.
- `print_h5repack_args` builds the HDF5 filter's `cd_values` for `h5repack`.
- A smoke test for `sz3`.

**Changes**
- The command-line tool was renamed from `sz` to `sz3`.
- The tools and examples moved from `examples/` to `tools/sz3/`.
- CMake 3.18 or newer is required.
- The HDF5 filter's `cd_values` now use the same layout as H5Z-SZ, the SZ2 filter.
- Data with more than four dimensions is rejected instead of being folded into four.
- When GSL is found, SZ3 links it, and the installed CMake package looks for GSL. [#27](https://github.com/szcompressor/SZ3/pull/27)

**Fixes**
- Fixes to OpenMP (also faster), `SZFastFrontend::size_est()`, the polynomial regression predictor and quantization-bin handling.
- Stride types are `size_t`. [#22](https://github.com/szcompressor/SZ3/pull/22)
- The H5Z-SZ3 header compiles as C, and the filter builds without warnings. [#21](https://github.com/szcompressor/SZ3/pull/21), [#25](https://github.com/szcompressor/SZ3/pull/25)

## 3.1.8 (2023-11-30)

Tag `v3.1.8`. The GitHub release has no notes and a source zip.

**New**
- Native Windows builds with Visual Studio. With MSVC the bundled Zstd is used by default. [#29](https://github.com/szcompressor/SZ3/pull/29)
- The HDF5 filter reads a configuration file named by the `SZ3_CONFIG_PATH` environment variable. `cd_values` still take precedence. [#45](https://github.com/szcompressor/SZ3/pull/45)
- A faster mode for 1D Lorenzo.
- An absolute error bound of 0 compresses losslessly with Zstd.
- MDZ accepts 3D input. [#41](https://github.com/szcompressor/SZ3/pull/41)
- The Python API finds the library in a default location when no path is given. [#32](https://github.com/szcompressor/SZ3/pull/32)

**Changes**
- The C++ namespace changed from `SZ` to `SZ3`, to avoid clashing with SZ2. [#44](https://github.com/szcompressor/SZ3/pull/44)
- `SZ3_DEBUG_TIMINGS` is off by default.
- OpenMP is optional: SZ3 links it only when it is found.
- The Python `decompress` returns data of the original dtype. [#43](https://github.com/szcompressor/SZ3/pull/43)

**Fixes**
- `Config` fields have default values. [#47](https://github.com/szcompressor/SZ3/pull/47)
- `SZFastFrontend::size_est()` no longer returns an undefined value. [#46](https://github.com/szcompressor/SZ3/pull/46)
- Builds with GCC 13 (missing `uint8_t`). [#38](https://github.com/szcompressor/SZ3/pull/38)
- The HDF5 filter links MPI when HDF5 is parallel. [#37](https://github.com/szcompressor/SZ3/pull/37)
- Portable 64-bit integer types, and an uninitialized `max_err`. [#30](https://github.com/szcompressor/SZ3/pull/30), [#28](https://github.com/szcompressor/SZ3/pull/28)
- Fixes for OpenMP and interpolation+Lorenzo sampling. The "OpenMP threads" message is gone.

## 3.2.0 (2024-08-16)

Tag `v3.2.0`, no GitHub release. The v3.2.1 release notes list this version's changes. [#60](https://github.com/szcompressor/SZ3/pull/60)

**New**
- The compressed data starts with a magic number and a data-format version. Decompression rejects data that SZ3 did not compress, or that has another data-format version, and names the version that wrote it.
- `SZ_compress(conf, data, cmpData, cmpCap)` compresses into a buffer the caller allocated and returns the compressed size.
- `ALGO_NOPRED`, which quantizes without prediction.
- Builds on 32-bit systems, such as wasm32, by using `std::unordered_map` there. [#56](https://github.com/szcompressor/SZ3/pull/56)
- Demos of custom pipelines in `tools/sz3/demo`.
- The README links a third-party Fortran API.

**Changes**
- The internal API was restructured into decomposition, quantizer, encoder and lossless modules (`SZGenericCompressor`, `InterpolationDecomposition`, …). Code written against the old frontend/predictor classes must be updated.
- `SZ_compress` and `SZ_decompress` throw `std::invalid_argument` when the output buffer is too small or the input is not SZ3 data.
- The HDF5 filter was rewritten. `cd_values` now hold the serialized SZ3 configuration, and the SZ2-style `cd_values` of 3.1.x are no longer understood.
- Data format 3.2.0. Data from 3.1.x cannot be decompressed.

**Fixes**
- An integer overflow in quantization.
- Fixes for OpenMP, interpolation and MDZ (1D input, and the first frame counted in the compression ratio).

## 3.2.1 (2024-10-03)

Tag `v3.2.1`. The GitHub release lists #56 and #60, which were first in 3.2.0.

**New**
- `ALGO_LOSSLESS`. SZ3 switches to Zstd only when the error bound is 0, when the output buffer is too small for lossy compression, or when Zstd gives the smaller output at a low compression ratio.
- The HDF5 filter has a README. `print_h5repack_args -c sz3.config` and `h5repack.sh` generate `h5repack` parameters from a configuration file.

**Changes**
- The stored configuration is smaller, and the HDF5 filter's `cd_values` are compact. Data format 3.2.1: data from 3.2.0 cannot be decompressed.
- `IntegerQuantizer` was renamed `LinearQuantizer`.
- The code is formatted with clang-format and builds with stricter warnings.

**Fixes**
- When the output buffer was too small, Zstd wrote a truncated stream without reporting an error. `Lossless_zstd` now checks the buffer size first.

## 3.3.0 (2025-08-06)

Tag `v3.3.0`. The GitHub release links [#94](https://github.com/szcompressor/SZ3/pull/94).

**New**
- The interpolation algorithms gain key QoZ 1 and 2 features, anchor points and level-wise error bounds, which improve speed and quality. New configuration keys: `InterpolationAnchorStride`, `InterpolationAlpha` and `InterpolationBeta`. The complete QoZ is on the `QoZ` branch. [#94](https://github.com/szcompressor/SZ3/pull/94)
- Lorenzo and regression run on a new blockwise iterator.
- `sz3 -p` prints the configuration stored in a compressed file, and `sz3 -v` also prints the data-format version.
- `SZ_compress_size_bound()` gives the output-buffer size `SZ_compress` needs.
- CMake options `BUILD_TESTING` (unit and integration tests under `tools/test`) and `SZ3_USE_SKA_HASH`.
- The C API exports its functions from a Windows DLL and adds `free_buf`.
- The README links the third-party Rust and numcodecs APIs.

**Changes**
- The INI parser is SZ3's own. Section names, keys and enum values are case-insensitive, and inih is no longer bundled.
- The polynomial regression, meta-Lorenzo and meta-regression predictors were removed, along with `SZIterateCompressor` and `LorenzoRegressionDecomposition` (replaced by `BlockwiseDecomposition`). The demos became `tools/sz3/sz3_customized_demo.cpp`.
- Dimensions of size 1 are dropped from `Config`.
- `ENABLE_GSL` was renamed `SZ3_ENABLE_GSL`. [#100](https://github.com/szcompressor/SZ3/pull/100)
- The API throws exceptions and writes to stderr, instead of calling `exit(0)`.
- Data format 3.3.0: data from 3.2.x cannot be decompressed.

**Fixes**
- NaN values are stored exactly instead of being replaced by predictions. [#91](https://github.com/szcompressor/SZ3/pull/91)
- Compressing small arrays no longer writes past the buffer. [#70](https://github.com/szcompressor/SZ3/pull/70)
- `SZ_decompress` takes the compressed data as `const`. [#73](https://github.com/szcompressor/SZ3/pull/73)
- MSVC builds again. [#74](https://github.com/szcompressor/SZ3/pull/74)
- A small AddressSanitizer fix. [#96](https://github.com/szcompressor/SZ3/pull/96)
- The Python API copies buffers faster. [#99](https://github.com/szcompressor/SZ3/pull/99)
- Fixes for the regression predictor, OpenMP and the HDF5 filter.

From December 2024 until 3.3.0, master reported version 3.2.2. That version was never tagged or released.

## 3.3.1 (2025-10-23)

No tag and no GitHub release. This is commit `f7f6989` ([#106](https://github.com/szcompressor/SZ3/pull/106)). pysz 1.0.2 was built from it.

**New**
- pysz 1.0, a new Python package (Cython) published on PyPI: `pip install pysz`. [#106](https://github.com/szcompressor/SZ3/pull/106)
- CMake option `BUILD_SZ3_BINARY` (default ON) for the `sz3` tool and the C API.

**Fixes**
- A sampling bug in interpolation tuning. [#104](https://github.com/szcompressor/SZ3/pull/104)
- OpenMP decompression failed when a thread's slice had fewer dimensions than the whole array. [#105](https://github.com/szcompressor/SZ3/pull/105)

The data format is still 3.3.0, so 3.3.1 and 3.3.0 read each other's data.

## 3.3.2 (2025-11-19)

Tag `v3.3.2`. The GitHub release lists #104–#116, including the 3.3.1 changes.

**New**
- `ALGO_BIOMD` and `ALGO_BIOMDXTC`, for molecular-dynamics data, plus a Huffman encoder with less storage overhead. [#115](https://github.com/szcompressor/SZ3/pull/115)
- SZ3Reader, a ParaView plugin that opens SZ3-compressed files, built with `-DBUILD_PARAVIEW_PLUGIN=ON`. [#112](https://github.com/szcompressor/SZ3/pull/112)
- The HDF5 filter builds on Windows with Visual Studio. CI builds and tests SZ3 with both Visual Studio and MinGW. [#107](https://github.com/szcompressor/SZ3/pull/107)
- `cdvalueHelper` converts between an SZ3 configuration and the filter's `cd_values`. [#109](https://github.com/szcompressor/SZ3/pull/109)
- pysz 1.0.3.

**Changes**
- The compressed layout is now a 16-byte header (magic number, data-format version, payload size), then the payload, then the configuration. Data format 3.3.2: data from 3.3.0 and 3.3.1 cannot be decompressed. [#116](https://github.com/szcompressor/SZ3/pull/116)
- Compressed data is little-endian on every host, including big-endian ones. [#109](https://github.com/szcompressor/SZ3/pull/109)

**Fixes**
- A bug in the compressed format; the PR does not describe it further. [#116](https://github.com/szcompressor/SZ3/pull/116)
- The HDF5 filter works on Windows. [#107](https://github.com/szcompressor/SZ3/pull/107)

## 3.4.0

<!-- Entry to be added when #163 (new Huffman encoder, data format 3.4.0) is merged. -->
