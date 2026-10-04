# SZ3 advanced guide

This guide is for choosing an algorithm and error bound, building and packaging SZ3, using its C++ API beyond the
defaults, and changing SZ3's code. To compress with the default settings, the [README](../README.md) is enough.

## Algorithms and error-bound modes

Set the algorithm with `Config::cmprAlgo`, or `CmprAlgo` in a configuration file.

| Algorithm | Use it for |
|---|---|
| `ALGO_INTERP_LORENZO` (default) | Most data. It tunes interpolation and Lorenzo prediction on a sample of the data and keeps whichever is better. |
| `ALGO_INTERP` | Interpolation with the parameters you set, without auto-tuning; for users who tune those parameters themselves. |
| `ALGO_LORENZO_REG` | Blockwise Lorenzo and regression prediction, the SZ2 algorithm. |
| `ALGO_NOPRED` | Quantization without prediction: a fast baseline. |
| `ALGO_LOSSLESS` | Zstd only. SZ3 also switches to it by itself when the error bound is 0, negative or not finite, or when Zstd alone gives a smaller result. |
| `ALGO_BIOMD`, `ALGO_BIOMDXTC` | Molecular-dynamics coordinates, `{frames, atoms, 3}`. `ALGO_BIOMDXTC` follows GROMACS's xtc and can, like xtc, round a coordinate slightly past the bound. See [molecular-dynamics.md](molecular-dynamics.md). |

Set the error-bound mode with `Config::errorBoundMode`, or `-M` on the command line.

| Mode | Bound |
|---|---|
| `EB_ABS` (`ABS`) | Every value is within `absErrorBound` of the original. |
| `EB_REL` (`REL`) | Within `relErrorBound` × (max − min) of the data. |
| `EB_ABS_AND_REL` / `EB_ABS_OR_REL` | The smaller / the larger of the two bounds above. |
| `EB_PSNR` (`PSNR`), `EB_L2NORM` (`NORM`) | A target PSNR or L2-norm error. SZ3 converts it to an absolute bound, which is the pointwise guarantee. |

## Building

Requirements: a C++17 compiler and CMake 3.19 or newer. Zstd, OpenMP and HDF5 are optional (see
[Dependencies](#dependencies)).

| Option | Default | Effect |
|---|---|---|
| `BUILD_SHARED_LIBS` | ON at top level; the parent's in a subproject | shared `libhdf5sz3` and `libSZ3c`; OFF builds static ones |
| `BUILD_SZ3_BINARY` | ON at top level, OFF in a subproject | `sz3`, `sz3_smoke_test`, the C API library `libSZ3c`, the filter's tools |
| `BUILD_H5Z_FILTER` | OFF | the HDF5 filter (needs HDF5) |
| `BUILD_MDZ` | OFF | `mdz`, for molecular dynamics of solid materials ([tools/mdz](../tools/mdz/README.md)) |
| `BUILD_PARAVIEW_PLUGIN` | OFF | the ParaView reader (needs ParaView) |
| `BUILD_TESTING` | OFF | the tests; only when SZ3 is the top-level project |
| `SZ3_INSTALL` | ON at top level, OFF in a subproject | install rules |
| `SZ3_USE_BUNDLED_ZSTD` | OFF (ON with MSVC) | use `tools/zstd` even when a system Zstd is found |
| `SZ3_DEBUG_TIMINGS` | OFF | print timings |
| `H5Z_SZ3_PLUGIN_INSTALL_DIR` | `lib/plugin` | where the plugin copy of the filter goes; empty for none |
| `CMAKE_DISABLE_FIND_PACKAGE_OpenMP` (CMake's own) | OFF | ON builds SZ3 without OpenMP |

On Linux with glibc 2.35 or newer, `GLIBC_TUNABLES=glibc.malloc.hugetlb=1` in the environment backs SZ3's large
buffers with huge pages, which makes compressing and decompressing large fields 10-20% faster.

### Dependencies

- **Zstd.** pkg-config first, else the vendored Zstd in `tools/zstd`, built as the private static target `sz3_zstd`
  with hidden symbols. Nothing is downloaded at configure time, so SZ3 builds offline. No SZ3 header includes
  `zstd.h` (`Lossless_zstd.hpp` declares the four functions it calls), so consumers need no Zstd include directory;
  keep it that way.
- **OpenMP** is optional and used when CMake finds it; `find_package(OpenMP)` is not `REQUIRED`.
- **HDF5** (C component) for the filter. Never add an MPI search: the filter calls no MPI function, and a parallel
  HDF5 brings its own MPI.
- Do not add a dependency for code nothing builds.

### What is installed

| Path | Content |
|---|---|
| `bin/` | `sz3`, `sz3_smoke_test`, `mdz`; DLLs on Windows |
| `include/SZ3/`, `include/SZ3c/`, `include/hdf5_sz3/` | the C++ headers with the generated `version.hpp`, the C API, the filter |
| `lib/` | `libSZ3c`, `libhdf5sz3` |
| `lib/plugin/` | a copy of the filter for `HDF5_PLUGIN_PATH` (shared builds only) |
| `lib/cmake/SZ3/` | the CMake package |

- `libhdf5sz3` and `libSZ3c` carry `VERSION` and a major.minor `SOVERSION`. Their signatures and `SZ3::Config`'s
  layout change only in a new major.minor.
- Both are built with hidden visibility; only `HDF5SZ3_EXPORT` / `SZ3C_API` functions are exported, so another SZ3
  copy in the process cannot take over their calls.
- Keep `INSTALL_RPATH_USE_LINK_PATH` on `hdf5sz3`: without a RUNPATH the installed filter cannot find a conda, Spack
  or module HDF5. On Windows the filter DLL is in `bin/`, which must be on `PATH`.
- A static `hdf5sz3` is position independent and defines `HDF5SZ3_STATIC` for its consumers. Its own functions
  are hidden in the shared library that links it; HDF5's `H5PLextern.h` makes `H5PLget_plugin_type` and
  `H5PLget_plugin_info` exported.
- Install rules use `$<TARGET_FILE:...>`. `bin/sz3_smoke_test` exits 0 on a working install; Spack runs it.

## Using SZ3 from CMake

### An installed SZ3

```cmake
find_package(SZ3 3.4 REQUIRED)                       # C++ API
target_link_libraries(app PRIVATE SZ3::SZ3)          # or SZ3::SZ3core
find_package(SZ3 3.4 REQUIRED COMPONENTS hdf5sz3)    # the HDF5 filter
target_link_libraries(app PRIVATE SZ3::hdf5sz3)
```

| Target | What it is |
|---|---|
| `SZ3::SZ3core` | The header-only C++17 library and Zstd. No optional dependencies. |
| `SZ3::SZ3` | `SZ3::SZ3core` plus OpenMP, when SZ3 was built with OpenMP and the consumer's compiler has it. Keep it carrying OpenMP: libpressio relies on it. |
| `SZ3::hdf5sz3` | The HDF5 filter library. Usable from C. Its consumers get no OpenMP flags. |

- Data compressed through `SZ3::SZ3core` or `SZ3::SZ3` decompresses through either.
- `SZ3` and `hdf5sz3` without the namespace are defined too, so `target_link_libraries(app PRIVATE hdf5sz3)` works
  against an installed and an in-tree SZ3 alike.
- `find_package(SZ3 <version>)` accepts the same major version. `libSZ3c` has no CMake target; link it by name.
- `SZ3Config.cmake` never ends the consumer's configure. Every search that can decline runs before any target is
  created; when one declines, it sets `SZ3_FOUND` false with a reason and leaves no targets behind, so a project can
  fall back to its own copy. No `REQUIRED` searches inside it.
- Never export build-machine paths or build-only targets: wrap such paths in `$<BUILD_INTERFACE:...>` and quote
  generator expressions that hold paths.
- If the consumer already has `HDF5::HDF5`, SZ3 does not search HDF5 again. OpenMP is searched per enabled language.

### As a subproject, or copied into another tree

```cmake
set(BUILD_H5Z_FILTER ON)   # only if you need the filter
FetchContent_Declare(SZ3 GIT_REPOSITORY https://github.com/szcompressor/SZ3.git GIT_TAG v3.4.0)
FetchContent_MakeAvailable(SZ3)
target_link_libraries(app PRIVATE SZ3::hdf5sz3)
```

- As a subproject SZ3 installs nothing, builds no tools, sets no cache defaults (`BUILD_SHARED_LIBS`,
  `CMAKE_BUILD_TYPE`), and builds neither its tests nor the filter's, even when the parent sets `BUILD_TESTING`. Any
  new top-level-only behaviour tests `SZ3_IS_TOP_LEVEL`.
- Depend on a release tag, never on a `master` commit.
- A project that copies SZ3 needs `CMakeLists.txt`, `include/`, `tools/zstd` (without
  `lib/decompress/huf_decompress_amd64.S`, which `ZSTD_DISABLE_ASM` leaves out) and, for the filter,
  `tools/H5Z-SZ3` without `test/` and `tools/`. Keep the license files: `XtcBasedEncoder.hpp` is LGPL-2.1 or later.
- The options of the parts left out still appear in the parent's cache, and turning one on then fails. Set them as
  normal variables before adding SZ3 (CMP0077), and mark the rest advanced:

  ```cmake
  set(SZ3_INSTALL OFF)
  set(BUILD_SZ3_BINARY OFF)
  set(BUILD_MDZ OFF)
  set(BUILD_PARAVIEW_PLUGIN OFF)
  # after adding SZ3:
  mark_as_advanced(SZ3_DEBUG_TIMINGS SZ3_USE_BUNDLED_ZSTD)
  ```
- Link `SZ3::SZ3core` (and `hdf5sz3`) unless you want SZ3's OpenMP; `SZ3::SZ3` adds `-fopenmp` and libgomp to the
  consumer. Put warning suppressions on `SZ3core`. To build the bundled SZ3 without OpenMP, set
  `CMAKE_DISABLE_FIND_PACKAGE_OpenMP` before adding it.

## The C++ API in depth

- **Output buffers.** To compress into your own buffer, size it with `SZ3::SZ_compress_size_bound`; `SZ_compress`
  refuses a smaller one. "Twice the input" is too small for small arrays and for integers.

  ```cpp
  std::vector<char> buf(SZ3::SZ_compress_size_bound<float>(conf));
  size_t cmpSize = SZ_compress(conf, data.data(), buf.data(), buf.size());
  ```
- **Decompressing into your own buffer.** `SZ_decompress(conf, cmpData, cmpSize, decData)` writes every value of the
  original data and cannot check the buffer's size. With `decData == nullptr` it allocates with `new[]`.
- **Bounds that are not positive and finite.** An error bound that is 0, negative or not finite compresses losslessly
  with every algorithm; a relative bound over data with Inf ends up there too.
- **Integer data** is compressed as floating point (float for 1- and 2-byte types, double for 4- and 8-byte types),
  then rounded and clamped back to the type, so every value is within `floor(bound)` of the original. 8-byte integers
  beyond ±2^53 are refused.
- **Recompression.** Data decompressed and compressed again can exceed the bound, except with `ALGO_BIOMDXTC`, with
  `ALGO_NOPRED` under an absolute bound, and with `ALGO_BIOMD` under an absolute bound. HDF5 recompresses a chunk
  that is written again.
- `ALGO_BIOMDXTC` is not strict, like GROMACS's xtc: a coordinate can come back up to 10% past the bound.
- `include/SZ3/api/sz.hpp` documents the rest of the API.

## OpenMP

OpenMP is used only when you ask for it: set `conf.openmp = true`, or `OpenMP = YES` in a configuration file.

- SZ3 splits the data along its slowest dimension (the first one in `SZ3::Config`, the last one on the `sz3` command
  line) into one part per OpenMP thread, at most one per index of that dimension, and compresses the parts separately.
  The number of parts is stored in the data.
- Decompression is correct with any number of threads: fewer than were used to compress, under `OMP_THREAD_LIMIT`, or
  inside your own parallel region.
- A build without OpenMP decompresses data compressed with OpenMP, one part after another, and ignores `conf.openmp`
  when compressing.
- SZ3 does not change your program's thread count and prints nothing. The HDF5 filter does not use OpenMP.
- `SZ3::SZ3` carries OpenMP to your target; `SZ3::SZ3core` does not.
- With MinGW, nested OpenMP regions crash in libgomp itself.

Rules for SZ3's own code:

- Never call `omp_set_num_threads`; use `num_threads(n)` on the region.
- Hand out chunks with `parallel for`, never by `omp_get_thread_num()`: a region may get fewer threads than asked.
- No exception may leave a parallel region. Store it per chunk and rethrow after the region.
- Index with `size_t`. Respect the stream's `openmp` flag on decompression.
- Test with `OMP_THREAD_LIMIT`, not `OMP_NUM_THREADS`; CI asserts that CMake found OpenMP. Skip nested-region tests
  on MinGW.

## Floating point and compiler flags

Builds with different floating-point options (FMA, fast math in Clang, icx and MSVC) decompress the same values.

- GCC's `-ffast-math` / `-Ofast` and x87 floating point (`FLT_EVAL_METHOD > 0`) stop the build with an error; on
  32-bit x86 build with `-msse2 -mfpmath=sse`.
- Clang, icx and MSVC fast modes are made precise inside SZ3's code with `#pragma float_control(precise, on)` in
  `SZImpl.hpp`, for code reached through `SZ3/api/sz.hpp`: include `sz.hpp` before any other SZ3 header.
- GCC's `-fassociative-math`, `-funsafe-math-optimizations` and `-ffinite-math-only` on their own are not detected and
  not supported.

Rules for SZ3's own code:

- The decompressor recomputes the compressor's reconstructed values, so both must round every step the same way.
  Every product that feeds a prediction or a reconstruction goes through `SZ3::nofma()` (`def.hpp`), which rounds it
  before the add; otherwise a build that fuses `a * b + c` (Apple Silicon, aarch64, `-march=native`, Spack's default)
  decodes other builds' data far past the bound.
- `nofma()` uses a `volatile`. Do not switch it to `__builtin_assoc_barrier`: GCC 13 and 14 drop it in vectorized
  loops.
- The compressor reconstructs through the decompressor's code (`LinearQuantizer::recover_pred`).
- Check the error in double: checked in `T`, float errors are rounded before the check and can pass the bound.
- Never test for NaN or Inf with floating-point comparisons in code that must hold under `-ffinite-math-only`; GCC
  removes them. Compare the bits.
- New codecs keep fusable floating-point arithmetic off the reconstruction path; `ALGO_BIOMD` predicts on an integer
  lattice for that reason.
- Verify with `tools/test/crossBuildDecode.sh` (two builds decode each other's streams) and `test_reconstruction`. CI
  builds a `-march=x86-64-v3` job and checks that it emits FMA.

## HDF5 filter (H5Z-SZ3, id 32024)

[tools/H5Z-SZ3/README.md](../tools/H5Z-SZ3/README.md) covers installing the filter, using it from `h5repack`, h5py
and C, registering it in an application, chunks and other filters, appending, and the rules for the filter's code.

## Data format

SZ3 can decompress data compressed by some earlier versions; [CHANGELOG.md](../CHANGELOG.md) says which.

Rules for SZ3's own code:

- `SZ3_DATA_VERSION` (CMakeLists.txt) is the version of the compressed format, separate from the program version. Bump
  it when, and only when, the bytes a build writes change. CI pins the bytes to the declared version in
  `tools/test/data_format_digest.txt`; if your change moves the digest, either the format changed (bump the version,
  update the digest and the CHANGELOG) or you broke something.
- `SZ_compress` writes `SZ3_MAGIC_NUMBER` and this build's data version, never the values in the caller's `Config`.
- Refuse data you cannot decode correctly, and say why.
- The same input compresses to the same bytes on every run and every build. Initialise every member that `save()`
  writes, zero buffers that go into the output, and never let hash-map order or ties decide anything that reaches the
  stream.
- Stored data is little-endian on every host. Use `SZ3::read` / `SZ3::write`, not `memcpy` of host words.
- `size_est()` is an upper bound on what `save()` writes for the current input. One buffer is sized from it and the
  writers get no capacity, so an under-report corrupts the heap silently.
- `Config::save()` has a one-byte length prefix; outgrowing it is an error.
- Changing a `concepts::` interface or removing a public header breaks out-of-tree code (libpressio, hdf5plugin,
  GROMACS). Check them and say so in the PR.
- The library never prints and never calls `exit()`. To tell the caller something, throw; in the HDF5 filter, refuse
  in `set_local` with the reason on the HDF5 error stack.

## Untrusted input

- The compressed stream is untrusted. Check every length, count or index read from it against the bytes that remain
  before using it to index, allocate or loop; `decode(..., remaining_length)` and `load(c, remaining_length)` charge
  what they consume.
- Lossless `decompress(src, srcLen, dst, dstCap)` takes a capacity.
- Bound allocations by the input, not by a count the input declares.
- Subtract budgets through a checked helper. Do not guard stream values with `assert`. Own buffers with RAII.
- Put checks where a symbol is produced, not on every bit, to keep decompression fast.

## ALGO_BIOMD and ALGO_BIOMDXTC

[molecular-dynamics.md](molecular-dynamics.md) is for users. For SZ3's code:

- `ALGO_BIOMD` takes `{frames, atoms, 3}` (a one-frame chunk becomes `{atoms, 3}`); other shapes throw
  `std::invalid_argument`, and a chunk of more values than an `int` counts throws.
- Coordinates go on an integer lattice spanning `|x| < B`, B the smallest power of two above 2^21 times the bound for
  float (2^27 for double). Coordinates past B, and NaN or Inf outside trailing fill frames, store the chunk losslessly.
  B depends on the bound only, so decompressed data compressed again gives the same values; keep it that way.
- BIOMD and `ALGO_LORENZO_REG` compress the caller's data without a copy because neither modifies it; a change that
  writes to the input restores the copy.
- `ALGO_BIOMDXTC` uses `XtcBasedEncoder.hpp`, LGPL-2.1 or later.
- `{1, 1024, 3}` is a 2D case, because `setDims` drops a dimension of 1: tests of the multi-frame path need two frames
  or more. Fill-frame tests need values off the quantization grid.

## Integrating SZ3 into an application

- Update the bundled SZ3 in the same commit as the code that adapts to it, so every commit builds and passes.
- Choose the algorithm per dataset: `ALGO_BIOMD` for coordinates only, `ALGO_NOPRED` for velocities and forces,
  lossless for the rest. Asking for BIOMD on another type is a programming error, not a file error.
- HDF5 filter tests: see the [filter README](../tools/H5Z-SZ3/README.md#writing-a-chunk-again-and-appending).
- Say in the integrating commit which SZ3 data the new version cannot read.

## Tests and CI

- A new test must fail on the code it guards; break the guarded code once to see it go red.
- Check determinism without relying on dirty heap memory: build the object over storage filled with two patterns and
  compare what `save()` writes.
- The ASan + UBSan job runs ctest and the CLI; it leaves the HDF5 filter out.
- `size_est()` is tested directly: save into a larger buffer and require the cursor to stay within what was declared.
- Installed-package checks build a separate project against the install (`consumeInstalledSZ3.sh`), including one
  without the filter; `filterAccessModes.sh` covers `h5repack`, `h5dump`, `h5ls` and the registration shapes. CI also
  builds libpressio and GROMACS against SZ3.
- Integration tests check compressed sizes against `compression_baseline.json` within 5% either way; timings are
  recorded, not gated.
- `SZ3::uint` is not the global `uint`; MSVC and MinGW have none.
- A failing tool in a pipe fails the check (`set -o pipefail`).
- Headers stay warning-clean under consumers' warning flags and compile on their own under libstdc++ and libc++:
  include what you use.
- Format with the repository's `.clang-format`, and leave lines you did not change alone.

## Commits, PRs and comments

- PRs are squash-merged; the commit title is the PR title plus `(#N)`, and the body says what changed in behaviour and
  why.
- PR titles state the outcome. A description says what went wrong (with a reproducer or a measurement), what changed,
  what it costs, whether the format changed, and how it was tested, including what was not run.
- Numbers, not adjectives: ratios, speed as a ratio of times with how it was measured, and the number of cases.
- The CHANGELOG says what changed for users; details belong in the PR. The README is for first-time users.
- Comments state the reason or constraint a reader needs, not the history of the change.
