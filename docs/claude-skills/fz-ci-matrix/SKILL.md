---
name: fz-ci-matrix
description: Use when changing how fz is built, installed, packaged or consumed — CMake, SZ3Config.cmake.in, the bundled Zstd, the HDF5 filter, or anything a downstream project links. Lists every way fz is reached in the wild and what each one has to be tested with. Triggers include "packaging", "find_package", "zstd", "HDF5 filter", "why did the consumer build break".
---

# How fz is actually reached, and what each way has to be tested with

fz is a header-only library with an HDF5 filter and a CLI. Almost every packaging defect it has
shipped came from a path nobody exercised: the repository's own tests reach SZ3 through
`add_subdirectory` or FetchContent, and **neither reads `SZ3Config.cmake`**, so `find_package(SZ3)`
was broken for at least two releases without a single test noticing.

The rule this file exists to enforce: **a way of reaching fz that nothing runs is a way of reaching
fz that is broken.**

## The axes

A "case" is one point in the product of these. They are not independent — some combinations are
impossible — but a change to packaging can break any of them individually.

| Axis | Values |
|---|---|
| How fz is obtained | `add_subdirectory` · FetchContent · **installed + `find_package(SZ3)`** · distro/Spack/conda package |
| Zstd | system via pkg-config · system via `find_package(zstd CONFIG)` · vendored (`SZ3_USE_BUNDLED_ZSTD=ON`) · **none present** · a version different from the build machine's |
| Prefix shape | one prefix per package (Debian) · one shared prefix (conda, Spack view, Homebrew) · relocated after install · **build tree deleted** |
| HDF5 | none · serial · **parallel (MPI)** · found as a CMake *module* vs as a *config* · 1.10.x, 1.14.x, 2.x |
| How the filter is reached | linked + `H5Z_SZ3_initialize()` · linked but never initialized, found via `HDF5_PLUGIN_PATH` · both · **plugin only, caller links no SZ3** |
| Who drives the filter | the C++ API · a linked C application · `h5repack` · `h5dump` · `h5ls` |
| Consumer language | C++ · **C** (an HDF5 filter is exactly what a C program links) |
| Library shape | static · shared |
| Threading | serial · OpenMP, and **write with N threads, read with M** |
| Platform | Linux/gcc · macOS/clang · Windows MSVC · Windows MinGW |
| Standard library | libstdc++ · libc++ |

The bolded values are the ones that have actually caught defects in this repository.

## The cases, and what each has caught

Each row is a thing that broke, found by running the case in the left column.

| Case | What it caught |
|---|---|
| Install, then `find_package(SZ3)` from a separate project | `SZ3Targets.cmake` named `PkgConfig::ZSTD`, an imported target that exists only inside the build that created it, so every consumer died at generate time |
| …linking by the plain name `hdf5sz3`, not `SZ3::hdf5sz3` | an installed SZ3 exported only namespaced names, so the line became `-lhdf5sz3` |
| …as a C project | `find_dependency(OpenMP)` without `COMPONENTS CXX` reported `SZ3_FOUND=1`, then failed at link |
| …as an optional `QUIET` probe that is allowed to decline | `find_package(... REQUIRED)` inside the config file called `FATAL_ERROR` and ended the consumer's configure, taking down `GMX_USE_SZ3=AUTO` |
| …after the build tree is deleted | the export carried the build machine's resolved `libzstd.so` path |
| …with a *different* HDF5 than the build machine's | `${HDF5_INCLUDE_DIRS}` was `PUBLIC` with no `$<BUILD_INTERFACE:>`, so the consumer compiled against the wrong `hdf5.h` |
| …on a shared prefix (conda, Spack, Homebrew) | exporting the Zstd include directory put a whole prefix, holding other packages' headers, ahead of the consumer's own |
| …where `zstd.h` is not on the default include path | the installed public header `#include`d `<zstd.h>` |
| Build with no network at configure time | the bundled Zstd was fetched from a URL; there was no offline path at all |
| Build with `tools/zstd/` deleted | Debian and Fedora forbid embedded code copies, so the vendored tree must be deletable, not merely unused |
| Install the bundled Zstd into a prefix that already has Zstd | it installed `libzstd.so` with an unversioned SONAME over the system one, plus four headers |
| `ldd` / `nm` the installed filter | no RUNPATH (`libhdf5.so.320 => not found`); and the bundled Zstd's symbols interposed on the consumer's own |
| Create a dataset creation property list *after* `H5Zregister` | `H5Zfilter_avail()` answers whether the filter is registered with the *library*, not whether it is on *this* property list, so the code called `H5Pmodify_filter` on a fresh list and faulted inside HDF5 |
| Read a property list that has no SZ3 filter | the return of `H5Pget_filter_by_id` was ignored, silently overwriting the caller's `Config` with zeros |
| `h5repack` with the plugin unreachable | exits **0** and quietly writes an *unfiltered* copy |
| `h5repack -f NONE` on a filtered file | exits **0** and loses the dataset outright |
| `h5dump`/`h5ls` on a filtered file with no plugin | header must still print; data must fail with a diagnostic naming the filter |
| A reader that links no SZ3 at all | the only way in is the plugin; this is what every third-party tool does |
| An application that links the filter and calls nothing | whether the link survives `--as-needed` — and the control that proves the check is not vacuous |
| OpenMP: write with N threads, read with M | `OMP_NUM_THREADS` **cannot** lower a `num_threads(n)` request, so a test that sets it on the reader tests nothing; `OMP_THREAD_LIMIT` is the one that works |
| Decompress a file written by an older fz | `Config::load` parsed foreign `cd_values` before anything established the format; its guards were the only defence, and `if (N > 4)` did not catch a negative `char` |
| Compress, then compare bytes against a recorded digest | v3.3.0 and v3.3.1 shipped different layouts under the same `SZ3_DATA_VERSION` and nothing caught it |
| Every header compiled on its own, under both standard libraries | a header that only works because another was included first |

## Who consumes fz, and how

Every one of these is a way fz can break that this repository's tests cannot see. When changing
packaging or the compressed format, this is the blast radius.

| Consumer | How it takes fz | Version it has | Notes |
|---|---|---|---|
| **hdf5plugin** (silx-kit, what `h5py` users get) | git subtree of the SZ3 source, built into its wheels | **3.1.7**, pinned at `4bbe9df` since 2022-12-08 and unchanged across 13 releases | The largest population by far. Their update PR (#359) has been a draft since 2025-09-24; the blocker is that they build `cd_values` in Python against the 9-integer layout. `tools/H5Z-SZ3/test/cdvalueHelper.py` in this repository answers it and they have not picked it up. Their docs already warn "Backward compatibility is currently not guaranteed". |
| **GROMACS** | links an installed SZ3 (`GMX_USE_SZ3` = EXTERNAL/AUTO/INTERNAL/OFF), plus a vendored copy at `src/external/SZ3-bio` | tracking current | H5MD lossy output is in the prototype branch only; 2026.2 release notes say "only lossless output is supported". The merge request that links an installed SZ3 ([gromacs!5560](https://gitlab.com/gromacs/gromacs/-/merge_requests/5560)) sat in Draft waiting on `find_package(SZ3)` working. `AUTO` **must** fall back to internal rather than fail, which is why `SZ3Config.cmake` may never `FATAL_ERROR`. |
| **MPTRAC** (Jülich Lagrangian transport model) | vendors `libs/SZ3-3.2.1.tar.bz2`, calls `SZ_compress`/`SZ_decompress` over raw `FILE*` | **3.2.1** | Writes its own `.sz3` container, **not** HDF5, so no filter id is involved. Self-consistent: writes and reads with the same vendored version. |
| **ClickHouse** | SZ3 column codec, native parts | **unknown — worth checking** | Codec merged 2026-06-29 (#108788). Open issue #111139: "CHECK TABLE reports a checksum mismatch for any part written with the lossy SZ3 codec". Another case of persisted SZ3 data that does not read back. |
| **Spack** | `sz3` package, filter opt-in | **3.2.0**, untouched since 2024-09-06 | Builds the same way a consumer does, so it met the `find_package` wall. ADIOS2 >= 2.12 defaults to `+sz3`. |
| **Arch Linux** | `extra/sz` | 3.3.2 | HDF5 filter **ON**. |
| **FreeBSD** | `science/sz3` | 3.3.2 | HDF5 filter **OFF**. |
| **Debian/Ubuntu** | SZ3 3.1.7 source inside `python-hdf5plugin` | 3.1.7 | Stripped (`HDF5PLUGIN_STRIP=all`), so nothing is compiled. No standalone SZ3 package. |
| **The HDF5 filter registry** | filter id **32024** | n/a | Maintained by disheng222, ayzk and robertu94. The registry entry documents **no `cd_values` layout at all** — unlike ZFP (32013) and Delta-Rice (32025), which give parameter tables and `h5repack --filter` examples. `community/sz3/README.md` in `HDFGroup/hdf5_plugins` is an unfilled template. The HDF Group does not build or ship SZ3. |

Not recorded here because the details were not confirmed: **libpressio**, the **sz3-rs** Rust crate,
and **conda-forge**. Each is known to package or wrap SZ3; the version each carries should be
established before the next format change.

### What this implies for a format change

- The dangerous direction is **a current fz writing a file that hdf5plugin 3.1.7 reads**, which is
  today's default install for anyone with both. Measured: `H5Dread` returns 0, the reader reports
  success, and every value is wrong. For files written by 3.2.0 through 3.3.1 the same reader hits
  `printf("Decompression Error: Unknown Datatype"); exit(0)` — **exit code zero**.
- `cd_values` carries **no magic and no version**; the payload carries both. That asymmetry is why
  the raw API refuses foreign data cleanly and the HDF5 path can hand back garbage.
- Consumers that vendor a copy and write and read with it (MPTRAC, ClickHouse) are self-consistent
  and only at risk when they upgrade.

## What GitHub Actions can and cannot run

It can run more than is there now. The limits worth knowing:

- **6 hours per job**, 35 days per workflow; public repositories get generous concurrency.
- Runners have **no GPU**, ~4 CPUs, ~14 GB RAM, ~14 GB free disk on `ubuntu-latest`.
- **MPI is fine**: `apt-get install -y libhdf5-openmpi-dev openmpi-bin` gives a parallel HDF5 in
  `/usr/lib/x86_64-linux-gnu/hdf5/openmpi`. Point `HDF5_ROOT` at it.
- **Several HDF5 versions are fine** via a matrix: apt for the distro one, conda-forge
  (`h5py/setup-conda` or micromamba) for 1.10.6 / 1.14.x / 2.x.
- **True offline is awkward.** A network namespace needs privileges the runner does not give.
  What works: point the fetch at a black-hole proxy, or simply assert that no `FetchContent`
  download step exists. Prefer the second — it is what is actually being claimed.
- **A different Zstd version** is easy: build one from source into a throwaway prefix.
- **GROMACS is feasible but expensive.** Its `h5md-test` target alone is ~20-30 min on a runner;
  a full `gmx` build is longer. Put it on a schedule or a label, not on every push.
- **libc++ on Linux** needs `apt-get install -y libc++-dev libc++abi-dev` and
  `-stdlib=libc++` with clang.

## The gap, as of this writing

What CI runs today, and what only ever ran on a developer's machine:

| Case | In CI? |
|---|---|
| Install + `find_package` + link + run | Linux only |
| `filter_access_modes.sh`, 31 assertions over h5repack/h5dump/h5ls and 4 application shapes | Linux only |
| Export carries no build-machine path; survives the build tree being deleted | Linux only |
| Packager deletes `tools/zstd/` | Linux only |
| Bundled Zstd: private, hidden symbols, no header installed | Linux only (MSVC defaults to it, so the Windows job exercises the build) |
| Format digest pinned to `SZ3_DATA_VERSION` | Linux only |
| OpenMP write-N read-M | Linux only |
| Linux vs macOS produce identical bytes | yes |
| Windows MSVC, Windows MinGW | build + ctest + round trip only |
| **Parallel HDF5 / MPI** | **no** |
| **More than one HDF5 version** | **no** |
| **`find_package(SZ3)` from a C project** | **no** |
| **A consumer on a shared prefix (conda/Homebrew)** | **no** |
| **Standalone-header compile, libstdc++ and libc++** | **no** — `tools/test/check_headers.py` exists and is run by hand |
| **Cross-version decompression against released fz** | **no** |
| **GROMACS as a consumer** | **no** |
| **static vs shared, `add_subdirectory` as a consumer** | **no** |

## How to add a case

1. Write it as an assertion that **names the path taken**, not just an exit code. `filter_access_modes.sh`
   is the model: `want <name> <expected substring> <file>`.
2. **Prove it can fail.** Break the thing it tests, watch the check go red, put it back. A check that
   has never failed is not known to be a check.
3. Two shell forms in this repository have produced checks that could not fail:
   - `set -e` is disabled inside a subshell on the left of `&&`
   - `set -e` exempts a `!`-inverted command, so `! grep ...` reads as a check and is not one
4. If the case needs something the runner lacks, say so in the step rather than skipping silently.
   A suite that reports more checks than it ran is the failure this whole file exists to catch.
