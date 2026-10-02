SZ3: A Modular Error-bounded Lossy Compression Framework for Scientific Datasets
=====

SZ3 compresses floating-point and integer arrays from simulations and instruments, and guarantees that every
decompressed value differs from the original by no more than the error bound you set. It is a header-only C++17
library. It can also be used from C, Python, the `sz3` command line, and HDF5 through a filter.

## Installation

Requirements:
* A C++17 compiler and CMake 3.19 or newer.
* Zstd (optional). If pkg-config does not find libzstd, the vendored copy of Zstd 1.5.6 in `tools/zstd` is built and
  linked statically as `libsz3_zstd`, private to SZ3.
* OpenMP (optional). SZ3 uses it when CMake finds it; `-DCMAKE_DISABLE_FIND_PACKAGE_OpenMP=ON` builds without it.
* HDF5 for the HDF5 filter, and ParaView for the ParaView plugin.

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>
cmake --build build -j
cmake --install build
```
The tools go to `<INSTALL_DIR>/bin` and the headers to `<INSTALL_DIR>/include`.

Build options, all ON/OFF switches passed to `cmake` as `-D<option>=ON` or `-D<option>=OFF`:

| Option | Default | Enables |
|---|---|---|
| `BUILD_SHARED_LIBS` | ON, the parent's when SZ3 is added with `add_subdirectory` or FetchContent | shared libraries; OFF builds static ones |
| `BUILD_SZ3_BINARY` | ON, OFF when SZ3 is added with `add_subdirectory` or FetchContent | the `sz3` executable, the C API library SZ3c and the H5Z-SZ3 tools |
| `BUILD_H5Z_FILTER` | OFF | the HDF5 filter H5Z-SZ3 (needs HDF5) |
| `BUILD_MDZ` | OFF | the compressor from the [MDZ paper](https://ieeexplore.ieee.org/document/9835212) (`tools/mdz`), for molecular dynamics of solid materials; for biomolecular trajectories use `ALGO_BIOMD` or `ALGO_BIOMDXTC` |
| `BUILD_PARAVIEW_PLUGIN` | OFF | the ParaView reader plugin (needs ParaView) |
| `BUILD_TESTING` | OFF, the parent's when SZ3 is added with `add_subdirectory` or FetchContent | the unit tests |
| `SZ3_USE_BUNDLED_ZSTD` | OFF (ON with MSVC) | Zstd from `tools/zstd` instead of the system one |
| `SZ3_DEBUG_TIMINGS` | OFF | debug timing output |
| `SZ3_INSTALL` | ON, OFF when SZ3 is added with `add_subdirectory` or FetchContent | the install rules |

On Linux with glibc 2.35 or newer, `GLIBC_TUNABLES=glibc.malloc.hugetlb=1` in the environment backs SZ3's large buffers with huge pages, which makes compressing and decompressing large fields 10-20% faster.

## Interfaces

| Interface | How to use / where | Maintained by |
|---|---|---|
| C++ | `#include <SZ3/api/sz.hpp>`; see [below](#c) | SZ3 |
| C | [tools/sz3c/include/sz3c.h](tools/sz3c/include/sz3c.h), library `SZ3c`; SZ2-compatible functions | SZ3 |
| Python | `pip install pysz`; [tools/pysz](tools/pysz/README.md) | SZ3 |
| Command line | `sz3`; see [below](#command-line) | SZ3 |
| HDF5 filter | H5Z-SZ3, filter ID 32024; [tools/H5Z-SZ3](tools/H5Z-SZ3/README.md) | SZ3 |
| ParaView | SZ3Reader plugin; [tools/paraview](tools/paraview/README.md) | SZ3 |
| Fortran | [ofmla/sz3_simple_example](https://github.com/ofmla/sz3_simple_example) | [Oscar Mojica](https://github.com/ofmla) |
| Rust | [sz3-rs](https://github.com/apertus-open-source-cinema/sz3-rs) | [Juniper Tyree](https://github.com/juntyr) and [Robin Heinemann](https://github.com/rroohhh) |
| Python numcodecs | [numcodecs-rs codecs/sz3](https://github.com/juntyr/numcodecs-rs/blob/main/codecs/sz3/) | [Juniper Tyree](https://github.com/juntyr) |

### Command line

```bash
# Compress 8x8x128 float data (dimensions fastest-varying first) with a relative error bound of 1e-3
sz3 -f -i tools/sz3/testfloat_8_8_128.dat -z test.sz -3 8 8 128 -M REL 1e-3
# Decompress to test.out; -a compares with the original (-i) and prints the maximum error and compression ratio
sz3 -f -z test.sz -o test.out -i tools/sz3/testfloat_8_8_128.dat -a
```
Run `sz3 -h` for every option; the common ones:

| Option | Meaning |
|---|---|
| `-f`, `-d`, `-I 32`, `-I 64` | data type: float, double, int32, int64 |
| `-i <file>` | original data, raw binary |
| `-z <file>` | compressed file: written when compressing, read when decompressing |
| `-o <file>` | decompressed data, raw binary (`-t` writes text) |
| `-1 nx`, `-2 nx ny`, `-3 nx ny nz`, `-4 nx ny nz np` | dimensions, fastest-varying first: `-3 nx ny nz` is `data[nz][ny][nx]` |
| `-M <mode> <bound>` | error-bound mode and bound: `ABS`, `REL`, `PSNR`, `NORM`; `ABS_AND_REL` and `ABS_OR_REL` take `-A <abs> -R <rel>` |
| `-c <file>` | configuration file, such as [tools/sz3/sz3.config](tools/sz3/sz3.config); choose the algorithm here |
| `-a` | after decompression, print error statistics (needs `-i`) |
| `-p` | print the configuration of the compressed data |
| `-v` | print the SZ3 version and the data-format version |

### C++

```cpp
#include <SZ3/api/sz.hpp>
#include <vector>

int main() {
    std::vector<float> data(100 * 200 * 300, 1.0f);
    SZ3::Config conf(100, 200, 300);  // 300 is the fastest-varying dimension
    conf.errorBoundMode = SZ3::EB_ABS;
    conf.absErrorBound = 1e-3;

    size_t cmpSize;
    char *cmpData = SZ_compress(conf, data.data(), cmpSize);

    SZ3::Config decConf;  // filled from the compressed data
    float *decData = SZ_decompress<float>(decConf, cmpData, cmpSize);

    delete[] cmpData;
    delete[] decData;
}
```

To use SZ3 in a CMake project, add `<INSTALL_DIR>` to `CMAKE_PREFIX_PATH`, call `find_package(SZ3)`, and link one of:
* `SZ3::SZ3core`: SZ3 with only the dependencies it cannot work without (Zstd).
* `SZ3::SZ3`: `SZ3::SZ3core` plus the optional dependencies SZ3 was built with (OpenMP, when the consumer's compiler supports it).
* `SZ3::hdf5sz3`: the HDF5 filter; see [tools/H5Z-SZ3/README.md](tools/H5Z-SZ3/README.md).

Data compressed with `SZ3::SZ3core` or `SZ3::SZ3` can be decompressed with either.

`include/SZ3/api/sz.hpp` documents the rest of the API.

## Algorithms and error-bound modes

Set the algorithm with `Config::cmprAlgo`, or `CmprAlgo` in a configuration file.

| Algorithm | Use it for |
|---|---|
| `ALGO_INTERP_LORENZO` (default) | Most data. It tunes interpolation and Lorenzo prediction on a sample of the data and keeps whichever is better. |
| `ALGO_INTERP` | Interpolation with the parameters you set, without auto-tuning; for users who tune those parameters themselves. |
| `ALGO_LORENZO_REG` | Blockwise Lorenzo and regression prediction, the SZ2 algorithm. |
| `ALGO_NOPRED` | Quantization without prediction: a fast baseline. |
| `ALGO_LOSSLESS` | Zstd only. SZ3 also switches to it by itself when the error bound is 0, or when Zstd alone gives a smaller result. |
| `ALGO_BIOMD`, `ALGO_BIOMDXTC` | Molecular-dynamics coordinates. `ALGO_BIOMD` takes `{frames, atoms, 3}` in nm with an absolute (or relative) bound and gives other input to `ALGO_LORENZO_REG`. `ALGO_BIOMDXTC` follows GROMACS's xtc and can, like xtc, round a coordinate slightly past the bound. |

Set the error-bound mode with `Config::errorBoundMode`, or `-M` on the command line.

| Mode | Bound |
|---|---|
| `EB_ABS` (`ABS`) | Every value is within `absErrorBound` of the original. |
| `EB_REL` (`REL`) | Within `relErrorBound` × (max − min) of the data. |
| `EB_ABS_AND_REL` / `EB_ABS_OR_REL` | The smaller / the larger of the two bounds above. |
| `EB_PSNR` (`PSNR`), `EB_L2NORM` (`NORM`) | A target PSNR or L2-norm error. SZ3 converts it to an absolute bound, which is the pointwise guarantee. |

## Data format and compatibility

* SZ3 can decompress data compressed by some earlier versions; see [CHANGELOG.md](CHANGELOG.md) for which versions.
* Builds with different floating-point options (FMA, fast math) decompress the same values, except data compressed by earlier versions built with FMA; see [#162](https://github.com/szcompressor/SZ3/pull/162) and [#166](https://github.com/szcompressor/SZ3/pull/166).

## Citing SZ3

[//]: # (**Kindly note**: If you mention SZ3 in your paper, the most appropriate citation is to include these three references &#40;**TBD22, ICDE21, Bigdata18**&#41; because they cover the design and implementation of the latest version of SZ.)
* QOZv2 (the enhanced interpolation-based algorithm): [High-performance Effective Scientific Error-bounded Lossy Compression with Auto-tuned Multi-component Interpolation](https://dl.acm.org/doi/10.1145/3639259).
* SZ3's interpolation-based algorithm: [Optimizing Error-Bounded Lossy Compression for Scientiﬁc Data by Dynamic Spline Interpolation](https://ieeexplore.ieee.org/document/9458791).
* The software engineering design of SZ3: [SZ3: A modular framework for composing prediction-based error-bounded lossy compressors](https://ieeexplore.ieee.org/abstract/document/9866018).

## Version history

See [CHANGELOG.md](CHANGELOG.md).

## License and contact

SZ3 is released under a BSD license; see [copyright-and-BSD-license.txt](copyright-and-BSD-license.txt).
`include/SZ3/encoder/XtcBasedEncoder.hpp` is based on GROMACS and licensed under the LGPL, version 2.1 or later.
The vendored Zstd in `tools/zstd` keeps its own license.

* Lead developer and maintainer: Kai Zhao
* Contributors: Robert Underwood, Xin Liang, Jinyang Liu, Sheng Di, and [everyone else on GitHub](https://github.com/szcompressor/SZ3/graphs/contributors)
* SZ project lead: Franck Cappello

Report bugs and ask questions in [GitHub issues](https://github.com/szcompressor/SZ3/issues).
