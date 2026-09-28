SZ3: A Modular Error-bounded Lossy Compression Framework for Scientific Datasets
=====

SZ3 compresses floating-point and integer arrays from simulations and instruments, and guarantees that every
decompressed value differs from the original by no more than the error bound you set. It is a header-only C++17
library. It can also be used from C, Python, the `sz3` command line, and HDF5 through a filter.

## Quick start

Python:
```bash
pip install pysz
```
```python
import numpy as np
from pysz import sz, szConfig

data = np.random.rand(100, 200).astype(np.float32)
compressed, ratio = sz.compress(data, szConfig())  # default: absolute error bound 1e-3
decompressed, _ = sz.decompress(compressed, np.float32, data.shape)
```
For other error bounds, see [tools/pysz/README.md](tools/pysz/README.md).

Command line, on the 8×8×128 float sample shipped with SZ3:
```bash
git clone https://github.com/szcompressor/SZ3.git && cd SZ3
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build -j
build/tools/sz3/sz3 -f -i tools/sz3/testfloat_8_8_128.dat -z test.sz -3 8 8 128 -M REL 1e-3
build/tools/sz3/sz3 -f -z test.sz -o test.out -i tools/sz3/testfloat_8_8_128.dat -a
```
The first `sz3` command compresses with a bound of 1e-3 × the data's value range. The second decompresses to
`test.out`, and `-a` prints the maximum error, PSNR and compression ratio.

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

Build options. Pass them to `cmake` as `-D<option>=ON` or `OFF`.

| Option | Default | Enables |
|---|---|---|
| `BUILD_SHARED_LIBS` | ON, the parent's when SZ3 is added with `add_subdirectory` or FetchContent | shared libraries; OFF builds static ones |
| `BUILD_SZ3_BINARY` | ON, OFF when SZ3 is added with `add_subdirectory` or FetchContent | the `sz3` executable, the C API library SZ3c and the H5Z-SZ3 tools |
| `BUILD_H5Z_FILTER` | OFF | the HDF5 filter H5Z-SZ3 (needs HDF5) |
| `BUILD_PARAVIEW_PLUGIN` | OFF | the ParaView reader plugin (needs ParaView) |
| `BUILD_TESTING` | OFF, the parent's when SZ3 is added with `add_subdirectory` or FetchContent | the unit tests |
| `SZ3_USE_BUNDLED_ZSTD` | OFF (ON with MSVC) | Zstd from `tools/zstd` instead of the system one |
| `SZ3_DEBUG_TIMINGS` | OFF | debug timing output |
| `SZ3_INSTALL` | ON, OFF when SZ3 is added with `add_subdirectory` or FetchContent | the install rules |

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

`sz3` compresses when it is given `-i`, `-z`, dimensions and an error bound (`-M` or `-c`). It decompresses when it is
given `-z` and `-o`. With `-i` and `-o` and no `-z`, it does both. The data type defaults to float. Run `sz3 -h` for
every option.

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

`include/SZ3/api/sz.hpp` documents the rest of the API: compressing into your own buffer, `SZ_compress_size_bound`,
and decompressing into a buffer you allocated.

## Algorithms and error-bound modes

Set the algorithm with `Config::cmprAlgo`, or `CmprAlgo` in a configuration file.

| Algorithm | Use it for |
|---|---|
| `ALGO_INTERP_LORENZO` (default) | Most data. It tunes interpolation and Lorenzo prediction on a sample of the data and keeps whichever is better. |
| `ALGO_INTERP` | Smooth data, when you want interpolation without the tuning step. |
| `ALGO_LORENZO_REG` | Blockwise Lorenzo and regression prediction, the SZ2 algorithm. |
| `ALGO_NOPRED` | Quantization without prediction: a fast baseline. |
| `ALGO_LOSSLESS` | Zstd only. SZ3 also switches to it by itself when the error bound is 0, or when Zstd alone gives a smaller result. |
| `ALGO_BIOMD`, `ALGO_BIOMDXTC` | Molecular-dynamics coordinates. `ALGO_BIOMDXTC` follows GROMACS's xtc and can, like xtc, round a coordinate slightly past the bound. |

Set the error-bound mode with `Config::errorBoundMode`, or `-M` on the command line.

| Mode | Bound |
|---|---|
| `EB_ABS` (`ABS`) | Every value is within `absErrorBound` of the original. |
| `EB_REL` (`REL`) | Within `relErrorBound` × (max − min) of the data. |
| `EB_ABS_AND_REL` / `EB_ABS_OR_REL` | The smaller / the larger of the two bounds above. |
| `EB_PSNR` (`PSNR`), `EB_L2NORM` (`NORM`) | A target PSNR or L2-norm error. SZ3 converts it to an absolute bound, which is the pointwise guarantee. |

## Data format and compatibility

Compressed data starts with a magic number and its data-format version, and is little-endian on every host. SZ3 3.4.0
writes data format 3.4.0. It reads 3.4.0 data and 3.3.2 data, except 3.3.2 data compressed with `ALGO_NOPRED` or `ALGO_BIOMD`.
Data from earlier versions has to be decompressed by the version that wrote it; SZ3 refuses it and, for data from 3.2.0 on,
names that version. `sz3 -v` prints the data-format version a build writes.

From 3.4.0 on, builds with and without fused multiply-add (FMA) instructions decode the same data to the same values.
Data compressed by an earlier build that used FMA (Apple Silicon, aarch64, or x86 built for a specific CPU, such as
with `-march=native`) should be decompressed by that build. Building SZ3 with `-ffast-math` is not supported, because
it lets the compiler reorder arithmetic that compression and decompression must repeat exactly.

[CHANGELOG.md](CHANGELOG.md) lists which version changed the format.

## Citing SZ3

[//]: # (**Kindly note**: If you mention SZ3 in your paper, the most appropriate citation is to include these three references &#40;**TBD22, ICDE21, Bigdata18**&#41; because they cover the design and implementation of the latest version of SZ.)
* QOZv2 (the enhanced interpolation-based algorithm): [High-performance Effective Scientific Error-bounded Lossy Compression with Auto-tuned Multi-component Interpolation](https://dl.acm.org/doi/10.1145/3639259).
* SZ3's interpolation-based algorithm: [Optimizing Error-Bounded Lossy Compression for Scientiﬁc Data by Dynamic Spline Interpolation](https://ieeexplore.ieee.org/document/9458791).
* The software engineering design of SZ3: [SZ3: A modular framework for composing prediction-based error-bounded lossy compressors](https://ieeexplore.ieee.org/abstract/document/9866018).

## Version history

See [CHANGELOG.md](CHANGELOG.md).

## License and contact

(C) 2016 by Mathematics and Computer Science (MCS), Argonne National Laboratory. SZ3 is released under a BSD license;
see [copyright-and-BSD-license.txt](copyright-and-BSD-license.txt).
`include/SZ3/encoder/XtcBasedEncoder.hpp` is based on GROMACS and licensed under the LGPL, version 2.1 or later.
The vendored Zstd in `tools/zstd` keeps its own license.

* Lead developer and maintainer: Kai Zhao
* Contributors: Arham Khan, Xin Liang, Robert Underwood, Jinyang Liu, and [everyone else on GitHub](https://github.com/szcompressor/SZ3/graphs/contributors)
* SZ project lead: Franck Cappello

Report bugs and ask questions in [GitHub issues](https://github.com/szcompressor/SZ3/issues).
