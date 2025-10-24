SZ3: A Modular Error-bounded Lossy Compression Framework for Scientific Datasets
=====
(C) 2016 by Mathematics and Computer Science (MCS), Argonne National Laboratory. See COPYRIGHT in the top-level directory.

* Major Authors: Sheng Di, Kai Zhao, Xin Liang, Jinyang Liu
* Supervisor: Franck Cappello
* Other Contributors: Robert Underwood, Sihuan Li, Ali M. Gok

## Installation

* mkdir build && cd build
* cmake -DCMAKE_INSTALL_PREFIX:PATH=[INSTALL_DIR] ..
* make
* make install

Then, you'll find all the executables in [INSTALL_DIR]/bin and header files in [INSTALL_DIR]/include


## How to run

#### SZ3 Executable
* You can use the executable 'tools/sz3/sz3' to do the compression/decompression.

#### SZ3 C++ API
* Located in 'include/SZ3/api/sz.hpp'. 
* Requiring a modern C++ compiler.  
* Different with SZ2 API.

#### SZ3 C API
* Located in 'tools/sz3c/include/sz3c.h'
* Compatible with SZ2 API

#### SZ3 Python API
* available via `pip install pysz`
* [Source code in 'tools/pysz'](https://github.com/szcompressor/SZ3/tree/master/tools/pysz)

#### H5Z-SZ3
* Located in 'tools/H5Z-SZ3'
* Please add "-DBUILD_H5Z_FILTER=ON" to enable this function for CMake.
* sz3ToHDF5 and HDF5ToSz3 are provided for testing.

#### Third-Party APIs

* [SZ3 Fortran API](https://github.com/ofmla/sz3_simple_example) (by [Oscar Mojica](https://github.com/ofmla))
* [SZ3 Rust API](https://github.com/juntyr/sz3-rs) (by [Juniper Tyree](https://github.com/juntyr))
* [SZ3 Numcodecs API](https://github.com/juntyr/numcodecs-rs/blob/main/codecs/sz3/) (by [Juniper Tyree](https://github.com/juntyr))


## Citations
[//]: # (**Kindly note**: If you mention SZ3 in your paper, the most appropriate citation is to include these three references &#40;**TBD22, ICDE21, Bigdata18**&#41; because they cover the design and implementation of the latest version of SZ.)
* QOZv2 (the enhanced interpolation-based algorithm): [High-performance Effective Scientific Error-bounded Lossy Compression with Auto-tuned Multi-component Interpolation](https://dl.acm.org/doi/10.1145/3639259).
* SZ3's interpolation-based algorithm: [Optimizing Error-Bounded Lossy Compression for Scientiﬁc Data by Dynamic Spline Interpolation](https://ieeexplore.ieee.org/document/9458791).
* The software engineering design of SZ3: [SZ3: A modular framework for composing prediction-based error-bounded lossy compressors](https://ieeexplore.ieee.org/abstract/document/9866018).


## Version history

Version New features

* SZ 3.0.0 SZ3 is the C++ version of SZ with a modular and composable design.
* SZ 3.0.1 Improve the build process.
* SZ 3.1.0 The default algorithm is now interpolation+Lorenzo.
* SZ 3.1.1 Add OpenMP support. Works for all algorithms. Please enable it using the config file. 
* SZ 3.1.2 Support configuration file (INI format). An example can be found in 'tools/sz3/sz3.config'.
* SZ 3.1.3 Support more error control mode: PSNR, L2Norm, ABS_AND_REL, ABS_OR_REL. Support INT32 and INT64 datatype.
* SZ 3.1.4 Support running on Windows natively with Visual Studio. Please use CMake to generate Visual Studio solution files.
* SZ 3.1.5 Support HDF5 by H5Z-SZ3. Please add "-DBUILD_H5Z_FILTER=ON" to enable this function for CMake.
* SZ 3.1.6 Support C API and Python API.
* SZ 3.1.7 Initial MDZ(https://github.com/szcompressor/SZ3/tree/master/tools/mdz) support.
* SZ 3.1.8 namespace changed from SZ to SZ3. H5Z-SZ3 supports configuration files now.
* SZ 3.2.0 API reconstructed for FZ. H5Z-SZ3 rewrite. Compression version checking.
* SZ 3.3.0 Add key QoZ v1 and v2 features to improve compression speed and data quality. The full QoZ is available from **a separate branch** (https://github.com/szcompressor/SZ3/tree/QoZ). 
* SZ 3.3.1: SZ3 Windows support for both Visual Studio and MinGW toolchains. pySZ v1 released and available via `pip install pysz`.

## 3rd party libraries/tools
* [Zstandard](https://facebook.github.io/zstd/) v1.4.5 will be fetched if libzstd can not be found by pkg-config.
* The source code of ska_hash is included in SZ3.