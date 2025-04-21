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

## 3rd party libraries/tools

* Zstandard (https://facebook.github.io/zstd/). Zstandard v1.4.5 is included and will be used if libzstd can not be found by
  pkg-config.

## Testing Examples

You can use the executable 'sz3' command to do the compression/decompression.

SZ3 simplifies command line arguments in the previous version. If you are a new user, please follow the instructions
given by the executable.


## Backward Compatibility with SZ2
For backward compatibility, most of the SZ2 command line parameters are supported in SZ3. **Exceptions are listed below**.
Scripts without the parameters below should work fine by replacing SZ2 with SZ3.

| Parameter | Explanation                     | SZ3 roadmap                              |
|-----------|---------------------------------|------------------------------------------|
| -c        | Config file                     | SZ3 has different config format with SZ2 |
| -p        | Print configuration info        | Will be supported soon                   |
| -T        | Tucker Tensor Decomposition     | Will be supported later                  |
| -P        | Point-wise relative error bound | Will be supported later                  |


## API

#### SZ3 C++ API
* Located in 'include/SZ3/api/sz.hpp'. 
* Requiring a modern C++ compiler.  
* Different with SZ2 API.

#### SZ3 C API
* Located in 'tools/sz3c/include/sz3c.h'
* Compatible with SZ2 API

#### Python API
* Located in 'tools/pysz/pysz.py'
* Test file provided ('tools/pysz/test.py')
* Compatible with both SZ3 and SZ2
* Requiring SZ2/3 dynamic library

#### Fortran API
* Special thanks to [Oscar Mojica](https://github.com/ofmla) for providing the Fortran API
* Visit [this Github repository](https://github.com/ofmla/sz3_simple_example) for details

#### H5Z-SZ3
* Located in 'tools/H5Z-SZ3'
* Please add "-DBUILD_H5Z_FILTER=ON" to enable this function for CMake.
* sz3ToHDF5 and HDF5ToSz3 are provided for testing.

[//]: # (* Use examples/print_h5repack_args.c to construct the cd_values parameters based on the specified error configuration.)
[//]: # ()
[//]: # (* Compression example: )
[//]: # (`h5repack -f UD=32024,0,5,0,981668463,0,0,0 -i ~/Data/CESM-ATM-tylor/1800x3600/CLDLOW_1_1800_3600.dat.h5 -o ~/Data/CESM-ATM-tylor/1800x3600/CLDLOW_1_1800_3600.dat.sz3.h5`)
[//]: # ()
[//]: # (* Decompression example:)
[//]: # (`h5repack -f NONE -i ~/Data/CESM-ATM-tylor/1800x3600/CLDLOW_1_1800_3600.dat.sz3.h5 -o ~/Data/CESM-ATM-tylor/1800x3600/CLDLOW_1_1800_3600.dat.sz3.out.h5`)
[//]: # ()
[//]: # (* Alternatively, the error bound information can also be given through sz3.config &#40; when there are no cd_values for h5repack&#41;. Example &#40;You need to put sz3.config in the current local directory so that it will read sz3.config to get error bounds&#41;:)
[//]: # (`h5repack -f UD=32024,0 -i ~/Data/CESM-ATM-tylor/1800x3600/CLDLOW_1_1800_3600.dat.h5 -o ~/Data/CESM-ATM-tylor/1800x3600/CLDLOW_1_1800_3600.dat.sz3.h5`)



## Version history

Version New features

* SZ 3.0.0 SZ3 is the C++ version of SZ with a modular and composable design.
* SZ 3.0.1 Improve the build process.
* SZ 3.1.0 The default algorithm is now interpolation+Lorenzo.
* SZ 3.1.1 Add OpenMP support. Works for all algorithms. Please enable it using the config file. 
* SZ 3.1.2 Support configuration file (INI format). An example can be found in 'tools/sz3/sz3.config'.
* SZ 3.1.3 Support more error control mode: PSNR, L2Norm, ABS_AND_REL, ABS_OR_REL. Support INT32 and INT64 datatype.
* SZ 3.1.4 Support running on Windows. Please refer to https://github.com/szcompressor/SZ3/issues/5#issuecomment-1094039224 for instructions.
* SZ 3.1.5 Support HDF5 by H5Z-SZ3. Please add "-DBUILD_H5Z_FILTER=ON" to enable this function for CMake.
* SZ 3.1.6 Support C API and Python API.
* SZ 3.1.7 Initial MDZ(https://github.com/szcompressor/SZ3/tree/master/tools/mdz) support.
* SZ 3.1.8 namespace changed from SZ to SZ3. H5Z-SZ3 supports configuration files now.
* SZ 3.2.0 API reconstructed for FZ. H5Z-SZ3 rewrite. Compression version checking.
* SZ 3.3.0 Add key QoZ v1 and v2 features to improve compression speed and data quality. The full QoZ is available from **a separate branch** (https://github.com/szcompressor/SZ3/tree/QoZ). 


## Citations
[//]: # (**Kindly note**: If you mention SZ3 in your paper, the most appropriate citation is to include these three references &#40;**TBD22, ICDE21, Bigdata18**&#41; because they cover the design and implementation of the latest version of SZ.)
* SZ3 with the original interpolation-based algorithm: [Optimizing Error-Bounded Lossy Compression for Scientiﬁc Data by Dynamic Spline Interpolation](https://ieeexplore.ieee.org/document/9458791).
* QOZv2 (the enhanced interpolation algorithm): [High-performance Effective Scientific Error-bounded Lossy Compression with Auto-tuned Multi-component Interpolation](https://dl.acm.org/doi/10.1145/3639259).
* The software engineering design of SZ3: [SZ3: A modular framework for composing prediction-based error-bounded lossy compressors](https://ieeexplore.ieee.org/abstract/document/9866018).
