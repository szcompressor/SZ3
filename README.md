SZ3: A Modular Error-bounded Lossy Compression Framework for Scientific Datasets
=====
(C) 2016 by Mathematics and Computer Science (MCS), Argonne National Laboratory. See COPYRIGHT in top-level directory.

* Major Authors: Sheng Di, Kai Zhao, Xin Liang
* Supervisor: Franck Cappello
* Other Contributors: Robert Underwood, Sihuan Li, Ali M. Gok

## Citations

**Kindly note**: If you mention SZ in your paper, the most appropriate citation is including these three references (**ICDE21, HPDC2020, Bigdata2018**), because they cover the design and implementation of the latest version of SZ.

* SZ3 Algorithm: Kai Zhao, Sheng Di, Maxim Dmitriev, Thierry-Laurent D. Tonellot, Zizhong Chen, and Franck
  Cappello. "[Optimizing Error-Bounded Lossy Compression for Scientiﬁc Data by Dynamic Spline Interpolation](https://ieeexplore.ieee.org/document/9458791)"
  , Proceeding of the 37th IEEE International Conference on Data Engineering (ICDE 21), Chania, Crete, Greece, Apr 19 -
  22, 2021.

* SZauto: Kai Zhao, Sheng Di, Xin Liang, Sihuan Li, Dingwen Tao, Zizhong Chen, and Franck
  Cappello. "[Significantly Improving Lossy Compression for HPC Datasets with Second-Order Prediction and Parameter Optimization](https://dl.acm.org/doi/10.1145/3369583.3392688)"
  , Proceedings of the 29th International Symposium on High-Performance Parallel and Distributed Computing (HPDC 20),
  Stockholm, Sweden, 2020. (code: https://github.com/szcompressor/SZauto/)

* SZ 2.0+: Xin Liang, Sheng Di, Dingwen Tao, Zizhong Chen, Franck
  Cappello, "[Error-Controlled Lossy Compression Optimized for High Compression Ratios of Scientific Datasets](https://ieeexplore.ieee.org/document/8622520)"
  , in IEEE International Conference on Big Data (Bigdata 2018), Seattle, WA, USA, 2018.

* SZ 1.4.0-1.4.13: Dingwen Tao, Sheng Di, Franck
  Cappello. "[Significantly Improving Lossy Compression for Scientific Data Sets Based on Multidimensional Prediction and Error-Controlled Quantization](https://ieeexplore.ieee.org/document/7967203)"
  , in IEEE International Parallel and Distributed Processing Symposium (IPDPS 2017), Orlando, Florida, USA, 2017.

* SZ 0.1-1.0: Sheng Di, Franck
  Cappello. "[Fast Error-bounded Lossy HPC Data Compression with SZ](https://ieeexplore.ieee.org/document/7516069)", in
  IEEE International Parallel and Distributed Processing Symposium (IPDPS 2016), Chicago, IL, USA, 2016.

* Point-wise relative error bound mode (i.e., PW_REL): Xin Liang, Sheng Di, Dingwen Tao, Zizhong Chen, Franck
  Cappello, "[An Efficient Transformation Scheme for Lossy Data Compression with Point-wise Relative Error Bound](https://ieeexplore.ieee.org/document/8514879)"
  , in IEEE International Conference on Clustering Computing (CLUSTER 2018), Belfast, UK, 2018. (Best Paper)

## 3rd party libraries/tools

* Zstandard (https://facebook.github.io/zstd/). Zstandard v1.3.5 is included and will be used if libzstd can be found by
  pkg-config.

## Installation

* mkdir build && cd build
* cmake -DCMAKE_INSTALL_PREFIX:PATH=[INSTALL_DIR] ..
* make
* make install

Then, you'll find all the executables in [INSTALL_DIR]/bin and header files in [INSTALL_DIR]/include

## Testing Examples

You can use the executable 'sz' command to do the compression/decompression.

SZ3 simplifies command line arguments in the previous version. If you are a new user, please follow the instructions
given by the executable.


## Backward Compatibility with SZ2
For backward compatibility, most of the SZ2 command line parameters are supported in SZ3. **Exceptions are listed below**.
Scripts without parameters below should work fine by replacing SZ2 with SZ3.

| Parameter | Explanation                     | SZ3 roadmap                              |
|-----------|---------------------------------|------------------------------------------|
| -c        | Config file                     | SZ3 has different config format with SZ2 |
| -p        | Print configuration info        | Will be supported soon                   |
| -T        | Tucker Tensor Decomposition     | Will be supported later                  |
| -P        | Point-wise relative error bound | Will be supported later                  |


## API

Please refer to 'include/SZ3/api/sz.hpp' for API and instructions. SZ3 API is different with SZ2.

## Version history

Version New features

* SZ 3.0.0 SZ3 is the C++ version of SZ with modular and composable design.
* SZ 3.0.1 Improve the build process.
* SZ 3.0.2 Support point-wise relative error bound mode.
* SZ 3.1.0 The default algorithm is now interpolation+Lorenzo.
* SZ 3.1.1 Add OpenMP support. Works for all algorithms.
* SZ 3.1.2 Support configuration file (INI format). Example can be found in 'test/sz.config'.
* SZ 3.1.3 Support more error control mode: PSNR, L2Norm, ABS_AND_REL, ABS_OR_REL. Support INT32 and INT64 datatype. 
