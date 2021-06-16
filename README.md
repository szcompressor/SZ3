MMD-SZ: A Modular Error-bounded Lossy Compressor Optimized for Material Molecular Dynamics
=====
(C) 2021 by Mathematics and Computer Science (MCS), Argonne National Laboratory.
See COPYRIGHT in top-level directory.

* Major Authors: Kai Zhao, Sheng Di, Danny Perez 
* Supervisor: Franck Cappello

## Citations
**Kindly note**: If you mention SZ in your paper, the most appropriate citation is including these four references (***HPDC2020, Bigdata2018, IPDPS2017 and IPDPS2016***), because they cover the whole design and implementation of the latest version of SZ.**

* SZauto: Kai Zhao, Sheng Di, Xin Liang, Sihuan Li, Dingwen Tao, Zizhong Chen, and Franck Cappello. "[Significantly Improving Lossy Compression for HPC Datasets with Second-Order Prediction and Parameter Optimization](https://dl.acm.org/doi/10.1145/3369583.3392688)", Proceedings of the 29th International Symposium on High-Performance Parallel and Distributed Computing (HPDC 20), Stockholm, Sweden, 2020. (code: https://github.com/szcompressor/SZauto/)

* SZ 2.0+: Xin Liang, Sheng Di, Dingwen Tao, Zizhong Chen, Franck Cappello, "[Error-Controlled Lossy Compression Optimized for High Compression Ratios of Scientific Datasets](https://ieeexplore.ieee.org/document/8622520)", in IEEE International Conference on Big Data (Bigdata 2018), Seattle, WA, USA, 2018.

* SZ 1.4.0-1.4.13: Dingwen Tao, Sheng Di, Franck Cappello, "[Significantly Improving Lossy Compression for Scientific Data Sets Based on Multidimensional Prediction and Error-Controlled Quantization](https://ieeexplore.ieee.org/document/7967203)", in IEEE International Parallel and Distributed Processing Symposium (IPDPS 2017), Orlando, Florida, USA, 2017.

* SZ 0.1-1.0: Sheng Di, Franck Cappello, "[Fast Error-bounded Lossy HPC Data Compression with SZ](https://ieeexplore.ieee.org/document/7516069)", in IEEE International Parallel and Distributed Processing Symposium (IPDPS 2016), Chicago, IL, USA, 2016.

* Point-wise relative error bound mode (i.e., PW_REL): Xin Liang, Sheng Di, Dingwen Tao, Zizhong Chen, Franck Cappello, "[An Efficient Transformation Scheme for Lossy Data Compression with Point-wise Relative Error Bound](https://ieeexplore.ieee.org/document/8514879)", in IEEE International Conference on Clustering Computing (CLUSTER 2018), Belfast, UK, 2018. (Best Paper)
## 3rd party libraries/tools
* Zstandard (https://facebook.github.io/zstd/). Zstandard v1.3.5 is included and will be used if libzstd can be found by pkg-config.
## Installation

* mkdir build && cd build
* cmake .. -DCMAKE_INSTALL_PREFIX:PATH=[INSTALL_DIR]
* make
* make install

Then, you'll find all the executables in [INSTALL_DIR]/bin and header files in [INSTALL_DIR]/include

## Testing Examples
build/test/sz_md datafile -2 dim1 dim2 -r reb buffer_size compressor
#### options:
* datafile: FP32 binary format. Contains single axis (X or Y or Z) only.
* dim1: number of timesteps
* dim2: number of atoms 
* reb: relative error bound, for example, 1E-3
* buffer_size: default 10
* compressor: 0: VQ, 1:VQT, 2:MT, default is ADP

#### examples:
* build/test/sz_md helium-mode-b-7852x1037/x.f32.dat -2 7852 1037 -r 1E-3 10
* build/test/sz_md helium-mode-b-7852x1037/x.f32.dat -2 7852 1037 -r 1E-3 10 0
* build/test/sz_md helium-mode-b-7852x1037/x.f32.dat -2 7852 1037 -r 1E-3 10 1
* build/test/sz_md helium-mode-b-7852x1037/x.f32.dat -2 7852 1037 -r 1E-3 10 2
