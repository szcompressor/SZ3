Toward Quantity-of-Interest Preserving Lossy Compression for Scientific Data
## Introduction

This is the source code of SZ-QoI-Error-Control. 
## Dependencies

Please Install the following dependencies before running the evaluation experiments:

* cmake>=3.13
* gcc>=8.0
* sz-qoi-error-control: (from https://github.com/szcompressor/SZ3/tree/qoi_error_control)

## 3rd party libraries/tools

* Zstd >= 1.3.5 (https://facebook.github.io/zstd/). Not mandatory to be mannually installed as Zstandard v1.4.5 is included and will be used if libzstd can not be found by
  pkg-config.

## Installation

* mkdir build && cd build
* cmake -DCMAKE_INSTALL_PREFIX:PATH=[INSTALL_DIR] ..
* make
* make install

Then, you'll find all the executables in [INSTALL_DIR]/bin and header files in [INSTALL_DIR]/include. A Cmake version >= 3.13.0 is needed and we recommend to use gcc version 8.x to compile the code. 

## Single compression/decompression testing examples

You can use the executable 'sz' command to do the compression/decompression. Just run "sz" command to check the instructions for its arguments.
Currently you need to add a configuration file to the argument line (-c) to ebable QoI error control. 

## Evaluation guides

Step 1: Download the dataset from the following links,then unzip them:

* Hurricane-ISABEL: https://g-8d6b0.fd635.8443.data.globus.org/ds131.2/Data-Reduction-Repo/raw-data/Hurricane-ISABEL/SDRBENCH-Hurricane-ISABEL-100x500x500_log.tar.gz
* NYX: https://g-8d6b0.fd635.8443.data.globus.org/ds131.2/Data-Reduction-Repo/raw-data/EXASKY/NYX/SDRBENCH-EXASKY-NYX-512x512x512_log.tar.gz
* SCALE-LETKF: https://g-8d6b0.fd635.8443.data.globus.org/ds131.2/Data-Reduction-Repo/raw-data/SCALE_LETKF/SDRBENCH-SCALE_98x1200x1200_log.tar.gz
* QMCPACK: https://g-8d6b0.fd635.8443.data.globus.org/ds131.2/Data-Reduction-Repo/raw-data/QMCPack/SDRBENCH-QMCPack.tar.gz
* Data download command: wget {DataLink} --no-check-certificate (DataLink: The link of data)
* Data unzip command: tar -zxvf {DataFile} (DataFile: the .tar.gz file of downloaded data)

Step 2: Run test on different QoIs

This is achieved with sz.config file. This file can be found under folder ./test. With different QoI, we need to tell the compressor what is the the target QoI as well as the target tolrence. 

For example, at the end of `test/sz.config`, modify the follwing lines to test $ x^2 $ whith its a tolrence of $1E-1$

`qoi=1` # index of QoI, see the table below for a full list of supported QoIs. 

`qoiEB=1E-1`  # This is a value-range based relative value.

To run a case, for example `Uf48.bin.f32` in Hurricane dataset.


`sz -f -c sz.config -i Uf48.bin.f32 -o test.out -3 500 500 100 -M REL 8E-1 -a`

In this command line, we set the input error bound as 8 times of the target QoI error bound (`qoiEB`). This relation was determined by a series of tuning based on real dataset tests based on the objective that "the resulting QoI error should be as close to the given target `qoiEB` as possible". You may change this value to someting else larger than your `qoiEB`. However, if the input error bound is set to be smaller than `qoiEB`, the "objective" can not be guranteed. 



Tabld of QoI index 

| Index | Description | Config |
| --- | ---- |----| 
| 1 | $x^2$: square of single value | `qoi=1` |
| 2 | $\log(x)$:logrithm of single value(base 2) | `qoi=2`
| 3 | $\frac{1}{k} \displaystyle\sum_{i=1}^{i=k} x^i $: regional average of $x_i^2$| `qoi=3`, `qoiRegionSize=k`|
| 4 | Isoline(specify the number of isolines/isosuraces. By defalt, only extract $1$ isovalue, the mean of data.)| `qoi=4`, `qoiIsoNum=k`|
| 5 | $x^2$ and $\log(x)$ | `qoi=5`|
| 6 | $x^2$ and isoline | `qoi=6`, `qoiIsoNum=k`|
| 7 | $\log(x)$ and isoline | `qoi=7`, `qoiIsoNum=k`|
| 8 | $x^2$, $\log(x)$ and isoline | `qoi=8`, `qoiIsoNum=k`|



