MMD-SZ(MDZ): A Modular Error-bounded Lossy Compressor Optimized for Molecular Dynamics
=====
(C) 2021 by Mathematics and Computer Science (MCS), Argonne National Laboratory.
See COPYRIGHT in top-level directory.

* Major Authors: Kai Zhao, Sheng Di, Danny Perez 
* Supervisor: Franck Cappello

## Citations

## Installation

* mkdir build && cd build
* cmake .. -DCMAKE_INSTALL_PREFIX:PATH=[INSTALL_DIR]
* make
* make install

Then, you'll find all the executables in [INSTALL_DIR]/tools/mdz and header files in [INSTALL_DIR]/include

## Testing Examples
build/test/mdz datafile -2 dim1 dim2 -r reb buffer_size cmpr_opt
#### options:
* datafile: FP32 binary format. Contains single axis (X or Y or Z) only.
* dim1: number of timesteps
* dim2: number of atoms 
* reb: relative error bound, for example, 1E-3
* buffer_size: default 10
* cmpr_opt: cmpr_opt<=0 means manually choose a compressor from 0: VQ, -1:VQT, -2:MT, -3: Lorenzo+Regression; cmpr_opt>0 controls the interval to automatically update the best compressor 

#### examples:
* mdz helium-mode-b-7852x1037/x.f32.dat -2 7852 1037 -r 1E-3 10
* mdz helium-mode-b-7852x1037/x.f32.dat -2 7852 1037 -r 1E-3 10 0
* mdz helium-mode-b-7852x1037/x.f32.dat -2 7852 1037 -r 1E-3 10 -1
* mdz helium-mode-b-7852x1037/x.f32.dat -2 7852 1037 -r 1E-3 10 -2
* mdz helium-mode-b-7852x1037/x.f32.dat -2 7852 1037 -r 1E-3 10 10
