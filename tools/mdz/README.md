MMD-SZ(MDZ): A Modular Error-bounded Lossy Compressor Optimized for Molecular Dynamics
=====
(C) 2021 by Mathematics and Computer Science (MCS), Argonne National Laboratory.
See COPYRIGHT in top-level directory.

* Major Authors: Kai Zhao, Sheng Di, Danny Perez 
* Supervisor: Franck Cappello

## Citations
* Kai Zhao, Sheng Di, Danny Perez, Zizhong Chen, and Franck Cappello. "[MDZ: An Efficient Error-bounded Lossy Compressor for Molecular Dynamics Simulations](https://ieeexplore.ieee.org/document/9835212)", Proceeding of the 38th IEEE International Conference on Data Engineering (ICDE 22), Kuala Lumpur, Malaysia, May 9 -
  12, 2022.
 
## Installation

Build SZ3 with the cmake option "-DBUILD_MDZ=ON"
`mdz` and `mdz_smoke_test` are installed in `[INSTALL_DIR]/bin` and the header files in `[INSTALL_DIR]/include`.

`mdz` compresses and decompresses in memory and prints the compression ratio, the times and the error; it does not write a compressed file.

## Testing Examples
mdz datafile -2 dim1 dim2 -r reb buffer_size compressor

mdz datafile -3 dim1 dim2 dim3 -r reb buffer_size compressor

#### options:
* datafile: FP32 binary format. Contains single axis (X or Y or Z) only.
* dim1: number of timesteps
* dim2: number of atoms
* dim3: number of atom dimensions (x,y,z,etc.)
* reb: relative error bound, for example, 1E-3; `-a <bound>` instead of `-r <reb>` takes an absolute error bound
* buffer_size (optional): default 10
* compressor (optional): -1: ADP (chooses one of the others), 0: VQ, 1: VQT, 2: MT, 3: Lorenzo+Regression, 4: TS

#### examples:
* mdz helium-mode-b-7852x1037/x.f32.dat -2 7852 1037 -r 1E-3
* mdz helium-mode-b-7852x1037/x.f32.dat -2 7852 1037 -r 1E-3 10
* mdz helium-mode-b-7852x1037/xyz.f32.dat -3 7852 1037 3 -r 1E-3 10
