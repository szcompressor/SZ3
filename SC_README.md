Prerequisites
===
* gcc >=9
* cmake >= 3.13

How to build
===
* mkdir -p build
* cd build
* cmake ..
* make

How to Run
===
Datasets can be download from:
https://drive.google.com/drive/folders/1DNG0junSINWMeDSQOaWIsHVYTbOdb9cf?usp=sharing

Results of the SC env script: 
https://www.cs.ucr.edu/~kzhao016/sc21_pap138_envs.txt

cd project root dir
build/test/sz_md datafile -2 dim1 dim2 -r reb snapshot_length compressor selection_iteration
options:
dim1, dim2: the dimensions are shown in the dataset folder name, for example, 'helium-mode-b-7852x1037'-> dim1 is  7852, dim2 is 1037
reb: relative error bound, 1E-3
snapshot_length: 10
compressor: 0: VQ, 1:VQT, 2:MT, 2: ADP 
selection_iteration: 50 for ADP and -1 for others 

examples:
VQ
build/test/sz_md helium-mode-b-7852x1037/x.f32.dat -2 7852 1037 -r 1E-3 10 0 -1
VQT
build/test/sz_md helium-mode-b-7852x1037/x.f32.dat -2 7852 1037 -r 1E-3 10 1 -1
MT
build/test/sz_md helium-mode-b-7852x1037/x.f32.dat -2 7852 1037 -r 1E-3 10 2 -1
ADP
build/test/sz_md helium-mode-b-7852x1037/x.f32.dat -2 7852 1037 -r 1E-3 10 2 50