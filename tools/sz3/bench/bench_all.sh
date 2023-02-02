#!/bin/bash

exefolder=$1

./bench.sh $exefolder ~/data/nyx-512x512x512/temperature.dat 3 512 512 512
./bench.sh $exefolder ~/data/nyx-512x512x512/velocity_x.dat 3 512 512 512
#./bench.sh $exefolder ~/data/miranda-256x384x384 velocityx.f32.dat  3 256 384 384
#./bench.sh $exefolder ~/data/hurricane-100x500x500 Uf48.bin.dat 3 100 500 500
#./bench.sh $exefolder ~/data/scale-98x1200x1200 QV-98x1200x1200.dat 3 98 1200 1200
#./bench.sh $exefolder ~/data/cesm-26x1800x3600 V_1_26_1800_3600.dat 3 26 1800 3600
#./bench.sh $exefolder ~/data/qmcpack-288x115x69x69 einspline_288_115_69_69.pre.f32.dat 3 33120 69 69
#
