#!/bin/bash

exefolder=$1
data=$2
dim=$3
r1=$4
r2=$5
r3=$6

echo $data $dim $r1 $r2 $r3

~/code/sz2/build/bin/sz -z -f -i $data -$dim $r1 $r2 $r3 -M ABS -A 1e-2

$exefolder/bench_iterator_sz2 $data $dim $r1 $r2 $r3
$exefolder/bench_iterator_sz2 $data $dim $r1 $r2 $r3


#$exefolder/bench_iterator_sz3 $data $dim $r1 $r2 $r3
#$exefolder/bench_iterator_sz3 $data $dim $r1 $r2 $r3

$exefolder/bench_iterator_hybrid1 $data $dim $r1 $r2 $r3
$exefolder/bench_iterator_hybrid1 $data $dim $r1 $r2 $r3
#
#$exefolder/bench_iterator_hybrid2 $data $dim $r1 $r2 $r3
#$exefolder/bench_iterator_hybrid2 $data $dim $r1 $r2 $r3
#
#$exefolder/bench_iterator_hybrid3 $data $dim $r1 $r2 $r3
#$exefolder/bench_iterator_hybrid3 $data $dim $r1 $r2 $r3
#
#$exefolder/bench_iterator_hybrid4 $data $dim $r1 $r2 $r3
#$exefolder/bench_iterator_hybrid4 $data $dim $r1 $r2 $r3
#
#$exefolder/bench_iterator_hybrid5 $data $dim $r1 $r2 $r3
#$exefolder/bench_iterator_hybrid5 $data $dim $r1 $r2 $r3
#
#$exefolder/bench_iterator_hybrid6 $data $dim $r1 $r2 $r3
#$exefolder/bench_iterator_hybrid6 $data $dim $r1 $r2 $r3

echo ''