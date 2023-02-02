#!/bin/bash

exefolder=$1
dataset=$2
dim=$3
r1=$4
r2=$5
r3=$6

echo $dataset $dim $r1 $r2 $r3

~/code/sz2/build/bin/sz -z -f -i $dataset -$dim $r1 $r2 $r3 -M ABS -A 1e-2

$exefolder/bench_iterator_sz2 $dataset $dim $r1 $r2 $r3
$exefolder/bench_iterator_sz2 $dataset $dim $r1 $r2 $r3

#$exefolder/bench_iterator_sz3 $dataset $dim $r1 $r2 $r3
#$exefolder/bench_iterator_sz3 $dataset $dim $r1 $r2 $r3

$exefolder/bench_iterator_hybrid1 $dataset $dim $r1 $r2 $r3
$exefolder/bench_iterator_hybrid1 $dataset $dim $r1 $r2 $r3

$exefolder/bench_iterator_hybrid2 $dataset $dim $r1 $r2 $r3
$exefolder/bench_iterator_hybrid2 $dataset $dim $r1 $r2 $r3

$exefolder/bench_iterator_hybrid3 $dataset $dim $r1 $r2 $r3
$exefolder/bench_iterator_hybrid3 $dataset $dim $r1 $r2 $r3

$exefolder/bench_iterator_hybrid4 $dataset $dim $r1 $r2 $r3
$exefolder/bench_iterator_hybrid4 $dataset $dim $r1 $r2 $r3

$exefolder/bench_iterator_hybrid5 $dataset $dim $r1 $r2 $r3
$exefolder/bench_iterator_hybrid5 $dataset $dim $r1 $r2 $r3

$exefolder/bench_iterator_hybrid6 $dataset $dim $r1 $r2 $r3
$exefolder/bench_iterator_hybrid6 $dataset $dim $r1 $r2 $r3

echo ''