#!/bin/bash

exefolder=$1
dataset=$2
dim=$r3
r1=$4
r2=$5
r3=$6

echo $dataset $dim $r1 $r2 $r3
$exefolder/bench_iterator_sz2 $dataset 3 $dim $r1 $r2 $r3
$exefolder/bench_iterator_sz2 $dataset 3 $dim $r1 $r2 $r3

$exefolder/bench_iterator_sz3 $dataset 3 $dim $r1 $r2 $r3
$exefolder/bench_iterator_sz3 $dataset 3 $dim $r1 $r2 $r3

$exefolder/bench_iterator_hybrid1 $dataset 3 $dim $r1 $r2 $r3
$exefolder/bench_iterator_hybrid1 $dataset 3 $dim $r1 $r2 $r3

$exefolder/bench_iterator_hybrid2 $dataset 3 $dim $r1 $r2 $r3
$exefolder/bench_iterator_hybrid2 $dataset 3 $dim $r1 $r2 $r3

$exefolder/bench_iterator_hybrid3 $dataset 3 $dim $r1 $r2 $r3
$exefolder/bench_iterator_hybrid3 $dataset 3 $dim $r1 $r2 $r3

$exefolder/bench_iterator_hybrid4 $dataset 3 $dim $r1 $r2 $r3
$exefolder/bench_iterator_hybrid4 $dataset 3 $dim $r1 $r2 $r3

$exefolder/bench_iterator_hybrid5 $dataset 3 $dim $r1 $r2 $r3
$exefolder/bench_iterator_hybrid5 $dataset 3 $dim $r1 $r2 $r3

$exefolder/bench_iterator_hybrid6 $dataset 3 $dim $r1 $r2 $r3
$exefolder/bench_iterator_hybrid6 $dataset 3 $dim $r1 $r2 $r3

echo ''