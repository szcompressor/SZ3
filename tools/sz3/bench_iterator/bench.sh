#!/bin/bash

exefolder=$1
data=$2
dim=$3
r1=$4
r2=$5
r3=$6

echo "======================"
echo "======================"
echo "======================"
echo $data $dim $r1 $r2 $r3


#$exefolder/bench_iterator_sz2 $data $dim $r1 $r2 $r3
#$exefolder/bench_iterator_sz2 $data $dim $r1 $r2 $r3

#$exefolder/bench_iterator_sz3 $data $dim $r1 $r2 $r3
#$exefolder/bench_iterator_sz3 $data $dim $r1 $r2 $r3

#$exefolder/bench_iterator_hybrid1 $data $dim $r1 $r2 $r3
#$exefolder/bench_iterator_hybrid1 $data $dim $r1 $r2 $r3
##
#$exefolder/bench_iterator_hybrid2 $data $dim $r1 $r2 $r3
#$exefolder/bench_iterator_hybrid2 $data $dim $r1 $r2 $r3
#
#$exefolder/bench_iterator_hybrid3 $data $dim $r1 $r2 $r3
#$exefolder/bench_iterator_hybrid3 $data $dim $r1 $r2 $r3

#$exefolder/bench_iterator_hybrid4 $data $dim $r1 $r2 $r3
#$exefolder/bench_iterator_hybrid4 $data $dim $r1 $r2 $r3

#$exefolder/bench_iterator_hybrid5 $data $dim $r1 $r2 $r3
#$exefolder/bench_iterator_hybrid5 $data $dim $r1 $r2 $r3
#
#$exefolder/bench_iterator_hybrid6 $data $dim $r1 $r2 $r3
#$exefolder/bench_iterator_hybrid6 $data $dim $r1 $r2 $r3


echo "============Lorenzo1=========="
$exefolder/bench_iterator_lorenzo $data $dim $r1 $r2 $r3
$exefolder/bench_iterator_lorenzo $data $dim $r1 $r2 $r3

echo "============Lorenzo4=========="
$exefolder/bench_iterator_lorenzo4 $data $dim $r1 $r2 $r3
$exefolder/bench_iterator_lorenzo4 $data $dim $r1 $r2 $r3

echo "============Lorenzo5=========="
$exefolder/bench_iterator_lorenzo5 $data $dim $r1 $r2 $r3
$exefolder/bench_iterator_lorenzo5 $data $dim $r1 $r2 $r3
#
#echo "===========SZ2=============="
#sz -z -f -i $data -$dim $r3 $r2 $r1 -M ABS -A 1e-2 -c ~/code/sz2/example/sz.config
#sz -x -f -i $data -s $data.sz -$dim $r3 $r2 $r1 -a
#
#echo "============SZ3============="
#sz3 -f -i $data -o /var/tmp/data -$dim $r3 $r2 $r1 -c ~/lorenzo.config -M ABS 1e-2 -a