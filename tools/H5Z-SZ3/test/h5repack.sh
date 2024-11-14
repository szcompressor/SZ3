#!/bin/bash

if [[ $# != 3 ]]
then
	echo "Usage: $0 [SZ3 config file][inputf_hdf5_ile] [output_hdf5_file]"
	echo "Example: $0 sz3.conf testfloat_8_8_128.h5 testfloat_8_8_128_sz.h5"
	exit
fi

configFile=$1
inputFile=$2
outputFile=$3
script_dir="$(dirname "$(readlink -f "$0")")"

repackArgs=$("${script_dir}/print_h5repack_args" -c "$configFile")
echo $repackArgs

echo "Executing h5repack-shared with the following parameters: $repackArgs $inputFile $outputFile"
h5repack-shared $repackArgs $inputFile $outputFile