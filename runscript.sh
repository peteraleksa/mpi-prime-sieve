#!/bin/bash

file=$1
inputsize=$2
maxnumproc=$3
testsize=$4

echo Program Name: $file
echo Input Size: $inputsize
echo Max Processors: $maxnumproc
echo Number of Tests: $testsize

proc=1
while [ $proc -lt $maxnumproc ]; do 
	echo Running $testsize tests on $proc processors...
	for (( j=1; j<=$testsize; j++ ))
		do
	   	echo Test $j on $proc processors:
    		mpirun -np $proc $file $size
		done
	$i++
	$proc=$((2**($i)))
done
