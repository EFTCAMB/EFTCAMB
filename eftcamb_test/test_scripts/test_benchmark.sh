#
# Script that runs the benchmark of camb with different parameters
#

#!/bin/bash

# go to the eftcamb root folder:
cd ..

# compile the eftcamb benchmarker:
make clean
make eftbenchmark

echo

# run benchmarks for all the parameter files:
for file in tests_EFT/parameters/*

	do
	./eftbenchmark $file | tee -a tests_EFT/results/Benchmark_Results/benchmark.log
	echo ; 

done

# clean up:
make clean
rm eftbenchmark