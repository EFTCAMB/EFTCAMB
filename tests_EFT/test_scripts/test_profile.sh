#
# Script that runs the profiler of camb for different choices of parameters
#

#!/bin/bash

# go to the eftcamb root folder:
cd ..

# compile the eftcamb profiler:
make clean
make eftprofile

# run the profiler:
for file in tests_EFT/parameters/*
	do
	
	./eftprofile $file
	gprof -b eftprofile > tests_EFT/results/Profile/$(basename "$file")
	echo ; echo ; echo ;

done

# clean up:
make clean
rm eftprofile
rm gmon.out