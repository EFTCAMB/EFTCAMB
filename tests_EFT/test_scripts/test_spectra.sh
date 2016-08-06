#
# Script that runs camb with all the parameters in the tests_EFT/parameters/ folder
#

#!/bin/bash

# go to the eftcamb root folder:
cd ..

# compile eftcamb:
make clean
make camb

# create all the spectra by cycling over all the parameter files:
for file in tests_EFT/parameters/*
	do
	
	./camb $file | tee -a tests_EFT/results/Spectra_results/spectra.log
	echo ; echo ; echo ;

done

# clean up:
make clean
rm camb

