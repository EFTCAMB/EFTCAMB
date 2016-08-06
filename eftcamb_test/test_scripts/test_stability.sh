#
# Script that runs the stability sampler for some models in EFTCAMB
#

#!/bin/bash

# go to the eftcamb root folder:
cd ..

# compile the stability sampler executable:
make clean
make eftstability

# run it:
./eftstability tests_EFT/results/Stability_Results

# clean up:
make clean
rm eftstability
