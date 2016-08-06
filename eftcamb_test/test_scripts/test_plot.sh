#
# Script that plots all the spectra obtained with the spectra script.
#

#!/bin/bash

TEST_FOLDER=./tests_EFT

# do the spectra plotting:

NPROCS=4  # number of parallel jobs, 4 is pretty fair.

# definitions:
PLOTTER=${TEST_FOLDER}/python/compare_Cls.py  # plotter
PATH_TO_DATA=${TEST_FOLDER}/results/Spectra_results  # data, spectra results
PATH_TO_RES=${TEST_FOLDER}/results/Spectra_Plots/  # path to results place

REFERENCE_MODEL=${TEST_FOLDER}/results/Spectra_results/1_GR  # model to use as a reference

# start the script:
echo 'Plotting of test results. This might take a while.'
echo 'Plotting spectra.'
echo ; echo ; echo ;
# go to the eftcamb root dir
cd ..
# define a counter to give the wait commands every NPROCS commands
var=0
# start the loop:
for file in ${PATH_TO_DATA}/*_params.ini
	do
		
	# get the file name:
	fname=$(basename $file)
	fname=${fname%_params.*}
	# increment the counter
	var=$((var+1))
	# call the plotter
	python $PLOTTER $REFERENCE_MODEL $PATH_TO_DATA/$fname $PATH_TO_RES '1_GR' $fname &
	# wait if reached NPROCS jobs running
	if [ $((var%NPROCS)) -eq 0 ] ; then 
		wait
	fi

done;
# wait for all jobs to complete:
wait

echo 'Spectra plotting done.'

exit 0
