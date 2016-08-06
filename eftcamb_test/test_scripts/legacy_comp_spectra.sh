#
# Script that compares the output of the new version of the code with the results present
# in the legacy output folder.
#

#!/bin/bash

TEST_FOLDER=./tests_EFT

NPROCS=4  # number of parallel jobs, 4 is pretty fair.

PLOTTER=${TEST_FOLDER}/python/compare_Cls.py
PATH_TO_DATA_1=${TEST_FOLDER}/results_legacy/Spectra_Legacy
PATH_TO_DATA_2=${TEST_FOLDER}/results/Spectra_results
PATH_TO_RES=${TEST_FOLDER}/results/Legacy_Spectra_Plot/

cd ..

echo 'Plotting legacy results results. This might take a while.'
echo 'Plotting spectra.'

# define a counter to give the wait commands every NPROCS commands
var=0
# start the loop:
for file in $PATH_TO_DATA_1/*_params.ini
	do
	
	# get the file name:
	fname=$(basename $file)
	fname=${fname%_params.*}
	# increment the counter
	var=$((var+1))
	# call the plotter
	python $PLOTTER $PATH_TO_DATA_1/$fname $PATH_TO_DATA_2/$fname $PATH_TO_RES ${fname}_LEGACY ${fname}_NEW &
	# wait if reached NPROCS jobs running
	if [ $((var%NPROCS)) -eq 0 ] ; then 
		wait
	fi
	
done;
# wait for all jobs to complete:
wait

echo 'Legacy plotting done.'

exit 0