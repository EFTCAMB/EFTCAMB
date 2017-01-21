#----------------------------------------------------------------------------------------
#
# This file is part of EFTCAMB.
#
# Copyright (C) 2013-2017 by the EFTCAMB authors
#
# The EFTCAMB code is free software;
# You can use it, redistribute it, and/or modify it under the terms
# of the GNU General Public License as published by the Free Software Foundation;
# either version 3 of the License, or (at your option) any later version.
# The full text of the license can be found in the file eftcamb/LICENSE at
# the top level of the EFTCAMB distribution.
#
#----------------------------------------------------------------------------------------

#
# Script that plots all the spectra obtained with the spectra script.
#
# Developed by: Marco Raveri (mraveri@sissa.it) for the EFTCAMB code
#

#!/bin/bash

# get the path of the script:
SCRIPT_PATH="`dirname \"$0\"`"                  # relative
SCRIPT_PATH="`( cd \"$SCRIPT_PATH\" && pwd )`"  # absolutized and normalized
if [ -z "$SCRIPT_PATH" ] ; then
  exit 1
fi

# import things:
source $SCRIPT_PATH/paths.sh
source $SCRIPT_PATH/colors.sh

# get the number of cores. Plotting is slow so we need to parallelize it
NPROCS=$(nproc --all)

# definitions:
PLOTTER=$TEST_PYTHON_DIR/compare_Cls.py    # plotter
PATH_TO_DATA=$RESULTS_SPECTRA_RAW          # data, spectra results
PATH_TO_RES=$RESULTS_SPECTRA_PLOTS/        # path to results place

REFERENCE_MODEL=$RESULTS_SPECTRA_RAW/1_GR  # model to use as a reference

# start the script:
printf "${Green}********************************************\n"
printf "EFTCAMB TEST: plotting test results\n"
printf "This might take a while.\n"
printf "********************************************\n${Color_Off}"

printf "\n"

# go to the CAMB directory:
cd $CAMB_DIR

# define a counter to give the wait commands every NPROCS commands
var=0

# start the loop:
for i in ${PATH_TO_DATA}/*_params.ini
	do
		
	filename=$(basename "$i")
	extension="${filename##*.}"
	filename="${filename%_params.*}"
	
	# increment the counter
	var=$((var+1))
	
	# call the plotter
	python $PLOTTER $REFERENCE_MODEL $PATH_TO_DATA/$filename $PATH_TO_RES '1_GR' $filename &
	# wait if reached NPROCS jobs running
	if [ $((var%NPROCS)) -eq 0 ] ; then 
		wait
	fi

done;
# wait for all jobs to complete:
wait

# feedback for the run:
printf "\n"
printf "Results in: %s \n" "$PATH_TO_RES"
printf "\n"

printf "${Green}********************************************\n"
printf "EFTCAMB TEST: done spectra plotting\n"
printf "********************************************\n${Color_Off}"

exit 0
