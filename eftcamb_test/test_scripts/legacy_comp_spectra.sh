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
# Script that compares the output of the new version of the code with the results present
# in the legacy output folder.
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

# path of things:
PLOTTER=$TEST_PYTHON_DIR/compare_Cls.py          # plotter
PATH_TO_DATA_1=$TEST_LEGACY_DIR/spectra_results  # legacy reference results
PATH_TO_DATA_2=$RESULTS_SPECTRA_RAW              # spectra results
PATH_TO_RES=$RESULTS_LEGACY_PLOT_DIR/            # plot results, if needed.

# decide script options:
NUMDIFF_OPTIONS="-V --absolute-tolerance=0.001 --relative-tolerance=0.001"

# start the script:
printf "${Green}********************************************\n"
printf "EFTCAMB TEST: comparing results with legacy\n"
printf "********************************************\n${Color_Off}"

printf "\n"

# test diagnostic:
FAILED_TEST=()
FAILURE_REASON=()
TOTAL_TEST=0

# start the loop:
for i in $PATH_TO_DATA_1/*_params.ini;
	do
		
	filename=$(basename "$i")
	extension="${filename##*.}"
	filename="${filename%_params.*}"
		
	printf "  Doing %s: " "$filename"
	TOTAL_TEST=$((TOTAL_TEST+1))	
	
	# test with numdiff:
	SUCCESS=true
	for file in $PATH_TO_DATA_1/$filename*;
	do 
		
		temp_filename=$(basename "$file")
		temp_extension="${temp_filename##*.}"
		log_filename="${temp_filename%.*}"
		
		# test if file exists:
		if [ -e "$PATH_TO_DATA_2/$temp_filename" ];
		then
			
			# file exists, run numdiff:
			numdiff $file $PATH_TO_DATA_2/$temp_filename $NUMDIFF_OPTIONS &> /dev/null
		
			if [ $? -ne 0 ]; then
				# if it is the first test failed go to new line and add failure reason:
				if [ "$SUCCESS" = true ]; then
					printf "\n"
					FAILURE_REASON+=('diff failure')
				fi
				# give feedback:
		    	printf "${BRed}   DIFF FAIL: %s\n${Color_Off}" "$temp_filename"
		    	SUCCESS=false
		    	# redo numdiff to print the fail log:
				numdiff $file $PATH_TO_DATA_2/$temp_filename $NUMDIFF_OPTIONS &> $PATH_TO_DATA_2/$log_filename.fail
			fi
		
		else
			# file does not exist:
			if [ "$SUCCESS" = true ]; then
					printf "\n"
					FAILURE_REASON+=('file not found')
			fi
			printf "${BRed}   FAIL FILE NOT FOUND: %s\n${Color_Off}" "$temp_filename"
			SUCCESS=false
		fi
		
	done;

	# if the test failed plot the results:
	
	if [ "$SUCCESS" = false ]; then
		
		
		if [ -z ${TRAVIS} ]; then 
			if [ "${FAILURE_REASON[-1]}" = "diff failure" ]; then 
				printf "${BRed}   Plotting results\n${Color_Off}"
				python $PLOTTER $PATH_TO_DATA_1/$filename $PATH_TO_DATA_2/$filename $PATH_TO_RES ${filename}_LEGACY ${filename}_NEW &> /dev/null
			fi
		fi		
		FAILED_TEST+=($filename)
	
	else
		printf "${BGreen} OK\n${Color_Off}"
	fi

done;

printf "\n"
printf "${Green}********************************************\n"
printf "EFTCAMB TEST: legacy report\n"
printf "********************************************\n${Color_Off}"


printf "\n Passed test: "$(( $TOTAL_TEST-${#FAILED_TEST[@]} ))" / "$TOTAL_TEST"\n"

if [ "${#FAILED_TEST[@]}" -gt "0" ]; then
  printf "\n Failed test:\n\n"
  for ind in `seq 1 ${#FAILED_TEST[@]}`
  do
  echo " " ${FAILED_TEST[$ind-1]} ": " ${FAILURE_REASON[$ind-1]}
  done;
  printf "\n"
else
  printf " All test successfully passed.\n\n"
fi

# return the test result:
if [ "${#FAILED_TEST[@]}" -gt "0" ]; then
	exit 1
else
	exit 0
fi
