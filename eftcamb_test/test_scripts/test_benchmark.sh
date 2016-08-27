#----------------------------------------------------------------------------------------
#
# This file is part of EFTCAMB.
#
# Copyright (C) 2013-2016 by the EFTCAMB authors
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
# Script that runs the benchmarker of EFTCAMB for all the parameters in the eftcamb_test/parameters/ folder
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

# start the script:
printf "${Green}********************************************\n"
printf "EFTCAMB TEST: running benchmarker\n"
printf "********************************************\n${Color_Off}"

# go to the CAMB directory:
cd $CAMB_DIR

printf "\nCompiling EFTCAMB benchmarker: "

make clean &> /dev/null
make eftcamb_apps ifortErr=1 &> $RESULTS_LOGS/benchmarker_log.txt

# check if compilation succeded:
if [ $? -eq 0 ]; then
    printf "${BGreen} OK\n${Color_Off}"
else
    printf "${BRed} FAIL\n${Color_Off}"
    exit 1
fi

printf "\n"

# define flag to check wether all went correct:
SUCCESS=true

# run all the parameters:

for i in $TEST_PARAMS_DIR/*.ini;
	do

	filename=$(basename "$i")
	extension="${filename##*.}"
	filename="${filename%.*}"

	printf "  Doing %s: " "$filename"
	
	./benchmark.x $i &> $RESULTS_BENCHMARK_DIR/$filename.bench
	
	# check if run succeded:
	if [ $? -eq 0 ]; then
		TIMING=$(tail -1 $RESULTS_BENCHMARK_DIR/$filename.bench)
	    printf " %s" "$TIMING"
	    printf "${BGreen} OK\n${Color_Off}"
	else
	    printf "${BRed} FAIL\n${Color_Off}"
	    SUCCESS=false
	fi

done;

printf "\n"


# cleaning up:
make clean        &> /dev/null
rm   eftbenchmark &> /dev/null

# feedback for the run:

printf "Results in: %s \n" "$RESULTS_BENCHMARK_DIR"
printf "\n"

if [ "$SUCCESS" = true ]; then
	    printf "${Green}********************************************\n"
		printf "EFTCAMB TEST: done benchmarking\n"
		printf "********************************************\n${Color_Off}"
else
	    printf "${Yellow}********************************************\n"
		printf "EFTCAMB TEST: done benchmarking\n"
		printf "WARNING: some parameters failed.\n"
		printf "********************************************\n${Color_Off}"
fi

# return status:

if [ "$SUCCESS" = true ]; then
	    exit 0
else
	    exit 1
fi
