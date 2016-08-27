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
# This file contains the paths of the EFTCAMB test suite.
#
# Developed by: Marco Raveri (mraveri@sissa.it) for the EFTCAMB code
#

# get the path of the script:
SCRIPT_PATH="`dirname \"$0\"`"                  # relative
SCRIPT_PATH="`( cd \"$SCRIPT_PATH\" && pwd )`"  # absolutized and normalized
if [ -z "$SCRIPT_PATH" ] ; then
  exit 1
fi

# get the path of the things in the test suite:

# test directory:
TEST_DIR=$SCRIPT_PATH/..

# camb directory:
CAMB_DIR=$TEST_DIR/..

# test directory:
TEST_PARAMS_DIR=$TEST_DIR/parameters
TEST_PYTHON_DIR=$TEST_DIR/python
TEST_RESULTS_DIR=$TEST_DIR/results
TEST_LEGACY_DIR=$TEST_DIR/eftcamb_legacy

# results directory:
RESULTS_BENCHMARK_DIR=$TEST_RESULTS_DIR/Benchmark_Results
RESULTS_LEGACY_PLOT_DIR=$TEST_RESULTS_DIR/Legacy_Spectra_Plot
RESULTS_LOGS=$TEST_RESULTS_DIR/Logs
RESULTS_PROFILE=$TEST_RESULTS_DIR/Profile
RESULTS_SPECTRA_PLOTS=$TEST_RESULTS_DIR/Spectra_Plots
RESULTS_SPECTRA_RAW=$TEST_RESULTS_DIR/Spectra_results
