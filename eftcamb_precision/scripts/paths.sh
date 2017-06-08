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
# This file contains the paths of the EFTCAMB example suite.
#
# Developed by: Marco Raveri (mraveri@sissa.it) for the EFTCAMB code
#

# get the path of the script:
SCRIPT_PATH="`dirname \"$0\"`"                  # relative
SCRIPT_PATH="`( cd \"$SCRIPT_PATH\" && pwd )`"  # absolutized and normalized
if [ -z "$SCRIPT_PATH" ] ; then
  exit 1
fi

# get the path of the things in the example suite:

# example directory:
EXAMPLE_DIR=$SCRIPT_PATH/..

# camb directory:
CAMB_DIR=$EXAMPLE_DIR/..

# example directory:
EXAMPLE_PARAMS_DIR=$EXAMPLE_DIR/parameters
EXAMPLE_PYTHON_DIR=$EXAMPLE_DIR/python
EXAMPLE_RESULTS_DIR=$EXAMPLE_DIR/results

# results directory:
RESULTS_BENCHMARK_DIR=$EXAMPLE_RESULTS_DIR/benchmark_results
RESULTS_LOGS=$EXAMPLE_RESULTS_DIR/logs
RESULTS_PROFILE=$EXAMPLE_RESULTS_DIR/profile
RESULTS_SPECTRA_PLOTS=$EXAMPLE_RESULTS_DIR/spectra_plots
RESULTS_SPECTRA_RAW=$EXAMPLE_RESULTS_DIR/spectra_results
