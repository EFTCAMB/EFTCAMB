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
# This file contains the make targets for the example suite
#

# structure of the folder:

EXAMPLE_DIR=$(abspath .)
SCRIPT_DIR=$(EXAMPLE_DIR)/scripts
RESULTS_DIR=$(EXAMPLE_DIR)/results
PARAMETER_DIR=$(EXAMPLE_DIR)/parameters
PYTHON_DIR=$(EXAMPLE_DIR)/python

# make general targets:

all: spectra benchmark profile spectra_plot

examples: spectra spectra_plot

default: examples

# results creation target:

spectra:
	@bash $(SCRIPT_DIR)/example_spectra.sh

benchmark:
	@bash $(SCRIPT_DIR)/example_benchmark.sh

profile:
	@bash $(SCRIPT_DIR)/example_profile.sh

# plot target:

spectra_plot:
	@bash $(SCRIPT_DIR)/example_plot.sh
	
# clean target:

clean:
	@rm -rf $(RESULTS_DIR)/spectra_results/*
	@rm -rf $(RESULTS_DIR)/spectra_plots/*
	@rm -rf $(RESULTS_DIR)/profile/*
	@rm -rf $(RESULTS_DIR)/benchmark_results/*
