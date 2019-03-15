#----------------------------------------------------------------------------------------
#
# This file is part of EFTCAMB.
#
# Copyright (C) 2013-2019 by the EFTCAMB authors
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
# This file contains the make targets for the test suite
#

# structure of the folder:

TEST_DIR=$(abspath .)
SCRIPT_DIR=$(TEST_DIR)/test_scripts
RESULTS_DIR=$(TEST_DIR)/results
PARAMETER_DIR=$(TEST_DIR)/parameters
PYTHON_DIR=$(TEST_DIR)/python
LEGACY_DIR=$(TEST_DIR)/eftcamb_legacy

# make general targets:

all: spectra benchmark profile compare_legacy spectra_plot

test: spectra compare_legacy

default: test

# results creation target:

spectra:
	@bash $(SCRIPT_DIR)/test_spectra.sh

benchmark:
	@bash $(SCRIPT_DIR)/test_benchmark.sh

profile:
	@bash $(SCRIPT_DIR)/test_profile.sh

# plot target:

spectra_plot:
	@bash $(SCRIPT_DIR)/test_plot.sh

# legacy target:

init_legacy:
	@bash $(SCRIPT_DIR)/legacy_init.sh

compare_legacy: init_legacy
	@bash $(SCRIPT_DIR)/legacy_comp_spectra.sh

# clean target:

clean:
	@rm -rf $(RESULTS_DIR)/spectra_results/*
	@rm -rf $(RESULTS_DIR)/spectra_plots/*
	@rm -rf $(RESULTS_DIR)/profile/*
	@rm -rf $(RESULTS_DIR)/legacy_spectra_plot/*
	@rm -rf $(RESULTS_DIR)/benchmark_results/*
