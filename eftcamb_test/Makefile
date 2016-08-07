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
# This file contains the make targets for the test suite
#

# structure of the folder:

TEST_DIR=$(abspath .)
SCRIPT_DIR=$(TEST_DIR)/test_scripts
RESULTS_DIR=$(TEST_DIR)/results
PARAMETER_DIR=$(TEST_DIR)/parameters
PYTHON_DIR=$(TEST_DIR)/python
LEGACY_DIR=$(TEST_DIR)/results_legacy

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

compare_legacy: spectra
	@bash $(SCRIPT_DIR)/legacy_comp_spectra.sh

create_legacy: spectra
	@rm -rf $(LEGACY_DIR)/*
	@mv -f  $(RESULTS_DIR)/Spectra_results/* $(LEGACY_DIR)/
	@rm -rf $(LEGACY_DIR)/*.log

# clean target:

clean:
	@rm -rf $(RESULTS_DIR)/Spectra_results/*
	@rm -rf $(RESULTS_DIR)/Spectra_Plots/*
	@rm -rf $(RESULTS_DIR)/Profile/*
	@rm -rf $(RESULTS_DIR)/Legacy_Spectra_Plot/*
	@rm -rf $(RESULTS_DIR)/Benchmark_Results/*
