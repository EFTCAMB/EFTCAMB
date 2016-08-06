script_folder=test_scripts/

all: spectra benchmark profile stability spectra_plot stability_plot legacy_plot

default: spectra spectra_plot legacy_plot

# results creation target:

spectra:
	@./$(script_folder)test_spectra.sh

benchmark:
	@./$(script_folder)test_benchmark.sh

profile:
	@./$(script_folder)test_profile.sh

stability:
	@./$(script_folder)test_stability.sh

# plot target:

spectra_plot:
	@./$(script_folder)test_plot.sh
	
stability_plot:
	@./$(script_folder)test_stability_plot.sh

legacy_plot:
	@./$(script_folder)legacy_comp_spectra.sh

# legacy creation target:

create_legacy: spectra legacy_plot
	@rm -rf results_legacy/Spectra_Legacy/*
	@mv -f  results/Spectra_results/* results_legacy/Spectra_Legacy/

# clean target:

clean:
	@rm -rf results/Spectra_results/*
	@rm -rf results/Spectra_Plots/*
	@rm -rf results/Stability_Results/*
	@rm -rf results/Stability_Plots/*
	@rm -rf results/Profile/*
	@rm -rf results/Legacy_Spectra_Plot/*
	@rm -rf results/Benchmark_Results/*
