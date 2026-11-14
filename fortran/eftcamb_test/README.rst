===================
EFTCAMB: test suite
===================

This folder contains the EFTCAMB test suite.

The purpose of this code is to test wether the EFTCAMB code is correctly compiled.
This set of scripts is also used by EFTCAMB developers to make sure that the code is working properly and that changes introduced do not compromise the results of the code or its stability.

1. EFTCAMB test requirements:
=============================

The test suite requires numdiff to properly work and is based on a set of bash scripts.

2. Running tests:
=================

To run the test suite issue::

	make

in this folder.
This will create a set of spectra with the EFTCAMB code for several models and compare them with the thrusted results.

If you issue::

  make all

Then the test suite will run also benchmarks and profiling for the same models and plot the comparison between these models and GR.

3. Adding a test:
=================

To add a test just add the parameter file to the parameters folder. Everything else will be executed automatically.

``parameters/base_params.ini`` holds every setting the tests share: the cosmology (self-consistent, not a fit to
data: one massive neutrino with omnuh2 = 0.01, N_eff = 3.044, BBN helium, r = 0.1), the outputs and the accuracy.
It holds no EFTCAMB parameter, so every EFTCAMB flag takes its code default unless a test file sets it, and it is not
a test itself (the test scripts skip it). ``1_GR.ini`` is the GR reference, with EFTflag = 0 set explicitly.
Every other file sets only its model keys, its perturbation turn-on, ``EFTCAMB_write_background = T`` (the background
is part of the legacy comparison) and, where needed, the model background grid.
Its header says what it tests, where the parameter point comes from, and the measured turn-on convergence and
deviation from GR. A new test file should follow the same rules:

- the turn-on is the latest decade at which TT, EE, TE/sqrt(TT EE) (l <= 1400), C_l^phiphi (l <= 400) and
  P(k <= 0.1/Mpc) change by at most 1e-4 against an earlier turn-on (1e-8, or later when the model background starts
  less than two decades before it), preferring a turn-on without EFTCAMB warnings;
- the default stability checks are on. Beyond-GLPV models (gamma6, Horava, 3+1) switch the generic ghost and gradient
  checks off, because the code has no conditions for them; 6_HDSK_TG_skipchecks tests the skip flags;
- the parameter point modifies gravity as much as the checks allow;
- every output is finite and the spectra are positive: compare_legacy first runs ``python/check_spectra.py`` on the
  new results, since numdiff accepts nan against nan and zeros against zeros.

4. Baseline results:
====================

The folder eftcamb_legacy contains a git repository that hosts the legacy results.
This is enclosed as a submodule of the main EFTCAMB repository.

Running the test for the first time will initialize it.

5. Makefile targets:
====================

The test suite Makefile has several targets::

  all: run all targets
  test: creates the spectra and tests them

  clean: removes all test results

  spectra: runs EFTCAMB to produce the spectra
  benchmark: runs the benchmarker to asses performances of the code
  compare_legacy: compares the spectra results with the legacy ones
  profile: runs the profiler
  spectra_plot: plots the comparison between all the spectra and GR
