===================
EFTCAMB: test suite
===================

This folder contains the EFTCAMB test suite.

The purpose of this code is to test wether the EFTCAMB code is correctly compiled.
This set of scripts is also used by EFTCAMB developers to make sure that the code is working properly and that changes introduced do not compromise the results of the code or its stability.

1. EFTCAMB test requirements:
===============================

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

4. Replacing baseline results:
==============================

By running the command::

	make create_legacy

the test script will create all the spectra and replace the thrusted ones with the new ones.
This will totally overwrite the previous results so this command is dangerous to use. Use it with care.

5. Makefile targets:
====================

The test suite Makefile has several targets:

  all: run all targets
  default: creates the spectra and tests them

  clean: removes all test results

  spectra: runs EFTCAMB to produce the spectra
  benchmark: runs the benchmarker to asses performances of the code
  compare_legacy: compares the spectra results with the legacy ones
  profile: runs the profiler
  spectra_plot: plots the comparison between all the spectra and GR
