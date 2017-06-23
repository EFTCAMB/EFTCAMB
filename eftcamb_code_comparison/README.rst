==============================
EFTCAMB: code comparison suite
==============================

This folder contains the EFTCAMB parameter files used in the EFTCAMB code comparison.

The purpose of this code is to run EFTCAMB for a set of thrusted parameters and precision.

1. EFTCAMB example requirements:
================================

The example suite requires numdiff to properly work and is based on a set of bash and python scripts.

2. Running examples:
====================

To run the example suite issue::

	make

in this folder.
This will create a set of spectra with the EFTCAMB code for several models and compare then plot the results.

If you issue::

  make all

Then the example suite will run also benchmarks and profiling for the same models.

3. Adding an example:
=====================

To add an example just add the parameter file to the parameters folder. Everything else will be executed automatically.

4. Makefile targets:
====================

The test suite Makefile has several targets::

  all: run all targets
  examples: creates the spectra and plots them

  clean: removes all results

  spectra: runs EFTCAMB to produce the spectra
  spectra_plot: plots the comparison between all the spectra and GR
  benchmark: runs the benchmarker to asses performances of the code
  profile: runs the profiler

5. Reusing the example suite:
=============================

The example suite can just be duplicated to create small example packages.
Everything will work automatically, the only thing to change is the output path for results in the parameters files.
