EFTCAMB
=======

[![Build Status](https://travis-ci.org/EFTCAMB/EFTCAMB.svg?branch=new_features)](https://travis-ci.org/EFTCAMB/EFTCAMB)

This folder contains the EFTCAMB code.

### 1. EFTCAMB Requirements:

Compiling EFTCAMB requires a modern fortran compiler capable of handeling F2008 features.

ifort
gcc/gfortran 6.0

### 2. Installation procedure:


### 3. Examples:


### 4. Documentation:


### 5. Citing this work:


### 6. Licence Information:

EFTCAMB is a modification of the CAMB code.
The code part of CAMB that is not modified is copyrighted by the CAMB authors and released under their licence.

For the part of code that constitutes EFTCAMB see the LICENSE file in ``eftcamb/LICENSE``.

### 7. Build system target:

### 8. EFTCAMB source files:

In the folder ``eftcamb`` all the source files for EFTCAMB are stored. 
In an effort to have small and readable files the the naming convention allows to have an 
intuition of the hierarchy of the code from alphabetical order of files.

For this reason we use the following convention for the prefixes:

* ``01_`` compile time utilities
* ``02_`` pure algorithms
* ``03_`` general parametrizations for 1D functions
* ``04_`` general parametrizations for 2D functions
* ``05_`` abstract implementation of EFT models
* ``06_`` implementation of pure EFT models
* ``07_`` implementation of alternative EFT parametrizations
* ``08_`` implementation of designer mapping EFT models
* ``09_`` implementation of full mapping EFT models
* ``10_`` general EFT algorithms (RGR, stability, init)