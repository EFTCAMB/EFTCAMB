EFTCAMB
=======

[![Build Status](https://travis-ci.org/EFTCAMB/EFTCAMB.svg?branch=new_features)](https://travis-ci.org/EFTCAMB/EFTCAMB)

This folder contains the EFTCAMB code.

### 1. EFTCAMB Requirements:

Compiling EFTCAMB requires a modern fortran compiler capable of handeling F2008 features.
These includes:

	ifort v>15.0 (?)
	gcc/gfortran v>6.0

To use other parts of the code, like the test or documentation parts other requirements have to be met:

A docker with all the required libraries is available at: https://hub.docker.com/r/eftcamb/eftbox/

### 2. Installation procedure:


### 3. Examples:


### 4. Documentation:


### 5. Citing this work:

If you use the EFTCAMB/EFTCosmoMC package, please refer the original CAMB/ CosmoMC paper and ours:

* *Effective Field Theory of Cosmic Acceleration: an implementation in CAMB*  
    Bin Hu, Marco Raveri, Noemi Frusciante, Alessandra Silvestri,  
    arXiv:1312.5742 [astro-ph.CO] Phys.Rev. D89 (2014) 103530

* *Effective Field Theory of Cosmic Acceleration: constraining dark energy with CMB data*  
    Marco Raveri, Bin Hu, Noemi Frusciante, Alessandra Silvestri,  
    arXiv:1405.1022 [astro-ph.CO] Phys.Rev. D90 (2014) 043513


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