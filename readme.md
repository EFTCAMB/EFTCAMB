EFTCAMB
=======

This folder contains the EFTCAMB code.

### 1. EFTCAMB Requirements:

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

* ``1_`` compile time utilities
* ``2_`` pure algorithms
* ``3_`` general parametrizations of functions
* ``4_`` abstract implementation of EFT models
* ``5_`` implementation of pure EFT models
* ``6_`` implementation of alternative EFT parametrizations
* ``7_`` implementation of designer mapping EFT models
* ``8_`` implementation of full mapping EFT models
* ``9_`` general EFT algorithms (RGR, stability, init)