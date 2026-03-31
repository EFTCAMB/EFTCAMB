EFTCAMB
=======

[![Build Status](https://travis-ci.org/EFTCAMB/EFTCAMB.svg?branch=new_features)](https://travis-ci.org/EFTCAMB/EFTCAMB)

This folder contains the EFTCAMB code.

### 1. A Direct Installation Procedure:

To get started with the EFTCAMB and its Python wrapper, you need to first download the entire package using Git:

 	git clone https://github.com/EFTCAMB/EFTCAMB.git

Next, it's recommended to create a new environment and install the necessary dependencies before compiling and using the Python wrapper:

    gfortran; lapack; numpy; python; sympy; scipy; packaging

For example, using conda, you can create a new environment and install the required packages with the following command:

    conda create -n EFTCAMB_env gfortran lapack numpy python sympy scipy packaging -c conda-forge

**Note for Mac Users:**
If you're using a Mac with an Apple Silicon chip (e.g., M1, M2), ensure you install the appropriate version of gfortran for your architecture (arm64), for example:

    conda create -n EFTCAMB_env gfortran_osx-arm64 ... -c conda-forge

You can verify that you're using the correct version of gfortran by executing:

    which gfortran

Finally, navigate to the package directory and install the package in editable mode:

    cd fortran
    make python

When you want to use the package by the python wrapper, import it where you compile, for example:

    camb_installation_path = './../../'
    camb_path = os.path.realpath(os.path.join(os.getcwd(),camb_installation_path))
    sys.path.insert(0,camb_path)
    print('camb path:', camb_path)
    import camb

### 2. Documentation:

The contents of EFTCAMB are listed in this chart:

![Chart](/find_your_model/EFTCAMB_STRUCTURE(Charts).png)


We provide a set of notes that contain all the details and formulas of the EFTCAMB implementation:

* *EFTCAMB/EFTCosmoMC: Numerical Notes v3.0*
    Bin Hu, Marco Raveri, Noemi Frusciante, Alessandra Silvestri, [arXiv:1405.3590 [astro-ph.CO]](http://arxiv.org/abs/1405.3590)

The EFTCAMB source files documentation is automatically built at any modification of the code and can
be found at [this link](https://eftcamb.github.io/EFTCAMB/).

Besides above documentation, there are several flowchart or markdown files helping you easily find your model flags and set parameters in “find_your_model" folder.

### 3. Examples and Usage with Cobaya:

The EFTCAMB distribution contains a folder called ``example`` that can be used to produce some example notebooks for brief instruction and sevral example input files of cobaya for sampling and statistical modelling.

To use cobaya, just follow the instruction at:

https://readthedocs.org/projects/cobaya/badge/?version=latest

The example folder is an example of how to select EFTCAMB flags and set parameters needed from the flags. You can copy the folder and replace the input files to produce single purpose packages to easily produce and plot results.

### 4. Citing this work:

If you use the EFTCAMB/EFTCosmoMC package, please refer the original CAMB/ CosmoMC paper and ours:

* *Effective Field Theory of Cosmic Acceleration: an implementation in CAMB*
    Bin Hu, Marco Raveri, Noemi Frusciante, Alessandra Silvestri,
    [arXiv:1312.5742 [astro-ph.CO]](http://arxiv.org/abs/1312.5742) [Phys.Rev. D89 (2014) 103530](http://journals.aps.org/prd/abstract/10.1103/PhysRevD.89.103530)


* *Effective Field Theory of Cosmic Acceleration: constraining dark energy with CMB data*
    Marco Raveri, Bin Hu, Noemi Frusciante, Alessandra Silvestri,
    [arXiv:1405.1022 [astro-ph.CO]](https://arxiv.org/abs/1405.1022) [Phys.Rev. D90 (2014) 043513](http://journals.aps.org/prd/abstract/10.1103/PhysRevD.90.043513)

This is the usual, fair way of giving credit to contributors to a
scientific result. In addition, it helps us justify our effort in
developing the EFTCAMB code as an academic undertaking.

### 5. Licence Information:

EFTCAMB is a modification of the CAMB code.
The code part of CAMB that is not modified is copyrighted by the CAMB authors and released under their licence.

For the part of code that constitutes EFTCAMB see the LICENSE file in ``eftcamb/LICENSE``.

### 6. Build system target:

In addition to CAMB makefile targets EFTCAMB comes with the additional:

* ``eftcamb``: to compile EFTCAMB;
* ``eftcamb_apps``: to compile EFTCAMB applications like the EFTCAMB benchmarker;
* ``eftcamb_dep``: to automatically sort out EFTCAMB source file dependencies;
* ``eftcamb_doc``: to build the EFTCAMB automatic documentation;
* ``intel_profile``: to compile EFTCAMB with the options that allow profiling with the VTUNE profiler;
* ``profile``: to compile EFTCAMB with the options that allow general profiling;

### 7. EFTCAMB source files:

In the folder ``eftcamb`` all the source files for EFTCAMB are stored.
In an effort to have small and readable files the the naming convention allows to have an
intuition of the hierarchy of the code from alphabetical order of files.

For this reason we use the following convention for the prefixes:

* ``01_`` compile time utilities;
* ``02_`` pure and general purpose algorithms;
* ``03_`` EFTCAMB cache containing the storage for all cosmological quantities of interest;
* ``04_`` general parametrizations for 1D functions;
* ``05_`` general parametrizations for 2D functions;
* ``06_`` abstract implementation of EFT models;
* ``07_`` implementation of pure EFT models;
* ``08_`` implementation of alternative EFT parametrizations;
* ``09_`` implementation of designer mapping EFT models;
* ``10_`` implementation of full mapping EFT models;
* ``11_`` general EFT algorithms (RGR, stability, init);

If you modify or add one or more files make sure to issue ``make eftcamb_dep`` before compiling the code to ensure that all dependencies are properly sorted out and built.