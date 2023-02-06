===================
CAMB
===================
:CAMB: Code for Anisotropies in the Microwave Background
:Author: Antony Lewis and Anthony Challinor
:Homepage: https://camb.info/

.. image:: https://img.shields.io/pypi/v/camb.svg?style=flat
        :target: https://pypi.python.org/pypi/camb/
.. image:: https://img.shields.io/conda/vn/conda-forge/camb.svg
   :target: https://anaconda.org/conda-forge/camb
.. image:: https://readthedocs.org/projects/camb/badge/?version=latest
   :target: https://camb.readthedocs.io/en/latest
.. image:: https://travis-ci.com/cmbant/camb.svg?branch=master
  :target: https://app.travis-ci.com/cmbant/camb/builds
.. image:: https://mybinder.org/badge_logo.svg
  :target: https://mybinder.org/v2/gh/cmbant/CAMB/HEAD?filepath=docs%2FCAMBdemo.ipynb

Description and installation
=============================

CAMB is a cosmology code for calculating cosmological observables, including
CMB, lensing, source count and 21cm angular power spectra, matter power spectra, transfer functions
and background evolution. The code is in Python, with numerical code implemented in fast modern Fortran.

See the `CAMB python example notebook <https://camb.readthedocs.io/en/latest/CAMBdemo.html>`_ for a
quick introduction to how to use the CAMB Python package.

=======
EFTCAMB
=======

[![Build Status](https://travis-ci.org/EFTCAMB/EFTCAMB.svg?branch=new_features)](https://travis-ci.org/EFTCAMB/EFTCAMB)

This folder contains the EFTCAMB code.

### 1. EFTCAMB Requirements:

Compiling EFTCAMB requires a modern fortran compiler capable of handeling F2008 features.
These includes:

	ifort v>15.0 (?)
	gcc/gfortran v>6.0

To use other parts of the code, like the test or documentation parts other requirements have to be met.
These include a fully fledged python installation. We warmly suggest to install a
bundled package like (https://www.continuum.io/downloads).

A docker with all the required libraries is available at [dockerhub](https://hub.docker.com/r/eftcamb/eftbox/).

### 2. Installation procedure:

For a standard non-editable installation use::

    pip install camb [--user]

The --user is optional and only required if you don't have write permission to your main python installation.
To install from source, clone from github using::

    git clone --recursive https://github.com/EFTCAMB/EFTCAMB.git

Then install using::

    pip install -e ./CAMB [--user]

You will need gfortran 6 or higher installed to compile. Binary files for Windows are also provided, so these are used instead if no
gfortran installation is found on Windows machines. If you have gfortran installed, "python setup.py make"
(and other standard setup commands) will build the Fortran library on all systems (including Windows without directly using a Makefile).

<<<<<<< HEAD
The python wrapper provides a module called "camb" documented in the Python `CAMB documentation <https://camb.readthedocs.io/en/latest/>`_.

=======
>>>>>>> 9a3393f (All the work)
After installation you can also run CAMB from the command line reading parameters from a .ini file, e.g.::

  camb inifiles/planck_2018.ini

To compile the Fortran command-line code run "make camb" in the fortran directory.

### 3. Examples:


### 4. Documentation:

We provide a set of notes that contain all the details and formulas of the EFTCAMB implementation:

* *EFTCAMB/EFTCosmoMC: Numerical Notes v2.0*
    Bin Hu, Marco Raveri, Noemi Frusciante, Alessandra Silvestri, [arXiv:1405.3590 [astro-ph.CO]](http://arxiv.org/abs/1405.3590)

The EFTCAMB source files documentation is automatically built at any modification of the code and can
be found at [this link](https://eftcamb.github.io/EFTCAMB/).

### 5. Citing this work:

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

### 6. Licence Information:

EFTCAMB is a modification of the CAMB code.
The code part of CAMB that is not modified is copyrighted by the CAMB authors and released under their licence.

For the part of code that constitutes EFTCAMB see the LICENSE file in ``eftcamb/LICENSE``.

### 7. Build system target:

### 8. EFTCAMB source files:

<<<<<<< HEAD
===================
=======
In the folder ``eftcamb`` all the source files for EFTCAMB are stored.
In an effort to have small and readable files the the naming convention allows to have an
intuition of the hierarchy of the code from alphabetical order of files.
>>>>>>> 9a3393f (All the work)

For this reason we use the following convention for the prefixes:

<<<<<<< HEAD
    <a href="https://www.sussex.ac.uk/astronomy/"><img src="https://cdn.cosmologist.info/antony/Sussex_white.svg" style="height:200px" height="200px"></a>
    <a href="https://erc.europa.eu/"><img src="https://cdn.cosmologist.info/antony/ERC_white.svg" style="height:200px" height="200px"></a>
    <a href="https://stfc.ukri.org/"><img src="https://cdn.cosmologist.info/antony/STFC_white.svg" style="height:200px" height="200px"></a>
=======
* ``01_`` compile time utilities
* ``02_`` pure algorithms
* ``03_`` EFTCAMB cache containing the storage for all cosmological quantities of interest
* ``04_`` general parametrizations for 1D functions
* ``05_`` general parametrizations for 2D functions
* ``06_`` abstract implementation of EFT models
* ``07_`` implementation of pure EFT models
* ``08_`` implementation of alternative EFT parametrizations
* ``09_`` implementation of designer mapping EFT models
* ``10_`` implementation of full mapping EFT models
* ``11_`` general EFT algorithms (RGR, stability, init)
>>>>>>> 9a3393f (All the work)
