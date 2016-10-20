#CAMB Makefile

#Set FISHER=Y to compile bispectrum fisher matrix code
FISHER=

#Will detect ifort/gfortran or edit for your compiler
ifortErr = $(shell which ifort >/dev/null; echo $$?)
ifeq "$(ifortErr)" "0"

#Intel compiler
# For OSX replace shared by dynamiclib
F90C     = ifort
# EFTCAMB MOD START: flags for the latest ifort
#FFLAGS = -openmp -fast -W0 -WB -fpp2 -vec_report0
#DEBUGFLAGS = -openmp -g -check all -check noarg_temp_created -traceback -fpp -fpe0
FFLAGS = -qopenmp -O3 -W0 -WB -fpp -qopt-report=0
DEBUGFLAGS = -qopenmp -mkl=parallel -fpp -g -qopt-report=0 -fp-stack-check -O0 -traceback -check all -check bounds -check uninit -check noarg_temp_created -DDEBUG
# EFTCAMB MOD END.
SFFLAGS = -shared -fpic
## This is flag is passed to the Fortran compiler allowing it to link C++ if required (not usually):
F90CRLINK = -cxxlib
MODOUT = -module $(OUTPUT_DIR)
SMODOUT = -module $(DLL_DIR)
ifneq ($(FISHER),)
FFLAGS += -mkl
endif

else
gfortErr = $(shell which gfortran >/dev/null; echo $$?)
ifeq "$(gfortErr)" "0"

#Gfortran compiler:
#The options here work in v4.6+. Python wrapper needs v4.9+.
F90C     = gfortran
SFFLAGS =  -shared -fPIC

FFLAGS =  -O3 -fopenmp -ffast-math -fmax-errors=4 -ffree-line-length-none -funroll-loops -cpp -ffpe-summary=none
DEBUGFLAGS = -cpp -g -fbounds-check -fbacktrace -ffree-line-length-none -fmax-errors=4 -ffpe-trap=invalid,overflow,zero -DDEBUG
MODOUT =  -J$(OUTPUT_DIR)
SMODOUT = -J$(DLL_DIR)

ifneq ($(shell uname -s),Darwin)
#native optimization does not work on Mac
FFLAGS+=-march=native
endif
endif
endif

IFLAG = -I

#G95 compiler
#F90C   = g95
#FFLAGS = -O2

#SGI, -mp toggles multi-processor. Use -O2 if -Ofast gives problems.
#F90C     = f90
#FFLAGS  = -Ofast -mp

#Digital/Compaq fortran, -omp toggles multi-processor
#F90C    = f90
#FFLAGS  = -omp -O4 -arch host -math_library fast -tune host -fpe1

#Absoft ProFortran, single processor:
#F90C     = f95
#FFLAGS = -O2 -cpu:athlon -s -lU77 -w -YEXT_NAMES="LCS" -YEXT_SFX="_"

#NAGF95, single processor:
#F90C     = f95
#FFLAGS = -DNAGF95 -O3

#PGF90
#F90C = pgf90
#FFLAGS = -O2 -DESCAPEBACKSLASH -Mpreprocess

#Sun V880
#F90C = mpf90
#FFLAGS =  -O4 -openmp -ftrap=%none -dalign

#Sun parallel enterprise:
#F90C     = f95
#FFLAGS =  -O2 -xarch=native64 -openmp -ftrap=%none
#try removing -openmp if get bus errors. -03, -04 etc are dodgy.

#IBM XL Fortran, multi-processor (run gmake)
#F90C     = xlf90_r
#FFLAGS  = -DESCAPEBACKSLASH -DIBMXL -qsmp=omp -qsuffix=f=f90:cpp=F90 -O3 -qstrict -qarch=pwr3 -qtune=pwr3

#Settings for building camb_fits
#Location of FITSIO and name of library
FITSDIR       ?= /usr/local/lib
FITSLIB       = cfitsio
#Location of HEALPIX for building camb_fits
HEALPIXDIR    ?= /usr/local/healpix

ifneq ($(FISHER),)
FFLAGS += -DFISHER
EXTCAMBFILES = Matrix_utils.o
else
EXTCAMBFILES =
endif

DEBUGFLAGS ?= FFLAGS
Debug: FFLAGS=$(DEBUGFLAGS)

include ./Makefile_main
