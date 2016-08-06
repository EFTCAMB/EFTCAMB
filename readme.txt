Thanks for downloading EFTCAMB! This is EFTCAMB_Oct15, which is compatible with CAMB_Feb15. 

In order to build EFTCAMB, please follow these 6 steps:

1. Download CAMB version Feb15

2. Add these new files into your CAMB folder:

EFT_def.90
EFT_designer.f90
EFT_functions.f90
EFT_Horndeski.f90
EFT_main.f90
EFT_stabilitySpace.f90
equations_EFT.f90


3. Replace the default files with the following files:

cmbmain.f90
inidriver.F90
lensing.f90
Makefile
Makefile_main
modules.f90
params.ini
reionization.f90
subroutines.f90

4. configure the Makefile according to your Fortran compiler.

5. $make

6. $./camb params.ini

7. if you want to perform basic explorations on the parameter space of some simple models you can use a built in stability sampler. To compile it just issue:
   $make eftstability

8. $./eftstability

--
EFTCAMB team
Oct/2015

