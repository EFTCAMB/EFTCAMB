!----------------------------------------------------------------------------------------
!
! This file is part of EFTCAMB.
!
! Copyright (C) 2013-2016 by the EFTCAMB authors
!
! The EFTCAMB code is free software;
! You can use it, redistribute it, and/or modify it under the terms
! of the GNU General Public License as published by the Free Software Foundation;
! either version 3 of the License, or (at your option) any later version.
! The full text of the license can be found in the file eftcamb/LICENSE at
! the top level of the EFTCAMB distribution.
!
!----------------------------------------------------------------------------------------

!> @file 01_EFT_def.f90
!! This file contains the EFTCAMB compile time options.


!----------------------------------------------------------------------------------------
!> This module contains the definitions of all the EFTCAMB compile time flags.

!> @author Bin Hu, Marco Raveri, Simone Peirone

module EFT_def

    use Precision

    implicit none

    ! EFT compile time flags:

    real(dl), parameter :: EFTturnonpiInitial = 1.d-2 !< Turn on pi field flag:
        !!    Sets the scale factor at which the code starts to evolve the pi field.
        !!    At times earlier than these the code evolves perturbations as in GR.
        !!    This number is used as a lower bound and is refined if the teory is very close to GR
        !!    by the EFTreturntoGR module.

    real(dl), parameter :: EFTtoGR = 1.d-8 !< Return to GR flag:
        !!    This is the threshold at which a theory is considered to be exactly GR.

    integer , parameter :: EFT_names_max_length       = 20    !< maximum length of names for EFT functions and parameters.
    integer , parameter :: EFT_names_latex_max_length = 40    !< maximum length of latex names for EFT functions and parameters.

    integer , parameter :: EFT_RGR_num_points   = 1000        !< number of points to sample logaritmically the time in the return to GR module.

#ifdef DEBUG
    logical , parameter :: DebugEFTCAMB = .true.              !< EFTCAMB debug flag.This will turn on printing of many things to aid debugging the code.
#else
    logical , parameter :: DebugEFTCAMB = .false.             !< EFTCAMB debug flag.This will turn on printing of many things to aid debugging the code.
#endif

end module EFT_def

!----------------------------------------------------------------------------------------
