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

!> @file 1_EFT_def.f90
!! This file contains the EFTCAMB compile time options.


!----------------------------------------------------------------------------------------
!> This module contains the definitions of all the EFTCAMB compile time flags.

!> @author Bin Hu, Marco Raveri

module EFTDef

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

    logical , parameter :: EarlyTimeStability = .true. !< Early time stability:
        !!    EFTCAMB checks the stability of the theory at all times making sure that the choosen model
        !!    is viable during radiation, matter and DE eras.
        !!    It is possible to enforce stability only at late times (from EFTturnonpiInitial to today)
        !!    by setting this flag to false.
        !!    This choice will however make the results dependent on what one calls late time,
        !!    i.e. the choice of EFTturnonpiInitial, so we advice to state it clearly when reporting results.

    logical , parameter :: DebugEFTCAMB = .false. !< EFTCAMB debug flag.This will turn on printing of many
        !! things to aid debugging the code.

    integer , parameter :: EFT_names_max_length = 20   !< maximum length of names for EFT functions and parameters.

end module EFTDef

!----------------------------------------------------------------------------------------
