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

!> @file 9_EFTCAMB_main.f90
!! This file contains the general EFTCAMB driver. We have a type that encapsulate all
!! EFTCAMB parameters and stuff and this is owned by CAMB to do all EFT calculations.


!----------------------------------------------------------------------------------------
!> This module contains the general EFTCAMB driver. We have a type that encapsulate all
!! EFTCAMB parameters and stuff and this is owned by CAMB to do all EFT calculations.

!> @author Bin Hu, Marco Raveri

module EFTCAMB_main

    use precision

    implicit none

    !----------------------------------------------------------------------------------------
    !> This is the main object for EFTCAMB. Get one of these and you can use all the stuff
    !! in EFTCAMB.
    type EFTCAMB

        ! EFTCAMB model selection flags:
        logical   :: EFTflag              !< Main EFTCAMB model selection flag. Decides one of the four modes to run EFTCAMB.
        logical   :: PureEFTmodel         !< Model selection flag for pure EFT models.
        logical   :: AltParEFTmodel       !< Model selection flag for alternative EFT parametrizations.
        logical   :: DesignerEFTmodel     !< Model selection flag for designer mapping EFT models.
        logical   :: FullMappingEFTmodel  !< Model selection flag for full mapping EFT models.

    contains

        ! utility functions:
        procedure :: EFTCAMB_init_from_file => read_EFTCAMB_flags_from_file  !< subroutine that initializes EFTCAMB from an INI file.

    end type EFTCAMB

    ! ---------------------------------------------------------------------------------------------

contains

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the i-th parameter name of the model
    subroutine read_EFTCAMB_flags_from_file( self )

        implicit none

        class(EFTCAMB)      :: self       !< the base class

    end subroutine read_EFTCAMB_flags_from_file

    ! ---------------------------------------------------------------------------------------------

end module EFTCAMB_main

!----------------------------------------------------------------------------------------
