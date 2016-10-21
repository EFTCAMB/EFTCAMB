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

!> @file 06p1_abstract_EFTCAMB_full_map.f90
!! This file contains the abstract definition of all the places where EFTCAMB interacts
!! with CAMB in case of a full mapping model. All models implementing a model in which
!! the cosmological background is computed from the values of the EFT functions
!! should inherit from this class.


!----------------------------------------------------------------------------------------
!> This module contains the abstract definition of all the places where EFTCAMB interacts
!! with CAMB in case of a full mapping model. All models implementing a model in which
!! the cosmological background is computed from the values of the EFT functions
!! should inherit from this class.

!> @author Bin Hu, Marco Raveri

module EFTCAMB_abstract_model_full

    use precision
    use IniFile
    use EFTCAMB_cache
    use EFTCAMB_abstract_model

    implicit none

    private

    public EFTCAMB_full_model

    !----------------------------------------------------------------------------------------
    !> This is the abstract type for EFTCAMB models.
    type, extends ( EFTCAMB_model ), abstract :: EFTCAMB_full_model

    contains

        procedure :: compute_dtauda   => EFTCAMBFullModelComputeDtauda    !< function that computes dtauda = 1/sqrt(a^2H^2).
        procedure :: compute_adotoa   => EFTCAMBFullModelComputeAdotoa    !< subroutine that computes adotoa = H.
        procedure :: compute_H_derivs => EFTCAMBFullModelComputeHubbleDer !< subroutine that computes the two derivatives wrt conformal time of H.

    end type EFTCAMB_full_model

contains

    ! ---------------------------------------------------------------------------------------------
    !> Function that computes dtauda = 1/sqrt(a^2H^2).
    function EFTCAMBFullModelComputeDtauda( self, a, eft_par_cache, eft_cache )

        implicit none

        class(EFTCAMB_full_model)                    :: self          !< the base class
        real(dl), intent(in)                         :: a             !< the input scale factor.
        type(EFTCAMB_parameter_cache), intent(inout) :: eft_par_cache !< the EFTCAMB parameter cache that contains all the physical parameters.
        type(EFTCAMB_timestep_cache ), intent(inout) :: eft_cache     !< the EFTCAMB timestep cache that contains all the physical values.

        real(dl) :: EFTCAMBFullModelComputeDtauda                     !< the output dtauda

        write(*,*) 'IW'
        stop

    end function EFTCAMBFullModelComputeDtauda

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes adotoa = H.
    subroutine EFTCAMBFullModelComputeAdotoa( self, a, eft_par_cache, eft_cache )

        implicit none

        class(EFTCAMB_full_model)                    :: self          !< the base class
        real(dl), intent(in)                         :: a             !< the input scale factor
        type(EFTCAMB_parameter_cache), intent(inout) :: eft_par_cache !< the EFTCAMB parameter cache that contains all the physical parameters.
        type(EFTCAMB_timestep_cache ), intent(inout) :: eft_cache     !< the EFTCAMB timestep cache that contains all the physical values.

        write(*,*) 'IW'
        stop

    end subroutine EFTCAMBFullModelComputeAdotoa

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes the two derivatives wrt conformal time of H.
    subroutine EFTCAMBFullModelComputeHubbleDer( self, a, eft_par_cache, eft_cache )

        implicit none

        class(EFTCAMB_full_model)                    :: self          !< the base class
        real(dl), intent(in)                         :: a             !< the input scale factor
        type(EFTCAMB_parameter_cache), intent(inout) :: eft_par_cache !< the EFTCAMB parameter cache that contains all the physical parameters.
        type(EFTCAMB_timestep_cache ), intent(inout) :: eft_cache     !< the EFTCAMB timestep cache that contains all the physical values.

        write(*,*) 'IW'
        stop

    end subroutine EFTCAMBFullModelComputeHubbleDer

    !----------------------------------------------------------------------------------------

end module EFTCAMB_abstract_model_full

!----------------------------------------------------------------------------------------
