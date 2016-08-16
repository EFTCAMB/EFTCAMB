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

!> @file 05p2_abstract_EFTCAMB_full.f90
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
    use EFTCAMB_abstract_model

    implicit none

    !----------------------------------------------------------------------------------------
    !> This is the abstract type for EFTCAMB models.
    type, extends ( EFTCAMB_model ), abstract :: EFTCAMB_full_model

    contains

        procedure :: compute_dtauda => EFTCAMBFullModelComputeDtauda !< function that computes dtauda = 1/sqrt(a^2H^2).
        procedure :: compute_adotoa => EFTCAMBFullModelComputeAdotoa !< subroutine that computes adotoa = H and its two derivatives wrt conformal time.

    end type EFTCAMB_full_model

contains

    ! ---------------------------------------------------------------------------------------------
    !> Function that computes dtauda = 1/sqrt(a^2H^2).
    function EFTCAMBFullModelComputeDtauda( self , a, grhoa2, &
        & grhok, grhov, &
        & grhoc, grhob, &
        & grhog, grhornomass )

        implicit none

        class(EFTCAMB_full_model) :: self                      !< the base class
        real(dl), intent(in)  :: a                             !< the input scale factor
        real(dl), intent(in)  :: grhoa2                        !< the input value of 8 \piG \rho_tot a^2
        real(dl), intent(in)  :: grhok                         !< the input value of curvature density
        real(dl), intent(in)  :: grhov                         !< the input value of DE density
        real(dl), intent(in)  :: grhoc                         !< the input value of CDM density
        real(dl), intent(in)  :: grhob                         !< the input value of Baryon density
        real(dl), intent(in)  :: grhog                         !< the input value of Radiation density
        real(dl), intent(in)  :: grhornomass                   !< the input value of massless neutrinos density
        real(dl)              :: EFTCAMBFullModelComputeDtauda !< the output dtauda

        EFTCAMBFullModelComputeDtauda = 0._dl

    end function EFTCAMBFullModelComputeDtauda

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes adotoa = H and its two derivatives wrt conformal time.
    subroutine EFTCAMBFullModelComputeAdotoa( self, a, adotoa, Hdot, Hdotdot )

        implicit none

        class(EFTCAMB_full_model) :: self   !< the base class
        real(dl), intent(in)  :: a          !< the input scale factor
        real(dl), intent(out) :: adotoa     !< the output value of H at the given scale factor
        real(dl), intent(out) :: Hdot       !< the output value of dH/dtau at the given scale factor
        real(dl), intent(out) :: Hdotdot    !< the output value of d^2H/dtau^2 at the given scale factor

    end subroutine EFTCAMBFullModelComputeAdotoa

    !----------------------------------------------------------------------------------------

end module EFTCAMB_abstract_model_full

!----------------------------------------------------------------------------------------
