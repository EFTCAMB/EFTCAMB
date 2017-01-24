!----------------------------------------------------------------------------------------
!
! This file is part of EFTCAMB.
!
! Copyright (C) 2013-2017 by the EFTCAMB authors
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

        real(dl) :: a2

        a2=a*a

        EFTCAMBFullModelComputeDtauda = 1._dl/sqrt( ( eft_cache%grhoa2/3._dl +2._dl/3._dl*eft_cache%EFTc*a2 + a2*eft_par_cache%h0_Mpc**2*eft_par_cache%omegav*a2*(1._dl +eft_cache%EFTLambda ) )&
                                      & /( 1._dl +eft_cache%EFTOmegaV +a*eft_cache%EFTOmegaP ) )

    end function EFTCAMBFullModelComputeDtauda

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes adotoa = H.
    subroutine EFTCAMBFullModelComputeAdotoa( self, a, eft_par_cache, eft_cache )

        implicit none

        class(EFTCAMB_full_model)                    :: self          !< the base class
        real(dl), intent(in)                         :: a             !< the input scale factor
        type(EFTCAMB_parameter_cache), intent(inout) :: eft_par_cache !< the EFTCAMB parameter cache that contains all the physical parameters.
        type(EFTCAMB_timestep_cache ), intent(inout) :: eft_cache     !< the EFTCAMB timestep cache that contains all the physical values.

        real(dl)    :: temp

        temp = 1.0_dl/(1.0_dl + eft_cache%EFTOmegaV+ a*eft_cache%EFTOmegaP)*(eft_cache%grhom_t + 2.0_dl*eft_cache%EFTc -eft_cache%EFTLambda )/3.0_dl
        eft_cache%adotoa = sqrt(temp)

    end subroutine EFTCAMBFullModelComputeAdotoa

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes the two derivatives wrt conformal time of H.
    subroutine EFTCAMBFullModelComputeHubbleDer( self, a, eft_par_cache, eft_cache )

        implicit none

        class(EFTCAMB_full_model)                    :: self          !< the base class
        real(dl), intent(in)                         :: a             !< the input scale factor
        type(EFTCAMB_parameter_cache), intent(inout) :: eft_par_cache !< the EFTCAMB parameter cache that contains all the physical parameters.
        type(EFTCAMB_timestep_cache ), intent(inout) :: eft_cache     !< the EFTCAMB timestep cache that contains all the physical values.

        real(dl)    :: a2

        a2=a*a
        eft_cache%Hdot = 1.0_dl/(1.0_dl + eft_cache%EFTOmegaV+ 0.5_dl*a*eft_cache%EFTOmegaP)*( -0.5_dl*( 1.0_dl +eft_cache%EFTOmegaV +2.0_dl*a*eft_cache%EFTOmegaP +a2*eft_cache%EFTOmegaPP )*eft_cache%adotoa**2 &
                      & -0.5_dl*(eft_cache%gpresm_t ) -0.5_dl*eft_cache%EFTLambda)
        eft_cache%Hdotdot = 1.0_dl/(1.0_dl + eft_cache%EFTOmegaV+ 0.5_dl*a*eft_cache%EFTOmegaP)*( -0.5_dl*a*eft_cache%adotoa**3*( 3.0_dl*eft_cache%EFTOmegaP +4.0_dl*a*eft_cache%EFTOmegaPP +a2*eft_cache%EFTOmegaPPP )&
                      & -eft_cache%adotoa*eft_cache%Hdot*( 1.0_dl +eft_cache%EFTOmegaV +3.5_dl*a*eft_cache%EFTOmegaP +1.5_dl*a2*eft_cache%EFTOmegaPP ) &
                      & -0.5_dl*( eft_cache%gpresdotm_t  +2.0_dl*eft_cache%adotoa*eft_cache%gpresm_t  +eft_cache%EFTLambdadot +2.0_dl*eft_cache%adotoa*eft_cache%EFTLambda ))

    end subroutine EFTCAMBFullModelComputeHubbleDer

    !----------------------------------------------------------------------------------------

end module EFTCAMB_abstract_model_full

!----------------------------------------------------------------------------------------
