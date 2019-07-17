!----------------------------------------------------------------------------------------
!
! This file is part of EFTCAMB.
!
! Copyright (C) 2013-2019 by the EFTCAMB authors
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
    use EFTCAMB_mixed_algorithms

    implicit none

    private

    public EFTCAMB_full_model

    !----------------------------------------------------------------------------------------
    !> This is the abstract type for EFTCAMB models.
    type, extends ( EFTCAMB_model ), abstract :: EFTCAMB_full_model

    contains

        procedure :: compute_dtauda            => EFTCAMBFullModelComputeDtauda        !< function that computes dtauda = 1/sqrt(a^2H^2).
        procedure :: compute_adotoa            => EFTCAMBFullModelComputeAdotoa        !< subroutine that computes adotoa = H.
        procedure :: compute_H_derivs          => EFTCAMBFullModelComputeHubbleDer     !< subroutine that computes the two derivatives wrt conformal time of H.
        procedure :: compute_additional_derivs => EFTCAMBFullModelComputeAdditionalDer !< subroutine that computes additional derivatives of the EFT functions, as needed by the stability calculation.

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

        real(dl) :: a2, temp

        a2   = a*a
        temp = 1.0_dl/(1.0_dl + eft_cache%EFTOmegaV+ a*eft_cache%EFTOmegaP)*( eft_cache%grhoa2 +2._dl*eft_cache%EFTc*a2 -eft_cache%EFTLambda*a2 )/3._dl

        EFTCAMBFullModelComputeDtauda = 1._dl/sqrt( temp )

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

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes additional derivatives of the EFT functions, as needed by the
    !> stability calculation.
    subroutine EFTCAMBFullModelComputeAdditionalDer( self, a, eft_par_cache, eft_cache )

        implicit none

        class(EFTCAMB_full_model)                    :: self          !< the base class
        real(dl), intent(in)                         :: a             !< the input scale factor.
        type(EFTCAMB_parameter_cache), intent(inout) :: eft_par_cache !< the EFTCAMB parameter cache that contains all the physical parameters.
        type(EFTCAMB_timestep_cache ), intent(inout) :: eft_cache     !< the EFTCAMB timestep cache that contains all the physical values.

        real(dl) :: dEFTcdot_da, err

        eft_cache%EFTGamma1PP = dfridr( dummy_helper_gamma_1, a, 0.03_dl*a, err )
        eft_cache%EFTGamma2PP = dfridr( dummy_helper_gamma_2, a, 0.03_dl*a, err )
        eft_cache%EFTGamma3PP = dfridr( dummy_helper_gamma_3, a, 0.03_dl*a, err )

        dEFTcdot_da = dfridr( dummy_helper_cP, a, 0.03_dl*a, err )

        eft_cache%EFTcdotdot = a*eft_cache%adotoa*dEFTcdot_da - 2._dl*eft_cache%adotoa*eft_cache%EFTcdot

        contains

            subroutine dummy_helper( temp_a, gamma1P, gamma2P, gamma3P, cP )

                real(dl) temp_a, gamma1P, gamma2P, gamma3P, cP

                type(EFTCAMB_timestep_cache)  :: temp_cache
                type(EFTCAMB_parameter_cache) :: temp_par_cache

                real(dl) :: a2
                real(dl) :: grhob_t, grhoc_t, grhor_t, grhog_t, grho_matter, gpres, grho, adotoa
                real(dl) :: EFT_grhonudot, EFT_gpinudot, adotdota, EFT_grhonu, EFT_gpinu, grhormass_t

                integer :: nu_i

                ! initialize:
                temp_cache     = eft_cache
                temp_par_cache = eft_par_cache
                ! initialize the cache:
                call temp_cache%initialize()
                ! prepare:
                a2 = temp_a*temp_a
                ! compute background densities of different species
                grhob_t = temp_par_cache%grhob/temp_a    ! 8\pi G_N \rho_b a^2: baryons background density
                grhoc_t = temp_par_cache%grhoc/temp_a    ! 8\pi G_N \rho_{cdm} a^2: cold dark matter background density
                grhor_t = temp_par_cache%grhornomass/a2  ! 8\pi G_N \rho_{\nu} a^2: massless neutrinos background density
                grhog_t = temp_par_cache%grhog/a2        ! 8\pi G_N \rho_{\gamma} a^2: radiation background density
                ! start computing background total pressure and total density:
                grho_matter  = grhob_t +grhoc_t
                gpres        = (grhog_t+grhor_t)/3._dl
                ! add radiation, massless neutrinos and Lambda to total background density:
                grho = grho_matter +grhor_t +grhog_t
                if ( temp_par_cache%Num_Nu_Massive /= 0 ) then
                    do nu_i = 1, temp_par_cache%Nu_mass_eigenstates
                        EFT_grhonu    = 0._dl
                        EFT_gpinu     = 0._dl
                        grhormass_t=temp_par_cache%grhormass(nu_i)/temp_a**2
                        call temp_par_cache%Nu_background(temp_a*temp_par_cache%nu_masses(nu_i),EFT_grhonu,EFT_gpinu)
                        temp_cache%grhonu_tot = temp_cache%grhonu_tot + grhormass_t*EFT_grhonu
                        temp_cache%gpinu_tot  = temp_cache%gpinu_tot  + grhormass_t*EFT_gpinu
                    end do
                end if
                grho  = grho  +temp_cache%grhonu_tot
                gpres = gpres +temp_cache%gpinu_tot
                ! start to fill the cache:
                temp_cache%a        = temp_a
                temp_cache%grhom_t  = grho
                temp_cache%gpresm_t = gpres
                temp_cache%grhob_t  = grhob_t
                temp_cache%grhoc_t  = grhoc_t
                temp_cache%grhor_t  = grhor_t
                temp_cache%grhog_t  = grhog_t
                ! compute the other things:
                ! background for full models. Here the expansion history is computed from the
                ! EFT functions. Hence compute them first and then compute the expansion history.
                call self%compute_background_EFT_functions( a, temp_par_cache , eft_cache )
                call self%compute_adotoa( a, temp_par_cache , eft_cache )
                ! store adotoa:
                adotoa   = temp_cache%adotoa
                ! compute massive neutrinos stuff:
                ! Massive neutrinos mod:
                temp_cache%grhonu_tot = 0._dl
                temp_cache%gpinu_tot  = 0._dl
                if ( temp_par_cache%Num_Nu_Massive /= 0 ) then
                    do nu_i = 1, temp_par_cache%Nu_mass_eigenstates
                        EFT_grhonu    = 0._dl
                        EFT_gpinu     = 0._dl
                        EFT_grhonudot = 0._dl
                        EFT_gpinudot  = 0._dl
                        grhormass_t=temp_par_cache%grhormass(nu_i)/a**2
                        call temp_par_cache%Nu_background(a*temp_par_cache%nu_masses(nu_i),EFT_grhonu,EFT_gpinu)
                        temp_cache%grhonu_tot = temp_cache%grhonu_tot + grhormass_t*EFT_grhonu
                        temp_cache%gpinu_tot  = temp_cache%gpinu_tot  + grhormass_t*EFT_gpinu
                        temp_cache%grhonudot_tot = temp_cache%grhonudot_tot + grhormass_t*(temp_par_cache%Nu_drho(a*temp_par_cache%nu_masses(nu_i) ,adotoa, EFT_grhonu)&
                            & -4._dl*adotoa*EFT_grhonu)
                        temp_cache%gpinudot_tot  = temp_cache%gpinudot_tot  + grhormass_t*(temp_par_cache%Nu_pidot(a*temp_par_cache%nu_masses(nu_i),adotoa, EFT_gpinu )&
                            & -4._dl*adotoa*EFT_gpinu)
                    end do
                end if
                ! compute pressure dot:
                temp_cache%gpresdotm_t = -4._dl*adotoa*( grhog_t+grhor_t )/3._dl +temp_cache%gpinudot_tot
                ! compute remaining quantities related to H:
                call self%compute_H_derivs( temp_a, temp_par_cache , temp_cache )
                ! store:
                adotdota = temp_cache%Hdot +temp_cache%adotoa**2
                ! compute backgrond EFT functions if model is designer:
                ! compute all other background stuff:
                call self%compute_rhoQPQ( temp_a, temp_par_cache , temp_cache )
                ! compute second order EFT functions:
                call self%compute_secondorder_EFT_functions( temp_a, temp_par_cache, temp_cache )
                ! copy out value:
                gamma1P = temp_cache%EFTGamma1P
                gamma2P = temp_cache%EFTGamma2P
                gamma3P = temp_cache%EFTGamma3P
                cP      = eft_cache%EFTcdot/(a*eft_cache%adotoa) + 2._dl*eft_cache%EFTc/a

            end subroutine dummy_helper

            ! Helper function to compute second derivative of gamma1
            function dummy_helper_gamma_1(a_temp)
              real(dl):: a_temp
              real(dl):: dummy_helper_gamma_1
              real(dl):: gamma1P, gamma2P, gamma3P, cP
              call dummy_helper(a_temp,gamma1P, gamma2P, gamma3P, cP)
              dummy_helper_gamma_1 = gamma1P
            end function dummy_helper_gamma_1

            ! Helper function to compute second derivative of gamma2
            function dummy_helper_gamma_2(a_temp)
              real(dl):: a_temp
              real(dl):: dummy_helper_gamma_2
              real(dl):: gamma1P, gamma2P, gamma3P, cP
              call dummy_helper(a_temp,gamma1P, gamma2P, gamma3P, cP)
              dummy_helper_gamma_2 = gamma2P
            end function dummy_helper_gamma_2

            ! Helper function to compute second derivative of gamma3
            function dummy_helper_gamma_3(a_temp)
              real(dl):: a_temp
              real(dl):: dummy_helper_gamma_3
              real(dl):: gamma1P, gamma2P, gamma3P, cP
              call dummy_helper(a_temp,gamma1P, gamma2P, gamma3P, cP)
              dummy_helper_gamma_3 = gamma3P
            end function dummy_helper_gamma_3

            ! Helper function to compute the derivative of cdot*a^2/m_0^2
            function dummy_helper_cP(a_temp)
              real(dl):: a_temp
              real(dl):: dummy_helper_cP
              real(dl):: gamma1P, gamma2P, gamma3P, cP
              call dummy_helper(a_temp,gamma1P, gamma2P, gamma3P, cP)
              dummy_helper_cP = cP
            end function dummy_helper_cP

    end subroutine EFTCAMBFullModelComputeAdditionalDer

! ---------------------------------------------------------------------------------------------

end module EFTCAMB_abstract_model_full

!----------------------------------------------------------------------------------------
