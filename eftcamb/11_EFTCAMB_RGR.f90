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

!> @file 11_EFTCAMB_RGR.f90
!! This file contains the RGR algorithm. This operates on general EFT models.


!----------------------------------------------------------------------------------------
!> This module contains the RGR algorithm. This operates on general EFT models.

!> @author Bin Hu, Marco Raveri

module EFTCAMB_ReturnToGR

    use precision
    use ModelParams
    use MassiveNu
    use EFTCAMB_abstract_model
    use EFTCAMB_abstract_model_full
    use EFTCAMB_abstract_model_designer
    use EFTDef

    implicit none

contains

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that prints feedback for the RGR module. Prints something only if
    !! EFTCAMB_feedback_level is greater than 1.
    subroutine EFTCAMBReturnToGR_feedback()

        implicit none

        real(dl) :: eft_functions(21)
        integer  :: i

        if ( CP%EFTCAMB%EFTCAMB_feedback_level > 1 ) then
            write(*,'(a,F12.10)') 'EFTCAMB Return to GR time: ', CP%EFTCAMB%EFTCAMB_turn_on_time
        end if

        if ( CP%EFTCAMB%EFTCAMB_feedback_level > 2 ) then

            call EFTCAMBReturnToGR_functions( CP%EFTCAMB%EFTCAMB_turn_on_time, eft_functions )
            write(*,'(a)') 'EFT functions at RGR time: '
            do i=1, 21
                write(*,*) eft_functions(i)
            end do

        end if

    end subroutine

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes the return to GR of a theory.
    subroutine EFTCAMBReturnToGR()

        implicit none

        real(dl) :: eft_functions_t1(21), eft_functions_t2(21)
        real(dl) :: a_initial, a_final, log_a_initial, log_a_final, log_a_used, a_used
        real(dl) :: last_a
        integer  :: i

        ! initial and final scale factors:
        a_initial = CP%EFTCAMB%EFTCAMB_turn_on_time
        a_final   = 1._dl
        ! convert to logarithm of the scale factor:
        log_a_initial = log10( a_initial )
        log_a_final   = log10( a_final   )
        ! compute RGR functions at initial time:
        log_a_used       = log_a_initial
        a_used           = 10._dl**log_a_used
        eft_functions_t1 = 0._dl
        eft_functions_t2 = 0._dl
        call EFTCAMBReturnToGR_functions( a_used, eft_functions_t2 )
        last_a           = a_used

        ! check if already above threshold:
        if ( any( eft_functions_t2 > 0._dl ) ) then
            CP%EFTCAMB%EFTCAMB_turn_on_time = last_a
            call EFTCAMBReturnToGR_feedback()
            return
        end if

        ! loop over the logarithm of the scale factor:
        do i=2, EFT_RGR_num_points
            ! initialize:
            log_a_used       = log_a_initial +real( i-1 )/real( EFT_RGR_num_points -1 )*( log_a_final -log_a_initial )
            a_used           = 10._dl**log_a_used
            eft_functions_t1 = 0._dl
            ! compute RGR functions:
            call EFTCAMBReturnToGR_functions( a_used, eft_functions_t1 )

            ! check if above threshold:
            if ( any( eft_functions_t1 > 0._dl ) ) then
                CP%EFTCAMB%EFTCAMB_turn_on_time = last_a
                call EFTCAMBReturnToGR_feedback()
                return
            end if

            ! check if the threshold has been crossed:
            eft_functions_t2 = eft_functions_t1*eft_functions_t2
            if ( any( eft_functions_t2 < 0._dl ) ) then
                CP%EFTCAMB%EFTCAMB_turn_on_time = last_a
                call EFTCAMBReturnToGR_feedback()
                return
            else
                ! prepare everything to proceed to the next step:
                last_a           = a_used
                eft_functions_t2 = eft_functions_t1
            end if
        end do

        call EFTCAMBReturnToGR_feedback()

    end subroutine EFTCAMBReturnToGR

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes the EFT functions as we need them for the RGR function.
    subroutine EFTCAMBReturnToGR_functions( a, eft_functions )

        implicit none

        real(dl), intent(in)  :: a                   !< scale factor at which EFT funcitons are computed.
        real(dl), intent(out) :: eft_functions(21)   !< vector containing the values of the EFT functions minus their GR value.

        type(EFTCAMB_timestep_cache) :: eft_cache
        real(dl) :: a2, adotoa, EFT_grhonu, EFT_gpinu, EFT_grhonudot, EFT_gpinudot, grhormass_t
        integer  :: nu_i

        ! prepare:
        a2 = a*a
        ! initialize the cache:
        call eft_cache%initialize()
        ! start filling:
        eft_cache%a = a
        ! compute background densities of different species
        eft_cache%grhob_t = grhob/a         ! 8\pi G_N \rho_b a^2: baryons background density
        eft_cache%grhoc_t = grhoc/a         ! 8\pi G_N \rho_{cdm} a^2: cold dark matter background density
        eft_cache%grhor_t = grhornomass/a2  ! 8\pi G_N \rho_{\nu} a^2: massless neutrinos background density
        eft_cache%grhog_t = grhog/a2        ! 8\pi G_N \rho_{\gamma} a^2: radiation background density
        ! Massive neutrinos terms:
        if ( CP%Num_Nu_Massive /= 0 ) then
            do nu_i = 1, CP%Nu_mass_eigenstates
                EFT_grhonu    = 0._dl
                EFT_gpinu     = 0._dl
                grhormass_t=grhormass(nu_i)/a**2
                call Nu_background(a*nu_masses(nu_i),EFT_grhonu,EFT_gpinu)
                eft_cache%grhonu_tot = eft_cache%grhonu_tot + grhormass_t*EFT_grhonu
                eft_cache%gpinu_tot  = eft_cache%gpinu_tot  + grhormass_t*EFT_gpinu
            end do
        end if
        ! assemble total densities and pressure:
        eft_cache%grhom_t  = eft_cache%grhob_t +eft_cache%grhoc_t +eft_cache%grhor_t +eft_cache%grhog_t +eft_cache%grhonu_tot
        eft_cache%gpresm_t = (+eft_cache%grhor_t +eft_cache%grhog_t)/3._dl +eft_cache%gpinu_tot
        ! compute the other things:
        select type ( model => CP%EFTCAMB%model )
            ! compute the background and the background EFT functions.
            class is ( EFTCAMB_full_model )
            ! background for full models. Here the expansion history is computed from the
            ! EFT functions. Hence compute them first and then compute the expansion history.
            call CP%EFTCAMB%model%compute_background_EFT_functions( a, CP%eft_par_cache , eft_cache )
            call CP%EFTCAMB%model%compute_adotoa( a, CP%eft_par_cache , eft_cache )
            class is ( EFTCAMB_designer_model )
            ! background for designer models. Here the expansion history is parametrized
            ! and does not depend on the EFT functions. Hence compute first the expansion history
            ! and then the EFT functions.
            call CP%EFTCAMB%model%compute_adotoa( a, CP%eft_par_cache , eft_cache )
        end select
        ! store adotoa:
        adotoa   = eft_cache%adotoa
        ! compute massive neutrinos stuff:
        ! Massive neutrinos mod:
        eft_cache%grhonu_tot = 0._dl
        eft_cache%gpinu_tot  = 0._dl
        if ( CP%Num_Nu_Massive /= 0 ) then
            do nu_i = 1, CP%Nu_mass_eigenstates
                EFT_grhonu    = 0._dl
                EFT_gpinu     = 0._dl
                EFT_grhonudot = 0._dl
                EFT_gpinudot  = 0._dl
                grhormass_t=grhormass(nu_i)/a**2
                call Nu_background(a*nu_masses(nu_i),EFT_grhonu,EFT_gpinu)
                eft_cache%grhonu_tot = eft_cache%grhonu_tot + grhormass_t*EFT_grhonu
                eft_cache%gpinu_tot  = eft_cache%gpinu_tot  + grhormass_t*EFT_gpinu
                eft_cache%grhonudot_tot = eft_cache%grhonudot_tot + grhormass_t*(Nu_drho(a*nu_masses(nu_i) ,adotoa, EFT_grhonu)&
                    & -4._dl*adotoa*EFT_grhonu)
                eft_cache%gpinudot_tot  = eft_cache%gpinudot_tot  + grhormass_t*(Nu_pidot(a*nu_masses(nu_i),adotoa, EFT_gpinu )&
                    & -4._dl*adotoa*EFT_gpinu)
            end do
        end if
        ! compute pressure dot:
        eft_cache%gpresdotm_t = -4._dl*adotoa*( eft_cache%grhog_t +eft_cache%grhor_t )/3._dl +eft_cache%gpinudot_tot
        ! compute remaining quantities related to H:
        call CP%EFTCAMB%model%compute_H_derivs( a, CP%eft_par_cache , eft_cache )
        ! compute backgrond EFT functions if model is designer:
        select type ( model => CP%EFTCAMB%model )
            class is ( EFTCAMB_designer_model )
            call CP%EFTCAMB%model%compute_background_EFT_functions( a, CP%eft_par_cache , eft_cache )
        end select
        ! compute all other background stuff:
        call CP%EFTCAMB%model%compute_rhoQPQ( a, CP%eft_par_cache , eft_cache )
        ! compute second order EFT functions:
        call CP%EFTCAMB%model%compute_secondorder_EFT_functions( a, CP%eft_par_cache , eft_cache )

        ! get the EFT functions:
        eft_functions( 1) = abs( eft_cache%EFTOmegaV           )
        eft_functions( 2) = abs( a*adotoa*eft_cache%EFTOmegaP  )
        eft_functions( 3) = 0._dl
        eft_functions( 4) = 0._dl
        eft_functions( 5) = abs( eft_cache%EFTc/a2             )
        eft_functions( 6) = abs( eft_cache%EFTLambda/a2 +CP%eft_par_cache%grhov )
        eft_functions( 7) = abs( eft_cache%EFTcdot/a2          )
        eft_functions( 8) = abs( eft_cache%EFTLambdadot/a2     )
        eft_functions( 9) = abs( eft_cache%EFTGamma1V )
        eft_functions(10) = abs( eft_cache%EFTGamma1P )
        eft_functions(11) = abs( eft_cache%EFTGamma2V )
        eft_functions(12) = abs( eft_cache%EFTGamma2P )
        eft_functions(13) = abs( eft_cache%EFTGamma3V )
        eft_functions(14) = abs( eft_cache%EFTGamma3P )
        eft_functions(15) = abs( eft_cache%EFTGamma4V )
        eft_functions(16) = abs( eft_cache%EFTGamma4P )
        eft_functions(17) = 0._dl
        eft_functions(18) = abs( eft_cache%EFTGamma5V )
        eft_functions(19) = abs( eft_cache%EFTGamma5P )
        eft_functions(20) = abs( eft_cache%EFTGamma6V )
        eft_functions(21) = abs( eft_cache%EFTGamma6P )

        ! subtract the GR threshold:
        eft_functions = eft_functions -EFTtoGR

    end subroutine EFTCAMBReturnToGR_functions

    !----------------------------------------------------------------------------------------

end module EFTCAMB_ReturnToGR

!----------------------------------------------------------------------------------------
