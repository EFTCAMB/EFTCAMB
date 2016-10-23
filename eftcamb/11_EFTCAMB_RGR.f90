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
    use EFT_def
    use EFTCAMB_cache
    use EFTCAMB_abstract_model
    use EFTCAMB_abstract_model_full
    use EFTCAMB_abstract_model_designer

    implicit none

    private

    public EFTCAMBReturnToGR_feedback, EFTCAMBReturnToGR, EFTCAMBReturnToGR_functions

contains

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that prints feedback for the RGR module. Prints something only if
    !! EFTCAMB_feedback_level is greater than 1.
    subroutine EFTCAMBReturnToGR_feedback( feedback_level, input_model, params_cache, RGR_time )

        implicit none

        integer                      , intent(in) :: feedback_level  !< feedback level for the RGR results reporter. 0=no feedback; 1=some feedback; 2=a lot of feedback.
        class(EFTCAMB_model)         , intent(in) :: input_model     !< the EFTCAMB model for which the code is computing the RGR time.
        type(EFTCAMB_parameter_cache), intent(in) :: params_cache    !< a EFTCAMB parameter cache containing cosmological parameters.
        real(dl)                     , intent(in) :: RGR_time        !< input RGR time. For better performances the code does not compute it in this function.

        real(dl) :: eft_functions(21)
        integer  :: i

        if ( feedback_level > 1 ) then
            write(*,'(a)') "***************************************************************"
            write(*,'(a,F12.10)') ' EFTCAMB Return to GR time: ', RGR_time
        end if

        if ( feedback_level > 2 ) then

            call EFTCAMBReturnToGR_functions( RGR_time, input_model, params_cache, eft_functions )

            write(*,'(a)')
            write(*,'(a)') ' EFT functions at RGR time: '

            if ( eft_functions(1 ) .gt. 0._dl ) write(*,'(a18,E13.4)') '   OmegaV       = ', eft_functions(1 ) +EFTtoGR
            if ( eft_functions(2 ) .gt. 0._dl ) write(*,'(a18,E13.4)') '   OmegaP       = ', eft_functions(2 ) +EFTtoGR
            if ( eft_functions(3 ) .gt. 0._dl ) write(*,'(a18,E13.4)') '   OmegaPP      = ', eft_functions(3 ) +EFTtoGR
            if ( eft_functions(4 ) .gt. 0._dl ) write(*,'(a18,E13.4)') '   OmegaPPP     = ', eft_functions(4 ) +EFTtoGR
            if ( eft_functions(5 ) .gt. 0._dl ) write(*,'(a18,E13.4)') '   EFTc         = ', eft_functions(5 ) +EFTtoGR
            if ( eft_functions(6 ) .gt. 0._dl ) write(*,'(a18,E13.4)') '   EFTLambda    = ', eft_functions(6 ) +EFTtoGR
            if ( eft_functions(7 ) .gt. 0._dl ) write(*,'(a18,E13.4)') '   EFTcdot      = ', eft_functions(7 ) +EFTtoGR
            if ( eft_functions(8 ) .gt. 0._dl ) write(*,'(a18,E13.4)') '   EFTLambdadot = ', eft_functions(8 ) +EFTtoGR
            if ( eft_functions(9 ) .gt. 0._dl ) write(*,'(a18,E13.4)') '   EFTGamma1V   = ', eft_functions(9 ) +EFTtoGR
            if ( eft_functions(10) .gt. 0._dl ) write(*,'(a18,E13.4)') '   EFTGamma1P   = ', eft_functions(10) +EFTtoGR
            if ( eft_functions(11) .gt. 0._dl ) write(*,'(a18,E13.4)') '   EFTGamma2V   = ', eft_functions(11) +EFTtoGR
            if ( eft_functions(12) .gt. 0._dl ) write(*,'(a18,E13.4)') '   EFTGamma2P   = ', eft_functions(12) +EFTtoGR
            if ( eft_functions(13) .gt. 0._dl ) write(*,'(a18,E13.4)') '   EFTGamma3V   = ', eft_functions(13) +EFTtoGR
            if ( eft_functions(14) .gt. 0._dl ) write(*,'(a18,E13.4)') '   EFTGamma3P   = ', eft_functions(14) +EFTtoGR
            if ( eft_functions(15) .gt. 0._dl ) write(*,'(a18,E13.4)') '   EFTGamma4V   = ', eft_functions(15) +EFTtoGR
            if ( eft_functions(16) .gt. 0._dl ) write(*,'(a18,E13.4)') '   EFTGamma4P   = ', eft_functions(16) +EFTtoGR
            if ( eft_functions(18) .gt. 0._dl ) write(*,'(a18,E13.4)') '   EFTGamma5V   = ', eft_functions(18) +EFTtoGR
            if ( eft_functions(17) .gt. 0._dl ) write(*,'(a18,E13.4)') '   EFTGamma4PP  = ', eft_functions(17) +EFTtoGR
            if ( eft_functions(19) .gt. 0._dl ) write(*,'(a18,E13.4)') '   EFTGamma5P   = ', eft_functions(19) +EFTtoGR
            if ( eft_functions(20) .gt. 0._dl ) write(*,'(a18,E13.4)') '   EFTGamma6V   = ', eft_functions(18) +EFTtoGR
            if ( eft_functions(21) .gt. 0._dl ) write(*,'(a18,E13.4)') '   EFTGamma6P   = ', eft_functions(19) +EFTtoGR

        end if

    end subroutine

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes the return to GR of a theory.
    subroutine EFTCAMBReturnToGR( input_model, params_cache, initial_time, RGR_time )

        implicit none

        class(EFTCAMB_model)         , intent(in)  :: input_model     !< the EFTCAMB model for which the code is computing the RGR time.
        type(EFTCAMB_parameter_cache), intent(in)  :: params_cache    !< a EFTCAMB parameter cache containing cosmological parameters.
        real(dl)                     , intent(in)  :: initial_time    !< initial scale factor at which the code starts to look for the RGR of the theory.
        real(dl)                     , intent(out) :: RGR_time        !< output value of the RGR time.

        real(dl) :: eft_functions(21)
        real(dl) :: a_initial, a_final, log_a_initial, log_a_final, log_a_used, a_used
        integer  :: i

        ! initialize:
        RGR_time = initial_time

        ! initial and final scale factors:
        a_initial = initial_time
        a_final   = 1._dl
        ! convert to logarithm of the scale factor:
        log_a_initial = log10( a_initial )
        log_a_final   = log10( a_final   )

        ! loop over the logarithm of the scale factor:
        do i=1, EFT_RGR_num_points

            ! initialize:
            log_a_used       = log_a_initial +real( i-1 )/real( EFT_RGR_num_points -1 )*( log_a_final -log_a_initial )
            a_used           = 10._dl**log_a_used
            eft_functions    = 0._dl
            ! compute RGR functions:
            call EFTCAMBReturnToGR_functions( a_used, input_model, params_cache, eft_functions )

            ! check if above threshold:
            if ( any( eft_functions > 0._dl ) ) then
                RGR_time = a_used
                return
            end if

        end do

        RGR_time = 1.1_dl*a_final

    end subroutine EFTCAMBReturnToGR

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes the EFT functions as we need them for the RGR function.
    subroutine EFTCAMBReturnToGR_functions( a, input_model, params_cache, eft_functions )

        implicit none

        real(dl)                     , intent(in)    :: a                  !< scale factor at which EFT funcitons are computed.
        class(EFTCAMB_model)         , intent(in)    :: input_model        !< the EFTCAMB model for which the code is computing the RGR time.
        type(EFTCAMB_parameter_cache)                :: params_cache       !< a EFTCAMB parameter cache containing cosmological parameters.
        real(dl)                     , intent(out)   :: eft_functions(21)  !< vector containing the values of the EFT functions minus their GR value.

        type(EFTCAMB_timestep_cache) :: eft_cache

        real(dl) :: a2, adotoa, grhonu, gpinu, grhonudot, gpinudot, grhormass_t
        integer  :: nu_i

        ! prepare:
        a2 = a*a
        ! initialize the cache:
        call eft_cache%initialize()
        ! start filling:
        eft_cache%a = a
        ! compute background densities of different species
        eft_cache%grhob_t = params_cache%grhob/a         ! 8\pi G_N \rho_b a^2: baryons background density
        eft_cache%grhoc_t = params_cache%grhoc/a         ! 8\pi G_N \rho_{cdm} a^2: cold dark matter background density
        eft_cache%grhor_t = params_cache%grhornomass/a2  ! 8\pi G_N \rho_{\nu} a^2: massless neutrinos background density
        eft_cache%grhog_t = params_cache%grhog/a2        ! 8\pi G_N \rho_{\gamma} a^2: radiation background density
        ! Massive neutrinos terms:
        if ( params_cache%Num_Nu_Massive /= 0 ) then
            do nu_i = 1, params_cache%Nu_mass_eigenstates
                grhonu    = 0._dl
                gpinu     = 0._dl
                grhormass_t = params_cache%grhormass(nu_i)/a**2
                call params_cache%Nu_background(a*params_cache%nu_masses(nu_i), grhonu, gpinu)
                eft_cache%grhonu_tot = eft_cache%grhonu_tot + grhormass_t*grhonu
                eft_cache%gpinu_tot  = eft_cache%gpinu_tot  + grhormass_t*gpinu
            end do
        end if
        ! assemble total densities and pressure:
        eft_cache%grhom_t  = eft_cache%grhob_t +eft_cache%grhoc_t +eft_cache%grhor_t +eft_cache%grhog_t +eft_cache%grhonu_tot
        eft_cache%gpresm_t = (+eft_cache%grhor_t +eft_cache%grhog_t)/3._dl +eft_cache%gpinu_tot
        ! compute the other things:
        select type ( model => input_model )
            ! compute the background and the background EFT functions.
            class is ( EFTCAMB_full_model )
            ! background for full models. Here the expansion history is computed from the
            ! EFT functions. Hence compute them first and then compute the expansion history.
            call input_model%compute_background_EFT_functions( a, params_cache , eft_cache )
            call input_model%compute_adotoa( a, params_cache , eft_cache )
            class is ( EFTCAMB_designer_model )
            ! background for designer models. Here the expansion history is parametrized
            ! and does not depend on the EFT functions. Hence compute first the expansion history
            ! and then the EFT functions.
            call input_model%compute_adotoa( a, params_cache , eft_cache )
        end select
        ! store adotoa:
        adotoa   = eft_cache%adotoa
        ! compute massive neutrinos stuff:
        ! Massive neutrinos mod:
        eft_cache%grhonu_tot = 0._dl
        eft_cache%gpinu_tot  = 0._dl
        if ( params_cache%Num_Nu_Massive /= 0 ) then
            do nu_i = 1, params_cache%Nu_mass_eigenstates
                grhonu    = 0._dl
                gpinu     = 0._dl
                grhonudot = 0._dl
                gpinudot  = 0._dl
                grhormass_t= params_cache%grhormass(nu_i)/a**2
                call params_cache%Nu_background(a*params_cache%nu_masses(nu_i),grhonu,gpinu)
                eft_cache%grhonu_tot = eft_cache%grhonu_tot + grhormass_t*grhonu
                eft_cache%gpinu_tot  = eft_cache%gpinu_tot  + grhormass_t*gpinu
                eft_cache%grhonudot_tot = eft_cache%grhonudot_tot + grhormass_t*( params_cache%Nu_drho(a*params_cache%nu_masses(nu_i), adotoa, grhonu)&
                    & -4._dl*adotoa*grhonu )
                eft_cache%gpinudot_tot  = eft_cache%gpinudot_tot  + grhormass_t*( params_cache%Nu_pidot(a*params_cache%nu_masses(nu_i),adotoa, gpinu )&
                    & -4._dl*adotoa*gpinu )
            end do
        end if
        ! compute pressure dot:
        eft_cache%gpresdotm_t = -4._dl*adotoa*( eft_cache%grhog_t +eft_cache%grhor_t )/3._dl +eft_cache%gpinudot_tot
        ! compute remaining quantities related to H:
        call input_model%compute_H_derivs( a, params_cache, eft_cache )
        ! compute backgrond EFT functions if model is designer:
        select type ( model => input_model )
            class is ( EFTCAMB_designer_model )
            call input_model%compute_background_EFT_functions( a, params_cache, eft_cache )
        end select
        ! compute all other background stuff:
        call input_model%compute_rhoQPQ( a, params_cache, eft_cache )
        ! compute second order EFT functions:
        call input_model%compute_secondorder_EFT_functions( a, params_cache, eft_cache )

        ! get the EFT functions:
        eft_functions( 1) = abs( eft_cache%EFTOmegaV           )
        eft_functions( 2) = abs( a*adotoa*eft_cache%EFTOmegaP  )
        eft_functions( 3) = 0._dl
        eft_functions( 4) = 0._dl
        eft_functions( 5) = abs( eft_cache%EFTc/a2             )
        eft_functions( 6) = abs( eft_cache%EFTLambda/a2 +params_cache%grhov )
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
