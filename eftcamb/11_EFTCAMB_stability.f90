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

!> @file 11_EFTCAMB_stability.f90
!! This file contains the stability detection algorithm of EFTCAMB.


!----------------------------------------------------------------------------------------
!> This module contains the stability detection algorithm of EFTCAMB.

!> @author Bin Hu, Marco Raveri, Simone Peirone

module EFTCAMB_stability

    use precision
    use EFT_def
    use EFTCAMB_cache
    use EFTCAMB_abstract_model
    use EFTCAMB_abstract_model_full
    use EFTCAMB_abstract_model_designer
    use EFTCAMB_main

    implicit none

    private

    public EFTCAMB_Stability_Check, EFTTestStability, EFTStability_cleanup, EFTStabilityComputation

    ! storage for some utility values that needs to be stored for the stability check.
    real(dl), save :: PastA1 = 0._dl
    real(dl), save :: PastAT = 0._dl

contains

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that tests the stability of a theory in a time range.
    !! To ensure best time coverage scans with three different strategies.
    subroutine EFTCAMB_Stability_Check( success, input_EFTCAMB, params_cache, astart, aend, k_max )

        implicit none

        logical                      , intent(out)   :: success         !< Output of the subroutine. Tells whether the model is found stable or not.
        class(EFTCAMB)               , intent(in)    :: input_EFTCAMB   !< the EFTCAMB object for which the code is computing stability.
        type(EFTCAMB_parameter_cache), intent(inout) :: params_cache    !< the EFTCAMB parameter cache that contains all the physical parameters.
        real(dl)                     , intent(in)    :: astart          !< Initial scale factor.
        real(dl)                     , intent(in)    :: aend            !< Final scale factor.
        real(dl)                     , intent(inout) :: k_max           !< the input maximum k mode at which stability is computed.

        ! parameters of the stability sampler:
        integer , parameter :: indMax            = 10000   ! Number of points sampled.
        real(dl), parameter :: LogSamplingScale  = -10._dl ! Where to start with the log sampling

        real(dl) :: Atest, y
        integer  :: ind

        ! 0) initial feedback:
        if ( input_EFTCAMB%EFTCAMB_feedback_level > 1 ) then
            write(*,'(a)') '***************************************************************'
            write(*,'(a)') ' EFTCAMB: checking stability of the theory'
        end if
        if ( input_EFTCAMB%EFTCAMB_feedback_level > 2 ) then
            write(*,'(a)')
        end if

        ! 1) stability code:
        success = .true.

        !    - linear sampling:
        call EFTStability_cleanup()
        do ind=1, indMax
            Atest = astart + REAL(ind-1)*(aend-astart)/REAL(indMax-1)
            success = EFTTestStability( Atest, k_max, input_EFTCAMB, params_cache )
            if ( .not. success ) then
                if ( input_EFTCAMB%EFTCAMB_feedback_level > 2 ) then
                    write(*,*)
                    write(*,'(a,E14.4)') '   Instability detected at a =', Atest
                end if
                return
            end if
        end do

        !    - log sampling close to astart:
        call EFTStability_cleanup()
        do ind=1, indMax
            y = LogSamplingScale + REAL(ind-1)*(0._dl-LogSamplingScale)/REAL(indMax-1)
            Atest = astart +(aend-astart)*10._dl**y
            success = EFTTestStability( Atest, k_max, input_EFTCAMB, params_cache )
            if ( .not. success ) then
                if ( input_EFTCAMB%EFTCAMB_feedback_level > 2 ) then
                    write(*,*)
                    write(*,'(a,E14.4)') '   Instability detected at a =', Atest
                end if
                return
            end if
        end do

        !    - log sampling close to aend:
        call EFTStability_cleanup()
        do ind=1, indMax
            Atest = aend +(astart-aend)*10._dl**y
            success = EFTTestStability( Atest, k_max, input_EFTCAMB, params_cache )
            if ( .not. success ) then
                if ( input_EFTCAMB%EFTCAMB_feedback_level > 2 ) then
                    write(*,*)
                    write(*,'(a,E14.4)') '   Instability detected at a =', Atest
                end if
                return
            end if
        end do

        ! 2) final feedback:
        if ( input_EFTCAMB%EFTCAMB_feedback_level > 1 ) then
            write(*,'(a)') ' EFTCAMB: theory stable'
        end if

    end subroutine EFTCAMB_Stability_Check

    ! ---------------------------------------------------------------------------------------------
    !> Function that fills the caches to check the stability of the theory.
    function EFTTestStability( a, k_max, input_EFTCAMB, params_cache )

        implicit none

        real(dl)                     , intent(in)    :: a                       !< the input scale factor.
        real(dl)                     , intent(inout) :: k_max                   !< the input maximum k mode at which stability is computed.
        class(EFTCAMB)               , intent(in)    :: input_EFTCAMB           !< the EFTCAMB object for which the code is computing stability.
        type(EFTCAMB_parameter_cache), intent(inout) :: params_cache            !< the EFTCAMB parameter cache that contains all the physical parameters.
        logical                                      :: EFTTestStability        !< Logical value returned by the function. If the model is stable this is True, otherwise False.

        ! Definitions of variables:
        type(EFTCAMB_timestep_cache ) :: eft_cache
        logical  :: EFT_HaveNan_parameter, EFT_HaveNan_timestep
        real(dl) :: EFT_instability_rate, tempk, temp1, temp2, temp3, temp4, temp5
        integer  :: ind_max, ind

        ! Stability check initialization
        EFTTestStability = .true.
        ! reset the time-step cache:
        call eft_cache%initialize()
        ! fill it:
        call EFTStabilityComputation( a, input_EFTCAMB%model, params_cache, eft_cache )
        ! protect against k_max too small:
        if ( k_max < 0.1_dl ) k_max = 0.1_dl

        ! check stability of the theory:

        ! 1) everything inside the parameter cache should not be a NaN:
        call params_cache%is_nan( EFT_HaveNan_parameter )
        if ( EFT_HaveNan_parameter ) then
            EFTTestStability = .false.
            if ( input_EFTCAMB%EFTCAMB_feedback_level > 0 ) write(*,'(a)') '   Model has Nan in the parameter cache'
            return
        end if
        ! 2) everything inside the time-step cache should not be a NaN:
        call eft_cache%is_nan( EFT_HaveNan_timestep )
        if ( EFT_HaveNan_timestep ) then
            EFTTestStability = .false.
            if ( input_EFTCAMB%EFTCAMB_feedback_level > 0 ) write(*,'(a)') '   Model has Nan in the timestep cache'
            return
        end if
        ! 3) enforce mathematical stability:
        if ( input_EFTCAMB%EFT_mathematical_stability ) then

            ! 1- the A coefficient should not change sign in time and in k, i.e. it shall not be zero.
            !    This is the strongest stability constraint since violating it would violate the mathematical
            !    consistency of the pi field equation.
            !    The first condition is A1/=0. Implemented by detecting sign changes in A1.
            if ( eft_cache%EFTpiA1*PastA1 < 0._dl ) then
                EFTTestStability = .false.
                if ( input_EFTCAMB%EFTCAMB_feedback_level > 0 ) write(*,'(a)') '   Mathematical instability: A is zero in time'
            end if
            PastA1 = eft_cache%EFTpiA1
            !    The second one is the condition on k.
            if ( (eft_cache%EFTpiA1 > 0 .and. eft_cache%EFTpiA1 + k_max**2*eft_cache%EFTpiA2 < 0) .or. &
                &(eft_cache%EFTpiA1 < 0 .and. eft_cache%EFTpiA1 + k_max**2*eft_cache%EFTpiA2 > 0) ) then
                EFTTestStability = .false.
                if ( input_EFTCAMB%EFTCAMB_feedback_level > 0 ) write(*,'(a)') '   Mathematical instability: A is zero in k'
            end if

            ! 2- the AT coefficient should not change sign in time, i.e. it shall not be zero.
            !    This is the second strongest stability constraint since violating it would
            !    violate the mathematical consistency of the tensor perturbation equation.
            !    Implemented by detecting sign changes in AT.
            if ( eft_cache%EFTAT*PastAT < 0._dl ) then
                EFTTestStability = .false.
                if ( input_EFTCAMB%EFTCAMB_feedback_level > 0 ) write(*,'(a)') '   Mathematical instability:  AT is zero in time'
            end if
            PastAT = eft_cache%EFTAT

            ! 3- we do not want (fast) growing exponential modes.
            !    This condition prevents the pi field from growing exponentially and destroying everything.
            !    Even though this condition is neither completely related to physics nor mathematics,
            !    violating it would completely mess up cosmological observables.

            !    This is the maximum allowed rate of instability. Units shall be Mpc^-1.
            EFT_instability_rate = 0._dl

            !    This condition needs to be tested in k. Sample in k.
            ind_max = 10

            do ind = 1, ind_max
                ! kmode to test. Linear sampling. Should suffice... (??)
                tempk = 0._dl + REAL(ind-1)*(k_max)/REAL(ind_max-1)
                ! vaule that discriminates between different cases:
                temp1 = (eft_cache%EFTpiB1 +eft_cache%EFTpiB2*tempk**2)
                temp2 = (eft_cache%EFTpiA1 +eft_cache%EFTpiA2*tempk**2)
                temp3 = temp1**2 -4._dl*temp2*(eft_cache%EFTpiC +eft_cache%EFTpiD1*tempk**2 + eft_cache%EFTpiD2*tempk**4)

                ! case 1:
                if ( temp3 > 0._dl .and. temp2 /= 0._dl ) then
                    temp4 = +0.5_dl*(-temp1 +sqrt(temp3))/temp2
                    temp5 = +0.5_dl*(-temp1 -sqrt(temp3))/temp2
                    if ( temp4>EFT_instability_rate .or. temp5>EFT_instability_rate ) then
                        EFTTestStability = .false.
                        if ( input_EFTCAMB%EFTCAMB_feedback_level > 0 ) then
                            write(*,'(a,E11.4)')       '   Mathematical instability: growing exponential at k =', tempk
                            write(*,'(a,E11.4,E11.4)') '      Rate of instability: ', temp4, temp5
                        end if
                        exit
                    end if
                ! case 2:
                else if ( temp2 /= 0._dl ) then
                    temp4 = -0.5_dl*temp1/temp2
                    if ( temp4>EFT_instability_rate ) then
                        EFTTestStability = .false.
                        if ( input_EFTCAMB%EFTCAMB_feedback_level > 0 ) then
                            write(*,'(a,E11.4)')       '   Mathematical instability: growing exponential at k =', tempk
                            write(*,'(a,E11.4,E11.4)') '      Rate of instability: ', temp4
                        end if
                        exit
                    end if
                end if

            end do

        end if
        ! 4) enforce model specific priors:
        if ( input_EFTCAMB%EFT_additional_priors ) then
            EFTTestStability = input_EFTCAMB%model%additional_model_stability( a, params_cache, eft_cache )
            if ( .not. EFTTestStability ) then
                if ( input_EFTCAMB%EFTCAMB_feedback_level > 0 ) write(*,'(a)') '   Model specific stability criteria are not met'
            end if
        end if
        ! 5) enforce physical viability:
        if ( input_EFTCAMB%EFT_physical_stability ) then

            ! the present conditions extend up to Horndeski. Enforce that:
            if ( (eft_cache%EFTGamma6V /= 0._dl) .or.      &
                & ( (eft_cache%EFTGamma3V + eft_cache%EFTGamma4V) /= 0._dl) ) then
                write(*,'(a)') '   EFTCAMB WARNING: stability for model beyond GLPV has not been worked out.'
                write(*,'(a)') '      It will be added in a future release.'
                write(*,'(a)') '      If you want to run this model disable EFT_physical_stability.'
                EFTTestStability = .false.
                return
            end if

            ! 1- Positive gravitational constant:
            if ( 1._dl +eft_cache%EFTOmegaV <= 0 ) then
                EFTTestStability = .false.
                if ( input_EFTCAMB%EFTCAMB_feedback_level > 0 ) write(*,'(a,E11.4)') '   Physical instability: negative gravitational constant = ', 1._dl +eft_cache%EFTOmegaV
            end if

            ! 2- Ghost condition:
            if ( eft_cache%EFT_kinetic < 0._dl ) then
                EFTTestStability = .false.
                if ( input_EFTCAMB%EFTCAMB_feedback_level > 0 ) write(*,'(a,E11.4)') '   Physical instability: ghost instability. Kinetic term = ', eft_cache%EFT_kinetic
            end if

            ! 3- Gradient instability:
            if ( eft_cache%EFT_gradient < 0._dl ) then
                EFTTestStability = .false.
                if ( input_EFTCAMB%EFTCAMB_feedback_level > 0 ) write(*,'(a,E11.4)') '   Physical instability: gradient instability. Gradient term = ', eft_cache%EFT_gradient
            end if

            ! 4- No tensor ghosts:
            if ( eft_cache%EFTAT < 0 ) then
                EFTTestStability = .false.
                if ( input_EFTCAMB%EFTCAMB_feedback_level > 0 ) write(*,'(a,E11.4)') '   Physical instability: tensor ghost instability. Tensor kinetic term = ', eft_cache%EFTAT
            end if

            ! 5- No tensor gradient:
            if ( eft_cache%EFTDT < 0 ) then
                EFTTestStability = .false.
                if ( input_EFTCAMB%EFTCAMB_feedback_level > 0 ) write(*,'(a,E11.4)') '   Physical instability: tensor gradient instability. Tensor gradient term = ', eft_cache%EFTDT
            end if

        end if

    end function EFTTestStability

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that restores the values stored in the module to the default.
    !! Needed for successive calls.
    subroutine EFTStability_cleanup()

        implicit none

        PastA1  = 0._dl
        PastAT  = 0._dl

    end subroutine EFTStability_cleanup

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that fills the caches to check the stability of the theory.
    subroutine EFTStabilityComputation( a, input_model, params_cache, eft_cache )

        implicit none

        real(dl), intent(in)                         :: a                       !< the input scale factor.
        class(EFTCAMB_model)         , intent(in)    :: input_model             !< the EFTCAMB model for which the code is computing the RGR time.
        type(EFTCAMB_parameter_cache), intent(inout) :: params_cache            !< the EFTCAMB parameter cache that contains all the physical parameters.
        type(EFTCAMB_timestep_cache ), intent(inout) :: eft_cache               !< the EFTCAMB timestep cache that contains all the physical values.

        ! Definitions of variables:
        real(dl) :: grhonu, gpinu, grhormass_t, grhonudot, gpinudot
        integer  :: nu_i, ind, ind_max

        ! start filling the cache:
        eft_cache%a = a
        ! compute background densities of different species
        eft_cache%grhob_t = params_cache%grhob/a         ! 8\pi G_N \rho_b a^2: baryons background density
        eft_cache%grhoc_t = params_cache%grhoc/a         ! 8\pi G_N \rho_{cdm} a^2: cold dark matter background density
        eft_cache%grhor_t = params_cache%grhornomass/a/a ! 8\pi G_N \rho_{\nu} a^2: massless neutrinos background density
        eft_cache%grhog_t = params_cache%grhog/a/a       ! 8\pi G_N \rho_{\gamma} a^2: radiation background density
        ! Massive neutrinos terms:
        if ( params_cache%Num_Nu_Massive /= 0 ) then
            do nu_i = 1, params_cache%Nu_mass_eigenstates
                grhonu      = 0._dl
                gpinu       = 0._dl
                grhormass_t = params_cache%grhormass(nu_i)/a**2
                call params_cache%Nu_background(a*params_cache%nu_masses(nu_i), grhonu, gpinu )
                eft_cache%grhonu_tot = eft_cache%grhonu_tot +grhormass_t*grhonu
                eft_cache%gpinu_tot  = eft_cache%gpinu_tot  +grhormass_t*gpinu
            end do
        end if
        ! assemble total densities and pressure:
        eft_cache%grhom_t  = eft_cache%grhob_t +eft_cache%grhoc_t +eft_cache%grhor_t +eft_cache%grhog_t +eft_cache%grhonu_tot
        eft_cache%gpresm_t = (+eft_cache%grhor_t +eft_cache%grhog_t)/3._dl +eft_cache%gpinu_tot
        ! compute the other things:
        select type ( model_temp => input_model )
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
        ! compute massive neutrinos stuff:
        ! Massive neutrinos mod:
        if ( params_cache%Num_Nu_Massive /= 0 ) then
            do nu_i = 1, params_cache%Nu_mass_eigenstates
                grhonu      = 0._dl
                gpinu       = 0._dl
                grhonudot   = 0._dl
                gpinudot    = 0._dl
                grhormass_t = params_cache%grhormass(nu_i)/a**2
                call params_cache%Nu_background( a*params_cache%nu_masses(nu_i), grhonu, gpinu )
                eft_cache%grhonu_tot = eft_cache%grhonu_tot + grhormass_t*grhonu
                eft_cache%gpinu_tot  = eft_cache%gpinu_tot  + grhormass_t*gpinu
                eft_cache%grhonudot_tot = eft_cache%grhonudot_tot + grhormass_t*( params_cache%Nu_drho(a*params_cache%nu_masses(nu_i) ,eft_cache%adotoa, grhonu)&
                    & -4._dl*eft_cache%adotoa*grhonu)
                eft_cache%gpinudot_tot  = eft_cache%gpinudot_tot  + grhormass_t*( params_cache%Nu_pidot(a*params_cache%nu_masses(nu_i),eft_cache%adotoa, gpinu )&
                    & -4._dl*eft_cache%adotoa*gpinu)
            end do
        end if
        ! compute pressure dot:
        eft_cache%gpresdotm_t = -4._dl*eft_cache%adotoa*( eft_cache%grhog_t +eft_cache%grhor_t )/3._dl +eft_cache%gpinudot_tot
        ! compute remaining quantities related to H:
        call input_model%compute_H_derivs( a, params_cache , eft_cache )
        ! compute backgrond EFT functions if model is designer:
        select type ( model_temp => input_model )
            class is ( EFTCAMB_designer_model )
            call input_model%compute_background_EFT_functions( a, params_cache , eft_cache )
        end select
        ! compute all other background stuff:
        call input_model%compute_rhoQPQ( a, params_cache , eft_cache )
        ! compute second order EFT functions:
        call input_model%compute_secondorder_EFT_functions( a, params_cache , eft_cache )
        ! Compute pi field equations factors:
        call input_model%compute_pi_factors( a, params_cache , eft_cache )
        ! Compute coefficients for the tensor propagation equation:
        call input_model%compute_tensor_factors( a, params_cache , eft_cache )
        ! Compute kinetic and gradient terms:
        call input_model%compute_stability_factors( a, params_cache , eft_cache )

    end subroutine EFTStabilityComputation

    !----------------------------------------------------------------------------------------

end module EFTCAMB_stability

!----------------------------------------------------------------------------------------
