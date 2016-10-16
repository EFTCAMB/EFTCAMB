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
!! This operates on general EFT models but needs CP to be set.


!----------------------------------------------------------------------------------------
!> This module contains the stability detection algorithm of EFTCAMB.
!! This operates on general EFT models but needs CP to be set.

!> @author Bin Hu, Marco Raveri, Simone Peirone

module EFTCAMB_stability

    use precision
    use ModelParams
    use MassiveNu
    use EFTCAMB_abstract_model
    use EFTCAMB_abstract_model_full
    use EFTCAMB_abstract_model_designer
    use EFTDef
    use EFTCAMB_cache
    use EFTCAMB_main

    implicit none

    ! storage for some utility values that needs to be stored for the stability check.
    real(dl), save :: PastA1 = 0._dl
    real(dl), save :: PastAT = 0._dl

contains

    ! ---------------------------------------------------------------------------------------------
    !> Function that fills the caches to check the stability of the theory.
    function EFTStabilityCheck( EFTCAMB_in, a, eft_par_cache )

        implicit none

        type(EFTCAMB)                                :: EFTCAMB_in              !< the input EFTCAMB.
        real(dl), intent(in)                         :: a                       !< the input scale factor.
        type(EFTCAMB_parameter_cache), intent(inout) :: eft_par_cache           !< the EFTCAMB parameter cache that contains all the physical parameters.
        logical                                      :: EFTStabilityCheck       !< Logical value returned by the function. If the model is stable this is True, otherwise False.

        ! Definitions of variables:
        type(EFTCAMB_timestep_cache ) :: eft_cache
        logical  :: EFT_HaveNan_parameter, EFT_HaveNan_timestep
        real(dl) :: EFT_instability_rate, kmax, tempk, temp1, temp2, temp3, temp4, temp5
        integer  :: ind_max, ind

        ! Stability check initialization
        EFTStabilityCheck = .true.
        ! reset the time-step cache:
        call eft_cache%initialize()
        ! fill it:
        call EFTStabilityComputation( EFTCAMB_in, a, eft_par_cache, eft_cache )
        ! Compute k_max_CAMB. This is slightly overestimated to be safe...
        ! This is just to be safe if doing only background stuff. We still require stability of linear scales.
        kmax = 0.1_dl
        !    This is if we want scalar Cls.
        if ( CP%WantCls ) kmax = CP%Max_eta_k/CP%tau0
        !    This is if we want also (or only) transfer functions.
        if ( CP%WantTransfer ) then
            kmax = CP%Max_eta_k/CP%tau0*exp((int(log(CP%Transfer%kmax/(CP%Max_eta_k/CP%tau0))*(3*AccuracyBoost))+2)/(3*AccuracyBoost))
        end if

        ! check stability of the theory:

        ! 1) everything inside the parameter cache should not be a NaN:
        call eft_par_cache%is_nan( EFT_HaveNan_parameter )
        if ( EFT_HaveNan_parameter ) then
            EFTStabilityCheck = .false.
            if ( EFTCAMB_in%EFTCAMB_feedback_level > 0 ) write(*,*) 'EFTCAMB: model has Nan in the parameter cache'
            return
        end if
        ! 2) everything inside the time-step cache should not be a NaN:
        call eft_cache%is_nan( EFT_HaveNan_timestep )
        if ( EFT_HaveNan_timestep ) then
            EFTStabilityCheck = .false.
            if ( EFTCAMB_in%EFTCAMB_feedback_level > 0 ) write(*,*) 'EFTCAMB: model has Nan in the timestep cache'
            return
        end if
        ! 3) enforce mathematical stability:
        if ( EFTCAMB_in%EFT_mathematical_stability ) then

            ! 1- the A coefficient should not change sign in time and in k, i.e. it shall not be zero.
            !    This is the strongest stability constraint since violating it would violate the mathematical
            !    consistency of the pi field equation.
            !    The first condition is A1/=0. Implemented by detecting sign changes in A1.
            if ( eft_cache%EFTpiA1*PastA1 < 0._dl ) then
                EFTStabilityCheck = .false.
                if ( EFTCAMB_in%EFTCAMB_feedback_level > 0 ) write(*,*) 'EFTCAMB: mathematical instability, A is zero in time.'
            end if
            PastA1 = eft_cache%EFTpiA1
            !    The second one is the condition on k.
            if ( (eft_cache%EFTpiA1 > 0 .and. eft_cache%EFTpiA1 + kmax**2*eft_cache%EFTpiA2 < 0) .or. &
                &(eft_cache%EFTpiA1 < 0 .and. eft_cache%EFTpiA1 + kmax**2*eft_cache%EFTpiA2 > 0) ) then
                EFTStabilityCheck = .false.
                if ( EFTCAMB_in%EFTCAMB_feedback_level > 0 ) write(*,*) 'EFTCAMB: mathematical instability, A is zero in k.'
            end if

            ! 2- the AT coefficient should not change sign in time, i.e. it shall not be zero.
            !    This is the second strongest stability constraint since violating it would
            !    violate the mathematical consistency of the tensor perturbation equation.
            !    Implemented by detecting sign changes in AT.
            if ( eft_cache%EFTAT*PastAT < 0._dl ) then
                EFTStabilityCheck = .false.
                if ( EFTCAMB_in%EFTCAMB_feedback_level > 0 ) write(*,*) 'EFTCAMB: mathematical instability, AT is zero in time.'
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
                tempk = 0._dl + REAL(ind-1)*(kmax)/REAL(ind_max-1)
                ! vaule that discriminates between different cases:
                temp1 = (eft_cache%EFTpiB1 +eft_cache%EFTpiB2*tempk**2)
                temp2 = (eft_cache%EFTpiA1 +eft_cache%EFTpiA2*tempk**2)
                temp3 = temp1**2 -4._dl*temp2*(eft_cache%EFTpiC +eft_cache%EFTpiD1*tempk**2 + eft_cache%EFTpiD2*tempk**4)

                ! case 1:
                if ( temp3 > 0._dl .and. temp2 /= 0._dl ) then
                    temp4 = +0.5_dl*(-temp1 +sqrt(temp3))/temp2
                    temp5 = +0.5_dl*(-temp1 -sqrt(temp3))/temp2
                    if ( temp4>EFT_instability_rate .or. temp5>EFT_instability_rate ) then
                        EFTStabilityCheck = .false.
                        if ( EFTCAMB_in%EFTCAMB_feedback_level > 0 ) write(*,*) 'EFTCAMB: mathematical instability. Growing exponential at k', tempk, temp4, temp5
                        exit
                    end if

                ! case 2:
                else if ( temp2 /= 0._dl ) then
                    temp4 = -0.5_dl*temp1/temp2
                    if ( temp4>EFT_instability_rate ) then
                        EFTStabilityCheck = .false.
                        if ( EFTCAMB_in%EFTCAMB_feedback_level > 0 ) write(*,*) 'EFTCAMB: mathematical instability. Growing exponential at k', tempk, temp4
                        exit
                    end if
                end if

            end do

        end if
        ! 4) enforce model specific priors:
        if ( EFTCAMB_in%EFT_AdditionalPriors ) then
            EFTStabilityCheck = EFTCAMB_in%model%additional_model_stability( a, eft_par_cache, eft_cache )
            if ( .not. EFTStabilityCheck ) then
                if ( EFTCAMB_in%EFTCAMB_feedback_level > 0 ) write(*,*) 'EFTCAMB: model specific stability criteria are not met.'
            end if
        end if
        ! 5) enforce physical viability:
        if ( EFTCAMB_in%EFT_physical_stability ) then

            ! the present conditions extend up to Horndeski. Enforce that:
            if ( (eft_cache%EFTGamma6V /= 0._dl) .or.      &
                & ( (eft_cache%EFTGamma3V + eft_cache%EFTGamma4V) /= 0._dl) ) then
                write(*,*) 'EFTCAMB WARNING: stability for model beyond GLPV has not been worked out.'
                write(*,*) 'It will be added in a future release.'
                write(*,*) 'If you want to run this model disable EFT_physical_stability.'
                EFTStabilityCheck = .false.
                return
            end if

            ! 1- Positive gravitational constant:
            if ( 1._dl +eft_cache%EFTOmegaV <= 0 ) then
                EFTStabilityCheck = .false.
                if ( EFTCAMB_in%EFTCAMB_feedback_level > 0 ) write(*,*) 'EFTCAMB: negative gravitational constant', 1._dl +eft_cache%EFTOmegaV
            end if

            ! 2- Ghost condition:
            if ( eft_cache%EFT_kinetic < 0._dl ) then
                EFTStabilityCheck = .false.
                if ( EFTCAMB_in%EFTCAMB_feedback_level > 0 ) write(*,*) 'EFTCAMB: ghost instability. Kinetic term: ', eft_cache%EFT_kinetic
            end if

            ! 3- Gradient instability:
            if ( eft_cache%EFT_gradient < 0._dl ) then
                EFTStabilityCheck = .false.
                if ( EFTCAMB_in%EFTCAMB_feedback_level > 0 ) write(*,*) 'EFTCAMB: gradient instability. Gradient term: ', eft_cache%EFT_gradient
            end if

            ! 4- No tensor ghosts:
            if ( eft_cache%EFTAT < 0 ) then
                EFTStabilityCheck = .false.
                if ( EFTCAMB_in%EFTCAMB_feedback_level > 0 ) write(*,*) 'EFTCAMB: tensor ghost instability'
            end if

            ! 5- No tensor gradient:
            if ( eft_cache%EFTDT < 0 ) then
                EFTStabilityCheck = .false.
                if ( EFTCAMB_in%EFTCAMB_feedback_level > 0 ) write(*,*) 'EFTCAMB: tensor gradient instability'
            end if

        end if

    end function

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
    subroutine EFTStabilityComputation( EFTCAMB_in, a, eft_par_cache, eft_cache )

        implicit none

        type(EFTCAMB)                                :: EFTCAMB_in              !< the EFTCAMB model.
        real(dl), intent(in)                         :: a                       !< the input scale factor.
        type(EFTCAMB_parameter_cache), intent(inout) :: eft_par_cache           !< the EFTCAMB parameter cache that contains all the physical parameters.
        type(EFTCAMB_timestep_cache ), intent(inout) :: eft_cache               !< the EFTCAMB timestep cache that contains all the physical values.

        ! Definitions of variables:
        real(dl) :: EFT_grhonu, EFT_gpinu, grhormass_t, EFT_grhonudot, EFT_gpinudot, kmax
        real(dl) :: temp1, temp2, temp3, temp4, temp5, tempk, PastA1, PastAT
        integer  :: nu_i, ind, ind_max

        ! start filling the cache:
        eft_cache%a = a
        ! compute background densities of different species
        eft_cache%grhob_t = eft_par_cache%grhob/a         ! 8\pi G_N \rho_b a^2: baryons background density
        eft_cache%grhoc_t = eft_par_cache%grhoc/a         ! 8\pi G_N \rho_{cdm} a^2: cold dark matter background density
        eft_cache%grhor_t = eft_par_cache%grhornomass/a/a ! 8\pi G_N \rho_{\nu} a^2: massless neutrinos background density
        eft_cache%grhog_t = eft_par_cache%grhog/a/a       ! 8\pi G_N \rho_{\gamma} a^2: radiation background density
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
        select type ( model_temp => EFTCAMB_in%model )
            ! compute the background and the background EFT functions.
            class is ( EFTCAMB_full_model )
            ! background for full models. Here the expansion history is computed from the
            ! EFT functions. Hence compute them first and then compute the expansion history.
            call EFTCAMB_in%model%compute_background_EFT_functions( a, eft_par_cache , eft_cache )
            call EFTCAMB_in%model%compute_adotoa( a, eft_par_cache , eft_cache )
            class is ( EFTCAMB_designer_model )
            ! background for designer models. Here the expansion history is parametrized
            ! and does not depend on the EFT functions. Hence compute first the expansion history
            ! and then the EFT functions.
            call EFTCAMB_in%model%compute_adotoa( a, eft_par_cache , eft_cache )
        end select
        ! compute massive neutrinos stuff:
        ! Massive neutrinos mod:
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
                eft_cache%grhonudot_tot = eft_cache%grhonudot_tot + grhormass_t*(Nu_drho(a*nu_masses(nu_i) ,eft_cache%adotoa, EFT_grhonu)&
                    & -4._dl*eft_cache%adotoa*EFT_grhonu)
                eft_cache%gpinudot_tot  = eft_cache%gpinudot_tot  + grhormass_t*(Nu_pidot(a*nu_masses(nu_i),eft_cache%adotoa, EFT_gpinu )&
                    & -4._dl*eft_cache%adotoa*EFT_gpinu)
            end do
        end if
        ! compute pressure dot:
        eft_cache%gpresdotm_t = -4._dl*eft_cache%adotoa*( eft_cache%grhog_t +eft_cache%grhor_t )/3._dl +eft_cache%gpinudot_tot
        ! compute remaining quantities related to H:
        call EFTCAMB_in%model%compute_H_derivs( a, eft_par_cache , eft_cache )
        ! compute backgrond EFT functions if model is designer:
        select type ( model_temp => EFTCAMB_in%model )
            class is ( EFTCAMB_designer_model )
            call EFTCAMB_in%model%compute_background_EFT_functions( a, eft_par_cache , eft_cache )
        end select
        ! compute all other background stuff:
        call EFTCAMB_in%model%compute_rhoQPQ( a, eft_par_cache , eft_cache )
        ! compute second order EFT functions:
        call EFTCAMB_in%model%compute_secondorder_EFT_functions( a, eft_par_cache , eft_cache )
        ! Compute pi field equations factors:
        call EFTCAMB_in%model%compute_pi_factors( a, eft_par_cache , eft_cache )
        ! Compute coefficients for the tensor propagation equation:
        call EFTCAMB_in%model%compute_tensor_factors( a, eft_par_cache , eft_cache )
        ! Compute kinetic and gradient terms:
        call EFTCAMB_in%model%compute_stability_factors( a, eft_par_cache , eft_cache )

    end subroutine EFTStabilityComputation

    !----------------------------------------------------------------------------------------

end module EFTCAMB_stability

!----------------------------------------------------------------------------------------
