! ----------------------------------------------------------------------------------------
!
! This file is part of EFTCAMB.
!
! Copyright (C) 2013-2023 by the EFTCAMB authors
!
! The EFTCAMB code is free software;
! You can use it, redistribute it, and/or modify it under the terms
! of the GNU General Public License as published by the Free Software Foundation;
! either version 3 of the License, or (at your option) any later version.
! The full text of the license can be found in the file eftcamb/LICENSE at
! the top level of the EFTCAMB distribution.
!
! ----------------------------------------------------------------------------------------

!> @file 09_EFTCAMB_RGR.f90
!! This file contains the RGR algorithm. This operates on general EFT models.


! ----------------------------------------------------------------------------------------
!> This module contains the RGR algorithm. This operates on general EFT models.

!> @author Bin Hu, Marco Raveri

    module EFTCAMB_ReturnToGR

    use precision
    use config, only: GlobalError, error_unsupported_params
    use EFT_def
    use EFTCAMB_cache
    use EFTCAMB_abstract_model
    use EFTCAMB_abstract_model_full
    use EFTCAMB_abstract_model_designer
    use MassiveNu

    implicit none

    private

    public EFTCAMBReturnToGR_feedback, EFTCAMBReturnToGR

    !< number of points to sample logaritmically the time in the return to GR module.
    integer, parameter :: EFT_RGR_num_points = 1000
    !< scale (Mpc^-1) at which the quasi-static modification of gravity is tested: scale-dependent models
    !! (e.g. designer f(R)) modify gravity at k >> aH by O(1) while their EFT functions are tiny. The linear
    !! scales k <= 1/Mpc are the target: a larger probe k turns near-GR f(R) on much earlier, where its pi equation
    !! is very stiff (4_Des_F_R_LCDM_2: k = 10/Mpc, threshold 1e-3 gives turn-on 0.016 and 2.8x the runtime).
    real(dl), parameter :: EFT_RGR_kmax = 1._dl
    !< the EFT perturbations are also switched on when the quasi-static mu or Sigma at EFT_RGR_kmax deviate from
    !! 1 by more than this. Tuned on 4_Des_F_R_LCDM_2 (turn-on 0.039) against an early turn-on (0.007): P(k <= 1/Mpc)
    !! 4e-4, P(k <= 10/Mpc) 3.3e-2, sigma8 5e-6, lensed C_l 3e-5, C_l^phiphi 1.1e-4, runtime 1.1x the pre-C3 code.
    real(dl), parameter :: EFT_RGR_qs_threshold = 5.e-4_dl

    contains

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that prints feedback for the RGR module. Prints something only if
    !! EFTCAMB_feedback_level is greater than 1.
    subroutine EFTCAMBReturnToGR_feedback(feedback_level, input_model, params_cache, initial_time, RGR_time, &
        RGR_threshold)

    implicit none

    !< feedback level for the RGR results reporter. 0=no feedback; 1=some feedback; 2=a lot of feedback.
    integer, intent(in) :: feedback_level
    !< the EFTCAMB model for which the code is computing the RGR time.
    class(EFTCAMB_model), intent(in) :: input_model
    !< a EFTCAMB parameter cache containing cosmological parameters.
    type(TEFTCAMB_parameter_cache), intent(in) :: params_cache
    !< initial scale factor at which the code starts to look for the RGR of the theory.
    real(dl), intent(in) :: initial_time
    !< input RGR time. For better performances the code does not compute it in this function.
    real(dl), intent(in) :: RGR_time
    !< input RGR threshold. The threshold for an EFT functions to have reached its GR value.
    real(dl), intent(in) :: RGR_threshold

    real(dl) :: eft_functions(29)
    integer  :: i

    if (feedback_level > 1) then
        write(*, '(a)') "***************************************************************"
        write(*, '(a,E13.4)') ' EFTCAMB Return to GR time:     ', RGR_time
        write(*, '(a,E13.4)') ' EFTCAMB Return to GR threshold:', RGR_threshold
        ! (the RGR scan samples 10**log10(a): equal up to roundoff)
        if (abs(RGR_time - initial_time) <= 1.e-10_dl*initial_time) then
            write(*, '(a)') ' Warning, RGR time and turn on time coincide.'
            write(*, '(a)') ' The model might be far from GR at initial time.'
        end if
        if (RGR_time > 1._dl) then
            write(*, '(a)') ' Warning, the model is GR within the threshold at all times:'
            write(*, '(a)') ' EFT perturbations are never switched on.'
        end if
    end if

    if (feedback_level > 2) then

        call EFTCAMBReturnToGR_functions(RGR_time, input_model, params_cache, eft_functions)

        write(*, '(a)')
        write(*, '(a)') ' EFT functions at RGR time: '

        if (eft_functions(1) - RGR_threshold > 0._dl) write(*, '(a18,E13.4)') '   OmegaV       = ', eft_functions(1)
        if (eft_functions(2) - RGR_threshold > 0._dl) write(*, '(a18,E13.4)') '   a OmegaP     = ', eft_functions(2)
        if (eft_functions(3) > EFT_RGR_qs_threshold) write(*, '(a18,E13.4)') '   |QS mu-1|    = ', eft_functions(3)
        if (eft_functions(4) - RGR_threshold > 0._dl) write(*, '(a18,E13.4)') '   OmegaPPP     = ', eft_functions(4)
        if (eft_functions(5) - RGR_threshold > 0._dl) write(*, '(a18,E13.4)') '   c/3H^2       = ', eft_functions(5)
        if (eft_functions(6) - RGR_threshold > 0._dl) write(*, '(a18,E13.4)') '   dLambda/3H^2 = ', eft_functions(6)
        if (eft_functions(7) - RGR_threshold > 0._dl) write(*, '(a18,E13.4)') '   cdot/3H^3    = ', eft_functions(7)
        if (eft_functions(8) - RGR_threshold > 0._dl) write(*, '(a18,E13.4)') '   Lambdadot/3H3= ', eft_functions(8)
        if (eft_functions(9) - RGR_threshold > 0._dl) write(*, '(a18,E13.4)') '   EFTGamma1V   = ', eft_functions(9)
        if (eft_functions(10) - RGR_threshold > 0._dl) write(*, '(a18,E13.4)') '   EFTGamma1P   = ', eft_functions(10)
        if (eft_functions(11) - RGR_threshold > 0._dl) write(*, '(a18,E13.4)') '   EFTGamma2V   = ', eft_functions(11)
        if (eft_functions(12) - RGR_threshold > 0._dl) write(*, '(a18,E13.4)') '   EFTGamma2P   = ', eft_functions(12)
        if (eft_functions(13) - RGR_threshold > 0._dl) write(*, '(a18,E13.4)') '   EFTGamma3V   = ', eft_functions(13)
        if (eft_functions(14) - RGR_threshold > 0._dl) write(*, '(a18,E13.4)') '   EFTGamma3P   = ', eft_functions(14)
        if (eft_functions(15) - RGR_threshold > 0._dl) write(*, '(a18,E13.4)') '   EFTGamma4V   = ', eft_functions(15)
        if (eft_functions(16) - RGR_threshold > 0._dl) write(*, '(a18,E13.4)') '   EFTGamma4P   = ', eft_functions(16)
        if (eft_functions(18) - RGR_threshold > 0._dl) write(*, '(a18,E13.4)') '   EFTGamma5V   = ', eft_functions(18)
        if (eft_functions(17) - RGR_threshold > 0._dl) write(*, '(a18,E13.4)') '   EFTGamma4PP  = ', eft_functions(17)
        if (eft_functions(19) - RGR_threshold > 0._dl) write(*, '(a18,E13.4)') '   EFTGamma5P   = ', eft_functions(19)
        if (eft_functions(20) - RGR_threshold > 0._dl) write(*, '(a18,E13.4)') '   EFTGamma6V   = ', eft_functions(20)
        if (eft_functions(21) - RGR_threshold > 0._dl) write(*, '(a18,E13.4)') '   EFTGamma6P   = ', eft_functions(21)
        if (eft_functions(22) - RGR_threshold > 0._dl) write(*, '(a18,E13.4)') '   OmegaPPPP   = ', eft_functions(22)
        if (eft_functions(23) - RGR_threshold > 0._dl) write(*, '(a18,E13.4)') '   EFTLambdadotdot   = ', &
            eft_functions(23)
        if (eft_functions(24) - RGR_threshold > 0._dl) write(*, '(a18,E13.4)') '   EFTGamma1PP   = ', eft_functions(24)
        if (eft_functions(25) - RGR_threshold > 0._dl) write(*, '(a18,E13.4)') '   EFTGamma2PP   = ', eft_functions(25)
        if (eft_functions(26) - RGR_threshold > 0._dl) write(*, '(a18,E13.4)') '   EFTGamma2PPP   = ', eft_functions(26)
        if (eft_functions(27) - RGR_threshold > 0._dl) write(*, '(a18,E13.4)') '   EFTGamma3PP   = ', eft_functions(27)
        if (eft_functions(28) - RGR_threshold > 0._dl) write(*, '(a18,E13.4)') '   EFTGamma3PPP   = ', eft_functions(28)
        if (eft_functions(29) - RGR_threshold > 0._dl) write(*, '(a18,E13.4)') '   EFTGamma3PPPP   = ', &
            eft_functions(29)

    end if

    end subroutine

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes the return to GR of a theory.
    subroutine EFTCAMBReturnToGR(input_model, params_cache, initial_time, RGR_time, RGR_threshold)

    implicit none

    !< the EFTCAMB model for which the code is computing the RGR time.
    class(EFTCAMB_model), intent(in) :: input_model
    !< a EFTCAMB parameter cache containing cosmological parameters.
    type(TEFTCAMB_parameter_cache), intent(in) :: params_cache
    !< initial scale factor at which the code starts to look for the RGR of the theory.
    real(dl), intent(in)  :: initial_time
    real(dl), intent(out) :: RGR_time        !< output value of the RGR time.
    !< input RGR threshold. The threshold for an EFT functions to have reached its GR value.
    real(dl), intent(in)  :: RGR_threshold

    real(dl) :: eft_functions(29)
    real(dl) :: a_initial, a_final, log_a_initial, log_a_final, log_a_used, a_used, qs_deviation
    integer  :: i

    ! initialize:
    RGR_time = initial_time

    ! the turn-on time of a model already found to be GR at all times (parameters of previous results):
    if (initial_time == EFT_RGR_never_on) return

    ! validate the turn-on time (it can be set directly from python): for a non-positive or NaN time the
    ! logarithmic scan is NaN and no EFT perturbations would ever be switched on; a time > 1 does the same.
    if (.not. (initial_time > 0._dl .and. initial_time <= 1._dl)) then
        call GlobalError('EFTCAMB: EFTCAMB_turn_on_time must be in (0,1]', error_unsupported_params)
        return
    end if

    ! initial and final scale factors:
    a_initial = initial_time
    a_final = 1._dl
    ! convert to logarithm of the scale factor:
    log_a_initial = log10(a_initial)
    log_a_final = log10(a_final)

    ! loop over the logarithm of the scale factor:
    do i = 1, EFT_RGR_num_points

        ! initialize:
        log_a_used = log_a_initial + real(i - 1)/real(EFT_RGR_num_points - 1)*(log_a_final - log_a_initial)
        a_used = 10._dl**log_a_used
        eft_functions = 0._dl
        ! compute RGR functions:
        call EFTCAMBReturnToGR_functions(a_used, input_model, params_cache, eft_functions)

        ! check if above threshold (eft_functions(3), the observable-based quasi-static criterion, has its own):
        qs_deviation = eft_functions(3)
        eft_functions(3) = 0._dl
        if (any(eft_functions - RGR_threshold > 0._dl) .or. qs_deviation > EFT_RGR_qs_threshold) then
            RGR_time = a_used
            return
        end if

    end do

    RGR_time = EFT_RGR_never_on

    end subroutine EFTCAMBReturnToGR

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes the EFT functions as we need them for the RGR function.
    subroutine EFTCAMBReturnToGR_functions(a, input_model, params_cache, eft_functions)

    implicit none

    real(dl), intent(in) :: a                  !< scale factor at which EFT funcitons are computed.
    !< the EFTCAMB model for which the code is computing the RGR time.
    class(EFTCAMB_model), intent(in) :: input_model
    !< a EFTCAMB parameter cache containing cosmological parameters.
    type(TEFTCAMB_parameter_cache)   :: params_cache
    !< vector containing the values of the EFT functions minus their GR value.
    real(dl), intent(out) :: eft_functions(29)

    type(TEFTCAMB_timestep_cache) :: eft_cache

    real(dl) :: a2, adotoa, grhonu, gpinu, grhonudot, gpinudot, grhormass_t, gpinudotdot, qs_mu, qs_sigma
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
    eft_cache%grhonu_tot = 0._dl
    eft_cache%gpinu_tot = 0._dl
    if (params_cache%Num_Nu_Massive /= 0) then
        do nu_i = 1, params_cache%Nu_mass_eigenstates
            grhonu = 0._dl
            gpinu = 0._dl
            grhormass_t = params_cache%grhormass(nu_i)/a**2
            call ThermalNuBackground%rho_P(a*params_cache%nu_masses(nu_i), grhonu, gpinu)
            eft_cache%grhonu_tot = eft_cache%grhonu_tot + grhormass_t*grhonu
            eft_cache%gpinu_tot = eft_cache%gpinu_tot + grhormass_t*gpinu
        end do
    end if
    ! assemble total densities and pressure:
    eft_cache%grhom_t = &
        eft_cache%grhob_t + eft_cache%grhoc_t + eft_cache%grhor_t + eft_cache%grhog_t + eft_cache%grhonu_tot
    eft_cache%gpresm_t = (+eft_cache%grhor_t + eft_cache%grhog_t)/3._dl + eft_cache%gpinu_tot
    ! compute the other things:
    select type (model => input_model)
        ! compute the background and the background EFT functions.
    class is (EFTCAMB_full_model)
        ! background for full models. Here the expansion history is computed from the
        ! EFT functions. Hence compute them first and then compute the expansion history.
        call input_model%compute_background_EFT_functions(a, params_cache, eft_cache)
        call input_model%compute_adotoa(a, params_cache, eft_cache)
    class is (EFTCAMB_designer_model)
        ! background for designer models. Here the expansion history is parametrized
        ! and does not depend on the EFT functions. Hence compute first the expansion history
        ! and then the EFT functions.
        call input_model%compute_adotoa(a, params_cache, eft_cache)
    end select
    ! store adotoa:
    adotoa = eft_cache%adotoa
    ! compute massive neutrinos stuff:
    ! Massive neutrinos mod:
    eft_cache%grhonu_tot = 0._dl
    eft_cache%gpinu_tot = 0._dl
    eft_cache%grhonudot_tot = 0._dl
    eft_cache%gpinudot_tot = 0._dl
    if (params_cache%Num_Nu_Massive /= 0) then
        do nu_i = 1, params_cache%Nu_mass_eigenstates
            grhonu = 0._dl
            gpinu = 0._dl
            grhonudot = 0._dl
            gpinudot = 0._dl
            gpinudotdot = 0._dl
            grhormass_t = params_cache%grhormass(nu_i)/a**2
            call ThermalNuBackground%rho_P(a*params_cache%nu_masses(nu_i), grhonu, gpinu)
            eft_cache%grhonu_tot = eft_cache%grhonu_tot + grhormass_t*grhonu
            eft_cache%gpinu_tot = eft_cache%gpinu_tot + grhormass_t*gpinu
            eft_cache%grhonudot_tot = eft_cache%grhonudot_tot + &
                grhormass_t*(ThermalNuBackground%drho(a*params_cache%nu_masses(nu_i), adotoa) - 4._dl*adotoa*grhonu)
            eft_cache%gpinudot_tot = eft_cache%gpinudot_tot + &
                grhormass_t*(ThermalNuBackground%pidot(a*params_cache%nu_masses(nu_i), adotoa, gpinu) - &
                4._dl*adotoa*gpinu)
            eft_cache%gpinudotdot_tot = 0._dl
        end do
    end if
    ! compute pressure dot:
    eft_cache%gpresdotm_t = -4._dl*adotoa*(eft_cache%grhog_t + eft_cache%grhor_t)/3._dl + eft_cache%gpinudot_tot
    ! compute remaining quantities related to H:
    call input_model%compute_H_derivs(a, params_cache, eft_cache)
    ! compute massive neutrinos stuff:
    ! Massive neutrinos mod:
    eft_cache%grhonu_tot = 0._dl
    eft_cache%gpinu_tot = 0._dl
    eft_cache%grhonudot_tot = 0._dl
    eft_cache%gpinudot_tot = 0._dl
    eft_cache%gpinudotdot_tot = 0._dl
    if (params_cache%Num_Nu_Massive /= 0) then
        do nu_i = 1, params_cache%Nu_mass_eigenstates
            grhonu = 0._dl
            gpinu = 0._dl
            grhonudot = 0._dl
            gpinudot = 0._dl
            gpinudotdot = 0._dl
            grhormass_t = params_cache%grhormass(nu_i)/a**2
            call ThermalNuBackground%rho_P(a*params_cache%nu_masses(nu_i), grhonu, gpinu)
            eft_cache%grhonu_tot = eft_cache%grhonu_tot + grhormass_t*grhonu
            eft_cache%gpinu_tot = eft_cache%gpinu_tot + grhormass_t*gpinu
            eft_cache%grhonudot_tot = eft_cache%grhonudot_tot + &
                grhormass_t*(ThermalNuBackground%drho(a*params_cache%nu_masses(nu_i), adotoa) - 4._dl*adotoa*grhonu)
            gpinudot = ThermalNuBackground%pidot(a*params_cache%nu_masses(nu_i), eft_cache%adotoa, gpinu)
            eft_cache%gpinudot_tot = eft_cache%gpinudot_tot + grhormass_t*(gpinudot &
                - 4._dl*adotoa*gpinu)
            eft_cache%gpinudotdot_tot = eft_cache%gpinudotdot_tot - &
                4._dl*adotoa*grhormass_t*(gpinudot - 4._dl*adotoa*gpinu) + &
                grhormass_t*(ThermalNuBackground%pidotdot(a*params_cache%nu_masses(nu_i), adotoa, eft_cache%Hdot, &
                gpinu, gpinudot) - 4._dl*eft_cache%Hdot*gpinu - 4._dl*adotoa*gpinudot)
        end do
    end if
    ! compute pressure ddot:
    eft_cache%gpresdotdotm_t = -(4._dl/3._dl)*eft_cache%Hdot*(eft_cache%grhog_t + eft_cache%grhor_t) + &
        (16._dl/3._dl)*(eft_cache%grhog_t + eft_cache%grhor_t)*adotoa**2 + eft_cache%gpinudotdot_tot
    ! compute remaining quantities related to H:
    call input_model%compute_H_derivs(a, params_cache, eft_cache)
    ! compute backgrond EFT functions if model is designer:
    select type (model => input_model)
    class is (EFTCAMB_designer_model)
        call input_model%compute_background_EFT_functions(a, params_cache, eft_cache)
    end select
    ! compute all other background stuff:
    call input_model%compute_rhoQPQ(a, params_cache, eft_cache)
    ! compute second order EFT functions:
    call input_model%compute_secondorder_EFT_functions(a, params_cache, eft_cache)

    ! quasi-static modification of the Poisson and lensing equations at k = EFT_RGR_kmax:
    call input_model%compute_QS_Mu_Sigma(a, EFT_RGR_kmax, params_cache, eft_cache)

    ! get the EFT functions, all dimensionless: the dimensionful c, Lambda (H^2 units, a^2-weighted in the cache)
    ! are normalised by 3 adotoa^2 = 3 a^2 H^2 and their conformal time derivatives by 3 adotoa^3.
    eft_functions(1) = abs(eft_cache%EFTOmegaV)
    eft_functions(2) = abs(a*eft_cache%EFTOmegaP)
    ! observable-based criterion (scale-dependent models): deviation of the quasi-static mu, Sigma from GR at
    ! k = EFT_RGR_kmax (NaN values, e.g. at poles of the quasi-static ratios, are ignored).
    qs_mu = abs(eft_cache%QS_mu_gen - 1._dl)
    qs_sigma = abs(eft_cache%QS_sigma_gen - 1._dl)
    if (qs_mu /= qs_mu) qs_mu = 0._dl
    if (qs_sigma /= qs_sigma) qs_sigma = 0._dl
    eft_functions(3) = max(qs_mu, qs_sigma)
    eft_functions(4) = 0._dl
    eft_functions(5) = abs(eft_cache%EFTc)/(3._dl*adotoa**2)
    eft_functions(6) = abs(eft_cache%EFTLambda + params_cache%grhov*a2)/(3._dl*adotoa**2)
    eft_functions(7) = abs(eft_cache%EFTcdot)/(3._dl*adotoa**3)
    eft_functions(8) = abs(eft_cache%EFTLambdadot)/(3._dl*adotoa**3)
    eft_functions(9) = abs(eft_cache%EFTGamma1V)
    eft_functions(10) = abs(eft_cache%EFTGamma1P)
    eft_functions(11) = abs(eft_cache%EFTGamma2V)
    eft_functions(12) = abs(eft_cache%EFTGamma2P)
    eft_functions(13) = abs(eft_cache%EFTGamma3V)
    eft_functions(14) = abs(eft_cache%EFTGamma3P)
    eft_functions(15) = abs(eft_cache%EFTGamma4V)
    eft_functions(16) = abs(eft_cache%EFTGamma4P)
    eft_functions(17) = 0._dl
    eft_functions(18) = abs(eft_cache%EFTGamma5V)
    eft_functions(19) = abs(eft_cache%EFTGamma5P)
    eft_functions(20) = abs(eft_cache%EFTGamma6V)
    eft_functions(21) = abs(eft_cache%EFTGamma6P)
    eft_functions(22) = 0._dl
    eft_functions(23) = 0._dl
    eft_functions(24) = 0._dl
    eft_functions(25) = 0._dl
    eft_functions(26) = 0._dl
    eft_functions(27) = 0._dl
    eft_functions(28) = 0._dl
    eft_functions(29) = 0._dl

    ! subtract the GR threshold:
    eft_functions = eft_functions

    end subroutine EFTCAMBReturnToGR_functions

    ! ----------------------------------------------------------------------------------------

    end module EFTCAMB_ReturnToGR

! ----------------------------------------------------------------------------------------
