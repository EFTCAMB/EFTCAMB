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

!> @file 09_EFTCAMB_stability.f90
!! This file contains the stability detection algorithm of EFTCAMB.


! ----------------------------------------------------------------------------------------
!> This module contains the stability detection algorithm of EFTCAMB.

!> @author Bin Hu, Marco Raveri, Simone Peirone

    module EFTCAMB_stability

    use precision
    use config, only: GlobalError
    use EFT_def
    use EFTCAMB_cache
    use EFTCAMB_abstract_model
    use EFTCAMB_abstract_model_full
    use EFTCAMB_abstract_model_designer
    use EFTCAMB_main
    use EFTCAMB_ReturnToGR
    use MassiveNu

    implicit none

    private

    public EFTCAMB_Stability_Check, EFTTestStability, EFTStability_cleanup, EFTStabilityComputation
    public EFTStabilityMassEigenvalues
    public EFTCAMB_compute_effective_w0wa

    ! storage for some utility values that needs to be stored for the stability check.
    real(dl), save :: PastA1 = 0._dl
    real(dl), save :: PastAT = 0._dl
    !< sign of the pi kinetic coefficient A (k -> 0) at the previous sampled scale factor after the turn-on.
    real(dl), save :: PastAsign = 0._dl
    !< the two previous sampled scale factors after the turn-on and the margin of A there (its value closest to zero
    !! in k, times its sign), to bracket a local minimum of A where it may touch zero between two samples.
    real(dl), save :: PastAa(2) = 0._dl, PastAm(2) = 0._dl
    !< c_s^2 above which the stability check warns that the pi equation is stiff: it oscillates at the rate c_s k and
    !! the run time grows about as c_s (test inis: 1.5x the GR time at c_s^2 = 5e2, 6x at 4e4, 15x at 2e5).
    real(dl), parameter :: EFT_stiff_cs2 = 1.e4_dl
    !< the first condition that failed in the last call of EFTTestStability (for the refusal message).
    character(len=120), save :: EFT_instability_reason = ''

    contains

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that tests the stability of a theory in a time range.
    !! To ensure best time coverage scans with three different strategies.
    subroutine EFTCAMB_Stability_Check(success, input_EFTCAMB, params_cache, astart, aend, k_max, reason)

    implicit none

    !< Output of the subroutine. Tells whether the model is found stable or not.
    logical, intent(out) :: success
    type(TEFTCAMB), intent(in) :: input_EFTCAMB   !< the EFTCAMB object for which the code is computing stability.
    !< the EFTCAMB parameter cache that contains all the physical parameters.
    type(TEFTCAMB_parameter_cache), intent(inout) :: params_cache
    real(dl), intent(in)    :: astart          !< Initial scale factor.
    real(dl), intent(in)    :: aend            !< Final scale factor.
    real(dl), intent(inout) :: k_max           !< the input maximum k mode at which stability is computed.
    !< on failure, the first condition that failed and the scale factor (for the refusal message).
    character(len=*), intent(out), optional :: reason

    ! parameters of the stability sampler:
    integer, parameter :: indMax = 1000    ! Number of points sampled.
    real(dl), parameter :: LogSamplingScale = -10._dl ! Where to start with the log sampling

    real(dl) :: Atest, y, cs2_max, a_cs2_max, a_fail
    integer  :: ind
    type(TEFTCAMB_timestep_cache) :: eft_cache

    ! 0) initial feedback:
    if (input_EFTCAMB%EFTCAMB_feedback_level > 1) then
        write(*, '(a)') '***************************************************************'
        write(*, '(a)') ' EFTCAMB: checking stability of the theory'
    end if
    if (input_EFTCAMB%EFTCAMB_feedback_level > 2) then
        write(*, '(a)')
    end if

    ! the mathematical mass condition is deprecated (EFTTestStability step 4): it rejects healthy models.
    if (input_EFTCAMB%EFT_mass_math_stability .and. input_EFTCAMB%EFTCAMB_feedback_level > 0) then
        write(*, '(a)') ' EFTCAMB WARNING: EFT_mass_math_stability is deprecated. It rejects any growing root of the'
        write(*, '(a)') '   k-sampled pi equation at frozen coefficients (threshold 0, from EFTCAMB_stability_time on)'
        write(*, '(a)') '   and refuses healthy models (11 of the 67 test models). Use EFT_mass_stability.'
    end if

    ! the mass stability check rejects models since the mass eigenvalues were corrected (step 7 of EFTTestStability);
    ! the former warning:
    ! ! the mass stability check is experimental and is ignored:
    ! if (input_EFTCAMB%EFT_mass_stability .and. input_EFTCAMB%EFTCAMB_feedback_level > 0) then
    !     write(*, '(a)') ' EFTCAMB WARNING: EFT_mass_stability is experimental and is ignored:'
    !     write(*, '(a)') '   the mass eigenvalues are computed but models are not rejected on them.'
    ! end if

    ! debug open cache files:
    if (DebugEFTCAMB) then
        call eft_cache%open_cache_files(input_EFTCAMB%outroot // 'stability_' )
    end if

    ! 1) stability code:
    success = .true.
    a_fail = 0._dl
    if (present(reason)) reason = ''
    ! the largest scalar sound speed after the perturbation turn-on (a stiff pi equation: reported below)
    cs2_max = 0._dl
    a_cs2_max = 0._dl

    ! - linear sampling:
    if (success) then
        call EFTStability_cleanup()
        do ind = 1, indMax
            Atest = astart + real(ind - 1)*(aend - astart)/real(indMax - 1)
            success = EFTTestStability(Atest, k_max, input_EFTCAMB, params_cache, eft_cache)
            if (success) call track_sound_speed(Atest)
            if (.not. success) then
                a_fail = Atest
                if (input_EFTCAMB%EFTCAMB_feedback_level > 2) then
                    write(*, *)
                    write(*, '(a,E14.4)') '   Instability detected at a =', Atest
                end if
                exit
            end if
        end do
    end if

    ! - log sampling close to astart:
    if (success) then
        call EFTStability_cleanup()
        do ind = 1, indMax
            y = LogSamplingScale + real(ind - 1)*(0._dl - LogSamplingScale)/real(indMax - 1)
            Atest = astart + (aend - astart)*10._dl**y
            success = EFTTestStability(Atest, k_max, input_EFTCAMB, params_cache, eft_cache)
            if (success) call track_sound_speed(Atest)
            if (.not. success) then
                a_fail = Atest
                if (input_EFTCAMB%EFTCAMB_feedback_level > 2) then
                    write(*, *)
                    write(*, '(a,E14.4)') '   Instability detected at a =', Atest
                end if
                exit
            end if
        end do
    end if

    ! - log sampling close to aend:
    if (success) then
        call EFTStability_cleanup()
        do ind = 1, indMax
            y = LogSamplingScale + real(ind - 1)*(0._dl - LogSamplingScale)/real(indMax - 1)
            Atest = aend + (astart - aend)*10._dl**y
            success = EFTTestStability(Atest, k_max, input_EFTCAMB, params_cache, eft_cache)
            if (success) call track_sound_speed(Atest)
            if (.not. success) then
                a_fail = Atest
                if (input_EFTCAMB%EFTCAMB_feedback_level > 2) then
                    write(*, *)
                    write(*, '(a,E14.4)') '   Instability detected at a =', Atest
                end if
                exit
            end if
        end do
    end if

    ! debug close cache files:
    if (DebugEFTCAMB) then
        call eft_cache%close_cache_files()
    end if

    ! the failed condition and where, for the refusal message:
    if (present(reason) .and. .not. success) then
        write(reason, '(a,a,es10.3)') trim(EFT_instability_reason), ' at a =', a_fail
    end if

    ! 2) final feedback:
    if (input_EFTCAMB%EFTCAMB_feedback_level > 1) then
        if (success) then
            write(*, '(a)') ' EFTCAMB: theory stable'
        else
            write(*, '(a)') ' EFTCAMB: theory unstable'
        end if
        if (success .and. cs2_max > 0._dl) write(*, '(a,es10.3,a,es10.3)') &
            '   largest c_s^2 after the turn-on:', cs2_max, ' at a =', a_cs2_max
    end if
    ! a large sound speed makes the pi equation stiff and the run slow (a diagnostic only: the model is kept)
    if (success .and. cs2_max > EFT_stiff_cs2 .and. input_EFTCAMB%EFTCAMB_feedback_level > 0) then
        write(*, '(a,es10.3,a,es10.3)') ' EFTCAMB WARNING: large scalar sound speed after the turn-on: c_s^2 =', &
            cs2_max, ' at a =', a_cs2_max
        write(*, '(a)') '   The pi equation is stiff (it oscillates at the rate c_s k) and the run can be slow.'
        write(*, '(a)') '   E.g. braiding with no kinetic term (alpha_K = 0) near the GR limit: c_s^2 ~ 1/alpha_B.'
    end if

    contains

    !> Subroutine that records the largest scalar sound speed squared, EFT_gradient/EFT_kinetic, after the turn-on.
    !! Only where both are positive beyond their roundoff floor and the model is within GLPV (where the ratio is
    !! c_s^2).
    subroutine track_sound_speed(a)

    real(dl), intent(in) :: a   !< the sampled scale factor (eft_cache holds the values there).

    real(dl) :: floor, kinetic_scale, gradient_scale
    real(dl), parameter :: domain_eps = 1.e-12_dl

    if (a < input_EFTCAMB%EFTCAMB_pert_turn_on) return
    if (abs(eft_cache%EFTGamma6V) > domain_eps*(1._dl + abs(eft_cache%EFTOmegaV)) .or. &
        abs(eft_cache%EFTGamma3V + eft_cache%EFTGamma4V) > &
        domain_eps*(1._dl + abs(eft_cache%EFTGamma3V) + abs(eft_cache%EFTGamma4V))) return
    call EFTCAMBStabilityScales(a, params_cache, eft_cache, kinetic_scale, gradient_scale)
    floor = max(input_EFTCAMB%EFTCAMB_stability_threshold, 1.e-12_dl)
    if (eft_cache%EFT_kinetic <= floor*kinetic_scale .or. eft_cache%EFT_gradient <= floor*gradient_scale) return
    if (eft_cache%EFT_gradient/eft_cache%EFT_kinetic > cs2_max) then
        cs2_max = eft_cache%EFT_gradient/eft_cache%EFT_kinetic
        a_cs2_max = a
    end if

    end subroutine track_sound_speed

    end subroutine EFTCAMB_Stability_Check

    ! ---------------------------------------------------------------------------------------------
    !> Function that fills the caches to check the stability of the theory.
    function EFTTestStability(a, k_max, input_EFTCAMB, params_cache, eft_cache)

    implicit none

    real(dl), intent(in)    :: a                       !< the input scale factor.
    real(dl), intent(inout) :: k_max                   !< the input maximum k mode at which stability is computed.
    !< the EFTCAMB object for which the code is computing stability.
    type(TEFTCAMB), intent(in) :: input_EFTCAMB
    !< the EFTCAMB parameter cache that contains all the physical parameters.
    type(TEFTCAMB_parameter_cache), intent(inout) :: params_cache
    !< the EFTCAMB timestep cache that contains all the physical values.
    type(TEFTCAMB_timestep_cache), intent(inout) :: eft_cache
    !< Logical value returned by the function. If the model is stable this is True, otherwise False.
    logical                                      :: EFTTestStability

    ! Definitions of variables:
    logical  :: EFT_HaveNan_parameter, EFT_HaveNan_timestep, mass_resolved, mass_unstable
    real(dl) :: EFT_instability_rate, tempk, temp1, temp2, temp3, temp4, temp5
    integer  :: ind_max, ind
    real(dl) :: test_dtauda, thr, kinetic_scale, gradient_scale
    real(dl), parameter :: domain_eps = 1.e-12_dl   ! roundoff floor of the GLPV domain test
    real(dl), external :: dtauda
    character(len=160) :: message
    real(dl) :: A_sign, A_margin, A_touch, a_zero, eom_G00, eom_G00_scale
    logical  :: beyond_GLPV, own_conditions, own_stable

    ! Stability check initialization
    EFTTestStability = .true.
    ! reset the time-step cache:
    call eft_cache%initialize()
    EFT_instability_reason = ''
    ! fill it:
    call EFTStabilityComputation(a, input_EFTCAMB, input_EFTCAMB%model, params_cache, eft_cache)
    ! protect against k_max too small:
    if (k_max < 0.1_dl) k_max = 0.1_dl
    ! relative numerical floor of the physical stability conditions: a condition passes if its value is
    ! >= -thr*(sum of the magnitudes of its terms), i.e. if it is non-negative up to cancellation roundoff.
    thr = max(input_EFTCAMB%EFTCAMB_stability_threshold, 0._dl)
    call EFTCAMBStabilityScales(a, params_cache, eft_cache, kinetic_scale, gradient_scale)

    ! check stability of the theory:

    ! 0) dtauda should be finite:
    test_dtauda = input_EFTCAMB%model%compute_dtauda(a, params_cache, eft_cache)
    if (test_dtauda > huge(test_dtauda) .or. IsNaN(test_dtauda)) then
        call fail('dtauda is not finite')
        if (input_EFTCAMB%EFTCAMB_feedback_level > 1) write(*, '(a)') '   Model dtauda is Nan'
        return
    end if

    ! 1) everything inside the parameter cache should not be a NaN:
    call params_cache%is_nan(EFT_HaveNan_parameter)
    if (EFT_HaveNan_parameter) then
        call fail('NaN in the parameter cache')
        if (input_EFTCAMB%EFTCAMB_feedback_level > 1) write(*, '(a)') '   Model has Nan in the parameter cache'
        return
    end if

    ! 2) everything inside the time-step cache should not be a NaN:
    call eft_cache%is_nan(EFT_HaveNan_timestep, input_EFTCAMB%EFTCAMB_feedback_level)
    if (EFT_HaveNan_timestep) then
        call fail('NaN in the time-step cache')
        if (input_EFTCAMB%EFTCAMB_feedback_level > 1) write(*, '(a)') '   Model has Nan in the timestep cache'
        return
    end if

    ! 2b) where the EFT perturbations are evolved and the model is not GR, the pi equation needs a kinetic term:
    !     A = A1 + k^2 A2 must not vanish, nor change sign in k (0 < k <= k_max) or in time. Without one pi is not
    !     dynamical and EFTCAMB does not solve its constraint. This is a limit of the solver, enforced whatever the
    !     stability flags. A negative A at every k and time is a ghost, which is left to the ghost checks.
    if (a >= input_EFTCAMB%EFTCAMB_pert_turn_on) then
        A_sign = 0._dl
        A_touch = 1._dl
        a_zero = a
        if (EFTPiHasKineticTerm(a, k_max, params_cache, eft_cache, max(thr, 1.e-12_dl), A_sign, A_margin)) then
            ! a sign change in time between two sampled scale factors:
            if (A_sign*PastAsign < 0._dl) A_sign = 0._dl
            PastAsign = A_sign
        end if
        ! A touching zero at an isolated time between two samples (e.g. w_DE = -1 at one time in quintessence, where
        ! A = c): the samples bracket a local minimum of A, which is refined. A minimum below 1e-10 of the bracketing
        ! values is a zero; near-touch models (1 + w_min >~ 1e-16) are kept.
        if (A_sign /= 0._dl) then
            if (PastAm(2) > 0._dl .and. PastAm(1) > PastAm(2) .and. A_margin > PastAm(2)) then
                A_touch = EFTPiKineticMinimum(min(PastAa(1), a), max(PastAa(1), a), A_sign, k_max, &
                    max(thr, 1.e-12_dl), input_EFTCAMB, params_cache, a_zero)/max(PastAm(1), A_margin)
                if (A_touch <= 1.e-10_dl) A_sign = 0._dl
            end if
            PastAa(1) = PastAa(2)
            PastAm(1) = PastAm(2)
            PastAa(2) = a
            PastAm(2) = A_margin
        else
            PastAm = 0._dl
        end if
        if (A_sign == 0._dl) then
            if (.not. EFTCAMBReturnToGR_is_GR(a_zero, input_EFTCAMB%model, params_cache, &
                input_EFTCAMB%EFTCAMB_GR_threshold)) then
                call fail('no pi kinetic term')
                if (A_touch <= 1.e-10_dl) then
                    write(message, '(a,es10.3,a)') 'EFTCAMB: the pi kinetic term vanishes (A = 0) at a =', a_zero, &
                        ' (e.g. w_DE touching -1): EFTCAMB cannot evolve through it'
                else
                    write(message, '(a,es10.3,a)') 'EFTCAMB: the pi field has no kinetic term (A = 0) at a =', a, &
                        ': EFTCAMB cannot evolve it, give it one (e.g. Gamma1, alpha_K)'
                end if
                call GlobalError(trim(message), 1)
                if (input_EFTCAMB%EFTCAMB_feedback_level > 1) write(*, '(a,es10.3)') &
                    '   No pi kinetic term: A1 + k^2 A2 vanishes or changes sign at a =', a_zero
                return
            end if
        end if
    end if

    ! 2c) where the EFT perturbations are evolved, the 00 Einstein constraint must determine h': its coefficient
    !     2 (1 + Omega) EFTeomG (2 in GR, 2 - alpha_B in the Bellini-Sawicki convention) must stay positive. At the
    !     alpha_B = 2 crossing (arXiv:2103.11195 App. A) z = h'/2k divides by zero and the integrator fails or stalls.
    !     Evolving through the crossing (h' from the trace equation) is not implemented: this is a limit of the
    !     solver, enforced whatever the stability flags.
    if (a >= input_EFTCAMB%EFTCAMB_pert_turn_on) then
        eom_G00 = 2._dl*(1._dl + eft_cache%EFTOmegaV) + a*eft_cache%EFTOmegaP + &
            a*params_cache%h0_Mpc*eft_cache%EFTGamma2V/eft_cache%adotoa + 3._dl*eft_cache%EFTGamma3V + &
            eft_cache%EFTGamma4V
        eom_G00_scale = 2._dl*abs(1._dl + eft_cache%EFTOmegaV) + abs(a*eft_cache%EFTOmegaP) + &
            abs(a*params_cache%h0_Mpc*eft_cache%EFTGamma2V/eft_cache%adotoa) + 3._dl*abs(eft_cache%EFTGamma3V) + &
            abs(eft_cache%EFTGamma4V)
        if (eom_G00 <= max(thr, 1.e-12_dl)*eom_G00_scale) then
            call fail('degenerate 00 Einstein constraint (alpha_B = 2)')
            write(message, '(a,es10.3,a)') 'EFTCAMB: the 00 Einstein constraint does not determine h'' at a =', a, &
                ' (alpha_B = 2 crossing): EFTCAMB cannot evolve through it'
            call GlobalError(trim(message), 1)
            if (input_EFTCAMB%EFTCAMB_feedback_level > 1) write(*, '(a,es11.4,a,es10.3)') &
                '   Degenerate 00 constraint: 2(1+Omega) EFTeomG =', eom_G00, ' at a =', a
            return
        end if
    end if

    ! 3) enforce ghost mathematical stability:
    if (input_EFTCAMB%EFT_ghost_math_stability) then

        ! 1- the A coefficient should not change sign in time and in k, i.e. it shall not be zero.
        !    This is the strongest stability constraint since violating it would violate the mathematical
        !    consistency of the pi field equation.
        !    The first condition is A1/=0. Implemented by detecting sign changes in A1.
        if (eft_cache%EFTpiA1*PastA1 < 0._dl) then
            call fail('mathematical ghost: A is zero in time')
            if (input_EFTCAMB%EFTCAMB_feedback_level > 1) write(*, &
                '(a)') '   Mathematical instability: A is zero in time'
        end if
        PastA1 = eft_cache%EFTpiA1
        ! The second one is the condition on k.
        if ((eft_cache%EFTpiA1 > 0 .and. eft_cache%EFTpiA1 + k_max**2*eft_cache%EFTpiA2 < 0) .or. &
            (eft_cache%EFTpiA1 < 0 .and. eft_cache%EFTpiA1 + k_max**2*eft_cache%EFTpiA2 > 0)) then
            call fail('mathematical ghost: A is zero in k')
            if (input_EFTCAMB%EFTCAMB_feedback_level > 1) write(*, '(a)') '   Mathematical instability: A is zero in k'
        end if

        ! 2- the AT coefficient should not change sign in time, i.e. it shall not be zero.
        !    This is the second strongest stability constraint since violating it would
        !    violate the mathematical consistency of the tensor perturbation equation.
        !    Implemented by detecting sign changes in AT.
        if (eft_cache%EFTAT*PastAT < 0._dl) then
            call fail('mathematical ghost: AT is zero in time')
            if (input_EFTCAMB%EFTCAMB_feedback_level > 1) write(*, &
                '(a)') '   Mathematical instability:  AT is zero in time'
        end if
        PastAT = eft_cache%EFTAT

    end if

    ! 4) enforce mass mathematical stability:
    if (input_EFTCAMB%EFT_mass_math_stability) then

        ! 1- we do not want (fast) growing exponential modes.
        !    This condition prevents the pi field from growing exponentially and destroying everything.
        !    Even though this condition is neither completely related to physics nor mathematics,
        !    violating it would completely mess up cosmological observables.

        ! This is the maximum allowed rate of instability. Units shall be Mpc^-1.
        EFT_instability_rate = 0._dl

        ! This condition needs to be tested in k. Sample in k.
        ind_max = 10

        do ind = 1, ind_max
            ! kmode to test. Linear sampling. Should suffice... (??)
            tempk = 0._dl + real(ind - 1)*(k_max)/real(ind_max - 1)
            ! vaule that discriminates between different cases:
            temp1 = (eft_cache%EFTpiB1 + eft_cache%EFTpiB2*tempk**2)
            temp2 = (eft_cache%EFTpiA1 + eft_cache%EFTpiA2*tempk**2)
            temp3 = temp1**2 - 4._dl*temp2*(eft_cache%EFTpiC + eft_cache%EFTpiD1*tempk**2 + eft_cache%EFTpiD2*tempk**4)

            ! case 1:
            if (temp3 > 0._dl .and. temp2 /= 0._dl) then
                temp4 = +0.5_dl*(-temp1 + sqrt(temp3))/temp2
                temp5 = +0.5_dl*(-temp1 - sqrt(temp3))/temp2
                if (temp4 > EFT_instability_rate .or. temp5 > EFT_instability_rate) then
                    call fail('growing exponential mode (EFT_mass_math_stability)')
                    if (input_EFTCAMB%EFTCAMB_feedback_level > 1) then
                        write(*, '(a,E11.4)')       '   Mathematical instability: growing exponential at k =', tempk
                        write(*, '(a,E11.4,E11.4)') '      Rate of instability: ', temp4, temp5
                    end if
                    exit
                end if
                ! case 2:
            else if (temp2 /= 0._dl) then
                temp4 = -0.5_dl*temp1/temp2
                if (temp4 > EFT_instability_rate) then
                    call fail('growing exponential mode (EFT_mass_math_stability)')
                    if (input_EFTCAMB%EFTCAMB_feedback_level > 1) then
                        write(*, '(a,E11.4)')       '   Mathematical instability: growing exponential at k =', tempk
                        write(*, '(a,E11.4,E11.4)') '      Rate of instability: ', temp4
                    end if
                    exit
                end if
            end if

        end do

    end if

    ! 5)-6) the generic scalar ghost and gradient conditions below extend up to GLPV. Beyond it (up to roundoff,
    ! domain_eps, not up to the physical floor thr: gamma6 and gamma3 + gamma4 enter the dispersion relation times
    ! (k/(aH))^2 ~ 1e9) a model with its own scalar conditions (Horava, 3+1) is tested with them instead, whatever
    ! EFT_additional_priors; any other model beyond GLPV is refused. The positive gravitational constant and the
    ! tensor conditions are generic.
    beyond_GLPV = abs(eft_cache%EFTGamma6V) > domain_eps*(1._dl + abs(eft_cache%EFTOmegaV)) .or. &
        abs(eft_cache%EFTGamma3V + eft_cache%EFTGamma4V) > &
        domain_eps*(1._dl + abs(eft_cache%EFTGamma3V) + abs(eft_cache%EFTGamma4V))
    own_conditions = .false.
    own_stable = .true.
    if (input_EFTCAMB%EFT_ghost_stability .or. input_EFTCAMB%EFT_gradient_stability) then
        call input_EFTCAMB%model%beyond_GLPV_stability(a, params_cache, eft_cache, input_EFTCAMB%EFT_ghost_stability, &
            input_EFTCAMB%EFT_gradient_stability, own_conditions, own_stable)
        if (own_conditions .and. .not. own_stable) then
            call fail('model specific ghost or gradient condition')
            if (input_EFTCAMB%EFTCAMB_feedback_level > 1) write(*, '(a,E11.4)') &
                '   Physical instability: the model specific ghost or gradient condition fails at a =', a
        end if
    end if

    ! 5) enforce ghost viability:
    if (input_EFTCAMB%EFT_ghost_stability) then

        if (beyond_GLPV .and. .not. own_conditions) then
            if (input_EFTCAMB%EFTCAMB_feedback_level > 0) then
                write(*, '(a)') '   EFTCAMB WARNING: stability for model beyond GLPV has not been worked out.'
                write(*, '(a)') '      It will be added in a future release.'
                write(*, '(a)') '      If you want to run this model disable EFT_ghost_stability.'
            end if
            call fail('beyond GLPV without ghost conditions (EFT_ghost_stability)')
            return
        end if

        ! 1- Positive gravitational constant:
        if (1._dl + eft_cache%EFTOmegaV <= -thr*(1._dl + abs(eft_cache%EFTOmegaV))) then
            call fail('negative gravitational constant')
            if (input_EFTCAMB%EFTCAMB_feedback_level > 1) write(*, &
                '(a,E11.4)') '   Physical instability: negative gravitational constant = ', 1._dl + eft_cache%EFTOmegaV
        end if

        ! 2- Ghost condition (up to GLPV):
        if (.not. beyond_GLPV .and. eft_cache%EFT_kinetic < -thr*kinetic_scale) then
            call fail('ghost instability')
            if (input_EFTCAMB%EFTCAMB_feedback_level > 1) write(*, &
                '(a,E11.4)') '   Physical instability: ghost instability. Kinetic term = ', eft_cache%EFT_kinetic
        end if

        ! 3- No tensor ghosts:
        if (eft_cache%EFTAT < -thr*(1._dl + abs(eft_cache%EFTOmegaV) + abs(eft_cache%EFTGamma4V))) then
            call fail('tensor ghost instability')
            if (input_EFTCAMB%EFTCAMB_feedback_level > 1) write(*, &
                '(a,E11.4)') '   Physical instability: tensor ghost instability. Tensor kinetic term = ', &
                eft_cache%EFTAT
        end if

    end if

    ! 6) enforce gradient viability:
    if (input_EFTCAMB%EFT_gradient_stability) then

        if (beyond_GLPV .and. .not. own_conditions) then
            if (input_EFTCAMB%EFTCAMB_feedback_level > 0) then
                write(*, '(a)') '   EFTCAMB WARNING: stability for model beyond GLPV has not been worked out.'
                write(*, '(a)') '      It will be added in a future release.'
                write(*, '(a)') '      If you want to run this model disable EFT_gradient_stability.'
            end if
            call fail('beyond GLPV without gradient conditions (EFT_gradient_stability)')
            return
        end if

        ! 1- Positive gravitational constant:
        if (1._dl + eft_cache%EFTOmegaV <= -thr*(1._dl + abs(eft_cache%EFTOmegaV))) then
            call fail('negative gravitational constant')
            if (input_EFTCAMB%EFTCAMB_feedback_level > 1) write(*, &
                '(a,E11.4)') '   Physical instability: negative gravitational constant = ', 1._dl + eft_cache%EFTOmegaV
        end if

        ! 2- Gradient instability (up to GLPV):
        if (.not. beyond_GLPV .and. eft_cache%EFT_gradient < -thr*gradient_scale) then
            call fail('gradient instability')
            if (input_EFTCAMB%EFTCAMB_feedback_level > 1) write(*, &
                '(a,E11.4)') '   Physical instability: gradient instability. Gradient term = ', eft_cache%EFT_gradient
        end if

        ! 3- No tensor gradient:
        if (eft_cache%EFTDT < -thr*(1._dl + abs(eft_cache%EFTOmegaV))) then
            call fail('tensor gradient instability')
            if (input_EFTCAMB%EFTCAMB_feedback_level > 1) write(*, &
                '(a,E11.4)') '   Physical instability: tensor gradient instability. Tensor gradient term = ', &
                eft_cache%EFTDT
        end if

    end if

    ! 7) mass viability (EFT_mass_stability, off by default): tachyonic modes of the scalar sector at k -> 0.
    !    EFTStabilityMassEigenvalues gives the masses EFT_mu1 (scalar) and EFT_mu2 (matter), cosmic time, with the
    !    method EFT_mass_stability_method. A mode grows as exp(int sqrt(-mu) dt): the model is rejected if
    !    mu < -(EFT_mass_stability_rate H)^2, i.e. if a mode grows faster than EFT_mass_stability_rate e-folds per
    !    Hubble time. Tested after the perturbation turn-on, where the masses are numerically resolved.
    !    Method 1 (arXiv:1609.03599 App. C, corrected) tests both eigenvalues. Method 2 tests the normal-form mass of
    !    pi, and only where pi itself also grows faster than the rate (largest root of A1 s^2 + B1 s + C = 0 at
    !    frozen coefficients): a fast change of the pi normalisation (hilltop quintessence, where pi = delta phi/phi'
    !    stays bounded while delta phi follows the rolling background) gives mu_1 < 0 without a growing pi;
    !    EFT_mu2 = -rho_dust/2 is reported, not tested.
    if (input_EFTCAMB%EFT_mass_stability .and. a > input_EFTCAMB%EFTCAMB_pert_turn_on) then
        call EFTStabilityMassEigenvalues(a, input_EFTCAMB, params_cache, eft_cache, mass_resolved, &
            floor=max(thr, 1.e-10_dl))
        if (mass_resolved) then
            if (input_EFTCAMB%EFT_mass_stability_method == 2) then
                mass_unstable = eft_cache%EFT_mu1 < -(input_EFTCAMB%EFT_mass_stability_rate*eft_cache%adotoa/a)**2 &
                    .and. EFTPiFrozenGrowthRate(eft_cache) > input_EFTCAMB%EFT_mass_stability_rate*eft_cache%adotoa
            else
                mass_unstable = min(eft_cache%EFT_mu1, eft_cache%EFT_mu2) < &
                    -(input_EFTCAMB%EFT_mass_stability_rate*eft_cache%adotoa/a)**2
            end if
            if (mass_unstable) then
                call fail('tachyonic mass (EFT_mass_stability)')
                if (input_EFTCAMB%EFTCAMB_feedback_level > 1) write(*, '(a,2E12.4)') &
                    '   Physical instability: tachyonic mass, mu_1/H^2, mu_2/H^2 = ', &
                    eft_cache%EFT_mu1*(a/eft_cache%adotoa)**2, eft_cache%EFT_mu2*(a/eft_cache%adotoa)**2
            end if
        end if
    end if
    ! the former check (experimental, it never rejected a model):
    ! ! 7) mass viability: EXPERIMENTAL, never rejects a model.
    ! !    The mass eigenvalues mu_1, mu_2 reject canonical quintessence (|mu|/H^2 ~ 1e9) and are dominated by
    ! !    roundoff near GR (audit findings F-042, F-063): until they are re-derived the check only reports.
    ! if (input_EFTCAMB%EFT_mass_stability .and. a > input_EFTCAMB%EFTCAMB_pert_turn_on) then
    !     if (input_EFTCAMB%EFTCAMB_feedback_level > 2) then
    !         if (eft_cache%EFT_mu1 < -(input_EFTCAMB%EFT_mass_stability_rate*eft_cache%adotoa/a)**2) &
    !             write(*, '(a,E11.4)') '   Experimental mass check (ignored): mu_1/H^2 = ', &
    !             eft_cache%EFT_mu1*(a/eft_cache%adotoa)**2
    !         if (eft_cache%EFT_mu2 < -(input_EFTCAMB%EFT_mass_stability_rate*eft_cache%adotoa/a)**2) &
    !             write(*, '(a,E11.4)') '   Experimental mass check (ignored): mu_2/H^2 = ', &
    !             eft_cache%EFT_mu2*(a/eft_cache%adotoa)**2
    !     end if
    ! end if

    ! 8) enforce model specific priors:
    if (input_EFTCAMB%EFT_additional_priors) then
        if (.not. input_EFTCAMB%model%additional_model_stability(a, params_cache, eft_cache)) then
            call fail('model specific prior (EFT_additional_priors)')
            if (input_EFTCAMB%EFTCAMB_feedback_level > 1) write(*, &
                '(a)') '   Model specific stability criteria are not met'
        end if
    end if

    ! 8b) optional phenomenological prior w_DE <= -1/3, with its own refusal message:
    if (input_EFTCAMB%EFT_wDE_prior) then
        if (.not. input_EFTCAMB%model%wDE_prior(a, params_cache, eft_cache)) then
            call fail('w_DE > -1/3 (EFT_wDE_prior)')
            write(message, '(a,es10.3,a)') 'EFTCAMB: w_DE > -1/3 at a =', a, ', refused by EFT_wDE_prior'
            call GlobalError(trim(message), 1)
            if (input_EFTCAMB%EFTCAMB_feedback_level > 1) write(*, '(a,es10.3)') &
                '   w_DE prior violated: w_DE > -1/3 at a =', a
        end if
    end if

    ! 9) enforce positivity bounds:
    ! works up to Horndeski
    if (input_EFTCAMB%EFT_positivity_bounds) then
        if ((eft_cache%posbound1 < 0._dl) .or. (eft_cache%posbound2 < 0._dl)) then
            call fail('positivity bounds')
            if (input_EFTCAMB%EFTCAMB_feedback_level > 1) write(*, '(a)') ' Positivity bounds are violated'
        end if
    end if

    ! 10) check the existence of a healthy Minkowski limit:
    ! works up to Horndeski
    if (input_EFTCAMB%EFT_minkowski_limit) then
        if ((-eft_cache%ghostpos <= 0) .or. (-eft_cache%tachpos < 0)) then
            call fail('no healthy Minkowski limit')
            if (input_EFTCAMB%EFTCAMB_feedback_level > 1) write(*, '(a)') 'Healthy Minkowski limit does not exist'
        end if
    end if


    contains

    !> Marks the theory unstable and records the first reason (it goes into the error message, with a).
    subroutine fail(reason)

    character(len=*), intent(in) :: reason

    EFTTestStability = .false.
    if (len_trim(EFT_instability_reason) == 0) EFT_instability_reason = reason

    end subroutine fail

    end function EFTTestStability

    ! ---------------------------------------------------------------------------------------------
    !> Function that tells whether the pi field equation has a kinetic term at every 0 < k <= k_max: A = A1 + k^2 A2
    !! does not vanish (relative to the magnitude of its terms, see compute_pi_factors) nor change sign in k.
    function EFTPiHasKineticTerm(a, k_max, params_cache, eft_cache, eps, A_sign, A_margin)

    implicit none

    real(dl), intent(in) :: a                     !< the input scale factor.
    real(dl), intent(in) :: k_max                 !< the maximum k mode of the stability check.
    type(TEFTCAMB_parameter_cache), intent(in) :: params_cache !< the EFTCAMB parameter cache.
    type(TEFTCAMB_timestep_cache), intent(in)  :: eft_cache    !< the EFTCAMB timestep cache, with the pi factors.
    real(dl), intent(in) :: eps                   !< relative floor below which A counts as zero.
    real(dl), intent(out) :: A_sign               !< the sign of A (at every k) when it has one, 0 otherwise.
    real(dl), intent(out) :: A_margin             !< min over 0 < k <= k_max of A_sign*A (0 without a kinetic term).
    logical :: EFTPiHasKineticTerm

    real(dl) :: A1, A2, scale1, A_max, sign0
    logical  :: A1_zero

    A1 = eft_cache%EFTpiA1
    A2 = eft_cache%EFTpiA2
    ! magnitude of the terms of A1: c, gamma1 and the mixing with the metric (Omega', gamma2):
    scale1 = abs(eft_cache%EFTc) + 2._dl*a**2*params_cache%h0_Mpc**2*abs(eft_cache%EFTGamma1V) + &
        abs(1.5_dl*a**2*(eft_cache%adotoa*eft_cache%EFTOmegaP + params_cache%h0_Mpc*eft_cache%EFTGamma2V)**2/ &
        (2._dl*(1._dl + eft_cache%EFTOmegaV) + 3._dl*eft_cache%EFTGamma3V + eft_cache%EFTGamma4V))
    ! sign of A for k -> 0 (A2 k^2 when A1 vanishes):
    A1_zero = .false.
    if (abs(A1) > eps*scale1) then
        sign0 = sign(1._dl, A1)
    else if (A2 /= 0._dl) then
        sign0 = sign(1._dl, A2)
        A1_zero = .true.
    else
        A_sign = 0._dl
        A_margin = 0._dl
        EFTPiHasKineticTerm = .false.
        return
    end if
    ! A is linear in k^2: it vanishes or changes sign in (0, k_max] only if it does so at k_max:
    A_max = A1 + k_max**2*A2
    EFTPiHasKineticTerm = abs(A_max) > eps*(scale1 + k_max**2*abs(A2)) .and. sign0*A_max > 0._dl
    A_sign = 0._dl
    A_margin = 0._dl
    if (EFTPiHasKineticTerm) then
        A_sign = sign0
        ! with A1 = 0 (up to the floor: roundoff of either sign) A = A2 k^2 vanishes only at k = 0, which is allowed:
        ! the margin is the value at k_max
        if (A1_zero) then
            A_margin = sign0*A_max
        else
            A_margin = min(sign0*A1, sign0*A_max)
        end if
    end if

    end function EFTPiHasKineticTerm

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the minimum over [a_lo, a_hi] (golden-section search) of the margin of the pi kinetic
    !! coefficient (as in EFTPiHasKineticTerm: 0 where A has no kinetic term or the other sign), and where it is
    !! reached.
    function EFTPiKineticMinimum(a_lo, a_hi, A_sign, k_max, eps, input_EFTCAMB, params_cache, a_min)

    implicit none

    real(dl), intent(in)       :: a_lo, a_hi    !< the scale factor interval.
    real(dl), intent(in)       :: A_sign        !< the sign of A in the interval.
    real(dl), intent(in)       :: k_max         !< the maximum k mode of the stability check.
    real(dl), intent(in)       :: eps           !< relative floor below which A counts as zero.
    type(TEFTCAMB), intent(in) :: input_EFTCAMB !< the EFTCAMB object.
    type(TEFTCAMB_parameter_cache), intent(inout) :: params_cache !< the EFTCAMB parameter cache.
    real(dl), intent(out)      :: a_min         !< the scale factor of the minimum.
    real(dl) :: EFTPiKineticMinimum

    real(dl), parameter :: golden = 0.6180339887498949_dl
    type(TEFTCAMB_timestep_cache) :: cache
    real(dl) :: lo, hi, x1, x2, f1, f2
    integer  :: i

    lo = a_lo
    hi = a_hi
    x1 = hi - golden*(hi - lo)
    x2 = lo + golden*(hi - lo)
    f1 = margin(x1)
    f2 = margin(x2)
    do i = 1, 200
        if (hi - lo <= 1.e-15_dl*hi) exit
        if (f1 < f2) then
            hi = x2
            x2 = x1
            f2 = f1
            x1 = hi - golden*(hi - lo)
            f1 = margin(x1)
        else
            lo = x1
            x1 = x2
            f1 = f2
            x2 = lo + golden*(hi - lo)
            f2 = margin(x2)
        end if
    end do
    if (f1 < f2) then
        a_min = x1
        EFTPiKineticMinimum = f1
    else
        a_min = x2
        EFTPiKineticMinimum = f2
    end if

    contains

    real(dl) function margin(x)
    real(dl), intent(in) :: x
    real(dl) :: x_sign, x_margin
    call cache%initialize()
    call EFTStabilityComputation(x, input_EFTCAMB, input_EFTCAMB%model, params_cache, cache)
    margin = 0._dl
    if (EFTPiHasKineticTerm(x, k_max, params_cache, cache, eps, x_sign, x_margin)) then
        if (x_sign == A_sign) margin = x_margin
    end if
    end function margin

    end function EFTPiKineticMinimum

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that restores the values stored in the module to the default.
    !! Needed for successive calls.
    subroutine EFTStability_cleanup

    implicit none

    PastA1 = 0._dl
    PastAT = 0._dl
    PastAsign = 0._dl
    PastAa = 0._dl
    PastAm = 0._dl

    end subroutine EFTStability_cleanup

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that fills the caches to check the stability of the theory.
    subroutine EFTStabilityComputation(a, input_EFTCAMB, input_model, params_cache, eft_cache)

    implicit none

    real(dl), intent(in)       :: a                       !< the input scale factor.
    !< the EFTCAMB model for which the code is computing the RGR time.
    type(TEFTCAMB), intent(in) :: input_EFTCAMB
    !< the EFTCAMB model for which the code is computing the RGR time.
    class(EFTCAMB_model), intent(in) :: input_model
    !< the EFTCAMB parameter cache that contains all the physical parameters.
    type(TEFTCAMB_parameter_cache), intent(inout) :: params_cache
    !< the EFTCAMB timestep cache that contains all the physical values.
    type(TEFTCAMB_timestep_cache), intent(inout) :: eft_cache

    ! Definitions of variables:
    real(dl) :: grhonu, gpinu, grhormass_t, grhonudot, gpinudot, gpinudotdot
    integer  :: nu_i, ind, ind_max

    ! start filling the cache:
    eft_cache%a = a
    ! compute background densities of different species
    eft_cache%grhob_t = params_cache%grhob/a         ! 8\pi G_N \rho_b a^2: baryons background density
    eft_cache%grhoc_t = params_cache%grhoc/a         ! 8\pi G_N \rho_{cdm} a^2: cold dark matter background density
    eft_cache%grhor_t = params_cache%grhornomass/a/a ! 8\pi G_N \rho_{\nu} a^2: massless neutrinos background density
    eft_cache%grhog_t = params_cache%grhog/a/a       ! 8\pi G_N \rho_{\gamma} a^2: radiation background density
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
    ! total matter density times a^4, used by compute_dtauda (as in dtauda of CAMB):
    eft_cache%grhoa2 = eft_cache%grhom_t*a**2
    ! compute the other things:
    select type (model_temp => input_model)
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
    ! compute massive neutrinos stuff:
    ! Massive neutrinos mod:
    eft_cache%gpinu_tot = 0._dl
    eft_cache%grhonu_tot = 0._dl
    eft_cache%grhonudot_tot = 0._dl
    eft_cache%gpinudot_tot = 0._dl
    if (params_cache%Num_Nu_Massive /= 0) then
        do nu_i = 1, params_cache%Nu_mass_eigenstates
            grhonu = 0._dl
            gpinu = 0._dl
            grhonudot = 0._dl
            gpinudot = 0._dl
            grhormass_t = params_cache%grhormass(nu_i)/a**2
            call ThermalNuBackground%rho_P(a*params_cache%nu_masses(nu_i), grhonu, gpinu)
            eft_cache%grhonu_tot = eft_cache%grhonu_tot + grhormass_t*grhonu
            eft_cache%gpinu_tot = eft_cache%gpinu_tot + grhormass_t*gpinu
            ! continuity equation, consistent with the rho and p fits (the drho fit is not, at 1e-4: U7)
            eft_cache%grhonudot_tot = eft_cache%grhonudot_tot - 3._dl*eft_cache%adotoa*grhormass_t*(grhonu + gpinu)
            eft_cache%gpinudot_tot = eft_cache%gpinudot_tot + &
                grhormass_t*(ThermalNuBackground%pidot(a*params_cache%nu_masses(nu_i), eft_cache%adotoa, gpinu) - &
                4._dl*eft_cache%adotoa*gpinu)
            eft_cache%gpinudotdot_tot = 0._dl
        end do
    end if
    ! compute pressure dot:
    eft_cache%gpresdotm_t = &
        -4._dl*eft_cache%adotoa*(eft_cache%grhog_t + eft_cache%grhor_t)/3._dl + eft_cache%gpinudot_tot
    ! compute remaining quantities related to H:
    call input_model%compute_H_derivs(a, params_cache, eft_cache)
    ! compute second derivative of massive neutrino pressure
    eft_cache%gpinu_tot = 0._dl
    eft_cache%grhonu_tot = 0._dl
    eft_cache%grhonudot_tot = 0._dl
    eft_cache%gpinudot_tot = 0._dl
    eft_cache%gpinudotdot_tot = 0._dl
    if (params_cache%Num_Nu_Massive /= 0) then
        do nu_i = 1, params_cache%Nu_mass_eigenstates
            grhonu = 0._dl
            gpinu = 0._dl
            grhonudot = 0._dl
            gpinudot = 0._dl
            grhormass_t = params_cache%grhormass(nu_i)/a**2
            call ThermalNuBackground%rho_P(a*params_cache%nu_masses(nu_i), grhonu, gpinu)
            eft_cache%grhonu_tot = eft_cache%grhonu_tot + grhormass_t*grhonu
            eft_cache%gpinu_tot = eft_cache%gpinu_tot + grhormass_t*gpinu
            ! continuity equation, consistent with the rho and p fits (the drho fit is not, at 1e-4: U7)
            eft_cache%grhonudot_tot = eft_cache%grhonudot_tot - 3._dl*eft_cache%adotoa*grhormass_t*(grhonu + gpinu)
            gpinudot = ThermalNuBackground%pidot(a*params_cache%nu_masses(nu_i), eft_cache%adotoa, gpinu)
            eft_cache%gpinudot_tot = eft_cache%gpinudot_tot + grhormass_t*(gpinudot &
                - 4._dl*eft_cache%adotoa*gpinu)
            eft_cache%gpinudotdot_tot = eft_cache%gpinudotdot_tot - &
                4._dl*eft_cache%adotoa*grhormass_t*(gpinudot - 4._dl*eft_cache%adotoa*gpinu) + &
                grhormass_t*(ThermalNuBackground%pidotdot(a*params_cache%nu_masses(nu_i), eft_cache%adotoa, &
                eft_cache%Hdot, gpinu, gpinudot) - 4._dl*eft_cache%Hdot*gpinu - 4._dl*eft_cache%adotoa*gpinudot)
        end do
    end if
    ! compute pressure ddot:
    eft_cache%gpresdotdotm_t = -(4._dl/3._dl)*eft_cache%Hdot*(eft_cache%grhog_t + eft_cache%grhor_t) + &
        (16._dl/3._dl)*(eft_cache%grhog_t + eft_cache%grhor_t)*eft_cache%adotoa**2 + eft_cache%gpinudotdot_tot
    ! compute third derivatives of H:
    call input_model%compute_H_derivs(a, params_cache, eft_cache)
    ! compute backgrond EFT functions if model is designer:
    select type (model_temp => input_model)
    class is (EFTCAMB_designer_model)
        call input_model%compute_background_EFT_functions(a, params_cache, eft_cache)
    end select
    ! compute all other background stuff:
    call input_model%compute_rhoQPQ(a, params_cache, eft_cache)
    ! compute second order EFT functions:
    call input_model%compute_secondorder_EFT_functions(a, params_cache, eft_cache)
    ! Compute pi field equations factors:
    call input_model%compute_pi_factors(a, params_cache, eft_cache)
    ! Compute coefficients for the tensor propagation equation:
    call input_model%compute_tensor_factors(a, params_cache, eft_cache)
    ! Compute kinetic, gradient and mass terms:
    call input_model%compute_stability_factors(a, params_cache, eft_cache)
    ! Higher derivatives of the EFT functions used by the mass eigenvalues of method 1, the positivity bounds and the
    ! Minkowski limit: compute them first (models that do not provide them analytically would otherwise use
    ! zeros), then recompute the Hubble derivatives so that Hdotdotdot includes EFTLambdadotdot and EFTOmegaPPPP.
    ! (The normal-form mass of method 2 needs only the pi factors.)
    if (input_EFTCAMB%EFT_positivity_bounds .or. input_EFTCAMB%EFT_minkowski_limit .or. &
        (input_EFTCAMB%EFT_mass_stability .and. input_EFTCAMB%EFT_mass_stability_method == 1 .and. &
        a > input_EFTCAMB%EFTCAMB_pert_turn_on)) then
        call input_model%compute_additional_derivs(a, params_cache, eft_cache)
        call input_model%compute_H_derivs(a, params_cache, eft_cache)
    end if
    ! Compute positivity bounds:
    call input_model%compute_positivity_bounds(a, params_cache, eft_cache)
    ! The mass eigenvalues are computed by the mass check (step 7 of EFTTestStability), which needs to know whether
    ! they are resolved; formerly here:
    ! ! Compute the mass eigenvalues:
    ! if (a > input_EFTCAMB%EFTCAMB_pert_turn_on .and. input_EFTCAMB%EFT_mass_stability) then
    !     call input_model%compute_mass_factors(a, params_cache, eft_cache)
    ! end if

    ! dump cache if in debug mode:
    if (DebugEFTCAMB) then
        call eft_cache%dump_cache_files()
    end if

    end subroutine EFTStabilityComputation

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes the k -> 0 masses EFT_mu1 (scalar) and EFT_mu2 (matter) at scale factor a, cosmic time
    !! (Mpc^-2), with the method EFT_mass_stability_method; eft_cache must be filled at a (EFTStabilityComputation,
    !! or the evolution). resolved = .false. where they are not numerically resolved (EFT_mu1, EFT_mu2 are then 0, or
    !! EFT_mu2 = -rho_dust/2 for method 2).
    !! Method 1: compute_mass_factors, the k -> 0 eigenvalues of arXiv:1609.03599 App. C (corrected).
    !! Method 2: the normal form of the k -> 0 pi equation A1 pi'' + B1 pi' + C pi = 0 (conformal time): with
    !!  p = EFTpiB1/EFTpiA1, q = EFTpiC/EFTpiA1 and pi = u exp(-int p/2) it reads u'' + (q - p^2/4 - p'/2) u = 0, and
    !!      EFT_mu1 = (q - p^2/4 - p'/2 + adotoa^2 + Hdot)/a^2,
    !!  the mass relative to a massless minimally coupled field (u'' - (a''/a) u = 0). EFT_mu1 does not change under
    !!  time-dependent rescalings of pi; it is V'' for canonical quintessence and the scalaron mass for f(R).
    !!  EFT_mu2 = -rho_dust/2 (the matter Jeans term). p' is a central difference of stability caches in ln a with
    !!  step mass_dlna (second-order backward where a exp(mass_dlna) > 1). Not resolved where EFTpiA1 is zero or below
    !!  floor times the sum of the magnitudes of its terms (no pi kinetic term at k -> 0: auxiliary pi, Horava, gamma6)
    !!  or changes sign or by more than |d ln EFTpiA1/d ln a| = mass_max_dlnA across the stencil (A1 -> 0, e.g. phi' = 0
    !!  for quintessence: there the truncation error of -p'/2 stays below 0.01 H^2 outside +-0.03 e-folds, where a
    !!  tachyon of rate 3 H grows by at most 20%). With recompute = .true. the pi factors at a are also taken from a
    !!  stability cache (for the evolution output: the same values as the stability check).
    subroutine EFTStabilityMassEigenvalues(a, input_EFTCAMB, params_cache, eft_cache, resolved, floor, recompute)

    implicit none

    real(dl), intent(in)       :: a              !< the input scale factor.
    type(TEFTCAMB), intent(in) :: input_EFTCAMB  !< the EFTCAMB object.
    !< the EFTCAMB parameter cache that contains all the physical parameters.
    type(TEFTCAMB_parameter_cache), intent(inout) :: params_cache
    !< the EFTCAMB timestep cache, filled at a. On output EFT_mu1, EFT_mu2 are set.
    type(TEFTCAMB_timestep_cache), intent(inout) :: eft_cache
    logical, intent(out)            :: resolved   !< whether EFT_mu1, EFT_mu2 are numerically resolved.
    real(dl), intent(in), optional  :: floor      !< relative numerical floor of EFTpiA1 (method 2; default 1e-10).
    logical, intent(in), optional   :: recompute  !< method 2: take the pi factors at a from a stability cache.

    real(dl), parameter :: mass_dlna = 1.e-4_dl      ! finite-difference step in ln a
    real(dl), parameter :: mass_max_dlnA = 63._dl    ! largest resolved |d ln EFTpiA1/d ln a|
    type(TEFTCAMB_timestep_cache) :: cache_0, cache_1, cache_2
    real(dl) :: a_1, a_2, p0, p1, p2, pdot, q, F1, A1_scale, dlnA, A1_floor
    logical  :: backward

    resolved = .false.
    if (input_EFTCAMB%EFT_mass_stability_method /= 2) then
        call input_EFTCAMB%model%compute_mass_factors(a, params_cache, eft_cache, resolved)
        return
    end if

    A1_floor = 1.e-10_dl
    if (present(floor)) A1_floor = floor
    ! the pi factors at a:
    call cache_0%initialize()
    if (present(recompute)) then
        if (recompute) then
            call EFTStabilityComputation(a, input_EFTCAMB, input_EFTCAMB%model, params_cache, cache_0)
        else
            cache_0 = eft_cache
        end if
    else
        cache_0 = eft_cache
    end if
    eft_cache%EFT_mu1 = 0._dl
    eft_cache%EFT_mu2 = -0.5_dl*(cache_0%grhob_t + cache_0%grhoc_t)/a**2
    if (cache_0%EFTpiA1 == 0._dl) return
    ! EFTpiA1 = a^2 K/(2 F1), K = F1 (2c + 4 H0^2 gamma1) + 3 (dOmega/dt + H0 gamma2)^2 (cosmic): numerical floor
    F1 = 2._dl*(1._dl + cache_0%EFTOmegaV) + 3._dl*cache_0%EFTGamma3V + cache_0%EFTGamma4V
    if (F1 == 0._dl) return
    A1_scale = (abs(F1)*(2._dl*abs(cache_0%EFTc) + 4._dl*a**2*params_cache%h0_Mpc**2*abs(cache_0%EFTGamma1V)) + &
        3._dl*a**2*(abs(cache_0%adotoa*cache_0%EFTOmegaP) + params_cache%h0_Mpc*abs(cache_0%EFTGamma2V))**2)/ &
        (2._dl*abs(F1))
    if (.not. abs(cache_0%EFTpiA1) > A1_floor*A1_scale) return
    ! pi factors at the neighbouring times:
    backward = a*exp(mass_dlna) > 1._dl
    a_1 = a*exp(-mass_dlna)
    if (backward) then
        a_2 = a*exp(-2._dl*mass_dlna)
    else
        a_2 = a*exp(mass_dlna)
    end if
    call cache_1%initialize()
    call EFTStabilityComputation(a_1, input_EFTCAMB, input_EFTCAMB%model, params_cache, cache_1)
    call cache_2%initialize()
    call EFTStabilityComputation(a_2, input_EFTCAMB, input_EFTCAMB%model, params_cache, cache_2)
    if (cache_1%EFTpiA1*cache_0%EFTpiA1 <= 0._dl .or. cache_2%EFTpiA1*cache_0%EFTpiA1 <= 0._dl) return
    if (backward) then
        dlnA = abs(log(cache_0%EFTpiA1/cache_2%EFTpiA1))/(2._dl*mass_dlna)
    else
        dlnA = abs(log(cache_2%EFTpiA1/cache_1%EFTpiA1))/(2._dl*mass_dlna)
    end if
    if (.not. dlnA <= mass_max_dlnA) return
    p0 = cache_0%EFTpiB1/cache_0%EFTpiA1
    p1 = cache_1%EFTpiB1/cache_1%EFTpiA1
    p2 = cache_2%EFTpiB1/cache_2%EFTpiA1
    if (backward) then
        pdot = cache_0%adotoa*(3._dl*p0 - 4._dl*p1 + p2)/(2._dl*mass_dlna)
    else
        pdot = cache_0%adotoa*(p2 - p1)/(2._dl*mass_dlna)
    end if
    q = cache_0%EFTpiC/cache_0%EFTpiA1
    eft_cache%EFT_mu1 = (q - 0.25_dl*p0**2 - 0.5_dl*pdot + cache_0%adotoa**2 + cache_0%Hdot)/a**2
    resolved = eft_cache%EFT_mu1 == eft_cache%EFT_mu1
    if (.not. resolved) eft_cache%EFT_mu1 = 0._dl

    end subroutine EFTStabilityMassEigenvalues

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the growth rate of pi at k -> 0 with frozen coefficients (conformal time, Mpc^-1): the
    !! largest real part of the roots of EFTpiA1 s^2 + EFTpiB1 s + EFTpiC = 0. Used by the mass condition of method 2
    !! to require that pi itself grows, not only its normal-form amplitude.
    function EFTPiFrozenGrowthRate(eft_cache)

    implicit none

    type(TEFTCAMB_timestep_cache), intent(in) :: eft_cache  !< the EFTCAMB timestep cache.
    real(dl) :: EFTPiFrozenGrowthRate
    real(dl) :: disc

    disc = eft_cache%EFTpiB1**2 - 4._dl*eft_cache%EFTpiA1*eft_cache%EFTpiC
    if (disc >= 0._dl) then
        EFTPiFrozenGrowthRate = max((-eft_cache%EFTpiB1 + sqrt(disc))/(2._dl*eft_cache%EFTpiA1), &
            (-eft_cache%EFTpiB1 - sqrt(disc))/(2._dl*eft_cache%EFTpiA1))
    else
        EFTPiFrozenGrowthRate = -eft_cache%EFTpiB1/(2._dl*eft_cache%EFTpiA1)
    end if

    end function EFTPiFrozenGrowthRate

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes the effective CPL equation of state (w0, wa) of the model background, as used by
    !! Halofit/HMcode. The effective dark energy is everything in the expansion history that is not matter,
    !! radiation or neutrinos: rho_DE = 3 H^2 - rho_m, p_DE = -(2 dH/dtau + H^2) - p_m (a^2-weighted, conformal H),
    !! which equals rho_Q, P_Q of compute_rhoQPQ when Omega = 0 and gives the wDE parametrization exactly for
    !! designer models. w(a) = p_DE/rho_DE is fitted by least squares with w0 + wa (1-a) on 51 points uniform in
    !! a in [0.5, 1] (z <= 1, where the non-linear corrections matter). Call after the background initialization.
    subroutine EFTCAMB_compute_effective_w0wa(input_EFTCAMB, params_cache)

    implicit none

    type(TEFTCAMB), intent(inout) :: input_EFTCAMB                  !< the EFTCAMB object
    type(TEFTCAMB_parameter_cache), intent(inout) :: params_cache   !< the EFTCAMB parameter cache

    integer, parameter :: num_points = 51
    real(dl), parameter :: a_min = 0.5_dl
    type(TEFTCAMB) :: temp_EFTCAMB
    type(TEFTCAMB_timestep_cache) :: eft_cache
    real(dl) :: a, x, w, rho_de, p_de, s0, s1, s2, sw, sxw, det
    integer  :: ind

    if (input_EFTCAMB%EFTflag == 0 .or. .not. allocated(input_EFTCAMB%model)) return
    input_EFTCAMB%model%background_w0wa_computed = .false.
    ! only the background is needed: no optional stability derivatives.
    temp_EFTCAMB%EFTflag = input_EFTCAMB%EFTflag
    temp_EFTCAMB%EFTCAMB_pert_turn_on = input_EFTCAMB%EFTCAMB_pert_turn_on
    temp_EFTCAMB%EFT_mass_stability = .false.
    temp_EFTCAMB%EFT_positivity_bounds = .false.
    temp_EFTCAMB%EFT_minkowski_limit = .false.
    s0 = 0._dl; s1 = 0._dl; s2 = 0._dl; sw = 0._dl; sxw = 0._dl
    do ind = 1, num_points
        a = a_min + real(ind - 1, dl)*(1._dl - a_min)/real(num_points - 1, dl)
        call eft_cache%initialize()
        call EFTStabilityComputation(a, temp_EFTCAMB, input_EFTCAMB%model, params_cache, eft_cache)
        rho_de = 3._dl*eft_cache%adotoa**2 - eft_cache%grhom_t
        p_de = -(2._dl*eft_cache%Hdot + eft_cache%adotoa**2) - eft_cache%gpresm_t
        if (.not. (rho_de > 0._dl)) cycle
        w = p_de/rho_de
        if (w /= w) cycle
        x = 1._dl - a
        s0 = s0 + 1._dl; s1 = s1 + x; s2 = s2 + x*x; sw = sw + w; sxw = sxw + x*w
    end do
    det = s0*s2 - s1*s1
    if (s0 < 2._dl .or. det <= 0._dl) return
    input_EFTCAMB%model%w0_background = (s2*sw - s1*sxw)/det
    input_EFTCAMB%model%wa_background = (s0*sxw - s1*sw)/det
    input_EFTCAMB%model%background_w0wa_computed = .true.

    end subroutine EFTCAMB_compute_effective_w0wa

    ! ----------------------------------------------------------------------------------------

    end module EFTCAMB_stability

! ----------------------------------------------------------------------------------------
