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

!> @file 06p1_abstract_EFTCAMB_full_map.f90
!! This file contains the abstract definition of all the places where EFTCAMB interacts
!! with CAMB in case of a full mapping model. All models implementing a model in which
!! the cosmological background is computed from the values of the EFT functions
!! should inherit from this class.


! ----------------------------------------------------------------------------------------
!> This module contains the abstract definition of all the places where EFTCAMB interacts
!! with CAMB in case of a full mapping model. All models implementing a model in which
!! the cosmological background is computed from the values of the EFT functions
!! should inherit from this class.

!> @author Bin Hu, Marco Raveri

    module EFTCAMB_abstract_model_full

    use precision
    use EFTCAMB_cache
    use EFTCAMB_abstract_model
    use EFTCAMB_mixed_algorithms
    use MassiveNu

    implicit none

    private

    public EFTCAMB_full_model

    ! ----------------------------------------------------------------------------------------
    !> This is the abstract type for EFTCAMB models.
    type, extends (EFTCAMB_model), abstract :: EFTCAMB_full_model

    contains

    !< function that computes dtauda = 1/sqrt(a^2H^2).
    procedure :: compute_dtauda => EFTCAMBFullModelComputeDtauda
    procedure :: compute_adotoa => EFTCAMBFullModelComputeAdotoa        !< subroutine that computes adotoa = H.
    !< subroutine that computes the two derivatives wrt conformal time of H.
    procedure :: compute_H_derivs => EFTCAMBFullModelComputeHubbleDer
    !< subroutine that computes additional derivatives of the EFT functions, as needed by the stability calculation.
    procedure :: compute_additional_derivs => EFTCAMBFullModelComputeAdditionalDer

    end type EFTCAMB_full_model

    contains

    ! ---------------------------------------------------------------------------------------------
    !> Function that computes dtauda = 1/sqrt(a^2H^2).
    function EFTCAMBFullModelComputeDtauda(self, a, eft_par_cache, eft_cache)

    implicit none

    class(EFTCAMB_full_model) :: self          !< the base class
    real(dl), intent(in)      :: a             !< the input scale factor.
    !< the EFTCAMB parameter cache that contains all the physical parameters.
    type(TEFTCAMB_parameter_cache), intent(inout) :: eft_par_cache
    !< the EFTCAMB timestep cache that contains all the physical values.
    type(TEFTCAMB_timestep_cache), intent(inout) :: eft_cache

    real(dl) :: EFTCAMBFullModelComputeDtauda                     !< the output dtauda

    real(dl) :: a2, temp

    a2 = a*a
    temp = 1.0_dl/(1.0_dl + eft_cache%EFTOmegaV + a*eft_cache%EFTOmegaP)* &
        (eft_cache%grhoa2 + 2._dl*eft_cache%EFTc*a2 - eft_cache%EFTLambda*a2)/3._dl

    EFTCAMBFullModelComputeDtauda = 1._dl/sqrt(temp)

    end function EFTCAMBFullModelComputeDtauda

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes adotoa = H.
    subroutine EFTCAMBFullModelComputeAdotoa(self, a, eft_par_cache, eft_cache)

    implicit none

    class(EFTCAMB_full_model) :: self          !< the base class
    real(dl), intent(in)      :: a             !< the input scale factor
    !< the EFTCAMB parameter cache that contains all the physical parameters.
    type(TEFTCAMB_parameter_cache), intent(inout) :: eft_par_cache
    !< the EFTCAMB timestep cache that contains all the physical values.
    type(TEFTCAMB_timestep_cache), intent(inout) :: eft_cache

    real(dl) :: temp

    temp = 1.0_dl/(1.0_dl + eft_cache%EFTOmegaV + a*eft_cache%EFTOmegaP)* &
        (eft_cache%grhom_t + 2.0_dl*eft_cache%EFTc - eft_cache%EFTLambda)/3.0_dl
    eft_cache%adotoa = sqrt(temp)

    end subroutine EFTCAMBFullModelComputeAdotoa

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes the two derivatives wrt conformal time of H.
    subroutine EFTCAMBFullModelComputeHubbleDer(self, a, eft_par_cache, eft_cache)

    implicit none

    class(EFTCAMB_full_model) :: self          !< the base class
    real(dl), intent(in)      :: a             !< the input scale factor
    !< the EFTCAMB parameter cache that contains all the physical parameters.
    type(TEFTCAMB_parameter_cache), intent(inout) :: eft_par_cache
    !< the EFTCAMB timestep cache that contains all the physical values.
    type(TEFTCAMB_timestep_cache), intent(inout) :: eft_cache

    real(dl) :: a2, denom, one_plus_omega, adotoa2, denomdot

    a2 = a*a
    adotoa2 = eft_cache%adotoa**2._dl

    one_plus_omega = 1.0_dl + eft_cache%EFTOmegaV

    eft_cache%Hdot = 1.0_dl/(one_plus_omega + 0.5_dl*a*eft_cache%EFTOmegaP)* &
        (-0.5_dl*(one_plus_omega + 2.0_dl*a*eft_cache%EFTOmegaP + a2*eft_cache%EFTOmegaPP)*eft_cache%adotoa**2 - &
        0.5_dl*(eft_cache%gpresm_t) - 0.5_dl*eft_cache%EFTLambda)

    eft_cache%Hdotdot = 1.0_dl/(one_plus_omega + 0.5_dl*a*eft_cache%EFTOmegaP)* &
        (-0.5_dl*a*eft_cache%adotoa**3* &
        (3.0_dl*eft_cache%EFTOmegaP + 4.0_dl*a*eft_cache%EFTOmegaPP + a2*eft_cache%EFTOmegaPPP) - &
        eft_cache%adotoa*eft_cache%Hdot* &
        (one_plus_omega + 3.5_dl*a*eft_cache%EFTOmegaP + 1.5_dl*a2*eft_cache%EFTOmegaPP) - &
        0.5_dl*(eft_cache%gpresdotm_t + 2.0_dl*eft_cache%adotoa*eft_cache%gpresm_t + eft_cache%EFTLambdadot + &
        2.0_dl*eft_cache%adotoa*eft_cache%EFTLambda))

    ! Hdotdotdot is the exact conformal time derivative of the Hdotdot expression above:
    ! Hdotdotdot = (d/dtau[numerator] - Hdotdot*d(denom)/dtau)/denom, with d/dtau = a H d/da on Omega and its
    ! derivatives and dX/dtau = Xdot + 2 H X for the a^2-weighted pressures and Lambda.
    ! It needs EFTLambdadotdot and EFTOmegaPPPP: full mapping models that do not compute them analytically get
    ! them from compute_additional_derivs, after which compute_H_derivs has to be called again.
    denom = one_plus_omega + 0.5_dl*a*eft_cache%EFTOmegaP
    denomdot = a*eft_cache%adotoa*(1.5_dl*eft_cache%EFTOmegaP + 0.5_dl*a*eft_cache%EFTOmegaPP)

    eft_cache%Hdotdotdot = 1._dl/denom*(-denomdot*eft_cache%Hdotdot &
        - a*adotoa2*(0.5_dl*adotoa2 + 1.5_dl*eft_cache%Hdot)* &
        (3._dl*eft_cache%EFTOmegaP + 4._dl*a*eft_cache%EFTOmegaPP + a2*eft_cache%EFTOmegaPPP) &
        - 0.5_dl*a2*adotoa2**2*(7._dl*eft_cache%EFTOmegaPP + 6._dl*a*eft_cache%EFTOmegaPPP + &
        a2*eft_cache%EFTOmegaPPPP) &
        - (eft_cache%Hdot**2 + eft_cache%adotoa*eft_cache%Hdotdot)* &
        (one_plus_omega + 3.5_dl*a*eft_cache%EFTOmegaP + 1.5_dl*a2*eft_cache%EFTOmegaPP) &
        - 0.5_dl*a*adotoa2*eft_cache%Hdot* &
        (9._dl*eft_cache%EFTOmegaP + 13._dl*a*eft_cache%EFTOmegaPP + 3._dl*a2*eft_cache%EFTOmegaPPP) &
        - 0.5_dl*(eft_cache%gpresdotdotm_t + 4._dl*eft_cache%adotoa*eft_cache%gpresdotm_t + &
        2._dl*eft_cache%Hdot*eft_cache%gpresm_t + 4._dl*adotoa2*eft_cache%gpresm_t) &
        - 0.5_dl*(eft_cache%EFTLambdadotdot + 4._dl*eft_cache%adotoa*eft_cache%EFTLambdadot + &
        2._dl*eft_cache%Hdot*eft_cache%EFTLambda + 4._dl*adotoa2*eft_cache%EFTLambda))

    end subroutine EFTCAMBFullModelComputeHubbleDer

! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes additional derivatives of the EFT functions, as needed by the
    !> stability calculation (mass eigenvalues, positivity bounds, Minkowski limit).
    !! The derivatives are finite differences (dfridr) of the lower derivatives evaluated, at the shifted scale
    !! factor, on a separate cache filled exactly as in the stability calculation. Values that the model already
    !! provides analytically (non-zero in the input cache) are kept. Derivatives of quantities that the model does
    !! not provide (e.g. EFTGamma3PPPP when EFTGamma3PPP is not computed) stay zero.
    subroutine EFTCAMBFullModelComputeAdditionalDer(self, a, eft_par_cache, eft_cache)

    implicit none

    class(EFTCAMB_full_model) :: self          !< the base class
    real(dl), intent(in)      :: a             !< the input scale factor.
    !< the EFTCAMB parameter cache that contains all the physical parameters.
    type(TEFTCAMB_parameter_cache), intent(inout) :: eft_par_cache
    !< the EFTCAMB timestep cache that contains all the physical values.
    type(TEFTCAMB_timestep_cache), intent(inout) :: eft_cache

    integer, parameter :: num_fields = 9, max_memo = 512
    real(dl) :: err, memo_a(max_memo), memo_v(num_fields, max_memo)
    integer  :: n_memo

    n_memo = 0

    if (eft_cache%EFTGamma1PP == 0._dl) eft_cache%EFTGamma1PP = dfridr(helper_1, a, 0.03_dl*a, err)
    if (eft_cache%EFTGamma2PP == 0._dl) eft_cache%EFTGamma2PP = dfridr(helper_2, a, 0.03_dl*a, err)
    if (eft_cache%EFTGamma3PP == 0._dl) eft_cache%EFTGamma3PP = dfridr(helper_3, a, 0.03_dl*a, err)
    if (eft_cache%EFTGamma2PPP == 0._dl) eft_cache%EFTGamma2PPP = dfridr(helper_4, a, 0.03_dl*a, err)
    if (eft_cache%EFTGamma3PPP == 0._dl) eft_cache%EFTGamma3PPP = dfridr(helper_5, a, 0.03_dl*a, err)
    if (eft_cache%EFTGamma3PPPP == 0._dl) eft_cache%EFTGamma3PPPP = dfridr(helper_6, a, 0.03_dl*a, err)
    if (eft_cache%EFTOmegaPPPP == 0._dl) eft_cache%EFTOmegaPPPP = dfridr(helper_7, a, 0.03_dl*a, err)
    ! the dotted EFT functions follow Xdot = dX/dtau - 2 H X (a^2-weighted), so Xdotdot = a H dXdot/da - 2 H Xdot:
    if (eft_cache%EFTcdotdot == 0._dl) eft_cache%EFTcdotdot = a*eft_cache%adotoa*dfridr(helper_8, a, 0.03_dl*a, err) &
        - 2._dl*eft_cache%adotoa*eft_cache%EFTcdot
    if (eft_cache%EFTLambdadotdot == 0._dl) eft_cache%EFTLambdadotdot = &
        a*eft_cache%adotoa*dfridr(helper_9, a, 0.03_dl*a, err) - 2._dl*eft_cache%adotoa*eft_cache%EFTLambdadot

    contains

    !> returns the i-th differentiated quantity at temp_a (memoized: all dfridr calls use the same points).
    function helper_value(temp_a, i)

    real(dl), intent(in) :: temp_a
    integer, intent(in)  :: i
    real(dl) :: helper_value

    type(TEFTCAMB_timestep_cache)  :: temp_cache
    type(TEFTCAMB_parameter_cache) :: temp_par_cache
    real(dl) :: a2, adotoa, EFT_grhonu, EFT_gpinu, EFT_gpinudot, grhormass_t
    integer  :: nu_i, j

    do j = 1, n_memo
        if (memo_a(j) == temp_a) then
            helper_value = memo_v(i, j)
            return
        end if
    end do

    temp_par_cache = eft_par_cache
    call temp_cache%initialize()
    a2 = temp_a*temp_a
    ! background densities of the different species at temp_a:
    temp_cache%a = temp_a
    temp_cache%grhob_t = temp_par_cache%grhob/temp_a
    temp_cache%grhoc_t = temp_par_cache%grhoc/temp_a
    temp_cache%grhor_t = temp_par_cache%grhornomass/a2
    temp_cache%grhog_t = temp_par_cache%grhog/a2
    temp_cache%grhonu_tot = 0._dl
    temp_cache%gpinu_tot = 0._dl
    if (temp_par_cache%Num_Nu_Massive /= 0) then
        do nu_i = 1, temp_par_cache%Nu_mass_eigenstates
            grhormass_t = temp_par_cache%grhormass(nu_i)/a2
            call ThermalNuBackground%rho_P(temp_a*temp_par_cache%nu_masses(nu_i), EFT_grhonu, EFT_gpinu)
            temp_cache%grhonu_tot = temp_cache%grhonu_tot + grhormass_t*EFT_grhonu
            temp_cache%gpinu_tot = temp_cache%gpinu_tot + grhormass_t*EFT_gpinu
        end do
    end if
    temp_cache%grhom_t = temp_cache%grhob_t + temp_cache%grhoc_t + temp_cache%grhor_t + temp_cache%grhog_t + &
        temp_cache%grhonu_tot
    temp_cache%gpresm_t = (temp_cache%grhor_t + temp_cache%grhog_t)/3._dl + temp_cache%gpinu_tot
    temp_cache%grhoa2 = temp_cache%grhom_t*a2
    ! full mapping: background EFT functions first, then the expansion history:
    call self%compute_background_EFT_functions(temp_a, temp_par_cache, temp_cache)
    call self%compute_adotoa(temp_a, temp_par_cache, temp_cache)
    adotoa = temp_cache%adotoa
    ! first derivatives of the massive neutrino quantities:
    temp_cache%grhonudot_tot = 0._dl
    temp_cache%gpinudot_tot = 0._dl
    if (temp_par_cache%Num_Nu_Massive /= 0) then
        do nu_i = 1, temp_par_cache%Nu_mass_eigenstates
            grhormass_t = temp_par_cache%grhormass(nu_i)/a2
            call ThermalNuBackground%rho_P(temp_a*temp_par_cache%nu_masses(nu_i), EFT_grhonu, EFT_gpinu)
            temp_cache%grhonudot_tot = temp_cache%grhonudot_tot + &
                grhormass_t*(ThermalNuBackground%drho(temp_a*temp_par_cache%nu_masses(nu_i), adotoa) - &
                4._dl*adotoa*EFT_grhonu)
            temp_cache%gpinudot_tot = temp_cache%gpinudot_tot + &
                grhormass_t*(ThermalNuBackground%pidot(temp_a*temp_par_cache%nu_masses(nu_i), adotoa, EFT_gpinu) - &
                4._dl*adotoa*EFT_gpinu)
        end do
    end if
    temp_cache%gpresdotm_t = -4._dl*adotoa*(temp_cache%grhog_t + temp_cache%grhor_t)/3._dl + temp_cache%gpinudot_tot
    call self%compute_H_derivs(temp_a, temp_par_cache, temp_cache)
    ! second derivative of the massive neutrino pressure (needs Hdot):
    temp_cache%gpinudotdot_tot = 0._dl
    if (temp_par_cache%Num_Nu_Massive /= 0) then
        do nu_i = 1, temp_par_cache%Nu_mass_eigenstates
            grhormass_t = temp_par_cache%grhormass(nu_i)/a2
            call ThermalNuBackground%rho_P(temp_a*temp_par_cache%nu_masses(nu_i), EFT_grhonu, EFT_gpinu)
            EFT_gpinudot = ThermalNuBackground%pidot(temp_a*temp_par_cache%nu_masses(nu_i), adotoa, EFT_gpinu)
            temp_cache%gpinudotdot_tot = temp_cache%gpinudotdot_tot - &
                4._dl*adotoa*grhormass_t*(EFT_gpinudot - 4._dl*adotoa*EFT_gpinu) + &
                grhormass_t*(ThermalNuBackground%pidotdot(temp_a*temp_par_cache%nu_masses(nu_i), adotoa, &
                temp_cache%Hdot, EFT_gpinu, EFT_gpinudot) - 4._dl*temp_cache%Hdot*EFT_gpinu - &
                4._dl*adotoa*EFT_gpinudot)
        end do
    end if
    temp_cache%gpresdotdotm_t = -(4._dl/3._dl)*temp_cache%Hdot*(temp_cache%grhog_t + temp_cache%grhor_t) + &
        (16._dl/3._dl)*(temp_cache%grhog_t + temp_cache%grhor_t)*adotoa**2 + temp_cache%gpinudotdot_tot
    call self%compute_H_derivs(temp_a, temp_par_cache, temp_cache)
    call self%compute_rhoQPQ(temp_a, temp_par_cache, temp_cache)
    call self%compute_secondorder_EFT_functions(temp_a, temp_par_cache, temp_cache)

    if (n_memo < max_memo) then
        n_memo = n_memo + 1
        j = n_memo
    else
        j = max_memo
    end if
    memo_a(j) = temp_a
    memo_v(:, j) = [temp_cache%EFTGamma1P, temp_cache%EFTGamma2P, temp_cache%EFTGamma3P, &
        temp_cache%EFTGamma2PP, temp_cache%EFTGamma3PP, temp_cache%EFTGamma3PPP, temp_cache%EFTOmegaPPP, &
        temp_cache%EFTcdot, temp_cache%EFTLambdadot]
    helper_value = memo_v(i, j)

    end function helper_value

    function helper_1(temp_a)
    real(dl) :: temp_a, helper_1
    helper_1 = helper_value(temp_a, 1)
    end function helper_1

    function helper_2(temp_a)
    real(dl) :: temp_a, helper_2
    helper_2 = helper_value(temp_a, 2)
    end function helper_2

    function helper_3(temp_a)
    real(dl) :: temp_a, helper_3
    helper_3 = helper_value(temp_a, 3)
    end function helper_3

    function helper_4(temp_a)
    real(dl) :: temp_a, helper_4
    helper_4 = helper_value(temp_a, 4)
    end function helper_4

    function helper_5(temp_a)
    real(dl) :: temp_a, helper_5
    helper_5 = helper_value(temp_a, 5)
    end function helper_5

    function helper_6(temp_a)
    real(dl) :: temp_a, helper_6
    helper_6 = helper_value(temp_a, 6)
    end function helper_6

    function helper_7(temp_a)
    real(dl) :: temp_a, helper_7
    helper_7 = helper_value(temp_a, 7)
    end function helper_7

    function helper_8(temp_a)
    real(dl) :: temp_a, helper_8
    helper_8 = helper_value(temp_a, 8)
    end function helper_8

    function helper_9(temp_a)
    real(dl) :: temp_a, helper_9
    helper_9 = helper_value(temp_a, 9)
    end function helper_9

    end subroutine EFTCAMBFullModelComputeAdditionalDer

    ! ---------------------------------------------------------------------------------------------

    end module EFTCAMB_abstract_model_full

! ----------------------------------------------------------------------------------------
