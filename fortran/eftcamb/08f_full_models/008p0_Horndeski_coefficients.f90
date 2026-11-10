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

!> @file 008p0_Horndeski_coefficients.f90
!! Coefficients of the Horndeski background equations: energy constraint, EFT alphas, background and designer
!! equations and their conformal-time derivatives. Generated from the symbolic derivation (arXiv:1105.5723) as C
!! code, then translated to Fortran keeping the operation order: there is no generator any more, edit by hand.

!> @author Gen Ye

    module EFTCAMB_Horndeski_coefficients

    use precision, only: dl
    use iso_c_binding, only: c_double

    implicit none

    private

    public horndeski_hubble_coefficients, horndeski_eft_alphas, horndeski_backeq_coefficients, &
        horndeski_backeq_coefficients_p, horndeski_designereq_coefficients, horndeski_designereq_coefficients_p

    contains

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes the energy constraint E3 H^3 + E2 H^2 + E1 H + E0 = 0
    subroutine horndeski_hubble_coefficients(a, dphi, rhotot, G2, G2p, G2x, G2pp, G2px, G2xx, G2ppp, G2ppx, G2pxx, &
        G2xxx, G3p, G3x, G3pp, G3px, G3xx, G3ppp, G3ppx, G3pxx, G3xxx, G4, G4p, G4x, G4pp, G4px, G4xx, G4ppp, G4ppx, &
        G4pxx, G4xxx, G4pppp, G4pppx, G4ppxx, G4pxxx, G4xxxx, G5, G5p, G5x, G5pp, G5px, G5xx, G5ppp, G5ppx, G5pxx, &
        G5xxx, G5pppp, G5pppx, G5ppxx, G5pxxx, G5xxxx, E0, E1, E2, E3) &
        bind(c, name='horndeski_hubble_coefficients')

    implicit none

    real(c_double), value :: a      !<- the scale factor
    real(c_double), value :: dphi   !<- conformal time derivative of the scalar field
    real(c_double), value :: rhotot !<- \f$ \rho_{matter} + p_{matter} \f$
    real(c_double), value :: G2     !<- Horndeski Lagrangian function \f$ G2 \f$
    real(c_double), value :: G2p    !<- Horndeski Lagrangian function \f$ G2_\phi \f$
    real(c_double), value :: G2x    !<- Horndeski Lagrangian function \f$ G2_X \f$
    real(c_double), value :: G2pp
    real(c_double), value :: G2px
    real(c_double), value :: G2xx
    real(c_double), value :: G2ppp
    real(c_double), value :: G2ppx
    real(c_double), value :: G2pxx
    real(c_double), value :: G2xxx
    real(c_double), value :: G3p
    real(c_double), value :: G3x
    real(c_double), value :: G3pp
    real(c_double), value :: G3px
    real(c_double), value :: G3xx
    real(c_double), value :: G3ppp
    real(c_double), value :: G3ppx
    real(c_double), value :: G3pxx
    real(c_double), value :: G3xxx
    real(c_double), value :: G4
    real(c_double), value :: G4p
    real(c_double), value :: G4x
    real(c_double), value :: G4pp
    real(c_double), value :: G4px
    real(c_double), value :: G4xx
    real(c_double), value :: G4ppp
    real(c_double), value :: G4ppx
    real(c_double), value :: G4pxx
    real(c_double), value :: G4xxx
    real(c_double), value :: G4pppp
    real(c_double), value :: G4pppx
    real(c_double), value :: G4ppxx
    real(c_double), value :: G4pxxx
    real(c_double), value :: G4xxxx
    real(c_double), value :: G5
    real(c_double), value :: G5p
    real(c_double), value :: G5x
    real(c_double), value :: G5pp
    real(c_double), value :: G5px
    real(c_double), value :: G5xx
    real(c_double), value :: G5ppp
    real(c_double), value :: G5ppx
    real(c_double), value :: G5pxx
    real(c_double), value :: G5xxx
    real(c_double), value :: G5pppp
    real(c_double), value :: G5pppx
    real(c_double), value :: G5ppxx
    real(c_double), value :: G5pxxx
    real(c_double), value :: G5xxxx
    real(c_double) :: E0
    real(c_double) :: E1
    real(c_double) :: E2
    real(c_double) :: E3

    real(c_double) :: X

    X = 0.5_dl * dphi * dphi / a / a
    E0 = G2 - rhotot + 2._dl * ((-G2x) + G3p) * X
    E1 = (6._dl * dphi * (G4p - G3x * X + 2._dl * G4px * X)) / a
    E2 = 6._dl * (G4 + X * ((-4._dl) * G4x + 3._dl * G5p - 4._dl * G4xx * X + 2._dl * G5px * X))
    E3 = ((-2._dl) * dphi * X * (5._dl * G5x + 2._dl * G5xx * X)) / a

    end subroutine horndeski_hubble_coefficients

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes the EFT functions alpha_i and a few derivatives
    subroutine horndeski_eft_alphas(a, H, dH, ddH, dphi, ddphi, dddphi, G2, G2p, G2x, G2pp, G2px, G2xx, G2ppp, G2ppx, &
        G2pxx, G2xxx, G3p, G3x, G3pp, G3px, G3xx, G3ppp, G3ppx, G3pxx, G3xxx, G4, G4p, G4x, G4pp, G4px, G4xx, G4ppp, &
        G4ppx, G4pxx, G4xxx, G4pppp, G4pppx, G4ppxx, G4pxxx, G4xxxx, G5, G5p, G5x, G5pp, G5px, G5xx, G5ppp, G5ppx, &
        G5pxx, G5xxx, G5pppp, G5pppx, G5ppxx, G5pxxx, G5xxxx, alphaM, dalphaM, alphaB, dalphaB, alphaK, dalphaK, &
        alphaT, dalphaT, M, dM, ddM) &
        bind(c, name='horndeski_eft_alphas')

    implicit none

    real(c_double), value :: a      !<- the scale factor
    real(c_double), value :: H      !<- the physcial Hubble parameter
    real(c_double), value :: dH     !<- conformal time derivative of the physical Hubble
    real(c_double), value :: ddH    !<- second conformal time derivative of the physical Hubble
    real(c_double), value :: dphi   !<- conformal time derivative of the scalar field
    real(c_double), value :: ddphi  !<- double conformal time derivative of the scalar field
    real(c_double), value :: dddphi !<- third conformal time derivative of the scalar field
    real(c_double), value :: G2     !<- Horndeski Lagrangian function \f$ G2 \f$
    real(c_double), value :: G2p    !<- Horndeski Lagrangian function \f$ G2_\phi \f$
    real(c_double), value :: G2x    !<- Horndeski Lagrangian function \f$ G2_X \f$
    real(c_double), value :: G2pp
    real(c_double), value :: G2px
    real(c_double), value :: G2xx
    real(c_double), value :: G2ppp
    real(c_double), value :: G2ppx
    real(c_double), value :: G2pxx
    real(c_double), value :: G2xxx
    real(c_double), value :: G3p
    real(c_double), value :: G3x
    real(c_double), value :: G3pp
    real(c_double), value :: G3px
    real(c_double), value :: G3xx
    real(c_double), value :: G3ppp
    real(c_double), value :: G3ppx
    real(c_double), value :: G3pxx
    real(c_double), value :: G3xxx
    real(c_double), value :: G4
    real(c_double), value :: G4p
    real(c_double), value :: G4x
    real(c_double), value :: G4pp
    real(c_double), value :: G4px
    real(c_double), value :: G4xx
    real(c_double), value :: G4ppp
    real(c_double), value :: G4ppx
    real(c_double), value :: G4pxx
    real(c_double), value :: G4xxx
    real(c_double), value :: G4pppp
    real(c_double), value :: G4pppx
    real(c_double), value :: G4ppxx
    real(c_double), value :: G4pxxx
    real(c_double), value :: G4xxxx
    real(c_double), value :: G5
    real(c_double), value :: G5p
    real(c_double), value :: G5x
    real(c_double), value :: G5pp
    real(c_double), value :: G5px
    real(c_double), value :: G5xx
    real(c_double), value :: G5ppp
    real(c_double), value :: G5ppx
    real(c_double), value :: G5pxx
    real(c_double), value :: G5xxx
    real(c_double), value :: G5pppp
    real(c_double), value :: G5pppx
    real(c_double), value :: G5ppxx
    real(c_double), value :: G5pxxx
    real(c_double), value :: G5xxxx
    real(c_double) :: alphaM          !<- The EFT function alpha_M
    real(c_double) :: dalphaM         !<- conformal time of the EFT function alpha_M
    real(c_double) :: alphaB          !<- The EFT function alpha_B
    real(c_double) :: dalphaB         !<- conformal time of the EFT function alpha_B
    real(c_double) :: alphaK          !<- The EFT function alpha_K
    real(c_double) :: dalphaK         !<- conformal time of the EFT function alpha_K
    real(c_double) :: alphaT          !<- The EFT function alpha_T
    real(c_double) :: dalphaT         !<- conformal time of the EFT function alpha_T
    real(c_double) :: M               !<- The effective planck mass squared
    real(c_double) :: dM              !<- conformal time of the effective planck mass squared
    real(c_double) :: ddM             !<- second conformal time of the effective planck mass squared

    real(c_double) :: X, M2, dM2, ddM2, aM, daM, aB, daB, aK, daK, aT, daT

    X = 0.5_dl * dphi * dphi / a / a
    M2 = 2._dl * (G4 - 2._dl * G4x * X + G5p * X - (dphi * G5x * H * X) / a)
    dM2 = (dphi * (2._dl * a**5 * G4p + 2._dl * a**4 * dphi * (G4x - G5p) * H - ddphi * dphi**3 * G5xx * H - a**2 * &
        dphi * (dH * dphi * G5x - 2._dl * dphi**2 * (G4xx - G5px) * H + 3._dl * ddphi * G5x * H) + a * (ddphi * &
        dphi**2 * ((-2._dl) * G4xx + G5px) + dphi**4 * G5xx * H**2) + a**3 * (2._dl * ddphi * ((-G4x) + G5p) + &
        dphi**2 * ((-2._dl) * G4px + G5pp + 3._dl * G5x * H**2)))) / a**5
    ddM2 = ((-(ddphi**2 * dphi**5 * G5xxx * H)) + a * ddphi * dphi**4 * (ddphi * ((-2._dl) * G4xxx + G5pxx) + 2._dl &
        * dphi**2 * G5xxx * H**2) - a**2 * dphi**3 * (7._dl * ddphi**2 * G5xx * H + dddphi * dphi * G5xx * H + &
        dphi**4 * G5xxx * H**3 + 2._dl * ddphi * dphi * (dH * G5xx + 2._dl * dphi * ((-G4xxx) + G5pxx) * H)) + 2._dl &
        * a**7 * (ddphi * G4p + dphi**2 * (G4pp + ((-G4x) + G5p) * H**2)) + 2._dl * a**6 * dphi * (dH * dphi * (G4x &
        - G5p) + 4._dl * ddphi * (G4x - G5p) * H + dphi**2 * H * (2._dl * G4px - 2._dl * G5pp - 3._dl * G5x * H**2)) &
        + a**5 * (2._dl * ddphi**2 * ((-G4x) + G5p) + 2._dl * dddphi * dphi * ((-G4x) + G5p) + 9._dl * dH * dphi**3 &
        * G5x * H + dphi**4 * ((-2._dl) * G4ppx + G5ppp + ((-8._dl) * G4xx + 11._dl * G5px) * H**2) + ddphi * &
        dphi**2 * ((-6._dl) * G4px + 5._dl * G5pp + 18._dl * G5x * H**2)) + a**3 * dphi**2 * (dddphi * dphi * &
        ((-2._dl) * G4xx + G5px) + ddphi**2 * ((-8._dl) * G4xx + 5._dl * G5px) + dphi**3 * H * (3._dl * dH * G5xx - &
        2._dl * dphi * G4xxx * H + 3._dl * dphi * G5pxx * H) + 2._dl * ddphi * dphi**2 * ((-2._dl) * G4pxx + G5ppx + &
        8._dl * G5xx * H**2)) - a**4 * dphi * (ddH * dphi**2 * G5x + dH * (dphi**3 * ((-2._dl) * G4xx + 3._dl * &
        G5px) + 6._dl * ddphi * dphi * G5x) + H * (ddphi * dphi**2 * ((-20._dl) * G4xx + 19._dl * G5px) + 6._dl * &
        ddphi**2 * G5x + 3._dl * dddphi * dphi * G5x + dphi**4 * ((-4._dl) * G4pxx + 3._dl * G5ppx + 7._dl * G5xx * &
        H**2)))) / a**7
    aM = (dphi * (2._dl * a**5 * G4p + 2._dl * a**4 * dphi * (G4x - G5p) * H - ddphi * dphi**3 * G5xx * H - a**2 * &
        dphi * (dH * dphi * G5x - 2._dl * dphi**2 * (G4xx - G5px) * H + 3._dl * ddphi * G5x * H) + a * (ddphi * &
        dphi**2 * ((-2._dl) * G4xx + G5px) + dphi**4 * G5xx * H**2) + a**3 * (2._dl * ddphi * ((-G4x) + G5p) + &
        dphi**2 * ((-2._dl) * G4px + G5pp + 3._dl * G5x * H**2)))) / (a**3 * H * (2._dl * a**3 * G4 + a * dphi**2 * &
        ((-2._dl) * G4x + G5p) - dphi**3 * G5x * H))
    daM = ((-4._dl) * a**11 * dphi * G4 * G4p * H**2 + ddphi**2 * dphi**8 * ((-G5xx**2) + G5x * G5xxx) * H**3 + a * &
        ddphi * dphi**7 * H**2 * (ddphi * (2._dl * G4xxx * G5x - G5pxx * G5x - 4._dl * G4xx * G5xx + 2._dl * G5px * &
        G5xx + 2._dl * G4x * G5xxx - G5p * G5xxx) + 2._dl * dphi**2 * (G5xx**2 - G5x * G5xxx) * H**2) + a**2 * &
        dphi**6 * H * (dddphi * dphi * G5x * G5xx * H**2 + dphi**4 * ((-G5xx**2) + G5x * G5xxx) * H**4 - ddphi * &
        dphi * H * (dH * G5x * G5xx + 2._dl * dphi * (2._dl * G4xxx * G5x - 2._dl * G5pxx * G5x - 4._dl * G4xx * &
        G5xx + 3._dl * G5px * G5xx + 2._dl * G4x * G5xxx - G5p * G5xxx) * H) + ddphi**2 * ((-4._dl) * G4xx**2 + &
        4._dl * G4x * G4xxx - 2._dl * G4xxx * G5p + 4._dl * G4xx * G5px - G5px**2 - 2._dl * G4x * G5pxx + G5p * &
        G5pxx + G5x * G5xx * H**2)) - 4._dl * a**10 * (dH * dphi * G4 * G4p - ddphi * G4 * G4p * H + dphi**2 * H * &
        (G4p**2 - G4 * (G4pp + 2._dl * ((-G4x) + G5p) * H**2))) + 2._dl * a**9 * dphi * H**2 * (10._dl * ddphi * G4 &
        * (G4x - G5p) + dphi**2 * ((-2._dl) * G4p * G4x + 3._dl * G4p * G5p + G4 * (6._dl * G4px - 5._dl * G5pp - &
        9._dl * G5x * H**2))) + a**6 * dphi**2 * (dH * ((-2._dl) * ddphi * dphi * (2._dl * G4x**2 - 2._dl * G4 * &
        G4xx - 3._dl * G4x * G5p + G5p**2 + G4 * G5px) + dphi**3 * ((-4._dl) * G4px * G4x + 2._dl * G4px * G5p + &
        2._dl * G4x * G5pp - G5p * G5pp - 10._dl * G4x * G5x * H**2 + 3._dl * G5p * G5x * H**2 + 4._dl * G4 * G5xx * &
        H**2)) + H * ((-2._dl) * ddphi**2 * (8._dl * G4 * G4xx - G4x * G5p + G5p**2 - 5._dl * G4 * G5px) + 2._dl * &
        dddphi * dphi * (2._dl * G4x**2 - 2._dl * G4 * G4xx - 3._dl * G4x * G5p + G5p**2 + G4 * G5px) + ddphi * &
        dphi**2 * ((-8._dl) * G4 * G4pxx + 4._dl * G4px * G4x + 8._dl * G4p * G4xx + 2._dl * G4px * G5p - 6._dl * &
        G4x * G5pp + G5p * G5pp + 4._dl * G4 * G5ppx - 4._dl * G4p * G5px - 28._dl * G4x * G5x * H**2 + 7._dl * G5p &
        * G5x * H**2 + 34._dl * G4 * G5xx * H**2) - dphi**4 * (4._dl * G4px**2 - 4._dl * G4ppx * G4x + 2._dl * G4ppx &
        * G5p - 4._dl * G4px * G5pp + G5pp**2 + 2._dl * G4x * G5ppp - G5p * G5ppp - 12._dl * G4x * G4xx * H**2 + &
        4._dl * G4 * G4xxx * H**2 + 2._dl * G4xx * G5p * H**2 + 18._dl * G4x * G5px * H**2 - 5._dl * G5p * G5px * &
        H**2 - 6._dl * G4 * G5pxx * H**2 - 6._dl * G4px * G5x * H**2 + G5pp * G5x * H**2 + 4._dl * G4p * G5xx * &
        H**2))) - a**3 * dphi**5 * H * (dddphi * dphi * ((-2._dl) * G4xx * G5x + G5px * G5x - 2._dl * G4x * G5xx + &
        G5p * G5xx) * H + ddphi**2 * (4._dl * G4xx * G5x - G5px * G5x - 10._dl * G4x * G5xx + 3._dl * G5p * G5xx + &
        2._dl * G4 * G5xxx) * H + dphi**4 * ((-2._dl) * G4xxx * G5x + 3._dl * G5pxx * G5x + 4._dl * G4xx * G5xx - &
        4._dl * G5px * G5xx - 2._dl * G4x * G5xxx + G5p * G5xxx) * H**3 + ddphi * dphi * (dH * (6._dl * G4xx * G5x - &
        3._dl * G5px * G5x - 2._dl * G4x * G5xx + G5p * G5xx) + dphi * H * ((-8._dl) * G4xx**2 + 8._dl * G4x * G4xxx &
        - 4._dl * G4xxx * G5p + 12._dl * G4xx * G5px - 4._dl * G5px**2 - 8._dl * G4x * G5pxx + 4._dl * G5p * G5pxx - &
        4._dl * G4pxx * G5x + 2._dl * G5ppx * G5x + 4._dl * G4px * G5xx - 2._dl * G5pp * G5xx + 5._dl * G5x * G5xx * &
        H**2))) + 2._dl * a**8 * (dH * dphi * (2._dl * ddphi * G4 * (G4x - G5p) + dphi**2 * (2._dl * G4 * G4px + &
        2._dl * G4p * G4x - G4p * G5p - G4 * G5pp + 7._dl * G4 * G5x * H**2)) + H * (2._dl * ddphi**2 * G4 * ((-G4x) &
        + G5p) + 2._dl * dddphi * dphi * G4 * ((-G4x) + G5p) + ddphi * dphi**2 * ((-6._dl) * G4 * G4px + 2._dl * G4p &
        * G4x - 3._dl * G4p * G5p + 5._dl * G4 * G5pp + 21._dl * G4 * G5x * H**2) + dphi**4 * (4._dl * G4p * G4px - &
        2._dl * G4pp * G4x + G4pp * G5p - 2._dl * G4p * G5pp + 2._dl * G4x**2 * H**2 - 2._dl * G4x * G5p * H**2 - &
        5._dl * G4p * G5x * H**2 + G4 * ((-2._dl) * G4ppx + G5ppp + ((-10._dl) * G4xx + 13._dl * G5px) * H**2)))) - &
        a**7 * dphi * ((-2._dl) * dH**2 * dphi**2 * G4 * G5x + 2._dl * dH * dphi * (3._dl * ddphi * G4 * G5x + &
        dphi**2 * (G4 * G5px - 3._dl * G4p * G5x)) * H + H * (2._dl * ddH * dphi**2 * G4 * G5x + H * (12._dl * &
        ddphi**2 * G4 * G5x + 6._dl * dddphi * dphi * G4 * G5x + 2._dl * ddphi * dphi**2 * (6._dl * G4x**2 - 22._dl &
        * G4 * G4xx - 7._dl * G4x * G5p + G5p**2 + 20._dl * G4 * G5px - 5._dl * G4p * G5x) + dphi**4 * ((-8._dl) * &
        G4 * G4pxx + 4._dl * G4px * G4x + 8._dl * G4p * G4xx + 2._dl * G4px * G5p - 6._dl * G4x * G5pp + G5p * G5pp &
        + 6._dl * G4 * G5ppx - 8._dl * G4p * G5px + 2._dl * G4pp * G5x - 10._dl * G4x * G5x * H**2 + G5p * G5x * &
        H**2 + 16._dl * G4 * G5xx * H**2)))) - a**5 * dphi**3 * (dH**2 * dphi**2 * (2._dl * G4x - G5p) * G5x + dH * &
        dphi * H * ((-3._dl) * ddphi * G5p * G5x + 2._dl * ddphi * G4 * G5xx + dphi**2 * ((-2._dl) * G4x * G5px + &
        G5p * G5px + G5x * (6._dl * G4px - 3._dl * G5pp + G5x * H**2))) + H * (ddH * dphi**2 * ((-2._dl) * G4x + &
        G5p) * G5x + H * ((-2._dl) * ddphi**2 * (G4x * G5x + 2._dl * G5p * G5x - 7._dl * G4 * G5xx) + dddphi * dphi &
        * ((-8._dl) * G4x * G5x + 5._dl * G5p * G5x + 2._dl * G4 * G5xx) + ddphi * dphi**2 * ((-8._dl) * G4 * G4xxx &
        - 6._dl * G4xx * G5p + 28._dl * G4x * (G4xx - G5px) + 8._dl * G5p * G5px + 8._dl * G4 * G5pxx + 6._dl * G4px &
        * G5x - G5pp * G5x - 4._dl * G4p * G5xx + 3._dl * G5x**2 * H**2) + dphi**4 * (8._dl * G4pxx * G4x - 8._dl * &
        G4px * G4xx - 4._dl * G4pxx * G5p + 4._dl * G4xx * G5pp - 6._dl * G4x * G5ppx + 3._dl * G5p * G5ppx + 8._dl &
        * G4px * G5px - 4._dl * G5pp * G5px - 2._dl * G4ppx * G5x + G5ppp * G5x + 2._dl * G4xx * G5x * H**2 + G5px * &
        G5x * H**2 - 12._dl * G4x * G5xx * H**2 + 4._dl * G5p * G5xx * H**2 + 2._dl * G4 * G5xxx * H**2)))) + a**4 * &
        dphi**4 * (ddphi**2 * H * (8._dl * G4x * G4xx - 4._dl * G4 * G4xxx - 6._dl * G4x * G5px + G5p * G5px + 2._dl &
        * G4 * G5pxx - 3._dl * G5x**2 * H**2) + ddphi * dphi * ((-(dH * (4._dl * G4x * G4xx - 2._dl * G4xx * G5p - &
        2._dl * G4x * G5px + G5p * G5px + 3._dl * G5x**2 * H**2))) + dphi * H * (8._dl * G4pxx * G4x - 8._dl * G4px &
        * G4xx - 4._dl * G4pxx * G5p + 4._dl * G4xx * G5pp - 4._dl * G4x * G5ppx + 2._dl * G5p * G5ppx + 4._dl * &
        G4px * G5px - 2._dl * G5pp * G5px + 2._dl * G4xx * G5x * H**2 + 2._dl * G5px * G5x * H**2 - 26._dl * G4x * &
        G5xx * H**2 + 9._dl * G5p * G5xx * H**2 + 4._dl * G4 * G5xxx * H**2)) + dphi * H * (dddphi * (4._dl * G4x * &
        G4xx - 2._dl * G4xx * G5p - 2._dl * G4x * G5px + G5p * G5px + 3._dl * G5x**2 * H**2) + dphi * ((-2._dl) * &
        dH**2 * G5x**2 + dH * dphi * (4._dl * G4xx * G5x - 3._dl * G5px * G5x - 4._dl * G4x * G5xx + 2._dl * G5p * &
        G5xx) * H + H * (ddH * G5x**2 + dphi**2 * H * ((-4._dl) * G4xx**2 + 4._dl * G4x * G4xxx - 2._dl * G4xxx * &
        G5p + 8._dl * G4xx * G5px - 4._dl * G5px**2 - 6._dl * G4x * G5pxx + 3._dl * G5p * G5pxx - 4._dl * G4pxx * &
        G5x + 3._dl * G5ppx * G5x + 4._dl * G4px * G5xx - 2._dl * G5pp * G5xx + 2._dl * G5x * G5xx * H**2)))))) / &
        (a**5 * H**2 * (2._dl * a**3 * G4 + a * dphi**2 * ((-2._dl) * G4x + G5p) - dphi**3 * G5x * H)**2)
    aB = (dphi * ((-2._dl) * a**4 * G4p + 4._dl * a**3 * dphi * (G4x - G5p) * H + 2._dl * a * dphi**3 * (2._dl * &
        G4xx - G5px) * H + dphi**4 * G5xx * H**2 + a**2 * dphi**2 * (G3x - 2._dl * G4px + 3._dl * G5x * H**2))) / &
        (a**2 * H * (2._dl * a**3 * G4 + a * dphi**2 * ((-2._dl) * G4x + G5p) - dphi**3 * G5x * H))
    daB = (4._dl * a**10 * dphi * G4 * G4p * H**2 + ddphi * dphi**9 * (G5xx**2 - G5x * G5xxx) * H**4 + a * dphi**8 * &
        H**3 * (ddphi * ((-4._dl) * G4xxx * G5x + 2._dl * G5pxx * G5x + 6._dl * G4xx * G5xx - 3._dl * G5px * G5xx - &
        2._dl * G4x * G5xxx + G5p * G5xxx) + dphi**2 * ((-G5xx**2) + G5x * G5xxx) * H**2) - 2._dl * a**8 * dphi * &
        H**2 * (8._dl * ddphi * G4 * ((-G4x) + G5p) + dphi**2 * (3._dl * G3x * G4 - 12._dl * G4 * G4px + 4._dl * G4p &
        * G4x - 3._dl * G4p * G5p + 4._dl * G4 * G5pp + 9._dl * G4 * G5x * H**2)) - 2._dl * a**7 * dphi**2 * (dH * &
        dphi * (G3x * G4 - 2._dl * G4 * G4px + 2._dl * G4p * G4x - G4p * G5p - 3._dl * G4 * G5x * H**2) - ddphi * H &
        * (3._dl * G3x * G4 - 8._dl * G4 * G4px + G4p * G5p + 9._dl * G4 * G5x * H**2) + dphi**2 * H * ((-(G3px * &
        G4)) + G3x * G4p + 2._dl * G4 * G4ppx - 2._dl * G4pp * G4x + G4pp * G5p - G4p * G5pp - 4._dl * G4x**2 * H**2 &
        + 20._dl * G4 * G4xx * H**2 + 4._dl * G4x * G5p * H**2 - 15._dl * G4 * G5px * H**2 + G4p * G5x * H**2)) + &
        a**2 * (dphi**9 * (4._dl * G4xxx * G5x - 3._dl * G5pxx * G5x - 6._dl * G4xx * G5xx + 4._dl * G5px * G5xx + &
        2._dl * G4x * G5xxx - G5p * G5xxx) * H**4 + ddphi * dphi**7 * H**2 * (8._dl * G4xx**2 - 8._dl * G4x * G4xxx &
        + 4._dl * G4xxx * G5p - 8._dl * G4xx * G5px + 2._dl * G5px**2 + 4._dl * G4x * G5pxx - 2._dl * G5p * G5pxx - &
        G3xx * G5x + 2._dl * G4pxx * G5x + G3x * G5xx - 2._dl * G4px * G5xx - 2._dl * G5x * G5xx * H**2)) + a**4 * &
        dphi**5 * H * (2._dl * dH * dphi * (G3x - 2._dl * G4px) * G5x + 2._dl * ddphi * ((-12._dl) * G4x * G4xx + &
        4._dl * G4 * G4xxx + 2._dl * G4xx * G5p + 8._dl * G4x * G5px - 2._dl * G5p * G5px - 2._dl * G4 * G5pxx + &
        G4px * G5x - G4p * G5xx) * H + dphi**2 * H * ((-2._dl) * G3x * G4xx + 12._dl * G4px * G4xx + G3xx * (2._dl * &
        G4x - G5p) + 6._dl * G4pxx * ((-2._dl) * G4x + G5p) - 4._dl * G4xx * G5pp + 4._dl * G4x * G5ppx - 2._dl * &
        G5p * G5ppx + 2._dl * G3x * G5px - 8._dl * G4px * G5px + 2._dl * G5pp * G5px - G3px * G5x + 2._dl * G4ppx * &
        G5x + 2._dl * G4xx * G5x * H**2 - 3._dl * G5px * G5x * H**2 + 10._dl * G4x * G5xx * H**2 - 2._dl * G5p * &
        G5xx * H**2 - 2._dl * G4 * G5xxx * H**2)) + a**5 * dphi**4 * (dphi**2 * H * (2._dl * G3x * G4px - 4._dl * &
        G4px**2 - 2._dl * G3px * G4x + 4._dl * G4ppx * G4x + G3px * G5p - 2._dl * G4ppx * G5p - G3x * G5pp + 2._dl * &
        G4px * G5pp + 24._dl * G4x * G4xx * H**2 - 8._dl * G4 * G4xxx * H**2 - 4._dl * G4xx * G5p * H**2 - 18._dl * &
        G4x * G5px * H**2 + 3._dl * G5p * G5px * H**2 + 6._dl * G4 * G5pxx * H**2 + G5pp * G5x * H**2) + dH * dphi * &
        (G3x * (2._dl * G4x - G5p) + G4px * ((-4._dl) * G4x + 2._dl * G5p) - (2._dl * G4x * G5x + G5p * G5x - 2._dl &
        * G4 * G5xx) * H**2) + ddphi * H * (2._dl * G3xx * G4 - 4._dl * G3x * G4x + 12._dl * G4px * G4x - 4._dl * &
        G4p * G4xx + G3x * G5p - 4._dl * G4px * G5p + 2._dl * G4p * G5px - 8._dl * G4x * G5x * H**2 - G5p * G5x * &
        H**2 - 4._dl * G4 * (G4pxx - 4._dl * G5xx * H**2))) - a**6 * dphi**3 * H * (4._dl * dH * dphi * G4p * G5x + &
        4._dl * ddphi * (2._dl * G4x**2 - 10._dl * G4 * G4xx - 2._dl * G4x * G5p + 6._dl * G4 * G5px + G4p * G5x) * &
        H + dphi**2 * H * (2._dl * G3xx * G4 - 4._dl * G3x * G4x + 12._dl * G4px * G4x + 4._dl * G4p * G4xx + G3x * &
        G5p - 4._dl * G4x * G5pp - 2._dl * G4pp * G5x - 8._dl * G4x * G5x * H**2 - G5p * G5x * H**2 + 4._dl * G4 * &
        ((-3._dl) * G4pxx + G5ppx + 4._dl * G5xx * H**2))) + a**3 * dphi**6 * H * (ddphi * (4._dl * G4pxx * G4x + &
        2._dl * G3x * G4xx - 4._dl * G4px * G4xx - 2._dl * G4pxx * G5p + G3xx * ((-2._dl) * G4x + G5p) - G3x * G5px &
        + 2._dl * G4px * G5px - 2._dl * G4xx * G5x * H**2 + 3._dl * G5px * G5x * H**2 - 10._dl * G4x * G5xx * H**2 + &
        2._dl * G5p * G5xx * H**2 + 2._dl * G4 * G5xxx * H**2) + dphi * H * (dH * (4._dl * G4xx * G5x - 2._dl * G5px &
        * G5x - 2._dl * G4x * G5xx + G5p * G5xx) + dphi * H * ((-8._dl) * G4xx**2 + 8._dl * G4x * G4xxx - 4._dl * &
        G4xxx * G5p + 12._dl * G4xx * G5px - 4._dl * G5px**2 - 6._dl * G4x * G5pxx + 3._dl * G5p * G5pxx + G3xx * &
        G5x - 6._dl * G4pxx * G5x + 2._dl * G5ppx * G5x - G3x * G5xx + 4._dl * G4px * G5xx - G5pp * G5xx + 2._dl * &
        G5x * G5xx * H**2))) + 4._dl * a**9 * (dH * dphi * G4 * G4p - H * (ddphi * G4 * G4p + dphi**2 * ((-G4p**2) + &
        G4 * (G4pp + 4._dl * (G4x - G5p) * H**2))))) / (a**4 * H**2 * (2._dl * a**3 * G4 + a * dphi**2 * ((-2._dl) * &
        G4x + G5p) - dphi**3 * G5x * H)**2)
    aK = (dphi**2 * (3._dl * a * dphi**4 * (2._dl * G4xxx - G5pxx) * H**2 + dphi**5 * G5xxx * H**3 + a**3 * dphi**2 &
        * (G2xx - G3px + 3._dl * (8._dl * G4xx - 5._dl * G5px) * H**2) + 6._dl * a**4 * dphi * H * (G3x - 3._dl * &
        G4px + G5x * H**2) + a**2 * dphi**3 * H * (3._dl * G3xx - 6._dl * G4pxx + 7._dl * G5xx * H**2) + a**5 * (G2x &
        - 2._dl * (G3p + 3._dl * ((-G4x) + G5p) * H**2)))) / (a**4 * H**2 * (2._dl * a**3 * G4 + a * dphi**2 * &
        ((-2._dl) * G4x + G5p) - dphi**3 * G5x * H))
    daK = (dphi * ((-2._dl) * a**2 * dH * dphi * (2._dl * a**3 * G4 + a * dphi**2 * ((-2._dl) * G4x + G5p) - dphi**3 &
        * G5x * H) * (3._dl * a * dphi**4 * (2._dl * G4xxx - G5pxx) * H**2 + dphi**5 * G5xxx * H**3 + a**3 * dphi**2 &
        * (G2xx - G3px + 3._dl * (8._dl * G4xx - 5._dl * G5px) * H**2) + 6._dl * a**4 * dphi * H * (G3x - 3._dl * &
        G4px + G5x * H**2) + a**2 * dphi**3 * H * (3._dl * G3xx - 6._dl * G4pxx + 7._dl * G5xx * H**2) + a**5 * (G2x &
        - 2._dl * (G3p + 3._dl * ((-G4x) + G5p) * H**2))) - dphi**2 * H * (3._dl * a * dphi**4 * (2._dl * G4xxx - &
        G5pxx) * H**2 + dphi**5 * G5xxx * H**3 + a**3 * dphi**2 * (G2xx - G3px + 3._dl * (8._dl * G4xx - 5._dl * &
        G5px) * H**2) + 6._dl * a**4 * dphi * H * (G3x - 3._dl * G4px + G5x * H**2) + a**2 * dphi**3 * H * (3._dl * &
        G3xx - 6._dl * G4pxx + 7._dl * G5xx * H**2) + a**5 * (G2x - 2._dl * (G3p + 3._dl * ((-G4x) + G5p) * H**2))) &
        * (2._dl * a**5 * G4p + 2._dl * a**4 * dphi * (G4x - G5p) * H - ddphi * dphi**3 * G5xx * H - a**2 * dphi * &
        (dH * dphi * G5x - 2._dl * dphi**2 * (G4xx - G5px) * H + 3._dl * ddphi * G5x * H) + a * (ddphi * dphi**2 * &
        ((-2._dl) * G4xx + G5px) + dphi**4 * G5xx * H**2) + a**3 * (2._dl * ddphi * ((-G4x) + G5p) + dphi**2 * &
        ((-2._dl) * G4px + G5pp + 3._dl * G5x * H**2))) + H * ((-2._dl) * a**3 * G4 + a * dphi**2 * (2._dl * G4x - &
        G5p) + dphi**3 * G5x * H) * ((-(ddphi * dphi**7 * G5xxxx * H**3)) + a * dphi**6 * H**2 * ((-6._dl) * ddphi * &
        G4xxxx + 3._dl * ddphi * G5pxxx + dphi**2 * G5xxxx * H**2) + 2._dl * a**8 * dphi * H * (G2x - 2._dl * (G3p + &
        3._dl * ((-G4x) + G5p) * H**2)) + a**4 * dphi**3 * (dphi**2 * H * (G2xxx - 4._dl * G3pxx + 6._dl * G4ppxx + &
        60._dl * G4xxx * H**2 - 40._dl * G5pxx * H**2) + ddphi * H * ((-21._dl) * G3xx + 48._dl * G4pxx - 41._dl * &
        G5xx * H**2) - 3._dl * dH * dphi * (G3xx - 2._dl * G4pxx + 7._dl * G5xx * H**2)) - a**2 * dphi**5 * H * &
        (dphi * H * (3._dl * dH * G5xxx - 6._dl * dphi * G4xxxx * H + 4._dl * dphi * G5pxxx * H) + ddphi * (3._dl * &
        G3xxx - 6._dl * G4pxxx + 14._dl * G5xxx * H**2)) + a**6 * dphi * ((-18._dl) * ddphi * H * (G3x - 3._dl * &
        G4px + G5x * H**2) - 6._dl * dH * dphi * (G3x - 3._dl * G4px + 3._dl * G5x * H**2) + dphi**2 * H * (5._dl * &
        G2xx + 6._dl * ((-2._dl) * G3px + 3._dl * G4ppx + 17._dl * G4xx * H**2 - 12._dl * G5px * H**2))) + a**7 * &
        ((-12._dl) * dH * dphi * (G4x - G5p) * H - 2._dl * ddphi * (G2x - 2._dl * (G3p + 3._dl * ((-G4x) + G5p) * &
        H**2)) + dphi**2 * ((-G2px) + 2._dl * G3pp + 6._dl * H**2 * (3._dl * G3x - 10._dl * G4px + G5pp + 3._dl * &
        G5x * H**2))) + a**5 * dphi**2 * ((-6._dl) * dH * dphi * (8._dl * G4xx - 5._dl * G5px) * H + ddphi * &
        ((-5._dl) * G2xx + 6._dl * (G3px + ((-17._dl) * G4xx + 11._dl * G5px) * H**2)) + dphi**2 * ((-G2pxx) + G3ppx &
        + H**2 * (21._dl * G3xx - 72._dl * G4pxx + 15._dl * G5ppx + 41._dl * G5xx * H**2))) + a**3 * dphi**4 * &
        (ddphi * ((-G2xxx) + G3pxx - 60._dl * G4xxx * H**2 + 33._dl * G5pxx * H**2) + dphi * H * (6._dl * dH * &
        ((-2._dl) * G4xxx + G5pxx) + dphi * H * (3._dl * G3xxx - 12._dl * G4pxxx + 3._dl * G5ppxx + 14._dl * G5xxx * &
        H**2)))))) / (4._dl * a**12 * H**3 * (G4 - (dphi**2 * (2._dl * a * G4x - a * G5p + dphi * G5x * H)) / (2._dl &
        * a**3))**2)
    aT = (dphi**2 * (2._dl * a**2 * (G4x - G5p) - ddphi * G5x + 2._dl * a * dphi * G5x * H)) / (a * (2._dl * a**3 * &
        G4 + a * dphi**2 * ((-2._dl) * G4x + G5p) - dphi**3 * G5x * H))
    daT = (dphi * (ddphi**2 * dphi**4 * ((-2._dl) * G4xx * G5x + G5px * G5x + 2._dl * G4x * G5xx - G5p * G5xx) + &
        8._dl * a**7 * dphi * G4 * ((-G4x) + G5p) * H + 4._dl * a**5 * dphi * (dH * dphi * G4 * G5x + 5._dl * ddphi &
        * G4 * G5x * H + dphi**2 * (G4x**2 - G4 * G4xx - G4x * G5p + 2._dl * G4 * G5px - G4p * G5x) * H) + 2._dl * &
        a**3 * dphi**3 * ((-(dH * dphi * G4x * G5x)) + dphi**2 * (G4xx * G5p - 2._dl * G4x * G5px + G4px * G5x) * H &
        + ddphi * ((-6._dl) * G4x * G5x + G5p * G5x + 3._dl * G4 * G5xx) * H) - a * dphi**3 * (ddphi**2 * G5x**2 * H &
        - dddphi * dphi * G5x**2 * H + ddphi * dphi * (dH * G5x**2 + dphi * ((-4._dl) * G4xx * G5x + G5px * G5x + &
        4._dl * G4x * G5xx - G5p * G5xx) * H)) - a**2 * dphi**2 * (dddphi * dphi * ((-2._dl) * G4x + G5p) * G5x + &
        ddphi**2 * ((-2._dl) * G4x * G5x + 2._dl * G4 * G5xx) + 2._dl * dphi**4 * (G4xx * G5x - G4x * G5xx) * H**2 + &
        ddphi * dphi**2 * (2._dl * G4xx * G5p - 4._dl * G4x * G5px + G5p * G5px + 2._dl * G4px * G5x - G5pp * G5x + &
        G5x**2 * H**2)) + 4._dl * a**6 * (2._dl * ddphi * G4 * (G4x - G5p) + dphi**2 * (G4p * ((-G4x) + G5p) + G4 * &
        (G4px - G5pp - 3._dl * G5x * H**2))) - 2._dl * a**4 * (2._dl * ddphi**2 * G4 * G5x + dddphi * dphi * G4 * &
        G5x + ddphi * dphi**2 * (2._dl * G4x**2 - 2._dl * G4 * G4xx - 2._dl * G4x * G5p + 3._dl * G4 * G5px - G4p * &
        G5x) + dphi**4 * (G4px * G5p + 2._dl * G4 * G5xx * H**2 - G4x * (G5pp + 3._dl * G5x * H**2))))) / (a**2 * &
        (2._dl * a**3 * G4 + a * dphi**2 * ((-2._dl) * G4x + G5p) - dphi**3 * G5x * H)**2)
    M = M2
    dM = dM2
    ddM = ddM2
    alphaM = aM
    dalphaM = daM
    alphaB = aB
    dalphaB = daB
    alphaK = aK
    dalphaK = daK
    alphaT = aT
    dalphaT = daT

    end subroutine horndeski_eft_alphas

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes the background EoM coefficients
    subroutine horndeski_backeq_coefficients(a, H, dphi, rptot, G2, G2p, G2x, G2pp, G2px, G2xx, G2ppp, G2ppx, G2pxx, &
        G2xxx, G3p, G3x, G3pp, G3px, G3xx, G3ppp, G3ppx, G3pxx, G3xxx, G4, G4p, G4x, G4pp, G4px, G4xx, G4ppp, G4ppx, &
        G4pxx, G4xxx, G4pppp, G4pppx, G4ppxx, G4pxxx, G4xxxx, G5, G5p, G5x, G5pp, G5px, G5xx, G5ppp, G5ppx, G5pxx, &
        G5xxx, G5pppp, G5pppx, G5ppxx, G5pxxx, G5xxxx, cp_1, ch_1, cs_1, cp_2, ch_2, cs_2) &
        bind(c, name='horndeski_backeq_coefficients')

    implicit none

    real(c_double), value :: a      !<- the scale factor
    real(c_double), value :: H      !<- the physcial Hubble parameter
    real(c_double), value :: dphi   !<- conformal time derivative of the scalar field
    real(c_double), value :: rptot  !<- \f$ \rho_{matter} + p_{matter} \f$
    real(c_double), value :: G2     !<- Horndeski Lagrangian function \f$ G2 \f$
    real(c_double), value :: G2p    !<- Horndeski Lagrangian function \f$ G2_\phi \f$
    real(c_double), value :: G2x    !<- Horndeski Lagrangian function \f$ G2_X \f$
    real(c_double), value :: G2pp
    real(c_double), value :: G2px
    real(c_double), value :: G2xx
    real(c_double), value :: G2ppp
    real(c_double), value :: G2ppx
    real(c_double), value :: G2pxx
    real(c_double), value :: G2xxx
    real(c_double), value :: G3p
    real(c_double), value :: G3x
    real(c_double), value :: G3pp
    real(c_double), value :: G3px
    real(c_double), value :: G3xx
    real(c_double), value :: G3ppp
    real(c_double), value :: G3ppx
    real(c_double), value :: G3pxx
    real(c_double), value :: G3xxx
    real(c_double), value :: G4
    real(c_double), value :: G4p
    real(c_double), value :: G4x
    real(c_double), value :: G4pp
    real(c_double), value :: G4px
    real(c_double), value :: G4xx
    real(c_double), value :: G4ppp
    real(c_double), value :: G4ppx
    real(c_double), value :: G4pxx
    real(c_double), value :: G4xxx
    real(c_double), value :: G4pppp
    real(c_double), value :: G4pppx
    real(c_double), value :: G4ppxx
    real(c_double), value :: G4pxxx
    real(c_double), value :: G4xxxx
    real(c_double), value :: G5
    real(c_double), value :: G5p
    real(c_double), value :: G5x
    real(c_double), value :: G5pp
    real(c_double), value :: G5px
    real(c_double), value :: G5xx
    real(c_double), value :: G5ppp
    real(c_double), value :: G5ppx
    real(c_double), value :: G5pxx
    real(c_double), value :: G5xxx
    real(c_double), value :: G5pppp
    real(c_double), value :: G5pppx
    real(c_double), value :: G5ppxx
    real(c_double), value :: G5pxxx
    real(c_double), value :: G5xxxx
    real(c_double) :: cp_1   !<- EoM cp_i * phi'' + ch_i * H' + cs_i = 0, i=1,2
    real(c_double) :: ch_1
    real(c_double) :: cs_1
    real(c_double) :: cp_2
    real(c_double) :: ch_2
    real(c_double) :: cs_2

    real(c_double) :: X

    X = 0.5_dl * dphi * dphi / a / a
    cp_1 = ((-2._dl) * a * G4p + 4._dl * dphi * H * (G4x - G5p + 2._dl * G4xx * X - G5px * X) + 2._dl * a * X * (G3x &
        - 2._dl * G4px + H**2 * (3._dl * G5x + 2._dl * G5xx * X))) / (3._dl * a**2)
    ch_1 = (4._dl * dphi * G5x * H * X) / (3._dl * a) - (4._dl * (G4 + ((-2._dl) * G4x + G5p) * X)) / 3._dl
    cs_1 = (4._dl * dphi * H * (G4p - X * (2._dl * G3x - 6._dl * G4px + G5pp + 3._dl * G5x * H**2 + 2._dl * G5xx * &
        H**2 * X)) - a * (rptot + 2._dl * X * (G2x - 2._dl * (G3p - G4pp + H**2 * ((-5._dl) * G4x + 5._dl * G5p - &
        10._dl * G4xx * X + 6._dl * G5px * X))))) / 3._dl
    cp_2 = (6._dl * dphi * (G3x - 3._dl * G4px + (G3xx - 2._dl * G4pxx) * X) + (a * (G2x - 2._dl * (G3p + ((-G2xx) + &
        G3px) * X))) / H + 6._dl * a * H * (G4x - G5p + X * (8._dl * G4xx - 5._dl * G5px + 4._dl * G4xxx * X - 2._dl &
        * G5pxx * X)) + 2._dl * dphi * H**2 * (3._dl * G5x + X * (7._dl * G5xx + 2._dl * G5xxx * X))) / a**3
    ch_2 = ((-6._dl) * G4p + 6._dl * G3x * X - 12._dl * G4px * X) / (a * H) + (dphi * (12._dl * (G4x - G5p) + 12._dl &
        * (2._dl * G4xx - G5px) * X)) / a**2 + (H * (18._dl * G5x * X + 12._dl * G5xx * X**2)) / a
    cs_2 = (2._dl * dphi * (G2x - 2._dl * G3p + 6._dl * G4x * H**2 - 6._dl * G5p * H**2 - G2xx * X + 4._dl * G3px * &
        X - 6._dl * G4ppx * X - 6._dl * G4xx * H**2 * X + 8._dl * G5px * H**2 * X - 12._dl * G4xxx * H**2 * X**2 + &
        8._dl * G5pxx * H**2 * X**2)) / a - (G2p + 12._dl * G4p * H**2 + 2._dl * X * ((-G2px) + G3pp + H**2 * &
        ((-3._dl) * G3x + 3._dl * G5pp - 3._dl * G5x * H**2 + 6._dl * G3xx * X - 24._dl * G4pxx * X + 6._dl * G5ppx &
        * X + 8._dl * G5xx * H**2 * X + 4._dl * G5xxx * H**2 * X**2))) / H

    end subroutine horndeski_backeq_coefficients

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes the conformal time derivative of background EoM coefficients
    subroutine horndeski_backeq_coefficients_p(a, H, dH, dphi, ddphi, rptot, drptot, G2, G2p, G2x, G2pp, G2px, G2xx, &
        G2ppp, G2ppx, G2pxx, G2xxx, G3p, G3x, G3pp, G3px, G3xx, G3ppp, G3ppx, G3pxx, G3xxx, G4, G4p, G4x, G4pp, G4px, &
        G4xx, G4ppp, G4ppx, G4pxx, G4xxx, G4pppp, G4pppx, G4ppxx, G4pxxx, G4xxxx, G5, G5p, G5x, G5pp, G5px, G5xx, &
        G5ppp, G5ppx, G5pxx, G5xxx, G5pppp, G5pppx, G5ppxx, G5pxxx, G5xxxx, cp_1p, ch_1p, cs_1p, cp_2p, ch_2p, &
        cs_2p) &
        bind(c, name='horndeski_backeq_coefficients_p')

    implicit none

    real(c_double), value :: a      !<- the scale factor
    real(c_double), value :: H      !<- the physcial Hubble parameter
    real(c_double), value :: dH     !<- conformal time derivative of the physical Hubble
    real(c_double), value :: dphi   !<- conformal time derivative of the scalar field
    real(c_double), value :: ddphi  !<- double conformal time derivative of the scalar field
    real(c_double), value :: rptot  !<- \f$ \rho_{matter} + p_{matter} \f$
    real(c_double), value :: drptot !<- conformal time derivative of \f$ \rho_{matter} + p_{matter} \f$
    real(c_double), value :: G2     !<- Horndeski Lagrangian function \f$ G2 \f$
    real(c_double), value :: G2p    !<- Horndeski Lagrangian function \f$ G2_\phi \f$
    real(c_double), value :: G2x    !<- Horndeski Lagrangian function \f$ G2_X \f$
    real(c_double), value :: G2pp
    real(c_double), value :: G2px
    real(c_double), value :: G2xx
    real(c_double), value :: G2ppp
    real(c_double), value :: G2ppx
    real(c_double), value :: G2pxx
    real(c_double), value :: G2xxx
    real(c_double), value :: G3p
    real(c_double), value :: G3x
    real(c_double), value :: G3pp
    real(c_double), value :: G3px
    real(c_double), value :: G3xx
    real(c_double), value :: G3ppp
    real(c_double), value :: G3ppx
    real(c_double), value :: G3pxx
    real(c_double), value :: G3xxx
    real(c_double), value :: G4
    real(c_double), value :: G4p
    real(c_double), value :: G4x
    real(c_double), value :: G4pp
    real(c_double), value :: G4px
    real(c_double), value :: G4xx
    real(c_double), value :: G4ppp
    real(c_double), value :: G4ppx
    real(c_double), value :: G4pxx
    real(c_double), value :: G4xxx
    real(c_double), value :: G4pppp
    real(c_double), value :: G4pppx
    real(c_double), value :: G4ppxx
    real(c_double), value :: G4pxxx
    real(c_double), value :: G4xxxx
    real(c_double), value :: G5
    real(c_double), value :: G5p
    real(c_double), value :: G5x
    real(c_double), value :: G5pp
    real(c_double), value :: G5px
    real(c_double), value :: G5xx
    real(c_double), value :: G5ppp
    real(c_double), value :: G5ppx
    real(c_double), value :: G5pxx
    real(c_double), value :: G5xxx
    real(c_double), value :: G5pppp
    real(c_double), value :: G5pppx
    real(c_double), value :: G5ppxx
    real(c_double), value :: G5pxxx
    real(c_double), value :: G5xxxx
    real(c_double) :: cp_1p          !<- EoM cp_i * phi'' + ch_i * H' + cs_i = 0, i=1,2
    real(c_double) :: ch_1p
    real(c_double) :: cs_1p
    real(c_double) :: cp_2p
    real(c_double) :: ch_2p
    real(c_double) :: cs_2p

    real(c_double) :: X

    X = 0.5_dl * dphi * dphi / a / a
    cp_1p = (2._dl * ddphi * dphi * (G3x - 3._dl * G4px + 3._dl * G5x * H**2 + G3xx * X - 2._dl * G4pxx * X + 7._dl &
        * G5xx * H**2 * X + 2._dl * G5xxx * H**2 * X**2) - 2._dl * a**2 * ((-2._dl) * dH * H * X * (3._dl * G5x + &
        2._dl * G5xx * X) + dphi * (G4pp + 4._dl * G4x * H**2 - 4._dl * G5p * H**2 - G3px * X + 2._dl * G4ppx * X + &
        20._dl * G4xx * H**2 * X - 15._dl * G5px * H**2 * X + 8._dl * G4xxx * H**2 * X**2 - 6._dl * G5pxx * H**2 * &
        X**2)) - 2._dl * a**3 * H * ((-G4p) + X * (3._dl * G3x - 12._dl * G4px + 4._dl * G5pp + 9._dl * G5x * H**2 + &
        2._dl * G3xx * X - 12._dl * G4pxx * X + 4._dl * G5ppx * X + 16._dl * G5xx * H**2 * X + 4._dl * G5xxx * H**2 &
        * X**2)) + 4._dl * a * (dH * dphi * (G4x - G5p + 2._dl * G4xx * X - G5px * X) + ddphi * H * (G4x - G5p + X * &
        (8._dl * G4xx - 5._dl * G5px + 4._dl * G4xxx * X - 2._dl * G5pxx * X)))) / (3._dl * a**3)
    ch_1p = ((-4._dl) * (2._dl * a**3 * H * X * (G4x - G5p + 2._dl * (G4xx - G5px) * X) + ddphi * dphi * ((-G4x) + &
        G5p + ((-2._dl) * G4xx + G5px) * X) - a * X * (dH * dphi * G5x + ddphi * H * (3._dl * G5x + 2._dl * G5xx * &
        X)) + a**2 * dphi * (G4p + X * ((-2._dl) * G4px + G5pp + H**2 * (3._dl * G5x + 2._dl * G5xx * X))))) / &
        (3._dl * a**2)
    cs_1p = ((-2._dl) * ddphi * dphi * (G2x - 2._dl * G3p + 2._dl * G4pp + 10._dl * G4x * H**2 - 10._dl * G5p * H**2 &
        + G2xx * X - 2._dl * G3px * X + 2._dl * G4ppx * X + 50._dl * G4xx * H**2 * X - 34._dl * G5px * H**2 * X + &
        20._dl * G4xxx * H**2 * X**2 - 12._dl * G5pxx * H**2 * X**2) + a**3 * H * ((-rptot) + 2._dl * X * (G2x + &
        2._dl * ((-G3p) + 3._dl * G4pp + 5._dl * G4x * H**2 - 5._dl * G5p * H**2 + G2xx * X - 6._dl * G3px * X + &
        14._dl * G4ppx * X - 2._dl * G5ppp * X + 40._dl * G4xx * H**2 * X - 34._dl * G5px * H**2 * X + 20._dl * &
        G4xxx * H**2 * X**2 - 16._dl * G5pxx * H**2 * X**2))) + a**2 * ((-drptot) - 2._dl * dphi * G2px * X - 8._dl &
        * dH * H * X * (5._dl * G4x - 5._dl * G5p + 10._dl * G4xx * X - 6._dl * G5px * X) + 4._dl * dphi * X * (G3pp &
        - (G4ppp + 22._dl * G4pxx * H**2 * X) + H**2 * (4._dl * G3x - 19._dl * G4px + 7._dl * G5pp + 6._dl * G5x * &
        H**2 + 4._dl * G3xx * X + 8._dl * G5ppx * X + 14._dl * G5xx * H**2 * X + 4._dl * G5xxx * H**2 * X**2))) + &
        4._dl * a * (dH * dphi * (G4p - X * (2._dl * G3x - 6._dl * G4px + G5pp + 9._dl * G5x * H**2 + 6._dl * G5xx * &
        H**2 * X)) + ddphi * H * (G4p - X * (6._dl * G3x - 20._dl * G4px + 3._dl * G5pp + 9._dl * G5x * H**2 + 4._dl &
        * G3xx * X - 12._dl * G4pxx * X + 2._dl * G5ppx * X + 16._dl * G5xx * H**2 * X + 4._dl * G5xxx * H**2 * &
        X**2)))) / (3._dl * a)
    cp_2p = (ddphi * dphi * H * (3._dl * G2xx - 4._dl * G3px + 54._dl * G4xx * H**2 - 36._dl * G5px * H**2 + 2._dl * &
        G2xxx * X - 2._dl * G3pxx * X + 96._dl * G4xxx * H**2 * X - 54._dl * G5pxx * H**2 * X + 24._dl * G4xxxx * &
        H**2 * X**2 - 12._dl * G5pxxx * H**2 * X**2) - 2._dl * a**3 * H**2 * (G2x - 2._dl * G3p + 6._dl * G4x * H**2 &
        - 6._dl * G5p * H**2 + 5._dl * G2xx * X - 12._dl * G3px * X + 18._dl * G4ppx * X + 102._dl * G4xx * H**2 * X &
        - 72._dl * G5px * H**2 * X + 2._dl * G2xxx * X**2 - 8._dl * G3pxx * X**2 + 12._dl * G4ppxx * X**2 + 120._dl &
        * G4xxx * H**2 * X**2 - 80._dl * G5pxx * H**2 * X**2 + 24._dl * G4xxxx * H**2 * X**3 - 16._dl * G5pxxx * &
        H**2 * X**3) + 2._dl * a * H**2 * (2._dl * dH * dphi * H * (3._dl * G5x + 7._dl * G5xx * X + 2._dl * G5xxx * &
        X**2) + ddphi * (3._dl * G3x - 9._dl * G4px + 3._dl * G5x * H**2 + 15._dl * G3xx * X - 36._dl * G4pxx * X + &
        27._dl * G5xx * H**2 * X + 6._dl * G3xxx * X**2 - 12._dl * G4pxxx * X**2 + 24._dl * G5xxx * H**2 * X**2 + &
        4._dl * G5xxxx * H**2 * X**3)) - a**2 * (dH * (G2x - 2._dl * (G3p + 3._dl * G4x * H**2 - 3._dl * G5p * H**2 &
        - G2xx * X + G3px * X + 24._dl * G4xx * H**2 * X - 15._dl * G5px * H**2 * X + 12._dl * G4xxx * H**2 * X**2 - &
        6._dl * G5pxx * H**2 * X**2)) + dphi * H * ((-G2px) + 2._dl * (G3pp + 9._dl * G3x * H**2 - 30._dl * G4px * &
        H**2 + 3._dl * G5pp * H**2 + 9._dl * G5x * H**4 - G2pxx * X + G3ppx * X + 21._dl * G3xx * H**2 * X - 72._dl &
        * G4pxx * H**2 * X + 15._dl * G5ppx * H**2 * X + 41._dl * G5xx * H**4 * X + 6._dl * G3xxx * H**2 * X**2 - &
        24._dl * G4pxxx * H**2 * X**2 + 6._dl * G5ppxx * H**2 * X**2 + 28._dl * G5xxx * H**4 * X**2 + 4._dl * G5xxxx &
        * H**4 * X**3)))) / (a**4 * H**2)
    ch_2p = 6._dl * (G4p + (ddphi * dphi * (G3x - 3._dl * G4px + 3._dl * G5x * H**2 + G3xx * X - 2._dl * G4pxx * X + &
        7._dl * G5xx * H**2 * X + 2._dl * G5xxx * H**2 * X**2)) / (a**3 * H) - X * (3._dl * G3x - 12._dl * G4px + &
        4._dl * G5pp + 9._dl * G5x * H**2 + 2._dl * G3xx * X - 12._dl * G4pxx * X + 4._dl * G5ppx * X + 16._dl * &
        G5xx * H**2 * X + 4._dl * G5xxx * H**2 * X**2) + (2._dl * ddphi * (G4x - G5p + X * (8._dl * G4xx - 5._dl * &
        G5px + 4._dl * G4xxx * X - 2._dl * G5pxx * X))) / a**2 + ((-(dphi * H * (G4pp + 4._dl * G4x * H**2 - 4._dl * &
        G5p * H**2 - G3px * X + 2._dl * G4ppx * X + 20._dl * G4xx * H**2 * X - 15._dl * G5px * H**2 * X + 8._dl * &
        G4xxx * H**2 * X**2 - 6._dl * G5pxx * H**2 * X**2))) + dH * (G4p + X * ((-G3x) + 2._dl * G4px + H**2 * &
        (3._dl * G5x + 2._dl * G5xx * X)))) / (a * H**2))
    cs_2p = 2._dl * a * X * (G2px - 2._dl * (G3pp + 3._dl * G3x * H**2 - 12._dl * G4px * H**2 + 3._dl * G5pp * H**2 &
        + 3._dl * G5x * H**4 + 2._dl * G2pxx * X - 5._dl * G3ppx * X + 6._dl * G4pppx * X - 9._dl * G3xx * H**2 * X &
        + 54._dl * G4pxx * H**2 * X - 23._dl * G5ppx * H**2 * X - 13._dl * G5xx * H**4 * X - 6._dl * G3xxx * H**2 * &
        X**2 + 36._dl * G4pxxx * H**2 * X**2 - 14._dl * G5ppxx * H**2 * X**2 - 20._dl * G5xxx * H**4 * X**2 - 4._dl &
        * G5xxxx * H**4 * X**3)) + (ddphi * dphi * (G2px - 2._dl * (G3pp - 3._dl * G3x * H**2 + 6._dl * G4px * H**2 &
        + 3._dl * G5pp * H**2 - 3._dl * G5x * H**4 - G2pxx * X + G3ppx * X + 9._dl * G3xx * H**2 * X - 48._dl * &
        G4pxx * H**2 * X + 15._dl * G5ppx * H**2 * X + 13._dl * G5xx * H**4 * X + 6._dl * G3xxx * H**2 * X**2 - &
        24._dl * G4pxxx * H**2 * X**2 + 6._dl * G5ppxx * H**2 * X**2 + 20._dl * G5xxx * H**4 * X**2 + 4._dl * G5xxxx &
        * H**4 * X**3))) / (a**2 * H) + (2._dl * (ddphi * (G2x - 2._dl * G3p + 6._dl * G4x * H**2 - 6._dl * G5p * &
        H**2 - G2xx * X + 8._dl * G3px * X - 18._dl * G4ppx * X - 6._dl * G4xx * H**2 * X + 12._dl * G5px * H**2 * X &
        - 2._dl * G2xxx * X**2 + 8._dl * G3pxx * X**2 - 12._dl * G4ppxx * X**2 - 72._dl * G4xxx * H**2 * X**2 + &
        56._dl * G5pxx * H**2 * X**2 - 24._dl * G4xxxx * H**2 * X**3 + 16._dl * G5pxxx * H**2 * X**3) + 4._dl * dH * &
        dphi * H * (3._dl * G4x - 3._dl * G5p + X * ((-3._dl) * G4xx + 4._dl * G5px - 6._dl * G4xxx * X + 4._dl * &
        G5pxx * X)))) / a + ((-(dphi * H * (G2pp + 2._dl * (G2x * H**2 - 2._dl * G3p * H**2 + 6._dl * G4pp * H**2 + &
        6._dl * G4x * H**4 - 6._dl * G5p * H**4 - G2ppx * X + G3ppp * X - G2xx * H**2 * X + 5._dl * G3px * H**2 * X &
        - 18._dl * G4ppx * H**2 * X + 3._dl * G5ppp * H**2 * X - 6._dl * G4xx * H**4 * X + 9._dl * G5px * H**4 * X - &
        2._dl * G2xxx * H**2 * X**2 + 14._dl * G3pxx * H**2 * X**2 - 36._dl * G4ppxx * H**2 * X**2 + 6._dl * G5pppx &
        * H**2 * X**2 - 72._dl * G4xxx * H**4 * X**2 + 64._dl * G5pxx * H**4 * X**2 - 24._dl * G4xxxx * H**4 * X**3 &
        + 20._dl * G5pxxx * H**4 * X**3)))) + dH * (G2p - 2._dl * (6._dl * G4p * H**2 + X * (G2px - G3pp + 3._dl * &
        H**2 * ((-G3x) + G5pp - 3._dl * G5x * H**2 + 2._dl * G3xx * X - 8._dl * G4pxx * X + 2._dl * G5ppx * X + &
        8._dl * G5xx * H**2 * X + 4._dl * G5xxx * H**2 * X**2))))) / H**2

    end subroutine horndeski_backeq_coefficients_p

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes the designer background EoM coefficients
    subroutine horndeski_designereq_coefficients(a, H, dH, ddH, dphi, rptot, drptot, G2, G2p, G2x, G2pp, G2px, G2xx, &
        G2ppp, G2ppx, G2pxx, G2xxx, G3p, G3x, G3pp, G3px, G3xx, G3ppp, G3ppx, G3pxx, G3xxx, G4, G4p, G4x, G4pp, G4px, &
        G4xx, G4ppp, G4ppx, G4pxx, G4xxx, G4pppp, G4pppx, G4ppxx, G4pxxx, G4xxxx, G5, G5p, G5x, G5pp, G5px, G5xx, &
        G5ppp, G5ppx, G5pxx, G5xxx, G5pppp, G5pppx, G5ppxx, G5pxxx, G5xxxx, d0, d1, d2) &
        bind(c, name='horndeski_designereq_coefficients')

    implicit none

    real(c_double), value :: a      !<- the scale factor
    real(c_double), value :: H      !<- the physcial Hubble parameter
    real(c_double), value :: dH     !<- conformal time derivative of the physical Hubble
    real(c_double), value :: ddH    !<- second conformal time derivative of the physical Hubble
    real(c_double), value :: dphi   !<- conformal time derivative of the scalar field
    real(c_double), value :: rptot  !<- \f$ \rho_{matter} + p_{matter} \f$
    real(c_double), value :: drptot !<- conformal time derivative of \f$ \rho_{matter} + p_{matter} \f$
    real(c_double), value :: G2     !<- Horndeski Lagrangian function \f$ G2 \f$
    real(c_double), value :: G2p    !<- Horndeski Lagrangian function \f$ G2_\phi \f$
    real(c_double), value :: G2x    !<- Horndeski Lagrangian function \f$ G2_X \f$
    real(c_double), value :: G2pp
    real(c_double), value :: G2px
    real(c_double), value :: G2xx
    real(c_double), value :: G2ppp
    real(c_double), value :: G2ppx
    real(c_double), value :: G2pxx
    real(c_double), value :: G2xxx
    real(c_double), value :: G3p
    real(c_double), value :: G3x
    real(c_double), value :: G3pp
    real(c_double), value :: G3px
    real(c_double), value :: G3xx
    real(c_double), value :: G3ppp
    real(c_double), value :: G3ppx
    real(c_double), value :: G3pxx
    real(c_double), value :: G3xxx
    real(c_double), value :: G4
    real(c_double), value :: G4p
    real(c_double), value :: G4x
    real(c_double), value :: G4pp
    real(c_double), value :: G4px
    real(c_double), value :: G4xx
    real(c_double), value :: G4ppp
    real(c_double), value :: G4ppx
    real(c_double), value :: G4pxx
    real(c_double), value :: G4xxx
    real(c_double), value :: G4pppp
    real(c_double), value :: G4pppx
    real(c_double), value :: G4ppxx
    real(c_double), value :: G4pxxx
    real(c_double), value :: G4xxxx
    real(c_double), value :: G5
    real(c_double), value :: G5p
    real(c_double), value :: G5x
    real(c_double), value :: G5pp
    real(c_double), value :: G5px
    real(c_double), value :: G5xx
    real(c_double), value :: G5ppp
    real(c_double), value :: G5ppx
    real(c_double), value :: G5pxx
    real(c_double), value :: G5xxx
    real(c_double), value :: G5pppp
    real(c_double), value :: G5pppx
    real(c_double), value :: G5ppxx
    real(c_double), value :: G5pxxx
    real(c_double), value :: G5xxxx
    real(c_double) :: d0            !<- EoM d2*ddphi^2 + d1*ddphi + d0 = 0
    real(c_double) :: d1
    real(c_double) :: d2

    real(c_double) :: X

    X = 0.5_dl * dphi * dphi / a / a
    d0 = (4._dl * dphi * G5x * (dH**2 + ddH * H) * X - 4._dl * a * (ddH * (G4 + ((-2._dl) * G4x + G5p) * X) + 2._dl &
        * dH * dphi * X * (G3x - 4._dl * G4px + G5pp + 6._dl * G5x * H**2 + 4._dl * G5xx * H**2 * X)) + a**3 * H * &
        ((-rptot) + 2._dl * X * (G2x + 2._dl * ((-G3p) + 3._dl * G4pp + 5._dl * G4x * H**2 - 5._dl * G5p * H**2 + &
        G2xx * X - 6._dl * G3px * X + 14._dl * G4ppx * X - 2._dl * G5ppp * X + 40._dl * G4xx * H**2 * X - 34._dl * &
        G5px * H**2 * X + 20._dl * G4xxx * H**2 * X**2 - 16._dl * G5pxx * H**2 * X**2))) + a**2 * ((-drptot) + 2._dl &
        * X * ((-(dphi * G2px)) - 8._dl * dH * H * (3._dl * G4x - 3._dl * G5p + 6._dl * G4xx * X - 4._dl * G5px * X) &
        + 2._dl * dphi * (G3pp - (G4ppp + 22._dl * G4pxx * H**2 * X) + H**2 * (4._dl * G3x - 19._dl * G4px + 7._dl * &
        G5pp + 6._dl * G5x * H**2 + 4._dl * G3xx * X + 8._dl * G5ppx * X + 14._dl * G5xx * H**2 * X + 4._dl * G5xxx &
        * H**2 * X**2))))) / (3._dl * a)
    d1 = ((-2._dl) * (4._dl * dH * dphi * ((-G4x) + G5p - 2._dl * G4xx * X + G5px * X) + a * ((-4._dl) * dH * H * X &
        * (3._dl * G5x + 2._dl * G5xx * X) + dphi * (G2x - 2._dl * G3p + 3._dl * G4pp + 14._dl * G4x * H**2 - 14._dl &
        * G5p * H**2 + G2xx * X - 3._dl * G3px * X + 4._dl * G4ppx * X + 70._dl * G4xx * H**2 * X - 49._dl * G5px * &
        H**2 * X + 28._dl * G4xxx * H**2 * X**2 - 18._dl * G5pxx * H**2 * X**2)) + a**2 * H * ((-3._dl) * G4p + X * &
        (15._dl * G3x - 52._dl * G4px + 10._dl * G5pp + 27._dl * G5x * H**2 + 10._dl * G3xx * X - 36._dl * G4pxx * X &
        + 8._dl * G5ppx * X + 48._dl * G5xx * H**2 * X + 12._dl * G5xxx * H**2 * X**2)))) / (3._dl * a**2)
    d2 = (2._dl * dphi * (G3x - 3._dl * G4px + 3._dl * G5x * H**2 + G3xx * X - 2._dl * G4pxx * X + 7._dl * G5xx * &
        H**2 * X + 2._dl * G5xxx * H**2 * X**2) + 4._dl * a * H * (G4x - G5p + X * (8._dl * G4xx - 5._dl * G5px + &
        4._dl * G4xxx * X - 2._dl * G5pxx * X))) / (3._dl * a**3)

    end subroutine horndeski_designereq_coefficients

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes the derived designer background EoM coefficients
    subroutine horndeski_designereq_coefficients_p(a, H, dH, ddH, dddH, dphi, ddphi, rptot, drptot, ddrptot, G2, G2p, &
        G2x, G2pp, G2px, G2xx, G2ppp, G2ppx, G2pxx, G2xxx, G3p, G3x, G3pp, G3px, G3xx, G3ppp, G3ppx, G3pxx, G3xxx, G4, &
        G4p, G4x, G4pp, G4px, G4xx, G4ppp, G4ppx, G4pxx, G4xxx, G4pppp, G4pppx, G4ppxx, G4pxxx, G4xxxx, G5, G5p, G5x, &
        G5pp, G5px, G5xx, G5ppp, G5ppx, G5pxx, G5xxx, G5pppp, G5pppx, G5ppxx, G5pxxx, G5xxxx, dp0, dp1) &
        bind(c, name='horndeski_designereq_coefficients_p')

    implicit none

    real(c_double), value :: a      !<- the scale factor
    real(c_double), value :: H      !<- the physcial Hubble parameter
    real(c_double), value :: dH     !<- conformal time derivative of the physical Hubble
    real(c_double), value :: ddH    !<- second conformal time derivative of the physical Hubble
    real(c_double), value :: dddH    !<- third conformal time derivative of the physical Hubble
    real(c_double), value :: dphi   !<- conformal time derivative of the scalar field
    real(c_double), value :: ddphi  !<- second conformal time derivative of the scalar field
    real(c_double), value :: rptot  !<- \f$ \rho_{matter} + p_{matter} \f$
    real(c_double), value :: drptot !<- conformal time derivative of \f$ \rho_{matter} + p_{matter} \f$
    real(c_double), value :: ddrptot !<- second conformal time derivative of \f$ \rho_{matter} + p_{matter} \f$
    real(c_double), value :: G2     !<- Horndeski Lagrangian function \f$ G2 \f$
    real(c_double), value :: G2p    !<- Horndeski Lagrangian function \f$ G2_\phi \f$
    real(c_double), value :: G2x    !<- Horndeski Lagrangian function \f$ G2_X \f$
    real(c_double), value :: G2pp
    real(c_double), value :: G2px
    real(c_double), value :: G2xx
    real(c_double), value :: G2ppp
    real(c_double), value :: G2ppx
    real(c_double), value :: G2pxx
    real(c_double), value :: G2xxx
    real(c_double), value :: G3p
    real(c_double), value :: G3x
    real(c_double), value :: G3pp
    real(c_double), value :: G3px
    real(c_double), value :: G3xx
    real(c_double), value :: G3ppp
    real(c_double), value :: G3ppx
    real(c_double), value :: G3pxx
    real(c_double), value :: G3xxx
    real(c_double), value :: G4
    real(c_double), value :: G4p
    real(c_double), value :: G4x
    real(c_double), value :: G4pp
    real(c_double), value :: G4px
    real(c_double), value :: G4xx
    real(c_double), value :: G4ppp
    real(c_double), value :: G4ppx
    real(c_double), value :: G4pxx
    real(c_double), value :: G4xxx
    real(c_double), value :: G4pppp
    real(c_double), value :: G4pppx
    real(c_double), value :: G4ppxx
    real(c_double), value :: G4pxxx
    real(c_double), value :: G4xxxx
    real(c_double), value :: G5
    real(c_double), value :: G5p
    real(c_double), value :: G5x
    real(c_double), value :: G5pp
    real(c_double), value :: G5px
    real(c_double), value :: G5xx
    real(c_double), value :: G5ppp
    real(c_double), value :: G5ppx
    real(c_double), value :: G5pxx
    real(c_double), value :: G5xxx
    real(c_double), value :: G5pppp
    real(c_double), value :: G5pppx
    real(c_double), value :: G5ppxx
    real(c_double), value :: G5pxxx
    real(c_double), value :: G5xxxx
    real(c_double) :: dp0            !<- EoM dp1*dddphi + dp1 = 0
    real(c_double) :: dp1

    real(c_double) :: X

    X = 0.5_dl * dphi * dphi / a / a
    dp0 = ((-2._dl) * dphi**9 * G5xxxx * H**5 + 4._dl * a**2 * ddphi**3 * dphi * H * (9._dl * G4xx - 6._dl * G5px + &
        X * (16._dl * G4xxx - 9._dl * G5pxx + 4._dl * G4xxxx * X - 2._dl * G5pxxx * X)) + 4._dl * a**4 * ddphi * &
        (3._dl * ddH * dphi * (G4x - G5p + 2._dl * G4xx * X - G5px * X) - ddphi * dphi * H * (9._dl * G3x - 33._dl * &
        G4px + 6._dl * G5pp + 18._dl * G5x * H**2 + 21._dl * G3xx * X - 78._dl * G4pxx * X + 18._dl * G5ppx * X + &
        82._dl * G5xx * H**2 * X + 6._dl * G3xxx * X**2 - 24._dl * G4pxxx * X**2 + 6._dl * G5ppxx * X**2 + 56._dl * &
        G5xxx * H**2 * X**2 + 8._dl * G5xxxx * H**2 * X**3) + 3._dl * ddphi * dH * (G4x - G5p + X * (8._dl * G4xx - &
        5._dl * G5px + 4._dl * G4xxx * X - 2._dl * G5pxx * X))) - 2._dl * a**9 * (2._dl * (G2ppx - 2._dl * G3ppp + &
        2._dl * G4pppp) * X**2 + H**2 * (rptot + 2._dl * X**2 * (3._dl * G2xx - 22._dl * G3px + 72._dl * G4ppx - &
        18._dl * G5ppp + 2._dl * G2xxx * X - 20._dl * G3pxx * X + 72._dl * G4ppxx * X - 20._dl * G5pppx * X)) + &
        20._dl * H**4 * X**2 * (9._dl * G4xx - 9._dl * G5px + 4._dl * X * (4._dl * G4xxx - 4._dl * G5pxx + G4xxxx * &
        X - G5pxxx * X))) + 2._dl * a**3 * ddphi**2 * (ddphi * (G3x - 3._dl * G4px + 3._dl * G5x * H**2 + 5._dl * &
        G3xx * X - 12._dl * G4pxx * X + 27._dl * G5xx * H**2 * X + 2._dl * G3xxx * X**2 - 4._dl * G4pxxx * X**2 + &
        24._dl * G5xxx * H**2 * X**2 + 4._dl * G5xxxx * H**2 * X**3) + 6._dl * dH * dphi * H * (3._dl * G5x + X * &
        (7._dl * G5xx + 2._dl * G5xxx * X))) - a**8 * (2._dl * drptot * H + dH * rptot - 4._dl * dphi * G2px * H * X &
        + 8._dl * dphi * H * X * (G3pp - 2._dl * G4ppp + 2._dl * G3x * H**2 - 12._dl * G4px * H**2 + 6._dl * G5pp * &
        H**2 + 3._dl * G5x * H**4 - G2pxx * X + 4._dl * G3ppx * X - 8._dl * G4pppx * X + G5pppp * X + 10._dl * G3xx &
        * H**2 * X - 72._dl * G4pxx * H**2 * X + 36._dl * G5ppx * H**2 * X + 27._dl * G5xx * H**4 * X + 4._dl * &
        G3xxx * H**2 * X**2 - 32._dl * G4pxxx * H**2 * X**2 + 16._dl * G5ppxx * H**2 * X**2 + 24._dl * G5xxx * H**4 &
        * X**2) - 2._dl * dH * X * (G2x + 2._dl * ((-G3p) + 3._dl * G4pp + 27._dl * G4x * H**2 - 27._dl * G5p * H**2 &
        + G2xx * X - 10._dl * G3px * X + 30._dl * G4ppx * X - 6._dl * G5ppp * X + 216._dl * G4xx * H**2 * X - &
        198._dl * G5px * H**2 * X + 108._dl * G4xxx * H**2 * X**2 - 96._dl * G5pxx * H**2 * X**2))) - 2._dl * a**5 * &
        ((-2._dl) * dphi * G5x * (3._dl * ddH * dH + dddH * H) * X + ddphi**2 * (G2x - 2._dl * G3p + 3._dl * G4pp + &
        18._dl * G4x * H**2 - 18._dl * G5p * H**2 + 5._dl * G2xx * X - 15._dl * G3px * X + 24._dl * G4ppx * X + &
        306._dl * G4xx * H**2 * X - 225._dl * G5px * H**2 * X + 2._dl * G2xxx * X**2 - 8._dl * G3pxx * X**2 + 12._dl &
        * G4ppxx * X**2 + 360._dl * G4xxx * H**2 * X**2 - 246._dl * G5pxx * H**2 * X**2 + 72._dl * G4xxxx * H**2 * &
        X**3 - 48._dl * G5pxxx * H**2 * X**3) - 6._dl * ddphi * (dH**2 * X * (3._dl * G5x + 2._dl * G5xx * X) + ddH &
        * H * X * (3._dl * G5x + 2._dl * G5xx * X) + dH * dphi * H * ((-10._dl) * G4x + 10._dl * G5p + X * &
        ((-50._dl) * G4xx + 37._dl * G5px - 20._dl * G4xxx * X + 14._dl * G5pxx * X)))) + a**7 * ((-ddrptot) + 2._dl &
        * X * ((-5._dl) * ddphi * G2px - 4._dl * ddH * H * (7._dl * G4x - 7._dl * G5p + 2._dl * (7._dl * G4xx - &
        5._dl * G5px) * X) - 12._dl * dH**2 * (2._dl * G4x - 2._dl * G5p + 4._dl * G4xx * X - 3._dl * G5px * X) + &
        4._dl * dH * dphi * H * (6._dl * G3x - 33._dl * G4px + 15._dl * G5pp + 24._dl * G5x * H**2 + 6._dl * G3xx * &
        X - 42._dl * G4pxx * X + 18._dl * G5ppx * X + 56._dl * G5xx * H**2 * X + 16._dl * G5xxx * H**2 * X**2) + &
        2._dl * ddphi * (5._dl * G3pp + 27._dl * G3x * H**2 - 126._dl * G4px * H**2 + 45._dl * G5pp * H**2 + 45._dl &
        * G5x * H**4 - 2._dl * G2pxx * X + 5._dl * G3ppx * X - 6._dl * G4pppx * X + 63._dl * G3xx * H**2 * X + &
        129._dl * G5ppx * H**2 * X + 205._dl * G5xx * H**4 * X + 18._dl * G3xxx * H**2 * X**2 - 108._dl * G4pxxx * &
        H**2 * X**2 + 42._dl * G5ppxx * H**2 * X**2 + 140._dl * G5xxx * H**4 * X**2 + 20._dl * G5xxxx * H**4 * X**3 &
        - 6._dl * (G4ppp + 57._dl * G4pxx * H**2 * X)))) - 2._dl * a**6 * (2._dl * dddH * (G4 + ((-2._dl) * G4x + &
        G5p) * X) + ddphi * ((-3._dl) * dH * G4p + 3._dl * dH * X * (9._dl * G3x - 36._dl * G4px + 10._dl * G5pp + &
        63._dl * G5x * H**2 + 6._dl * G3xx * X - 28._dl * G4pxx * X + 8._dl * G5ppx * X + 112._dl * G5xx * H**2 * X &
        + 28._dl * G5xxx * H**2 * X**2) - 2._dl * dphi * H * (G2x - 2._dl * G3p + 6._dl * G4pp + 12._dl * G4x * H**2 &
        - 12._dl * G5p * H**2 + 5._dl * G2xx * X - 27._dl * G3px * X + 66._dl * G4ppx * X - 9._dl * G5ppp * X + &
        204._dl * G4xx * H**2 * X - 174._dl * G5px * H**2 * X + 2._dl * G2xxx * X**2 - 14._dl * G3pxx * X**2 + &
        36._dl * G4ppxx * X**2 - 6._dl * G5pppx * X**2 + 240._dl * G4xxx * H**2 * X**2 - 200._dl * G5pxx * H**2 * &
        X**2 + 48._dl * G4xxxx * H**2 * X**3 - 40._dl * G5pxxx * H**2 * X**3)) + 2._dl * dphi * (9._dl * dH**2 * H * &
        X * (3._dl * G5x + 2._dl * G5xx * X) + ddH * (G4p + X * (2._dl * G3x - 10._dl * G4px + 3._dl * G5pp + 15._dl &
        * G5x * H**2 + 10._dl * G5xx * H**2 * X))))) / (3._dl * a**6)
    dp1 = (6._dl * ddphi * dphi * (G3x - 3._dl * G4px + 3._dl * G5x * H**2 + G3xx * X - 2._dl * G4pxx * X + 7._dl * &
        G5xx * H**2 * X + 2._dl * G5xxx * H**2 * X**2) - 2._dl * a**2 * ((-6._dl) * dH * H * X * (3._dl * G5x + &
        2._dl * G5xx * X) + dphi * (G2x - 2._dl * G3p + 4._dl * G4pp + 18._dl * G4x * H**2 - 18._dl * G5p * H**2 + &
        G2xx * X - 4._dl * G3px * X + 6._dl * G4ppx * X + 90._dl * G4xx * H**2 * X - 64._dl * G5px * H**2 * X + &
        36._dl * G4xxx * H**2 * X**2 - 24._dl * G5pxx * H**2 * X**2)) - 4._dl * a**3 * H * ((-2._dl) * G4p + X * &
        (9._dl * G3x - 32._dl * G4px + 7._dl * G5pp + 18._dl * G5x * H**2 + 6._dl * G3xx * X - 24._dl * G4pxx * X + &
        6._dl * G5ppx * X + 32._dl * G5xx * H**2 * X + 8._dl * G5xxx * H**2 * X**2)) + 12._dl * a * (dH * dphi * &
        (G4x - G5p + 2._dl * G4xx * X - G5px * X) + ddphi * H * (G4x - G5p + X * (8._dl * G4xx - 5._dl * G5px + &
        4._dl * G4xxx * X - 2._dl * G5pxx * X)))) / (3._dl * a**3)

    end subroutine horndeski_designereq_coefficients_p

    ! ---------------------------------------------------------------------------------------------

    end module EFTCAMB_Horndeski_coefficients
