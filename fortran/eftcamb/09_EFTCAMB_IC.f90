!----------------------------------------------------------------------------------------
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
!----------------------------------------------------------------------------------------

!> @file 09_EFTCAMB_IC.f90
!! This file contains the perturbation initial conditions calculator.


!----------------------------------------------------------------------------------------
!> This module contains the perturbation initial conditions calculator.

!> @author Bin Hu, Marco Raveri

submodule (GaugeInterface) EFTCAMB_IC

use precision
use model
use EFTCAMB_mixed_algorithms
use MassiveNu

implicit none

contains

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that initializes EFTCAMB perturbations initial conditions.
    module subroutine EFTCAMBInitialConditions( y, EV, tau )

        implicit none

        type(EvolutionVars) EV  !< CAMB evolutionary variables.
        real(dl) :: y(EV%nvar)  !< status vector with the values of the perturbations.
        real(dl) :: tau         !< conformal time at which the initial conditions are desired.

        real(dl) :: a, a2, k, k2, etak, grhob_t, grhoc_t, grhor_t, grhog_t, grhov_t, gpres, grho_matter
        real(dl) :: dgrho_matter, dgq, dgpnu, grho, adotoa, EFT_grhonu, EFT_gpinu, EFT_grhonudot
        real(dl) :: EFT_gpinudot, grhormass_t, adotdota, cothxor, cs2, dopacity, dgrho, z, dz, EFT_gpinudotdot
        real(dl) :: clxr, qr, pir, clxg, qg, pig, opacity, vb, clxc, clxb, EFTpiEdot
        real(dl) :: Omega0, Gamma10, Gamma20, Gamma30
        real(dl) :: wnu_arr(max_nu)

        integer  :: nu_i

        k    = EV%k_buf
        k2   = EV%k2_buf

        !  Get background scale factor, sound speed and ionisation fraction.
        if (EV%TightCoupling) then
            call EV%ThermoData%Values(tau,a,cs2,opacity,dopacity)
        else
            call EV%ThermoData%Values(tau,a,cs2,opacity)
        end if

        ! copy variables:
        etak = y(ix_etak) ! conformal time
        clxc = y(ix_clxc) ! cdm density perturbations
        clxb = y(ix_clxb) ! baryons density perturbations
        vb   = y(ix_vb) ! baryons velocity perturbations
        ! process:
        a2   = a*a
        ! compute background densities of different species
        grhob_t = State%grhob/a         ! 8\pi G_N \rho_b a^2: baryons background density
        grhoc_t = State%grhoc/a         ! 8\pi G_N \rho_{cdm} a^2: cold dark matter background density
        grhor_t = State%grhornomass/a2  ! 8\pi G_N \rho_{\nu} a^2: massless neutrinos background density
        grhog_t = State%grhog/a2        ! 8\pi G_N \rho_{\gamma} a^2: radiation background density
        grhov_t = 0._dl

        ! start computing background total pressure and total density:
        gpres        = 0._dl
        grho_matter  = grhob_t +grhoc_t
        ! start computing total density and velocity perturbations:
        dgrho_matter = grhob_t*clxb +grhoc_t*clxc ! 8\pi G_N \sum(rho_i*clx_i) a^2   : total density perturbation
        dgq          = grhob_t*vb                 ! 8\pi G_N \sum(rho_i+p_i)*v_i a^2 : total velocity perturbation
        ! add massive neutrinos to both background and perturbations:

        dgpnu = 0._dl
        if (State%CP%Num_Nu_Massive > 0) then
            call MassiveNuVars(EV,y,a,grho_matter,gpres,dgrho_matter,dgq, wnu_arr, dgp=dgpnu)
        end if

        ! add radiation, massless neutrinos and Lambda to total background density:
        grho = grho_matter +grhor_t +grhog_t +grhov_t

        ! compute gpres: add radiation and massless neutrinos to massive neutrinos
        gpres = gpres + (grhog_t+grhor_t)/3._dl
        ! initialize the cache:
        call EV%eft_cache%initialize()
        ! start to fill the cache:
        EV%eft_cache%a           = a
        EV%eft_cache%tau         = tau
        EV%eft_cache%k           = k
        EV%eft_cache%grhom_t     = grho
        EV%eft_cache%gpresm_t    = gpres
        EV%eft_cache%grhob_t     = grhob_t
        EV%eft_cache%grhoc_t     = grhoc_t
        EV%eft_cache%grhor_t     = grhor_t
        EV%eft_cache%grhog_t     = grhog_t
        ! compute the other things:
        select type ( model => CP%EFTCAMB%model )
            ! compute the background and the background EFT functions.
            class is ( EFTCAMB_full_model )
            ! background for full models. Here the expansion history is computed from the
            ! EFT functions. Hence compute them first and then compute the expansion history.
            call CP%EFTCAMB%model%compute_background_EFT_functions( a, CP%eft_par_cache , EV%eft_cache )
            call CP%EFTCAMB%model%compute_adotoa( a, CP%eft_par_cache , EV%eft_cache )
            class is ( EFTCAMB_designer_model )
            ! background for designer models. Here the expansion history is parametrized
            ! and does not depend on the EFT functions. Hence compute first the expansion history
            ! and then the EFT functions.
            call CP%EFTCAMB%model%compute_adotoa( a, CP%eft_par_cache , EV%eft_cache )
        end select
        ! store adotoa:
        adotoa   = EV%eft_cache%adotoa
        ! compute massive neutrinos stuff:
        ! Massive neutrinos mod:
        EV%eft_cache%grhonu_tot    = 0._dl
        EV%eft_cache%gpinu_tot     = 0._dl
        EV%eft_cache%grhonudot_tot = 0._dl
        EV%eft_cache%gpinudot_tot  = 0._dl
        if ( CP%Num_Nu_Massive /= 0 ) then
            do nu_i = 1, CP%Nu_mass_eigenstates
                EFT_grhonu    = 0._dl
                EFT_gpinu     = 0._dl
                EFT_grhonudot = 0._dl
                EFT_gpinudot  = 0._dl
                EFT_gpinudotdot = 0._dl
                grhormass_t=State%grhormass(nu_i)/a**2
                call ThermalNuBack%rho_P(a*State%nu_masses(nu_i),EFT_grhonu,EFT_gpinu)
                EV%eft_cache%grhonu_tot = EV%eft_cache%grhonu_tot + grhormass_t*EFT_grhonu
                EV%eft_cache%gpinu_tot  = EV%eft_cache%gpinu_tot  + grhormass_t*EFT_gpinu
                EV%eft_cache%grhonudot_tot = EV%eft_cache%grhonudot_tot + grhormass_t*(ThermalNuBack%drho(a*State%nu_masses(nu_i) ,adotoa)&
                    & -4._dl*adotoa*EFT_grhonu)
                EV%eft_cache%gpinudot_tot  = EV%eft_cache%gpinudot_tot  + grhormass_t*(ThermalNuBack%pidot(a*State%nu_masses(nu_i),adotoa, EFT_gpinu )&
                    & -4._dl*adotoa*EFT_gpinu)
                EV%eft_cache%gpinudotdot_tot = 0._dl
            end do
        end if
        ! compute pressure dot:
        EV%eft_cache%gpresdotm_t = -4._dl*adotoa*( grhog_t+grhor_t )/3._dl +EV%eft_cache%gpinudot_tot
        ! compute remaining quantities related to H:
        call CP%EFTCAMB%model%compute_H_derivs( a, CP%eft_par_cache , EV%eft_cache )
        ! compute second derivative of massive neutrino pressure
        EV%eft_cache%grhonu_tot    = 0._dl
        EV%eft_cache%gpinu_tot     = 0._dl
        EV%eft_cache%grhonudot_tot = 0._dl
        EV%eft_cache%gpinudot_tot  = 0._dl
        EV%eft_cache%gpinudotdot_tot = 0._dl
        if ( CP%Num_Nu_Massive /= 0 ) then
            do nu_i = 1, CP%Nu_mass_eigenstates
                EFT_grhonu    = 0._dl
                EFT_gpinu     = 0._dl
                EFT_grhonudot = 0._dl
                EFT_gpinudot  = 0._dl
                grhormass_t=State%grhormass(nu_i)/a**2
                call ThermalNuBack%rho_P(a*State%nu_masses(nu_i),EFT_grhonu,EFT_gpinu)
                EV%eft_cache%grhonu_tot = EV%eft_cache%grhonu_tot + grhormass_t*EFT_grhonu
                EV%eft_cache%gpinu_tot  = EV%eft_cache%gpinu_tot  + grhormass_t*EFT_gpinu
                EFT_gpinudot = ThermalNuBack%pidot(a*State%nu_masses(nu_i),adotoa,EFT_gpinu)
                EV%eft_cache%gpinudot_tot  = EV%eft_cache%gpinudot_tot  + grhormass_t*(EFT_gpinudot&
                    & -4._dl*adotoa*EFT_gpinu)
                EV%eft_cache%grhonudot_tot = EV%eft_cache%grhonudot_tot + grhormass_t*(ThermalNuBack%drho(a*State%nu_masses(nu_i) ,adotoa)&
                    & -4._dl*adotoa*EFT_grhonu)
                EV%eft_cache%gpinudotdot_tot = EV%eft_cache%gpinudotdot_tot - 4._dl*adotoa*grhormass_t*&
                (EFT_gpinudot-4._dl*adotoa*EFT_gpinu) &
                + grhormass_t*(ThermalNuBack%pidotdot(a*State%nu_masses(nu_i),adotoa,EV%eft_cache%Hdot,EFT_gpinu,EFT_gpinudot)-4._dl*EV%eft_cache%Hdot*EFT_gpinu &
                & -4._dl*adotoa*EFT_gpinudot)
            end do
        end if
        ! compute pressure ddot:
        EV%eft_cache%gpresdotdotm_t = -(4._dl/3._dl)*EV%eft_cache%Hdot*(grhog_t+grhor_t)+(16._dl/3._dl)*(grhog_t+grhor_t)*adotoa**2 + EV%eft_cache%gpinudotdot_tot
        ! compute remaining quantities related to H:
        call CP%EFTCAMB%model%compute_H_derivs( a, CP%eft_par_cache , EV%eft_cache )
        ! store:
        adotdota = EV%eft_cache%Hdot +EV%eft_cache%adotoa**2
        ! compute backgrond EFT functions if model is designer:
        select type ( model => CP%EFTCAMB%model )
            class is ( EFTCAMB_designer_model )
            call CP%EFTCAMB%model%compute_background_EFT_functions( a, CP%eft_par_cache , EV%eft_cache )
        end select
        ! compute all other background stuff:
        call CP%EFTCAMB%model%compute_rhoQPQ( a, CP%eft_par_cache , EV%eft_cache )
        ! compute second order EFT functions:
        call CP%EFTCAMB%model%compute_secondorder_EFT_functions( a, CP%eft_par_cache , EV%eft_cache )
        !
        cothxor=1._dl/tau

        ! define total density perturbations:
        dgrho = dgrho_matter

        ! compute z and zdot in GR:
        z  = (0.5_dl*dgrho/k + etak)/adotoa
        dz = -2._dl*adotoa*z + etak

        if ( EV%no_nu_multpoles ) then
            !RSA approximation of arXiv:1104.2933, dropping opactity terms in the velocity
            !Approximate total density variables with just matter terms
            clxr = -4*dz/k
            qr   = -4._dl/3*z
        else
            !  Massless neutrinos
            clxr = y(EV%r_ix)
            qr   = y(EV%r_ix+1)
        endif

        if ( EV%no_phot_multpoles ) then
            if ( .not. EV%no_nu_multpoles ) then
                clxg = -4*dz/k -4/k*opacity*( vb +z )
                qg   = -4._dl/3*z
            else
                clxg = clxr -4/k*opacity*( vb +z )
                qg   = qr
            end if
        else
            !  Photons
            clxg = y(EV%g_ix)
            qg   = y(EV%g_ix+1)
        end if

        ! add radiation to total density and velocity perturbations:
        dgrho = dgrho +grhog_t*clxg +grhor_t*clxr ! 8\pi G_N \sum(rho_i*clx_i) a^2   : total density perturbation
        dgq   = dgq   +grhog_t*qg   +grhor_t*qr   ! 8\pi G_N \sum(rho_i+p_i)*v_i a^2 : total velocity perturbation

        ! recompute z and zdot in GR:
        z     = (0.5_dl*dgrho/k + etak)/adotoa
        dz    = -2._dl*adotoa*z + etak

        ! complete the cache:
        EV%eft_cache%z     = z
        EV%eft_cache%dz    = dz
        EV%eft_cache%clxg  = clxg
        EV%eft_cache%clxr  = clxr
        EV%eft_cache%dgpnu = dgpnu
        EV%eft_cache%dgrho = dgrho
        EV%eft_cache%dgq   = dgq

        ! compute pi field equation factors:
        call CP%EFTCAMB%model%compute_pi_factors( a, CP%eft_par_cache , EV%eft_cache )

        Omega0  = EV%eft_cache%EFTOmegaV
        Gamma10 = EV%eft_cache%EFTGamma1V
        Gamma20 = EV%eft_cache%EFTGamma2V
        Gamma30 = EV%eft_cache%EFTGamma3V

        if ( EV%eft_cache%EFTGamma1V /= 0._dl .or. EV%eft_cache%EFTGamma2V /= 0._dl .or. &
            &EV%eft_cache%EFTGamma3V /= 0._dl .or. EV%eft_cache%EFTGamma4V /= 0._dl .or. &
            &EV%eft_cache%EFTGamma5V /= 0._dl .or. EV%eft_cache%EFTGamma6V /= 0._dl ) then

            EFTpiEdot = 0._dl

        else

            EFTpiEdot = ( EV%eft_cache%EFTcdot +0.75_dl*(2._dl*EV%eft_cache%adotoa*EV%eft_cache%Hdot*a*a*EV%eft_cache%EFTOmegaP**2/(1._dl+EV%eft_cache%EFTOmegaV)&
                & +2._dl*EV%eft_cache%adotoa**3*a**3*EV%eft_cache%EFTOmegaP*EV%eft_cache%EFTOmegaPP/(1._dl + EV%eft_cache%EFTOmegaV)-a**3*EV%eft_cache%adotoa**3*EV%eft_cache%EFTOmegaP**3/(1._dl+EV%eft_cache%EFTOmegaV)**2))*k*EV%eft_cache%z &
                & +(EV%eft_cache%EFTc+0.75_dl*EV%eft_cache%adotoa**2*(a*a*EV%eft_cache%EFTOmegaP**2)/(1._dl+EV%eft_cache%EFTOmegaV))*k*EV%eft_cache%dz +2._dl*EV%eft_cache%adotoa*EV%eft_cache%EFTpiE&
                & +(-EV%eft_cache%dgrho +(EV%eft_cache%grhog_t*EV%eft_cache%clxg+EV%eft_cache%grhor_t*EV%eft_cache%clxr))/4._dl*(+a*EV%eft_cache%adotoa**2*EV%eft_cache%EFTOmegaP +a*EV%eft_cache%Hdot*EV%eft_cache%EFTOmegaP+a*a*EV%eft_cache%adotoa**2*EV%eft_cache%EFTOmegaPP &
                & -a*a*EV%eft_cache%adotoa**2*EV%eft_cache%EFTOmegaP**2/(1._dl+EV%eft_cache%EFTOmegaV))/(1._dl+EV%eft_cache%EFTOmegaV)&
                & -EV%eft_cache%adotoa/4._dl*a*EV%eft_cache%EFTOmegaP/(1._dl+EV%eft_cache%EFTOmegaV)*(EV%eft_cache%grhob_t*(-k*(EV%eft_cache%z +vb)-3._dl*EV%eft_cache%adotoa*clxb) + EV%eft_cache%grhoc_t*(-k*EV%eft_cache%z -3._dl*EV%eft_cache%adotoa*clxc))

        end if

        if ( CP%EFTCAMB%EFTCAMB_evolve_delta_phi ) then
            ! for now use zero IC for delta phi
            y(EV%w_ix)   = 0._dl
            y(EV%w_ix+1) = 0._dl
            ! this is h_prime = 2kz
            if ( CP%EFTCAMB%EFTCAMB_evolve_metric_h ) y(EV%w_ix+2) = 2._dl * k * z
        else
            select case (CP%EFTCAMB%EFT_IC_type)
            ! --- CASE 1: Standard EFTCAMB / GR ---
            case (1)
                if ( EV%eft_cache%EFTpiC +k2*EV%eft_cache%EFTpiD1 +k2**2*EV%eft_cache%EFTpiD2 /= 0._dl ) then
                    y(EV%w_ix)   = -CP%eft_par_cache%h0_Mpc*EV%eft_cache%EFTpiE/( EV%eft_cache%EFTpiC +k2*EV%eft_cache%EFTpiD1 +k2**2*EV%eft_cache%EFTpiD2 )
                    y(EV%w_ix+1) = CP%eft_par_cache%h0_Mpc*(EV%eft_cache%EFTpiE/( EV%eft_cache%EFTpiC +k2*EV%eft_cache%EFTpiD1 +k2**2*EV%eft_cache%EFTpiD2  )**2*a*EV%eft_cache%adotoa*(EFTpiCdotFunction(a,k)+k*k*EFTpiDdotFunction(a,k))&
                    &  -EFTpiEdot/(EV%eft_cache%EFTpiC +k2*EV%eft_cache%EFTpiD1 +k2**2*EV%eft_cache%EFTpiD2 ) )
                else
                    y(EV%w_ix)   = 0._dl
                    y(EV%w_ix+1) = 0._dl
                end if
            case (2)
                ! --- CASE 2: Constant Omega ---
                y(EV%w_ix)   = CP%eft_par_cache%h0_Mpc*(1._dl/4._dl)*((1._dl+EV%eft_cache%EFTOmegaV)/(5._dl+3._dl*EV%eft_cache%EFTOmegaV))*k2*(tau)**3
                y(EV%w_ix+1) =  CP%eft_par_cache%h0_Mpc*(3._dl/4._dl)*((1._dl+EV%eft_cache%EFTOmegaV)/(5._dl+3._dl*EV%eft_cache%EFTOmegaV))*k2*(tau)**2

            case (3)
                ! --- CASE 3: Constant Gamma + Constant Omega ---
                
                ! -- Pi Field (y) --
                y(EV%w_ix) =  ((k2 * EV%eft_cache%tau**3 * CP%eft_par_cache%h0_Mpc * &
                    (3._dl*Gamma30**2 + 2._dl*Omega0*(1._dl+Omega0) + Gamma30*(4._dl+5._dl*Omega0))) / &
                    (8._dl*(3._dl*Gamma30*(1._dl+Omega0) + Omega0*(5._dl+3._dl*Omega0)))) - &
                    ((9._dl * k2 * EV%eft_cache%tau**4 * (CP%eft_par_cache%h0_Mpc)**2 * Gamma20 * (Gamma30+Omega0) * &
                    (5._dl+3._dl*Gamma30+3._dl*Omega0) * &
                    (4._dl*Gamma30 + 3._dl*Gamma30**2 + 2._dl*Omega0 + 5._dl*Gamma30*Omega0 + 2._dl*Omega0**2)) / &
                    (160._dl*(1._dl+Gamma30+Omega0) * &
                    (3._dl*Gamma30+5._dl*Omega0+3._dl*Gamma30*Omega0+3._dl*Omega0**2) * &
                    (3._dl*Gamma30+8._dl*Omega0+6._dl*Gamma30*Omega0+6._dl*Omega0**2)))

                ! -- Pi Velocity (y+1) --
                y(EV%w_ix+1) =  3._dl*((k2 * EV%eft_cache%tau**2 * CP%eft_par_cache%h0_Mpc * &
                    (3._dl*Gamma30**2 + 2._dl*Omega0*(1._dl+Omega0) + Gamma30*(4._dl+5._dl*Omega0))) / &
                    (8._dl*(3._dl*Gamma30*(1._dl+Omega0) + Omega0*(5._dl+3._dl*Omega0)))) - &
                    4._dl*((9._dl * k2 * EV%eft_cache%tau**3 * (CP%eft_par_cache%h0_Mpc)**2 * Gamma20 * (Gamma30+Omega0) * &
                    (5._dl+3._dl*Gamma30+3._dl*Omega0) * &
                    (4._dl*Gamma30 + 3._dl*Gamma30**2 + 2._dl*Omega0 + 5._dl*Gamma30*Omega0 + 2._dl*Omega0**2)) / &
                    (160._dl*(1._dl+Gamma30+Omega0) * &
                    (3._dl*Gamma30+5._dl*Omega0+3._dl*Gamma30*Omega0+3._dl*Omega0**2) * &
                    (3._dl*Gamma30+8._dl*Omega0+6._dl*Gamma30*Omega0+6._dl*Omega0**2)))
                
            case default
                y(EV%w_ix)   = 0._dl
                y(EV%w_ix+1) = 0._dl
            end select
        end if

        ! reset the cache:
        call EV%eft_cache%initialize()

    end subroutine EFTCAMBInitialConditions

    !----------------------------------------------------------------------------------------
    !> Function that computes the pi field coefficient C at a given time and scale.
    function EFTpiCfunction( a, k )

        implicit none

        real(dl), intent(in) :: a               !< input scale factor
        real(dl), intent(in) :: k               !< input scale
        real(dl)             :: EFTpiCfunction  !< return value of the function

        type(TEFTCAMB_timestep_cache) :: eft_cache
        integer  :: nu_i
        real(dl) :: grhormass_t, EFT_gpinudot, EFT_grhonudot, EFT_gpinu, EFT_grhonu
        real(dl) :: adotoa, adotdota, grho, grhob_t, grhoc_t, grhor_t, grhog_t, gpres, grho_matter, a2

        ! initialize the cache:
        call eft_cache%initialize()

        ! prepare:
        a2 = a*a

        ! compute background densities of different species
        grhob_t = State%grhob/a         ! 8\pi G_N \rho_b a^2: baryons background density
        grhoc_t = State%grhoc/a         ! 8\pi G_N \rho_{cdm} a^2: cold dark matter background density
        grhor_t = State%grhornomass/a2  ! 8\pi G_N \rho_{\nu} a^2: massless neutrinos background density
        grhog_t = State%grhog/a2        ! 8\pi G_N \rho_{\gamma} a^2: radiation background density

        ! start computing background total pressure and total density:
        grho_matter  = grhob_t +grhoc_t
        gpres        = (grhog_t+grhor_t)/3._dl

        ! add radiation, massless neutrinos and Lambda to total background density:
        grho = grho_matter +grhor_t +grhog_t
        ! Massive neutrinos mod
        eft_cache%grhonu_tot = 0._dl
        eft_cache%gpinu_tot = 0._dl
        if ( CP%Num_Nu_Massive /= 0 ) then
            do nu_i = 1, CP%Nu_mass_eigenstates
                EFT_grhonu    = 0._dl
                EFT_gpinu     = 0._dl
                grhormass_t=State%grhormass(nu_i)/a**2
                call ThermalNuBack%rho_P(a*State%nu_masses(nu_i),EFT_grhonu,EFT_gpinu)
                eft_cache%grhonu_tot = eft_cache%grhonu_tot + grhormass_t*EFT_grhonu
                eft_cache%gpinu_tot  = eft_cache%gpinu_tot  + grhormass_t*EFT_gpinu
            end do
        end if

        grho  = grho  +eft_cache%grhonu_tot
        gpres = gpres +eft_cache%gpinu_tot

        ! start to fill the cache:
        eft_cache%a           = a
        eft_cache%k           = k
        eft_cache%grhom_t     = grho
        eft_cache%gpresm_t    = gpres
        eft_cache%grhob_t     = grhob_t
        eft_cache%grhoc_t     = grhoc_t
        eft_cache%grhor_t     = grhor_t
        eft_cache%grhog_t     = grhog_t
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
        eft_cache%grhonudot_tot = 0._dl
        eft_cache%gpinudot_tot = 0._dl
        if ( CP%Num_Nu_Massive /= 0 ) then
            do nu_i = 1, CP%Nu_mass_eigenstates
                EFT_grhonu    = 0._dl
                EFT_gpinu     = 0._dl
                EFT_grhonudot = 0._dl
                EFT_gpinudot  = 0._dl
                grhormass_t=State%grhormass(nu_i)/a**2
                call ThermalNuBack%rho_P(a*State%nu_masses(nu_i),EFT_grhonu,EFT_gpinu)
                eft_cache%grhonu_tot = eft_cache%grhonu_tot + grhormass_t*EFT_grhonu
                eft_cache%gpinu_tot  = eft_cache%gpinu_tot  + grhormass_t*EFT_gpinu
                eft_cache%grhonudot_tot = eft_cache%grhonudot_tot + grhormass_t*(ThermalNuBack%drho(a*State%nu_masses(nu_i) ,adotoa)&
                    & -4._dl*adotoa*EFT_grhonu)
                eft_cache%gpinudot_tot  = eft_cache%gpinudot_tot  + grhormass_t*(ThermalNuBack%pidot(a*State%nu_masses(nu_i),adotoa, EFT_gpinu )&
                    & -4._dl*adotoa*EFT_gpinu)
            end do
        end if
        ! compute pressure dot:
        eft_cache%gpresdotm_t = -4._dl*adotoa*( grhog_t+grhor_t )/3._dl +eft_cache%gpinudot_tot
        ! compute remaining quantities related to H:
        call CP%EFTCAMB%model%compute_H_derivs( a, CP%eft_par_cache , eft_cache )
        ! store:
        adotdota = eft_cache%Hdot +eft_cache%adotoa**2
        ! compute backgrond EFT functions if model is designer:
        select type ( model => CP%EFTCAMB%model )
            class is ( EFTCAMB_designer_model )
            call CP%EFTCAMB%model%compute_background_EFT_functions( a, CP%eft_par_cache , eft_cache )
        end select
        ! compute all other background stuff:
        call CP%EFTCAMB%model%compute_rhoQPQ( a, CP%eft_par_cache , eft_cache )
        ! compute second order EFT functions:
        call CP%EFTCAMB%model%compute_secondorder_EFT_functions( a, CP%eft_par_cache , eft_cache )

        ! compute Einstein equations factors:
        call CP%EFTCAMB%model%compute_pi_factors( a, CP%eft_par_cache , eft_cache )

        EFTpiCfunction = eft_cache%EFTpiC

        return

    end function EFTpiCfunction

    !----------------------------------------------------------------------------------------
    !> Function that computes the pi field coefficient D at a given time and scale.
    function EFTpiDfunction(a,k)

        implicit none

        real(dl), intent(in) :: a               !< input value of the scale factor
        real(dl), intent(in) :: k               !< input scale
        real(dl)             :: EFTpiDfunction  !< return value of the function

        type(TEFTCAMB_timestep_cache) :: eft_cache

        integer  :: nu_i
        real(dl) :: grhormass_t, EFT_gpinudot, EFT_grhonudot, EFT_gpinu, EFT_grhonu
        real(dl) :: adotoa, adotdota, grho, grhob_t, grhoc_t, grhor_t, grhog_t, gpres, grho_matter, a2

        ! initialize the cache:
        call eft_cache%initialize()

        ! prepare:
        a2 = a*a

        ! compute background densities of different species
        grhob_t = State%grhob/a         ! 8\pi G_N \rho_b a^2: baryons background density
        grhoc_t = State%grhoc/a         ! 8\pi G_N \rho_{cdm} a^2: cold dark matter background density
        grhor_t = State%grhornomass/a2  ! 8\pi G_N \rho_{\nu} a^2: massless neutrinos background density
        grhog_t = State%grhog/a2        ! 8\pi G_N \rho_{\gamma} a^2: radiation background density

        ! start computing background total pressure and total density:
        grho_matter  = grhob_t +grhoc_t
        gpres        = (grhog_t+grhor_t)/3._dl

        ! add radiation, massless neutrinos and Lambda to total background density:
        grho = grho_matter +grhor_t +grhog_t
        ! Massive neutrinos mod:
        eft_cache%grhonu_tot = 0._dl
        eft_cache%gpinu_tot = 0._dl
        if ( CP%Num_Nu_Massive /= 0 ) then
            do nu_i = 1, CP%Nu_mass_eigenstates
                EFT_grhonu    = 0._dl
                EFT_gpinu     = 0._dl
                grhormass_t=State%grhormass(nu_i)/a**2
                call ThermalNuBack%rho_P(a*State%nu_masses(nu_i),EFT_grhonu,EFT_gpinu)
                eft_cache%grhonu_tot = eft_cache%grhonu_tot + grhormass_t*EFT_grhonu
                eft_cache%gpinu_tot  = eft_cache%gpinu_tot  + grhormass_t*EFT_gpinu
            end do
        end if

        grho  = grho  +eft_cache%grhonu_tot
        gpres = gpres +eft_cache%gpinu_tot

        ! start to fill the cache:
        eft_cache%a           = a
        eft_cache%k           = k
        eft_cache%grhom_t     = grho
        eft_cache%gpresm_t    = gpres
        eft_cache%grhob_t     = grhob_t
        eft_cache%grhoc_t     = grhoc_t
        eft_cache%grhor_t     = grhor_t
        eft_cache%grhog_t     = grhog_t
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
        eft_cache%grhonudot_tot = 0._dl
        eft_cache%gpinudot_tot = 0._dl
        if ( CP%Num_Nu_Massive /= 0 ) then
            do nu_i = 1, CP%Nu_mass_eigenstates
                EFT_grhonu    = 0._dl
                EFT_gpinu     = 0._dl
                EFT_grhonudot = 0._dl
                EFT_gpinudot  = 0._dl
                grhormass_t=State%grhormass(nu_i)/a**2
                call ThermalNuBack%rho_P(a*State%nu_masses(nu_i),EFT_grhonu,EFT_gpinu)
                eft_cache%grhonu_tot = eft_cache%grhonu_tot + grhormass_t*EFT_grhonu
                eft_cache%gpinu_tot  = eft_cache%gpinu_tot  + grhormass_t*EFT_gpinu
                eft_cache%grhonudot_tot = eft_cache%grhonudot_tot + grhormass_t*(ThermalNuBack%drho(a*State%nu_masses(nu_i) ,adotoa)&
                    & -4._dl*adotoa*EFT_grhonu)
                eft_cache%gpinudot_tot  = eft_cache%gpinudot_tot  + grhormass_t*(ThermalNuBack%pidot(a*State%nu_masses(nu_i),adotoa, EFT_gpinu )&
                    & -4._dl*adotoa*EFT_gpinu)
            end do
        end if
        ! compute pressure dot:
        eft_cache%gpresdotm_t = -4._dl*adotoa*( grhog_t+grhor_t )/3._dl +eft_cache%gpinudot_tot
        ! compute remaining quantities related to H:
        call CP%EFTCAMB%model%compute_H_derivs( a, CP%eft_par_cache , eft_cache )
        ! store:
        adotdota = eft_cache%Hdot +eft_cache%adotoa**2
        ! compute backgrond EFT functions if model is designer:
        select type ( model => CP%EFTCAMB%model )
            class is ( EFTCAMB_designer_model )
            call CP%EFTCAMB%model%compute_background_EFT_functions( a, CP%eft_par_cache , eft_cache )
        end select
        ! compute all other background stuff:
        call CP%EFTCAMB%model%compute_rhoQPQ( a, CP%eft_par_cache , eft_cache )
        ! compute second order EFT functions:
        call CP%EFTCAMB%model%compute_secondorder_EFT_functions( a, CP%eft_par_cache , eft_cache )

        ! compute Einstein equations factors:
        call CP%EFTCAMB%model%compute_pi_factors( a, CP%eft_par_cache , eft_cache )

        EFTpiDfunction = eft_cache%EFTpiD1 + k*k*eft_cache%EFTpiD1

        return

    end function EFTpiDfunction

    !----------------------------------------------------------------------------------------
    !> Function that computes the numerical derivative of the pi field coefficient C
    !! at a given time and scale.
    function EFTpiCdotFunction(a,k)

        implicit none

        real(dl), intent(in) :: a                  !< input value of the scale factor
        real(dl), intent(in) :: k                  !< input scale
        real(dl)             :: EFTpiCdotFunction  !< return value of the function

        real(dl) err

        EFTpiCdotfunction = dfridr( dummy_helper, a, 0.03_dl*a, err )

        return

    contains

        function dummy_helper(temp_a)
            real(dl) temp_a,dummy_helper
            dummy_helper = EFTpiCfunction(temp_a,k)
        end function dummy_helper

    end function EFTpiCdotFunction

    !----------------------------------------------------------------------------------------
    !> Function that computes the numerical derivative of the pi field coefficient D
    !! at a given time and scale.
    function EFTpiDdotFunction(a,k)

        implicit none

        real(dl), intent(in) :: a                  !< input value of the scale factor
        real(dl), intent(in) :: k                  !< input scale
        real(dl)             :: EFTpiDdotFunction  !< return value of the function


        real(dl) err

        EFTpiDdotfunction = dfridr( dummy_helper, a, 0.03_dl*a, err )

        return

    contains

        function dummy_helper(temp_a)
            real(dl) temp_a,dummy_helper
            dummy_helper = EFTpiDfunction(temp_a,k)
        end function dummy_helper

    end function EFTpiDdotFunction

    !----------------------------------------------------------------------------------------

end submodule EFTCAMB_IC

!----------------------------------------------------------------------------------------
