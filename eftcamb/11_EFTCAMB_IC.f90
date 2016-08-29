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

!> @file 11_EFTCAMB_IC.f90
!! This file contains the perturbation initial conditions calculator.


!----------------------------------------------------------------------------------------
!> This module contains the perturbation initial conditions calculator.

!> @author Bin Hu, Marco Raveri

submodule (GaugeInterface) EFTCAMB_IC

use precision

implicit none

contains

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that initializes EFTCAMB perturbations initial conditions.
    module subroutine EFTCAMBInitialConditions( y, EV, tau )

        use ThermoData

        implicit none

        type(EvolutionVars) EV  !< CAMB evolutionary variables.
        real(dl) :: y(EV%nvar)  !< status vector with the values of the perturbations.
        real(dl) :: tau         !< conformal time at which the initial conditions are desired.

        real(dl) :: a, a2, k, k2, etak, grhob_t, grhoc_t, grhor_t, grhog_t, grhov_t, gpres, grho_matter
        real(dl) :: dgrho_matter, dgq, dgpnu, grho, adotoa, EFT_grhonu, EFT_gpinu, EFT_grhonudot
        real(dl) :: EFT_gpinudot, grhormass_t, adotdota, cothxor, cs2, dopacity, dgrho, z, dz
        real(dl) :: clxr, qr, pir, clxg, qg, pig, opacity, vb, clxc, clxb, EFTpiEdot
        real(dl) :: wnu_arr(max_nu)

        integer  :: nu_i

        k    = EV%k_buf
        k2   = EV%k2_buf
        ! copy variables:
        a    = y(1) ! scale factor
        etak = y(2) ! conformal time
        clxc = y(3) ! cdm density perturbations
        clxb = y(4) ! baryons density perturbations
        vb   = y(5) ! baryons velocity perturbations
        ! process:
        a2   = a*a
        ! compute background densities of different species
        grhob_t = grhob/a         ! 8\pi G_N \rho_b a^2: baryons background density
        grhoc_t = grhoc/a         ! 8\pi G_N \rho_{cdm} a^2: cold dark matter background density
        grhor_t = grhornomass/a2  ! 8\pi G_N \rho_{\nu} a^2: massless neutrinos background density
        grhog_t = grhog/a2        ! 8\pi G_N \rho_{\gamma} a^2: radiation background density
        grhov_t = 0._dl

        ! start computing background total pressure and total density:
        gpres        = 0._dl
        grho_matter  = grhob_t +grhoc_t
        ! start computing total density and velocity perturbations:
        dgrho_matter = grhob_t*clxb +grhoc_t*clxc ! 8\pi G_N \sum(rho_i*clx_i) a^2   : total density perturbation
        dgq          = grhob_t*vb                 ! 8\pi G_N \sum(rho_i+p_i)*v_i a^2 : total velocity perturbation
        ! add massive neutrinos to both background and perturbations:

        dgpnu = 0._dl
        if (CP%Num_Nu_Massive > 0) then
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
        if ( CP%Num_Nu_Massive /= 0 ) then
            do nu_i = 1, CP%Nu_mass_eigenstates
                EFT_grhonu    = 0._dl
                EFT_gpinu     = 0._dl
                EFT_grhonudot = 0._dl
                EFT_gpinudot  = 0._dl
                grhormass_t=grhormass(nu_i)/a**2
                call Nu_background(a*nu_masses(nu_i),EFT_grhonu,EFT_gpinu)
                EV%eft_cache%grhonu_tot = EV%eft_cache%grhonu_tot + grhormass_t*EFT_grhonu
                EV%eft_cache%gpinu_tot  = EV%eft_cache%gpinu_tot  + grhormass_t*EFT_gpinu
                EV%eft_cache%grhonudot_tot = EV%eft_cache%grhonudot_tot + grhormass_t*(Nu_drho(a*nu_masses(nu_i) ,adotoa, EFT_grhonu)&
                    & -4._dl*adotoa*EFT_grhonu)
                EV%eft_cache%gpinudot_tot  = EV%eft_cache%gpinudot_tot  + grhormass_t*(Nu_pidot(a*nu_masses(nu_i),adotoa, EFT_gpinu )&
                    & -4._dl*adotoa*EFT_gpinu)
            end do
        end if
        ! compute pressure dot:
        EV%eft_cache%gpresdotm_t = -4._dl*adotoa*( grhog_t+grhor_t )/3._dl +EV%eft_cache%gpinudot_tot
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

        !  Get sound speed and ionisation fraction.
        if ( EV%TightCoupling ) then
            call thermo( tau, cs2, opacity, dopacity )
        else
            call thermo( tau, cs2, opacity )
        end if

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

        if ( EV%eft_cache%EFTpiC +k2*EV%eft_cache%EFTpiD1 +k2**2*EV%eft_cache%EFTpiD2 /= 0._dl ) then
            y(EV%w_ix)   = -CP%eft_par_cache%h0_Mpc*EV%eft_cache%EFTpiE/( EV%eft_cache%EFTpiC +k2*EV%eft_cache%EFTpiD1 +k2**2*EV%eft_cache%EFTpiD2 )
            y(EV%w_ix+1) = CP%eft_par_cache%h0_Mpc*(EV%eft_cache%EFTpiE/( EV%eft_cache%EFTpiC +k2*EV%eft_cache%EFTpiD1 +k2**2*EV%eft_cache%EFTpiD2  )**2*a*EV%eft_cache%adotoa*(EFTpiCdotFunction(a,k)+k*k*EFTpiDdotFunction(a,k))&
                &  -EFTpiEdot/(EV%eft_cache%EFTpiC +k2*EV%eft_cache%EFTpiD1 +k2**2*EV%eft_cache%EFTpiD2 ) )
        else
            y(EV%w_ix)   = 0._dl
            y(EV%w_ix+1) = 0._dl
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

        type(EFTCAMB_timestep_cache) :: eft_cache
        integer  :: nu_i
        real(dl) :: grhormass_t, EFT_gpinudot, EFT_grhonudot, EFT_gpinu, EFT_grhonu
        real(dl) :: adotoa, adotdota, grho, grhob_t, grhoc_t, grhor_t, grhog_t, gpres, grho_matter, a2

        ! initialize the cache:
        call eft_cache%initialize()

        ! prepare:
        a2 = a*a

        ! compute background densities of different species
        grhob_t = grhob/a         ! 8\pi G_N \rho_b a^2: baryons background density
        grhoc_t = grhoc/a         ! 8\pi G_N \rho_{cdm} a^2: cold dark matter background density
        grhor_t = grhornomass/a2  ! 8\pi G_N \rho_{\nu} a^2: massless neutrinos background density
        grhog_t = grhog/a2        ! 8\pi G_N \rho_{\gamma} a^2: radiation background density

        ! start computing background total pressure and total density:
        grho_matter  = grhob_t +grhoc_t
        gpres        = (grhog_t+grhor_t)/3._dl

        ! add radiation, massless neutrinos and Lambda to total background density:
        grho = grho_matter +grhor_t +grhog_t

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

        type(EFTCAMB_timestep_cache) :: eft_cache

        integer  :: nu_i
        real(dl) :: grhormass_t, EFT_gpinudot, EFT_grhonudot, EFT_gpinu, EFT_grhonu
        real(dl) :: adotoa, adotdota, grho, grhob_t, grhoc_t, grhor_t, grhog_t, gpres, grho_matter, a2

        ! initialize the cache:
        call eft_cache%initialize()

        ! prepare:
        a2 = a*a

        ! compute background densities of different species
        grhob_t = grhob/a         ! 8\pi G_N \rho_b a^2: baryons background density
        grhoc_t = grhoc/a         ! 8\pi G_N \rho_{cdm} a^2: cold dark matter background density
        grhor_t = grhornomass/a2  ! 8\pi G_N \rho_{\nu} a^2: massless neutrinos background density
        grhog_t = grhog/a2        ! 8\pi G_N \rho_{\gamma} a^2: radiation background density

        ! start computing background total pressure and total density:
        grho_matter  = grhob_t +grhoc_t
        gpres        = (grhog_t+grhor_t)/3._dl

        ! add radiation, massless neutrinos and Lambda to total background density:
        grho = grho_matter +grhor_t +grhog_t

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

        EFTpiCdotfunction = dfridr( EFTpiCfunction, a, k, 0.03_dl*a, err )

        return

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

        EFTpiDdotfunction = dfridr( EFTpiDfunction, a, k, 0.03_dl*a, err )

        return

    end function EFTpiDdotFunction

    !----------------------------------------------------------------------------------------
    !> Algorithm to compute the numerical derivative.
    !! Modified to accept two function arguments.
    function dfridr( func, x, k, h, err )

        implicit none

        integer,parameter :: ntab = 100
        real(dl) dfridr,err,h,x,k,func
        external func

        real(dl), parameter :: CON=1.4_dl       ! decrease of the stepsize.
        real(dl), parameter :: CON2=CON*CON
        real(dl), parameter :: BIG=1.d+30
        real(dl), parameter :: SAFE=2._dl

        integer i,j
        real(dl) errt, fac, hh
        real(dl), dimension(ntab,ntab) :: a

        if (h.eq.0._dl) h = 1.d-8

        hh=h
        a(1,1)=(func(x+hh,k)-func(x-hh,k))/(2.0_dl*hh)
        err=BIG

        do 12 i=2,NTAB
            hh=hh/CON
            a(1,i)=(func(x+hh,k)-func(x-hh,k))/(2.0_dl*hh)
            fac=CON2
            do 11 j=2,i
                a(j,i)=(a(j-1,i)*fac-a(j-1,i-1))/(fac-1._dl)
                fac=CON2*fac
                errt=max(abs(a(j,i)-a(j-1,i)),abs(a(j,i)-a(j-1,i-1)))
                if (errt.le.err) then
                    err=errt
                    dfridr=a(j,i)
                endif
11          continue
            if(abs(a(i,i)-a(i-1,i-1)).ge.SAFE*err)return
12      continue
        return
    end function dfridr

    !----------------------------------------------------------------------------------------

end submodule EFTCAMB_IC

!----------------------------------------------------------------------------------------
