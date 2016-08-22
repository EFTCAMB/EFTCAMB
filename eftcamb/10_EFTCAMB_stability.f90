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

!> @file 10_EFTCAMB_stability.f90
!! This file contains the stability detection algorithm of EFTCAMB.
!! This operates on general EFT models but needs CP to be set.


!----------------------------------------------------------------------------------------
!> This module contains the stability detection algorithm of EFTCAMB.
!! This operates on general EFT models but needs CP to be set.

!> @author Bin Hu, Marco Raveri

module EFTCAMB_stability

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
    !> Function that computes if the stability requirements are fullfilled in a given time.
    function EFTStabilityComputation(a)

        implicit none

        real(dl), intent(in) :: a                       !< Input scale factor
        logical              :: EFTStabilityComputation !< Logical value returned by the function. If the model is stable this is True, otherwise False.


        ! 1) Definitions of variables:
        logical :: EFT_HaveNan
        real(dl) :: a2,k,k2,grhob_t,grhor_t,grhoc_t,grhog_t,grhov_t, grho, gpres, EFT_H0
        real(dl) :: adotoa, adotdota, Hdot, Hdotdot
        real(dl) :: EFTOmegaV, EFTOmegaP,EFTOmegaPP,EFTOmegaPPP, EFTc, EFTLambda, EFTcdot, EFTLambdadot
        real(dl) :: EFTGamma1V, EFTGamma1P, EFTGamma2V, EFTGamma2P, EFTGamma3V, EFTGamma3P
        real(dl) :: EFTGamma4V, EFTGamma4P, EFTGamma5V, EFTGamma5P, EFTGamma6V, EFTGamma6P
        real(dl) :: EFTgrhoq,EFTgpresq,EFTgrhodotq,EFTgpresdotq
        real(dl) :: EFTpiA1, EFTpiA2, EFTpiB1, EFTpiB2, EFTpiC, EFTpiD, EFTpiD1, EFTpiD2, EFTAT
        real(dl) :: EFTtemp_H0,  EFTtemp_Omegac, EFTtemp_Omegab, EFTtemp_Omegav, EFTtemp_TCMB, EFTtemp_nu_massless_degeneracy
        real(dl) :: EFTtemp_grhom, EFTtemp_grhog, EFTtemp_grhor
        real(dl) :: EFT_grhonu, EFT_gpinu, grhormass_t
        real(dl) :: EFT_grhonu_tot, EFT_gpinu_tot, EFT_gpinudot_tot, grho_matter, gpres_matter
        real(dl) :: kmax
        real(dl) :: temp1, temp2, temp3, temp4, temp5, tempk
        integer  :: nu_i, ind, ind_max
        real(dl) :: EFT_instability_rate
        real(dl) :: EFT_W0, EFT_W1, EFT_W2, EFT_W3, EFT_W6, EFT_W2P, EFT_W3P, EFT_W6P  ! New stability
        real(dl) :: EFT_kinetic, EFT_gradient                                          ! New stability

!        ! 1) Stability check initialization
!        EFTStabilityComputation = .true.
!
!        ! prepare:
!        a2 = a*a
!        ! initialize the cache:
!        call eft_cache%initialize()
!        ! start filling:
!        eft_cache%a = a
!        ! compute background densities of different species
!        eft_cache%grhob_t = grhob/a         ! 8\pi G_N \rho_b a^2: baryons background density
!        eft_cache%grhoc_t = grhoc/a         ! 8\pi G_N \rho_{cdm} a^2: cold dark matter background density
!        eft_cache%grhor_t = grhornomass/a2  ! 8\pi G_N \rho_{\nu} a^2: massless neutrinos background density
!        eft_cache%grhog_t = grhog/a2        ! 8\pi G_N \rho_{\gamma} a^2: radiation background density
!        ! Massive neutrinos terms:
!        if ( CP%Num_Nu_Massive /= 0 ) then
!            do nu_i = 1, CP%Nu_mass_eigenstates
!                EFT_grhonu    = 0._dl
!                EFT_gpinu     = 0._dl
!                grhormass_t=grhormass(nu_i)/a**2
!                call Nu_background(a*nu_masses(nu_i),EFT_grhonu,EFT_gpinu)
!                eft_cache%grhonu_tot = eft_cache%grhonu_tot + grhormass_t*EFT_grhonu
!                eft_cache%gpinu_tot  = eft_cache%gpinu_tot  + grhormass_t*EFT_gpinu
!            end do
!        end if
!        ! assemble total densities and pressure:
!        eft_cache%grhom_t  = eft_cache%grhob_t +eft_cache%grhoc_t +eft_cache%grhor_t +eft_cache%grhog_t +eft_cache%grhonu_tot
!        eft_cache%gpresm_t = (+eft_cache%grhor_t +eft_cache%grhog_t)/3._dl +eft_cache%gpinu_tot
!        ! compute the other things:
!        select type ( model => CP%EFTCAMB%model )
!            ! compute the background and the background EFT functions.
!            class is ( EFTCAMB_full_model )
!            ! background for full models. Here the expansion history is computed from the
!            ! EFT functions. Hence compute them first and then compute the expansion history.
!            call CP%EFTCAMB%model%compute_background_EFT_functions( a, CP%eft_par_cache , eft_cache )
!            call CP%EFTCAMB%model%compute_adotoa( a, CP%eft_par_cache , eft_cache )
!            class is ( EFTCAMB_designer_model )
!            ! background for designer models. Here the expansion history is parametrized
!            ! and does not depend on the EFT functions. Hence compute first the expansion history
!            ! and then the EFT functions.
!            call CP%EFTCAMB%model%compute_adotoa( a, CP%eft_par_cache , eft_cache )
!        end select
!        ! store adotoa:
!        adotoa   = eft_cache%adotoa
!        ! compute massive neutrinos stuff:
!        ! Massive neutrinos mod:
!        if ( CP%Num_Nu_Massive /= 0 ) then
!            do nu_i = 1, CP%Nu_mass_eigenstates
!                EFT_grhonu    = 0._dl
!                EFT_gpinu     = 0._dl
!                EFT_grhonudot = 0._dl
!                EFT_gpinudot  = 0._dl
!                grhormass_t=grhormass(nu_i)/a**2
!                call Nu_background(a*nu_masses(nu_i),EFT_grhonu,EFT_gpinu)
!                eft_cache%grhonu_tot = eft_cache%grhonu_tot + grhormass_t*EFT_grhonu
!                eft_cache%gpinu_tot  = eft_cache%gpinu_tot  + grhormass_t*EFT_gpinu
!                eft_cache%grhonudot_tot = eft_cache%grhonudot_tot + grhormass_t*(Nu_drho(a*nu_masses(nu_i) ,adotoa, EFT_grhonu)&
!                    & -4._dl*adotoa*EFT_grhonu)
!                eft_cache%gpinudot_tot  = eft_cache%gpinudot_tot  + grhormass_t*(Nu_pidot(a*nu_masses(nu_i),adotoa, EFT_gpinu )&
!                    & -4._dl*adotoa*EFT_gpinu)
!            end do
!        end if
!        ! compute pressure dot:
!        eft_cache%gpresdotm_t = -4._dl*adotoa*( eft_cache%grhog_t +eft_cache%grhor_t )/3._dl +eft_cache%gpinudot_tot
!        ! compute remaining quantities related to H:
!        call CP%EFTCAMB%model%compute_H_derivs( a, CP%eft_par_cache , eft_cache )
!        ! compute backgrond EFT functions if model is designer:
!        select type ( model => CP%EFTCAMB%model )
!            class is ( EFTCAMB_designer_model )
!            call CP%EFTCAMB%model%compute_background_EFT_functions( a, CP%eft_par_cache , eft_cache )
!        end select
!        ! compute all other background stuff:
!        call CP%EFTCAMB%model%compute_rhoQPQ( a, CP%eft_par_cache , eft_cache )
!        ! compute second order EFT functions:
!        call CP%EFTCAMB%model%compute_secondorder_EFT_functions( a, CP%eft_par_cache , eft_cache )



!        ! 6) Pi field equation coefficients:
!        k = 0._dl
!        k2 = 0._dl
!        ! Coefficient A of tensors:
!        EFTAT = 1._dl + EFTOmegaV - EFTGamma4V
!        ! First part of the A coefficient:
!        EFTpiA1 = EFTc +2._dl*a2*EFT_H0**2*EFTGamma1V +1.5_dl*a2*(adotoa*EFTOmegaP+EFT_H0*EFTGamma2V)**2&
!            &/(2._dl*(1+EFTOmegaV)+EFTGamma3V+EFTGamma4V)
!        ! Second part of the A coefficient, the k^2 term:
!        EFTpiA2 = +4._dl*EFTGamma6V
!        ! First part of the B coefficient:
!        EFTpiB1 = EFTcdot +4._dl*adotoa*EFTc +8._dl*a2*adotoa*EFT_H0**2*(EFTGamma1V+ 0.25_dl*a*EFTGamma1P)&
!            & -a*(adotoa*EFTOmegaP+EFT_H0*EFTGamma2V)/(4._dl*(1._dl+EFTOmegaV)+6._dl*EFTGamma3V +2._dl*EFTGamma4V)*&
!            &(-3._dl*(EFTgrhoQ + EFTgpresQ) -3._dl*a*adotoa**2*EFTOmegaP*(4._dl +Hdot/(adotoa**2)) -3._dl*a2*adotoa**2*EFTOmegaPP&
!            & -3._dl*a*adotoa*EFT_H0*(4._dl*EFTGamma2V + a*EFTGamma2P) -(9._dl*EFTGamma3V -3._dl*EFTGamma4V)*&
!            &(Hdot-adotoa**2))&
!            & +1._dl/(1._dl+EFTOmegaV+2._dl*EFTGamma5V)*(a*adotoa*EFTOmegaP+2._dl*adotoa*(EFTGamma5V + EFTGamma5P)&
!            & -(1._dl+EFTOmegaV)*(a*adotoa*EFTOmegaP+a*EFT_H0*EFTGamma2V)/(2._dl*(1._dl+EFTOmegaV) +3._dl*EFTGamma3V +EFTGamma4V))*&
!            &(-EFTc +1.5_dl*a*adotoa**2*EFTOmegaP -2._dl*a2*EFT_H0*EFTGamma1V +1.5_dl*a*adotoa*EFT_H0*EFTGamma2V)
!        ! Second part of the B coefficient, the k^2 term:
!        k2      = 1._dl
!        EFTpiB2 = +4._dl*k2*adotoa*(2._dl*EFTGamma6V+ a*EFTGamma6P) +a*k2*(EFTGamma4V&
!            & +2._dl*EFTGamma5V)/(2._dl*(1._dl+EFTOmegaV)-2._dl*EFTGamma4V)*(adotoa*EFTOmegaP+EFT_H0*EFTGamma2V)&
!            & -a*(adotoa*EFTOmegaP+EFT_H0*EFTGamma2V)/(4._dl*(1._dl+EFTOmegaV)+6._dl*EFTGamma3V +2._dl*EFTGamma4V)*&
!            &( +k2*(3._dl*EFTGamma3V -EFTGamma4V +4._dl*EFTGamma5V))&
!            & +1._dl/(1._dl+EFTOmegaV+2._dl*EFTGamma5V)*(a*adotoa*EFTOmegaP+2._dl*adotoa*(EFTGamma5V + EFTGamma5P)&
!            & -(1._dl+EFTOmegaV)*(a*adotoa*EFTOmegaP+a*EFT_H0*EFTGamma2V)/(2._dl*(1._dl+EFTOmegaV) +3._dl*EFTGamma3V +EFTGamma4V))*&
!            &(-4._dl*EFTGamma6V*k2)
!        ! The C coefficient:
!        EFTpiC = +adotoa*EFTcdot +(6._dl*adotoa**2-2._dl*Hdot)*EFTc +1.5_dl*a*adotoa*EFTOmegaP*(Hdotdot-2._dl*adotoa**3) &
!            & +6._dl*a2*adotoa**2*EFT_H0**2*EFTGamma1V +2._dl*a2*Hdot*EFT_H0**2*EFTGamma1V &
!            & +2._dl*a2*a*adotoa**2*EFT_H0**2*EFTGamma1P +1.5_dl*(Hdot-adotoa**2)**2*(EFTGamma4V +3._dl*EFTGamma3V )&
!            & +4.5_dl*adotoa*EFT_H0*a*(Hdot-adotoa**2)*(EFTGamma2V + a*EFTGamma2P/3._dl )&
!            & +0.5_dl*a*EFT_H0*EFTGamma2V*(3._dl*Hdotdot -12._dl*Hdot*adotoa +6._dl*adotoa**3) &
!            & -a*(adotoa*EFTOmegaP+EFT_H0*EFTGamma2V)/(4._dl*(1._dl+EFTOmegaV)+6._dl*EFTGamma3V +2._dl*EFTGamma4V)*&
!            &(-3._dl*EFTgpresdotQ -3._dl*adotoa*(EFTgrhoQ +EFTgpresQ)-3._dl*a*adotoa**3*(a*EFTOmegaPP +6._dl*EFTOmegaP) &
!            & -6._dl*a*adotoa*Hdot*EFTOmegaP +3._dl*(Hdotdot-2._dl*adotoa*Hdot)*(EFTGamma4V +3._dl*EFTGamma3V)&
!            & +6._dl*adotoa*(Hdot-adotoa**2)*(3._dl*EFTGamma3V +1.5_dl*a*EFTGamma3P +EFTGamma4V + 0.5_dl*a*EFTGamma4P)&
!            & -3._dl*a*EFT_H0*(3._dl*adotoa**2*EFTGamma2V +Hdot*EFTGamma2V +a*adotoa**2*EFTGamma2P))&
!            & +1._dl/(1._dl+EFTOmegaV+2._dl*EFTGamma5V)*(a*adotoa*EFTOmegaP+2._dl*adotoa*(EFTGamma5V +a*EFTGamma5P)&
!            & -(1._dl+EFTOmegaV)*(a*adotoa*EFTOmegaP+a*EFT_H0*EFTGamma2V)/(2._dl*(1._dl+EFTOmegaV)+3._dl*EFTGamma3V +EFTGamma4V))*&
!            &(-0.5*EFTgrhodotQ -adotoa*EFTc +1.5_dl*a*adotoa*EFTOmegaP*(3._dl*adotoa**2-Hdot) -2._dl*a2*adotoa*EFT_H0**2*EFTGamma1V&
!            & -1.5_dl*a*EFT_H0*EFTGamma2V*(Hdot-2._dl*adotoa**2) -3._dl*adotoa*(Hdot-adotoa**2)*(1.5_dl*EFTGamma3V +0.5_dl*EFTGamma4V))
!        ! First part of the D coefficient, the k^2 term:
!        EFTpiD1 = EFTc -0.5_dl*a*adotoa*EFT_H0*(EFTGamma2V +a*EFTGamma2P) +(adotoa**2-Hdot)*(3._dl*EFTGamma3V +EFTGamma4V)&
!            & +4._dl*(Hdot*EFTGamma6V + adotoa**2*EFTGamma6V + a*adotoa**2*EFTGamma6P)&
!            & +2._dl*(Hdot*EFTGamma5V +a*adotoa**2*EFTGamma5P)&
!            & -a*(adotoa*EFTOmegaP+EFT_H0*EFTGamma2V)/(4._dl*(1._dl+EFTOmegaV)+6._dl*EFTGamma3V +2._dl*EFTGamma4V)*&
!            &(-2._dl*a*adotoa*EFTOmegaP +4._dl*adotoa*EFTGamma5V -2._dl*adotoa*(3._dl*EFTGamma3V +1.5_dl*a*EFTGamma3P &
!            & +EFTGamma4V +0.5_dl*a*EFTGamma4P))&
!            & +1._dl/(1._dl+EFTOmegaV+2._dl*EFTGamma5V)*(a*adotoa*EFTOmegaP+2._dl*adotoa*(EFTGamma5V +a*EFTGamma5P)&
!            & -(1._dl+EFTOmegaV)*(a*adotoa*EFTOmegaP +a*EFT_H0*EFTGamma2V)/(2._dl*(1._dl+EFTOmegaV)+3._dl*EFTGamma3V +EFTGamma4V))*&
!            &(+0.5_dl*a*adotoa*EFTOmegaP -2._dl*adotoa*EFTGamma5V +0.5_dl*a*EFT_H0*EFTGamma2V +1.5_dl*adotoa*EFTGamma3V&
!            & +0.5_dl*adotoa*EFTGamma4V -4._dl*adotoa*EFTGamma6V)&
!            & +(EFTGamma4V +2._dl*EFTGamma5V)/(2._dl*(1._dl+EFTOmegaV) -2._dl*EFTGamma4V)*(EFTgrhoQ +EFTgpresQ +a*adotoa**2*EFTOmegaP&
!            & -EFTGamma4V*(Hdot-adotoa**2) +a*adotoa*EFT_H0*EFTGamma2V +3._dl*EFTGamma3V*(adotoa**2-Hdot))
!        ! Second part of the D coefficient, the k^4 term:
!        EFTpiD2 = +(+0.5_dl*EFTGamma3V +0.5_dl*EFTGamma4V &
!            & +(EFTGamma4V +2._dl*EFTGamma5V)/(2._dl*(1._dl+EFTOmegaV) -2._dl*EFTGamma4V)*(EFTGamma3V +EFTGamma4V))
!
!        if ( CP%EFTflag==4 ) then ! Full mapping EFT overwrites some of these quantities:
!            if ( CP%FullMappingEFTmodel==1 ) then ! Horava gravity
!                EFTpiA1 = 0._dl
!                EFTpiA2 = 0.5_dl*CP%Horava_eta
!                EFTpiB1 = 0._dl
!                EFTpiB2 = adotoa*CP%Horava_eta
!                EFTpiC  = 0._dl
!                EFTpiD1 = 0.5_dl*(3._dl*CP%Horava_lambda-2._dl*CP%Horava_xi)*(adotoa**2-Hdot) + 0.5_dl*CP%Horava_eta*(adotoa**2+Hdot)
!                EFTpiD2 = 0.5_dl*CP%Horava_lambda*(1._dl+CP%Horava_xi)
!            end if
!        end if
!
!        EFT_kinetic  = 9._dl*( 1._dl +EFTOmegaV -EFTGamma4V )*( 4._dl*EFTc*( 1._dl +EFTOmegaV -EFTGamma4V ) &
!            & +3._dl*adotoa**2*EFTOmegaP**2*a2 + &
!            & a2*EFT_H0*( EFT_H0*( 3._dl*EFTGamma2V**2 +8._dl*EFTGamma1V*( 1._dl +EFTOmegaV -EFTGamma4V ) +6._dl*adotoa*EFTGamma2V*EFTOmegaP ) ) )
!        ! Gradient term:
!        EFT_gradient = 9._dl*(8._dl*a*adotoa**2*EFTGamma5P - 16._dl*adotoa**2*EFTGamma5V**2 + 16._dl*EFTc*EFTGamma5V**2 - 2._dl*a**2*adotoa**2*EFTGamma4P*EFTOmegaP + 4._dl*a**2*adotoa**2*EFTGamma5P*EFTOmegaP -&
!            &4._dl*a*adotoa**2*EFTGamma5V*EFTOmegaP - 4._dl*a**2*adotoa**2*EFTGamma4P*EFTGamma5V*EFTOmegaP - 8._dl*a*adotoa**2*EFTGamma5V**2*EFTOmegaP + 3._dl*a**2*adotoa**2*EFTOmegaP**2 +&
!            &4._dl*a**2*adotoa**2*EFTGamma5V*EFTOmegaP**2 + 16._dl*a*adotoa**2*EFTGamma5P*EFTOmegaV + 16._dl*EFTc*EFTGamma5V*EFTOmegaV - 16._dl*adotoa**2*EFTGamma5V**2*EFTOmegaV -&
!            &2._dl*a**2*adotoa**2*EFTGamma4P*EFTOmegaP*EFTOmegaV + 4._dl*a**2*adotoa**2*EFTGamma5P*EFTOmegaP*EFTOmegaV - 4._dl*a*adotoa**2*EFTGamma5V*EFTOmegaP*EFTOmegaV +&
!            &3._dl*a**2*adotoa**2*EFTOmegaP**2*EFTOmegaV + 8._dl*a*adotoa**2*EFTGamma5P*EFTOmegaV**2 - a**2*EFTGamma2V**2*EFT_H0**2*(1 + EFTOmegaV) +&
!            &4._dl*a**2*adotoa**2*EFTGamma5V*EFTOmegaPP*(1._dl + 2._dl*EFTGamma5V + EFTOmegaV) + 4._dl*EFTc*(4._dl*EFTGamma5V + (1._dl + EFTOmegaV)**2) -&
!            &2._dl*a*adotoa*EFT_H0*(a*EFTGamma2P*(1._dl - EFTGamma4V + EFTOmegaV)*(1._dl + 2._dl*EFTGamma5V + EFTOmegaV) +&
!            &EFTGamma2V*(EFTGamma4V*(-1._dl + 2._dl*a*EFTGamma5P + 2._dl*EFTGamma5V + a*EFTOmegaP - EFTOmegaV) +&
!            &(1._dl + EFTOmegaV)*(1._dl + a*(EFTGamma4P - 2._dl*EFTGamma5P - EFTOmegaP) + EFTOmegaV) - 2._dl*EFTGamma5V*(1._dl - a*EFTGamma4P + a*EFTOmegaP + EFTOmegaV))) +&
!            &8._dl*EFTGamma5V*Hdot + 16._dl*EFTGamma5V**2*Hdot + 4._dl*a*EFTGamma5V*EFTOmegaP*Hdot + 8._dl*a*EFTGamma5V**2*EFTOmegaP*Hdot + 16._dl*EFTGamma5V*EFTOmegaV*Hdot +&
!            &16._dl*EFTGamma5V**2*EFTOmegaV*Hdot + 4._dl*a*EFTGamma5V*EFTOmegaP*EFTOmegaV*Hdot + 8._dl*EFTGamma5V*EFTOmegaV**2*Hdot +&
!            &4._dl*EFTGamma4V**2*(adotoa**2*(1._dl + 2._dl*a*EFTGamma5P + 4._dl*EFTGamma5V + a*EFTOmegaP + EFTOmegaV) - (1._dl + 2._dl*EFTGamma5V + EFTOmegaV)*Hdot) +&
!            &2._dl*EFTGamma4V*(adotoa**2*(-(a**2*EFTOmegaP**2) + a**2*EFTOmegaPP*(1._dl + 2._dl*EFTGamma5V + EFTOmegaV) -&
!            &4._dl*(1._dl + EFTOmegaV)*(1._dl + 2._dl*a*EFTGamma5P + 4._dl*EFTGamma5V + EFTOmegaV) - a*EFTOmegaP*(3._dl + 2._dl*a*EFTGamma5P + 2._dl*EFTGamma5V + 3._dl*EFTOmegaV)) +&
!            &(1._dl + 2._dl*EFTGamma5V + EFTOmegaV)*(4._dl + a*EFTOmegaP + 4._dl*EFTOmegaV)*Hdot))
!
!        ! 7) Compute k_max_CAMB. This is slightly overestimated to be safe...
!
!        !    This is just to be safe if doing only background stuff. We still require stability of linear scales.
!        kmax = 0.1_dl
!        !    This is if we want scalar Cls.
!        if (CP%WantCls) kmax = CP%Max_eta_k/CP%tau0
!        !    This is if we want also (or only) transfer functions.
!        if (CP%WantTransfer) then
!            kmax = CP%Max_eta_k/CP%tau0*exp((int(log(CP%Transfer%kmax/(CP%Max_eta_k/CP%tau0))*(3*AccuracyBoost))+2)/(3*AccuracyBoost))
!        end if
!
!        ! All the coefficients should not be Nan. This can happen for strange values of the parameters
!        ! for which a division by zero may occur.
!
!        EFT_HaveNan = IsNaN(EFTAT).or.IsNaN(EFTpiA1).or.IsNaN(EFTpiA2).or.IsNaN(EFTpiB1).or.IsNaN(EFTpiB2)&
!            & .or.IsNaN(EFTpiC).or.IsNaN(EFTpiD1).or.IsNaN(EFTpiD2).or.IsNaN(EFT_kinetic).or.IsNaN(EFT_gradient)
!
!        if (EFT_HaveNan) then
!            EFTStabilityComputation = .false.
!            if (Feedbacklevel > 0) write(*,*) 'EFTCAMB: model has Nan.'
!            return
!        end if
!
!        ! Mathematical stability:
!        if (CP%EFT_mathematical_stability) then
!
!            ! 1- the A coefficient should not change sign in time and in k, i.e. it shall not be zero.
!            !    This is the stronger stability constraint since violating it would violate the mathematical
!            !    consistency of the pi field equation.
!
!            !    The first condition is A1/=0. Implemented by detecting sign changes in A1.
!            if ( EFTpiA1*PastA1 < 0._dl ) then
!                EFTStabilityComputation = .false.
!                if (Feedbacklevel > 0) write(*,*) 'EFTCAMB: mathematical instability, A is zero in time.'
!            end if
!            PastA1 = EFTpiA1
!            !    The second one is the condition on k.
!            if ( (EFTpiA1 > 0 .and. EFTpiA1 + kmax**2*EFTpiA2 < 0) .or. &
!                &(EFTpiA1 < 0 .and. EFTpiA1 + kmax**2*EFTpiA2 > 0) ) then
!                EFTStabilityComputation = .false.
!                if (Feedbacklevel > 0) write(*,*) 'EFTCAMB: mathematical instability, A is zero in k.'
!            end if
!
!            ! 2- the AT coefficient should not change sign in time, i.e. it shall not be zero.
!            !    This is the second stronger stability constraint since violating it would
!            !    violate the mathematical consistency of the tensor perturbation equation.
!            !    Implemented by detecting sign changes in AT.
!            if ( EFTAT*PastAT < 0._dl ) then
!                EFTStabilityComputation = .false.
!                if (Feedbacklevel > 0) write(*,*) 'EFTCAMB: mathematical instability, AT is zero in time.'
!            end if
!            PastAT = EFTAT
!
!            ! 3- we do not want (fast) growing exponential modes.
!            !    This condition prevents the pi field from growing exponentially and destroying everything.
!            !    Even though this condition is not completely related to physics not mathematics, violating it
!            !    would completely screw up cosmological observables.
!
!            !    This is the maximum allowed rate of instability. Units shall be Mpc^-1.
!            EFT_instability_rate = 0._dl
!
!            !    This condition needs to be tested in k. Sample in k.
!
!            ind_max = 10
!
!            do ind = 1, ind_max
!                ! kmode to test. Linear sampling. Should suffice... (??)
!                tempk = 0._dl + REAL(ind-1)*(kmax)/REAL(ind_max-1)
!                ! vaule that discriminates between different cases:
!                temp1 = (EFTpiB1 +EFTpiB2*tempk**2)
!                temp2 = (EFTpiA1 +EFTpiA2*tempk**2)
!                temp3 = temp1**2 -4._dl*temp2*(EFTpiC +EFTpiD1*tempk**2 + EFTpiD2*tempk**4)
!
!                !    case 1:
!                if (temp3 > 0._dl .and. temp2 /= 0._dl) then
!                    temp4 = +0.5_dl*(-temp1 +sqrt(temp3))/temp2
!                    temp5 = +0.5_dl*(-temp1 -sqrt(temp3))/temp2
!                    if (temp4>EFT_instability_rate .or. temp5>EFT_instability_rate) then
!                        EFTStabilityComputation = .false.
!                        if (Feedbacklevel > 0) write(*,*) 'EFTCAMB: mathematical instability. Growing exponential at k', tempk, temp4, temp5
!                        exit
!                    end if
!                !    case 2:
!                else if ( temp2 /= 0._dl ) then
!                    temp4 = -0.5_dl*temp1/temp2
!                    if (temp4>EFT_instability_rate) then
!                        EFTStabilityComputation = .false.
!                        if (Feedbacklevel > 0) write(*,*) 'EFTCAMB: mathematical instability. Growing exponential at k', tempk, temp4
!                        exit
!                    end if
!                end if
!
!            end do
!
!        end if
!
!        ! Additional priors:
!        if (CP%EFTAdditionalPriors) then
!            if (EFTw(a,0)>-1._dl/3._dl) EFTStabilityComputation = .false.
!        end if
!
!        ! Minkowsky prior: some theories have known stability properties on Minkowsky background:
!        if (CP%MinkowskyPriors) then
!            if (CP%EFTflag==4) then
!                if (CP%FullMappingEFTmodel==1) then ! Horava gravity
!
!                    if ( CP%Horava_lambda > -2._dl/3._dl .and. CP%Horava_lambda < 0._dl ) then
!                        EFTStabilityComputation = .false.
!                        if (Feedbacklevel > 0) write(*,*) 'EFTCAMB: Instability on Minkowsky backgound'
!                    end if
!
!                    if ( CP%Horava_eta < 0._dl .or. CP%Horava_eta > 2._dl*CP%Horava_xi +2._dl ) then
!                        EFTStabilityComputation = .false.
!                        if (Feedbacklevel > 0) write(*,*) 'EFTCAMB: Instability on Minkowsky backgound'
!                    end if
!
!                end if
!            end if
!        end if
!
!
!
!        ! Physical viability:
!        if (CP%EFT_physical_stability) then
!
!            if ( .not. EFT_old_stability .and. &
!                & CP%EFTFlag /= 4 .and. ( &
!                & (EFTGamma6V /= 0._dl) .or.      &
!                & ((EFTGamma3V + EFTGamma4V) /= 0._dl) ) ) then
!
!                write(*,*) 'EFTCAMB WARNING: stability for model beyond GLPV has not been worked out.'
!                write(*,*) 'It will be added in a future release.'
!                write(*,*) 'If you want to run this model disable EFT_physical_stability.'
!
!                EFTStabilityComputation = .false.
!                return
!            end if
!
!            ! 1- Positive gravitational constant:
!            if (1._dl +EFTOmegaV <= 0) then
!                EFTStabilityComputation = .false.
!                if (Feedbacklevel > 0) write(*,*) 'EFTCAMB: negative gravitational constant', 1._dl +EFTOmegaV
!            end if
!
!            ! 2- Old ghost and gradient conditions:
!            if ( EFT_old_stability .or. &
!                & (EFTGamma6V /= 0._dl) .or. &
!                & ((EFTGamma3V + EFTGamma4V) /= 0._dl) ) then
!
!                ! Ghost instability:
!                if (EFTpiA1 < 0 .or. ( EFTpiA1 + kmax**2*EFTpiA2 < 0)) then
!                    EFTStabilityComputation = .false.
!                    if (Feedbacklevel > 0) write(*,*) 'EFTCAMB: ghost instability', EFTpiA1, EFTpiA1 + kmax**2*EFTpiA2
!                end if
!                ! Gradient instability 1:
!                if (EFTpiD1 < 0) then
!                    EFTStabilityComputation = .false.
!                    if (Feedbacklevel > 0) write(*,*) 'EFTCAMB: gradient instability k^2', EFTpiD1
!                end if
!                ! Gradient instability 2:
!                if (EFTpiD2 < 0) then
!                    EFTStabilityComputation = .false.
!                    if (Feedbacklevel > 0) write(*,*) 'EFTCAMB: gradient instability k^4', EFTpiD2
!                end if
!
!            else
!                ! New ghost and gradient conditions:
!                ! ghost condition:
!                if ( EFT_kinetic < 0._dl ) then
!                    EFTStabilityComputation = .false.
!                    if (Feedbacklevel > 0) write(*,*) 'EFTCAMB new stability: ghost instability. Kinetic: ', EFT_kinetic
!                end if
!                ! gradient instability:
!                if ( EFT_gradient < 0._dl ) then
!                    EFTStabilityComputation = .false.
!                    if (Feedbacklevel > 0) write(*,*) 'EFTCAMB new stability: gradient instability. Gradient: ', EFT_gradient
!                end if
!
!            end if
!
!            !5- Positive effective mass of pi:
!            if (EFTpiC < 0.and.EFTpiMassPrior) then
!                EFTStabilityComputation = .false.
!                if (Feedbacklevel > 0) write(*,*) 'EFTCAMB: negative mass'
!            end if
!            ! 6- No tensor ghosts:
!            if (EFTAT < 0) then
!                EFTStabilityComputation = .false.
!                if (Feedbacklevel > 0) write(*,*) 'EFTCAMB: tensor ghost instability'
!            end if
!            ! 7- Sub-luminal propagation:
!            if (EFTlightspeedPrior) then
!                if (EFTpiA2==0.and.EFTpiD2==0.and.(EFTpiD1/EFTpiA1)>1.001_dl) then
!                    EFTStabilityComputation = .false.
!                    if (Feedbacklevel > 0) write(*,*) 'EFTCAMB: tachion perturbations'
!                else if (EFTpiA2/=0.and.EFTpiD2/=0.and.(EFTpiD2/EFTpiA2)>1.001_dl) then
!                    EFTStabilityComputation = .false.
!                    if (Feedbacklevel > 0) write(*,*) 'EFTCAMB: tachion perturbations'
!                end if
!            end if
!            ! 8- Every theory has it's own peculiarities...
!            ! 1) F(R): for this model it is easy to show that the positive mass condition requires that OmegaPrime
!            !    should be positive. We add this test for this models as it is numerically easier to check.
!            if (CP%EFTflag==2.and.CP%DesignerEFTmodel==1) then
!                if (EFTOmega(0.11_dl*EFTturnonpiInitial,1)*EFTOmegaP<0) EFTStabilityComputation = .false.
!            end if
!        end if
!
!        return

    end function EFTStabilityComputation

    !----------------------------------------------------------------------------------------

end module EFTCAMB_stability

!----------------------------------------------------------------------------------------
