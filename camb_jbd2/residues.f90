   !!!ENSURE FLATNESS TODAY!!!
    current_omegav = omegade(des_nstep)
    try = 0
    do while ( (abs(current_omegav - CP%omegav) .gt. precflat) .or. (abs(phibdf_final - CP%EFTphi_0) .gt. precflat) )


    !ADJUSTING THE WEIGHT OF THE POTENTIAL/LAMBDA
    if ((current_omegav - CP%omegav) .lt. precflat) then

    arglamda = arglamda + precflat*10.0

    call EvolveBackground(phibdi_final,phibdf_final,arglamda)

    else if ((current_omegav-CP%omegav) .gt. precflat) then

    arglamda = arglamda - precflat*10.0

    call EvolveBackground(phibdi_final,phibdf_final,arglamda)

    end if


    !TRYING NOT TO GET TOO FAR AWAY FROM PHI_0 TODAY
    if ((phibdf_final - CP%EFTphi_0) .lt. precflat) then

    phibdi_final = phibdi_final + precflat

    call EvolveBackground(phibdi_final,phibdf_final,arglamda)

    else if ((phibdf_final - CP%EFTphi_0) .gt. precflat) then

    phibdi_final = phibdi_final - precflat

    call EvolveBackground(phibdi_final,phibdf_final,arglamda)

    end if

    current_omegav = omegade(des_nstep)
    print*, current_omegav, CP%omegav
    print*, phibdf_final, CP%EFTphi_0
    end do







function Hdotbd(a,phi,phiprime,om,omr,H0,arglamda) !d H/dt
    implicit none
    real(dl) :: a,phi,phiprime,om,omr,H0,Hdotbd,hubble_aux,arglamda
    real(dl) :: rhonu, presnu, rhonu_tot, presnu_tot, grhormass_t, EFT_E_nu, EFT_EP_nu, EFT_pres_nu
    real(dl) :: presnudot_tot, adotoa, Hdot, presnudotdot_tot, presnudot, presnudotdot
    integer  :: nu_i
    real(dl) :: OmegaMassiveNu_EFT

    hubble_aux = Hubblebd(a,phi,phiprime,om,omr,H0,arglamda)

    Hdotbd = 0.5*(-3.0*hubble_aux**2.0 - H0**2.0*Eomr(a,omr)/phi - &
    0.5*CP%EFTwbds*hubble_aux**2.0*phiprime**2.0/phi**2.0 - 2.0*hubble_aux**2.0*phiprime/phi - &
    phidotdot(phi,phiprime,a,om,omr,H0,arglamda)/phi + H0**2.0*potbd(phi,0,arglamda)/phi)

end function Hdotbd

function Hdotdotbd(a,phi,phiprime,om,omr,H0,arglamda) !d^2 H/dt^2
    implicit none
    real(dl) :: a,phi,phiprime,om,omr,H0,Hdotdotbd,hubble_aux,hubbledot_aux,hubbledtau_aux,arglamda
        ! Massive neutrinos variables
        real(dl) :: rhonu, presnu, rhonu_tot, presnu_tot, grhormass_t, EFT_E_nu, EFT_EP_nu, EFT_EPP_nu, EFT_E3P_nu, EFT_pres_nu
        real(dl) :: presnudot_tot, adotoa, Hdot, presnudotdot_tot, presnudot, presnudotdot
        integer  :: nu_i
        real(dl) :: OmegaMassiveNu_EFT
	real(dl) :: Efun, H0_EFT

    H0_EFT = sqrt(grhom/3.0)
    hubble_aux = Hubblebd(a,phi,phiprime,om,omr,H0,arglamda)
    hubbledot_aux = Hdotbd(a,phi,phiprime,om,omr,H0,arglamda)
    hubbledtau_aux = Hdtaubd(a,phi,phiprime,om,omr,H0,arglamda)
    EFun = hubble_aux**2.0/H0**2.0

    Hdotdotbd = 0.5*(-6.0*(hubbledot_aux)*hubble_aux + 4.0*hubble_aux*H0**2.0*Eomr(a,omr)/phi + &
    hubble_aux*phiprime*H0**2.0*Eomr(a,omr)/phi**2.0 - CP%EFTwbds*phidotdot(phi,phiprime,a,om,omr,H0,arglamda)*hubble_aux*phiprime/phi**2.0 &
    + CP%EFTwbds*hubble_aux**3.0*phiprime**3.0/phi**3.0 - (2.0/phi)*(hubbledot_aux*hubble_aux*phiprime + &
    hubble_aux*phidotdot(phi,phiprime,a,om,omr,H0,arglamda)) + 2.0*hubble_aux**3.0*phiprime**2.0/phi**2.0 + &
    phidotdot(phi,phiprime,a,om,omr,H0,arglamda)*hubble_aux*phiprime/phi**2.0 + &
    H0**2.0*hubble_aux*potbd(phi,1,arglamda)*phiprime/phi - H0**2.0*hubble_aux*potbd(phi,0,arglamda)*phiprime/phi**2.0)

end function Hdotdotbd

function Hdtaubd(a,phi,phiprime,om,omr,H0,arglamda) !d H/ d tau
    implicit none
    real(dl) :: a,phi,phiprime,om,omr,H0,hubble_aux,hubbledot_aux,Hdtaubd,arglamda

    hubble_aux = Hubblebd(a,phi,phiprime,om,omr,H0,arglamda)
    hubbledot_aux = Hdotbd(a,phi,phiprime,om,omr,H0,arglamda)

    Hdtaubd = a**2.0*(hubbledot_aux + hubble_aux**2.0)

end function Hdtaubd    

function Hdtaudtaubd(a,phi,phiprime,om,omr,H0,arglamda) !d^2 H/ dtau^2
    implicit none
    real(dl) :: a,phi,phiprime,om,omr,H0,hubble_aux,hubbledot_aux,hubbledotdot_aux,Hdtaudtaubd,arglamda

    hubble_aux = Hubblebd(a,phi,phiprime,om,omr,H0,arglamda)
    hubbledot_aux = Hdotbd(a,phi,phiprime,om,omr,H0,arglamda)
    hubbledotdot_aux = Hdotdotbd(a,phi,phiprime,om,omr,H0,arglamda)

    Hdtaudtaubd = a**3.0*(hubbledotdot_aux + 2.0*hubble_aux**3.0 + 4.0*hubble_aux*hubbledot_aux)

end function Hdtaudtaubd

function dedensitybd(a,phi,phiprime,om,omr,H0,arglamda)
    implicit none
    real(dl) :: a,phi,phiprime,om,omr,H0,hubble_aux,dedensitybd,arglamda

    hubble_aux = Hubblebd(a,phi,phiprime,om,omr,H0,arglamda)

    dedensitybd =  0.5*CP%EFTwbds*phiprime**2.0*hubble_aux**2.0/phi - 3.0*hubble_aux**2.0*phiprime + H0**2.0*potbd(phi,0,arglamda)



function phidotdot(phi,phiprime,a,om,omr,H0,arglamda) !in physical time
    implicit none
    real(dl) :: phi, phiprime, a, om, omr, H0, hubble_aux, phidotdot, arglamda

    hubble_aux = Hubblebd(a,phi,phiprime,om,omr,H0,arglamda)

    phidotdot = -3.0*hubble_aux**2.0*phiprime + 3.0*H0**2.0*Eom(a,om)/(3.0+2.0*CP%EFTwbds) + H0**2.0*(4.0*potbd(phi,0,arglamda) - 2.0*phi*potbd(phi,1,arglamda))/(3.0+2.0*CP%EFTwbds) 

end function phidotdot


