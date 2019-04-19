! -------------------------------------------------------------------------------------------------
!
!   EFTCAMB
!
!   Developed and implemented by:
!       Bin Hu (hu@lorentz.leidenuniv.nl), Marco Raveri (mraveri@sissa.it)
!       Noemi Frusciante (fruscian@iap.fr), Alessandra Silvestri (silvestri@lorentz.leidenuniv.nl)
!
!
!   EFTCAMB module for Horndeski models.
!
!   For more informations about the methods contained in this file see
!   the documentation: arXiv:1405.3590
!
! -------------------------------------------------------------------------------------------------

! Definitions of the RPH functions later used by the code.
module RPHfunctions
    use EFTDef

contains

    ! ----------------------------

    function RPH_PlanckMass(a,deriv)
        !RPH function Planck mass.
        implicit none
        real(dl)              :: RPH_PlanckMass
        real(dl) , intent(in) :: a      ! scale factor
        integer  , intent(in) :: deriv  ! order of derivative wrt the scale factor

        select case (CP%RPHmassPmodel)
            case (0) ! operator neglected
                RPH_PlanckMass = 0._dl
            case (1) ! constant model
                if (deriv==0) then
                    RPH_PlanckMass = CP%RPHmassP0
                else
                    RPH_PlanckMass = 0._dl
                end if
            case (2) ! power law model
                if (deriv==0) then
                    RPH_PlanckMass= CP%RPHmassP0*a**CP%RPHmassPexp
                else if (deriv==1) then
                    RPH_PlanckMass= CP%RPHmassPexp*CP%RPHmassP0*a**(CP%RPHmassPexp-1.0_dl)
                else if (deriv==2) then
                    RPH_PlanckMass= (CP%RPHmassPexp-1.0_dl)*CP%RPHmassPexp*CP%RPHmassP0*&
                               &a**(CP%RPHmassPexp-2.0_dl)
                else if (deriv==3) then
                    RPH_PlanckMass= (CP%RPHmassPexp-2.0_dl)*(CP%RPHmassPexp-1.0_dl)*&
                                &CP%RPHmassPexp*CP%RPHmassP0*a**(CP%RPHmassPexp-3.0_dl)
                end if
            case (3) ! user defined
                stop 'You have to define your model for RPH_PlanckMass first! Go to the file EFT_Horndeski.f90.'
                !if (deriv==0) then
                !else if (deriv==1) then
                !else if (deriv==2) then
                !else if (deriv==3) then
                !end if
            case default ! not implemented
                stop 'Wrong model for RPH PlanckMass.'
        end select

    end function RPH_PlanckMass

    ! ----------------------------

    function RPH_Kineticity(a,deriv)
        !RPH function \alpha_K.
        implicit none
        real(dl)              :: RPH_Kineticity
        real(dl) , intent(in) :: a      ! scale factor
        integer  , intent(in) :: deriv  ! order of derivative wrt the scale factor

        select case (CP%RPHkineticitymodel)
            case (0) ! operator neglected
                RPH_Kineticity = 0._dl
            case (1) ! constant model
                if (deriv==0) then
                    RPH_Kineticity = CP%RPHkineticity0
                else
                    RPH_Kineticity = 0._dl
                end if
            case (2) ! power law model
                if (deriv==0) then
                    RPH_Kineticity= CP%RPHkineticity0*a**CP%RPHkineticityexp
                else if (deriv==1) then
                    RPH_Kineticity= CP%RPHkineticityexp*CP%RPHkineticity0*a**(CP%RPHkineticityexp-1.0_dl)
                else if (deriv==2) then
                    RPH_Kineticity= (CP%RPHkineticityexp-1.0_dl)*CP%RPHkineticityexp*CP%RPHkineticity0*&
                                  &a**(CP%RPHkineticityexp-2.0_dl)
                else if (deriv==3) then
                    RPH_Kineticity= (CP%RPHkineticityexp-2.0_dl)*(CP%RPHkineticityexp-1.0_dl)*&
                                  &CP%RPHkineticityexp*CP%RPHkineticity0*a**(CP%RPHkineticityexp-3.0_dl)
                end if
            case (3) ! user defined
                stop 'You have to define your model for RPH_Kineticity first! Go to the file EFT_Horndeski.f90.'
                !if (deriv==0) then
                !else if (deriv==1) then
                !else if (deriv==2) then
                !else if (deriv==3) then
                !end if
            case default ! not implemented
                stop 'Wrong model for RPH Kineticity.'
        end select

    end function RPH_Kineticity

    ! ----------------------------

    function RPH_Braiding(a,deriv)
        !RPH function \alpha_B.
        implicit none
        real(dl)              :: RPH_Braiding
        real(dl) , intent(in) :: a      ! scale factor
        integer  , intent(in) :: deriv  ! order of derivative wrt the scale factor

        select case (CP%RPHbraidingmodel)
            case (0) ! operator neglected
                RPH_Braiding = 0._dl
            case (1) ! constant model
                if (deriv==0) then
                    RPH_Braiding = CP%RPHbraiding0
                else
                    RPH_Braiding = 0._dl
                end if
            case (2) ! power law model
                if (deriv==0) then
                    RPH_Braiding= CP%RPHbraiding0*a**CP%RPHbraidingexp
                else if (deriv==1) then
                    RPH_Braiding= CP%RPHbraidingexp*CP%RPHbraiding0*a**(CP%RPHbraidingexp-1.0_dl)
                else if (deriv==2) then
                    RPH_Braiding= (CP%RPHbraidingexp-1.0_dl)*CP%RPHbraidingexp*CP%RPHbraiding0*&
                                  &a**(CP%RPHbraidingexp-2.0_dl)
                else if (deriv==3) then
                    RPH_Braiding= (CP%RPHbraidingexp-2.0_dl)*(CP%RPHbraidingexp-1.0_dl)*&
                                  &CP%RPHbraidingexp*CP%RPHbraiding0*a**(CP%RPHbraidingexp-3.0_dl)
                end if
            case (3) ! user defined
                stop 'You have to define your model for RPH_Braiding first! Go to the file EFT_Horndeski.f90.'
                !if (deriv==0) then
                !else if (deriv==1) then
                !else if (deriv==2) then
                !else if (deriv==3) then
                !end if
            case default ! not implemented
                stop 'Wrong model for RPH Braiding.'
        end select

    end function RPH_Braiding

    ! ----------------------------

    function RPH_Tensor(a,deriv)
        !RPH function \alpha_T.
        implicit none
        real(dl)              :: RPH_Tensor
        real(dl) , intent(in) :: a      ! scale factor
        integer  , intent(in) :: deriv  ! order of derivative wrt the scale factor

        select case (CP%RPHtensormodel)
            case (0) ! operator neglected
                RPH_Tensor = 0._dl
            case (1) ! constant model
                if (deriv==0) then
                    RPH_Tensor = CP%RPHtensor0
                else
                    RPH_Tensor = 0._dl
                end if
            case (2) ! power law model
                if (deriv==0) then
                    RPH_Tensor= CP%RPHtensor0*a**CP%RPHtensorexp
                else if (deriv==1) then
                    RPH_Tensor= CP%RPHtensorexp*CP%RPHtensor0*a**(CP%RPHtensorexp-1.0_dl)
                else if (deriv==2) then
                    RPH_Tensor= (CP%RPHtensorexp-1.0_dl)*CP%RPHtensorexp*CP%RPHtensor0*&
                                  &a**(CP%RPHtensorexp-2.0_dl)
                else if (deriv==3) then
                    RPH_Tensor= (CP%RPHtensorexp-2.0_dl)*(CP%RPHtensorexp-1.0_dl)*&
                                  &CP%RPHtensorexp*CP%RPHtensor0*a**(CP%RPHtensorexp-3.0_dl)
                end if
            case (3) ! user defined
                stop 'You have to define your model for RPH_Tensor first! Go to the file EFT_Horndeski.f90.'
                !if (deriv==0) then
                !else if (deriv==1) then
                !else if (deriv==2) then
                !else if (deriv==3) then
                !end if
            case default ! not implemented
                stop 'Wrong model for RPH Tensor.'
        end select

    end function RPH_Tensor

end module RPHfunctions

! -------------------------------------------------------------------------------------------------

! Defines all the functions used to map the RPH approach to EFT.
module RPH_to_EFT
    use Precision
    use ModelParams
    use EFTDef
    use RPHfunctions
    use EFTdeEOS

    implicit none

contains

    function RPH_Hubble(a)
        !Function that returns conformal H as a function of a
        use MassiveNu

        implicit none

        real(dl) RPH_Hubble
        real(dl), intent(in) :: a
        real(dl) a2, grho, rhonu
        integer nu_i

        a2=a**2

        grho=grhob/a +grhoc/a +grhornomass/a2 +grhog/a2 +grhov*EFTw(a,3)

        if (CP%Num_Nu_massive /= 0) then
            !Get massive neutrino density relative to massless
            do nu_i = 1, CP%nu_mass_eigenstates
                call Nu_rho(a*nu_masses(nu_i),rhonu)
                grho=grho+rhonu*grhormass(nu_i)/a2
            end do
        end if

        RPH_Hubble=sqrt(grho/3)

    end function RPH_Hubble

    ! ----------------------------

    function RPH_HubbleDot(a)
        !Function that returns time derivative of conformal H as a function of a
        use MassiveNu

        implicit none

        real(dl) RPH_HubbleDot
        real(dl), intent(in) :: a
        real(dl) a2, grho, gpres, rhonu, presnu
        integer nu_i

        a2=a*a

        gpres = 0._dl ; grho = 0._dl ;

        grho  = +grhob/a +grhoc/a +grhornomass/a2 +grhog/a2 +grhov*EFTw(a,3)
        gpres = +(grhog/a2 +grhornomass/a2)/3._dl +EFTw(a,0)*grhov*EFTw(a,3)

        if (CP%Num_Nu_massive /= 0) then
            !Get massive neutrino density relative to massless
            do nu_i = 1, CP%nu_mass_eigenstates
                rhonu = 0._dl ; presnu = 0._dl
                call Nu_background(a*nu_masses(nu_i), rhonu, presnu)
                grho  = grho  +rhonu*grhormass(nu_i)/a2
                gpres = gpres +presnu*grhormass(nu_i)/a2
            end do
        end if

        RPH_HubbleDot = -0.5_dl*(gpres +grho/3)

    end function RPH_HubbleDot

    ! ----------------------------

    function RPH_EFT_Omega(a, deriv)
        !Function that returns the function EFT c for designer approach.
        implicit none

        integer , intent(in)  :: deriv
        real(dl), intent(in)  :: a
        real(dl) :: RPH_EFT_Omega, RPH_PM_V

        select case (deriv)
            case (0)
                RPH_PM_V      = RPH_PlanckMass(a,0)
                RPH_EFT_Omega = RPH_PM_V +RPH_Tensor(a,0)*(1._dl +RPH_PM_V)
            case (1)
                RPH_EFT_Omega = RPH_PlanckMass(a,0)*RPH_Tensor(a,1) &
                              &+(1._dl +RPH_Tensor(a,0))*RPH_PlanckMass(a,1)
            case (2)
                RPH_EFT_Omega = RPH_PlanckMass(a,0)*RPH_Tensor(a,2) &
                              & +2._dl*RPH_PlanckMass(a,1)*RPH_Tensor(a,1) &
                              & +(1._dl +RPH_Tensor(a,0))*RPH_PlanckMass(a,2)
            case (3)
                RPH_EFT_Omega = RPH_PlanckMass(a,0)*RPH_Tensor(a,3) &
                              & +3._dl*RPH_PlanckMass(a,2)*RPH_Tensor(a,1) &
                              & +3._dl*RPH_PlanckMass(a,1)*RPH_Tensor(a,2) &
                              & +(1._dl +RPH_Tensor(a,0))*RPH_PlanckMass(a,3)
        end select

    end function RPH_EFT_Omega

    ! ----------------------------

    function RPH_EFTc(a)
        !Function that returns the function EFT c for designer approach.
        use MassiveNu

        implicit none

        real(dl) a, RPH_EFTc
        real(dl) a2, grhob_t,grhoc_t,grhor_t,grhog_t,grhov_t,grho, gpres, EFT_H0
        real(dl) adotoa, adotdota, Hdot
        real(dl) EFTOmegaV, EFTOmegaP, EFTOmegaPP
        real(dl) EFT_grhonu, EFT_gpinu, EFT_grhonu_tot, EFT_gpinu_tot, grhormass_t
        integer  nu_i

        a2 = a*a

        grhob_t=grhob/a
        grhoc_t=grhoc/a
        grhor_t=grhornomass/a2
        grhog_t=grhog/a2
        grhov_t=grhov*EFTw(a,3)

        grho=grhob_t+grhoc_t+grhor_t+grhog_t+grhov_t
        gpres= (grhog_t+grhor_t)/3._dl +EFTw(a,0)*grhov_t

        EFT_gpinu_tot = 0._dl
        EFT_grhonu_tot = 0._dl
        if (CP%Num_Nu_Massive /= 0) then
            do nu_i = 1, CP%Nu_mass_eigenstates
                EFT_grhonu    = 0._dl
                EFT_gpinu     = 0._dl
                grhormass_t=grhormass(nu_i)/a2
                call Nu_background(a*nu_masses(nu_i),EFT_grhonu,EFT_gpinu)
                EFT_grhonu_tot = EFT_grhonu_tot + grhormass_t*EFT_grhonu
                EFT_gpinu_tot  = EFT_gpinu_tot + grhormass_t*EFT_gpinu
            end do
        end if

        grho  = grho  + EFT_grhonu_tot
        gpres = gpres + EFT_gpinu_tot
        ! 3) Hubble and derivatives:
        adotoa=sqrt(grho/3)
        adotdota=(adotoa*adotoa-gpres)/2._dl
        Hdot =adotdota-adotoa**2._dl
        ! 4) EFT functions:
        EFTOmegaV  = RPH_EFT_Omega(a,0)
        EFTOmegaP  = RPH_EFT_Omega(a,1)
        EFTOmegaPP = RPH_EFT_Omega(a,2)

        RPH_EFTc = (adotoa*adotoa - Hdot)*(EFTOmegaV + a*EFTOmegaP*0.5_dl) &
             & - 0.5_dl*a2*adotoa*adotoa*EFTOmegaPP&
             & + 0.5_dl*grhov_t*(1._dl+EFTw(a,0))

    end function RPH_EFTc

    ! ----------------------------

    function RPH_EFTcdot(a)
        !Function that returns the time derivative of the EFT function c for designer approach.
        use MassiveNu

        implicit none

        real(dl) a, RPH_EFTcdot
        real(dl) a2, grhob_t,grhoc_t,grhor_t,grhog_t,grhov_t,grho, gpres, EFT_H0
        real(dl) adotoa, adotdota, Hdot,Hdotdot
        real(dl) EFTOmegaV, EFTOmegaP, EFTOmegaPP,EFTOmegaPPP
        real(dl) EFTc, EFTcdot, EFTLambda, EFTLambdadot
        real(dl) EFT_grhonu, EFT_gpinu, EFT_grhonu_tot, EFT_gpinu_tot, EFT_gpinudot_tot, grhormass_t
        integer  nu_i

        a2 = a*a
        grhob_t=grhob/a
        grhoc_t=grhoc/a
        grhor_t=grhornomass/a2
        grhog_t=grhog/a2
        grhov_t=grhov*EFTw(a,3)
        grho=grhob_t+grhoc_t+grhor_t+grhog_t+grhov_t
        gpres= (grhog_t+grhor_t)/3._dl +EFTw(a,0)*grhov_t
        EFT_H0 = (CP%h0/c_EFT)*1000._dl

        EFT_grhonu_tot = 0._dl
        EFT_gpinu_tot = 0._dl
        if (CP%Num_Nu_Massive /= 0) then
            do nu_i = 1, CP%Nu_mass_eigenstates
                EFT_grhonu    = 0._dl
                EFT_gpinu     = 0._dl
                grhormass_t=grhormass(nu_i)/a2
                call Nu_background(a*nu_masses(nu_i),EFT_grhonu,EFT_gpinu)
                EFT_grhonu_tot = EFT_grhonu_tot + grhormass_t*EFT_grhonu
                EFT_gpinu_tot  = EFT_gpinu_tot + grhormass_t*EFT_gpinu
            end do
        end if

        grho  = grho  + EFT_grhonu_tot
        gpres = gpres + EFT_gpinu_tot
        ! 3) Hubble and derivatives:
        adotoa=sqrt(grho/3)
        adotdota=(adotoa*adotoa-gpres)/2._dl
        Hdot =adotdota-adotoa**2._dl

        EFT_gpinu_tot = 0._dl
        EFT_grhonu_tot = 0._dl
        EFT_gpinudot_tot = 0._dl

        if (CP%Num_Nu_Massive /= 0) then
            do nu_i = 1, CP%Nu_mass_eigenstates
                EFT_grhonu    = 0._dl
                EFT_gpinu     = 0._dl
                grhormass_t=grhormass(nu_i)/a**2
                call Nu_background(a*nu_masses(nu_i),EFT_grhonu,EFT_gpinu)
                EFT_grhonu_tot = EFT_grhonu_tot + grhormass_t*EFT_grhonu
                EFT_gpinu_tot  = EFT_gpinu_tot + grhormass_t*EFT_gpinu
                EFT_gpinudot_tot  = EFT_gpinudot_tot + grhormass_t*(Nu_pidot(a*nu_masses(nu_i),adotoa, EFT_gpinu)&
                    & -4._dl*adotoa*EFT_gpinu)
            end do
        end if

        Hdotdot = 2._dl*adotoa*Hdot &
            & + 0.5_dl*adotoa*(grhob_t + grhoc_t + 8._dl*(grhog_t+grhor_t)/3._dl)&
            & + 0.5_dl*adotoa*grhov_t*((1._dl+EFTw(a,0))*(1._dl+3._dl*EFTw(a,0))-a*EFTw(a,1))&
            & + adotoa/6._dl*EFT_grhonu_tot -0.5_dl*adotoa*EFT_gpinu_tot -0.5_dl*EFT_gpinudot_tot

        ! 4) EFT functions:
        EFTOmegaV   = RPH_EFT_Omega(a,0)
        EFTOmegaP   = RPH_EFT_Omega(a,1)
        EFTOmegaPP  = RPH_EFT_Omega(a,2)
        EFTOmegaPPP = RPH_EFT_Omega(a,3)

        RPH_EFTcdot = -EFTOmegaV*(Hdotdot-4._dl*adotoa*Hdot+2._dl*adotoa*adotoa*adotoa) &
                & + 0.5_dl*a*EFTOmegaP*(-Hdotdot+adotoa*Hdot+adotoa*adotoa*adotoa)&
                & +0.5_dl*a2*adotoa*EFTOmegaPP*(adotoa*adotoa-3._dl*Hdot)&
                & -0.5_dl*a*a2*adotoa*adotoa*adotoa*EFTOmegaPPP&
                & +0.5_dl*adotoa*grhov_t*(-3._dl*(1._dl+EFTw(a,0))**2 + a*EFTw(a,1))

    end function RPH_EFTcdot

    ! ----------------------------

    function RPH_EFT_Gamma1(a, deriv)
        !Function that returns the EFT function gamma1 for RPH.
        implicit none

        integer , intent(in)  :: deriv
        real(dl), intent (in) :: a
        real(dl) :: RPH_EFT_Gamma1
        real(dl) :: EFT_H0, Hubble_temp, PM_temp

        EFT_H0 = (CP%h0/c_EFT)*1000._dl

        select case (deriv)
            case (0)
                RPH_EFT_Gamma1 = 0.25_dl*(RPH_Kineticity(a,0)*(1._dl+ RPH_PlanckMass(a,0))*RPH_Hubble(a)**2 &
                               & -2._dl*RPH_EFTc(a) )/(EFT_H0*EFT_H0*a*a)
            case (1)
                Hubble_temp = RPH_Hubble(a)
                PM_temp     = RPH_PlanckMass(a,0)
                RPH_EFT_Gamma1 = - 0.5_dl*(RPH_Kineticity(a,0)*(1._dl+ PM_temp)*Hubble_temp**2 &
                               & -2._dl*RPH_EFTc(a) )/(EFT_H0*EFT_H0*a*a*a) &
                               & +0.25_dl*( RPH_Kineticity(a,1)*(1._dl+ PM_temp)*Hubble_temp**2 &
                               & +RPH_Kineticity(a,0)*RPH_PlanckMass(a,1)*Hubble_temp**2 &
                               & +2._dl*RPH_Kineticity(a,0)*(1._dl+ PM_temp)*RPH_HubbleDot(a)/a &
                               & -4._dl*RPH_EFTc(a) -2._dl*RPH_EFTcdot(a)/a/Hubble_temp )/(EFT_H0*EFT_H0*a*a)
        end select

    end function RPH_EFT_Gamma1

    ! ----------------------------

    function RPH_EFT_Gamma2(a, deriv)
        !Function that returns the EFT function gamma2 for RPH.
        implicit none

        integer , intent(in)  :: deriv
        real(dl), intent (in) :: a
        real(dl) :: RPH_EFT_Gamma2
        real(dl) :: EFT_H0, Hubble_temp, PM_temp

        EFT_H0 = (CP%h0/c_EFT)*1000._dl

        select case (deriv)
            case (0)
                RPH_EFT_Gamma2 = ( +2._dl*RPH_Braiding(a,0)*(1._dl+ RPH_PlanckMass(a,0)) &
                               & -a*RPH_EFT_Omega(a,1) )*RPH_Hubble(a)/(EFT_H0*a)
            case (1)
                Hubble_temp = RPH_Hubble(a)
                PM_temp     = RPH_PlanckMass(a,0)
                RPH_EFT_Gamma2 = -0.5_dl*( +2._dl*RPH_Braiding(a,0)*(1._dl+ PM_temp) &
                               & -a*RPH_EFT_Omega(a,1) )*Hubble_temp/(EFT_H0*a) &
                               & -( -2._dl*(1._dl+ PM_temp)*( RPH_Braiding(a,1)*Hubble_temp**2 &
                               & + RPH_Braiding(a,0)*RPH_HubbleDot(a)/a) &
                               & - 2._dl*RPH_Braiding(a,0)*Hubble_temp**2*RPH_PlanckMass(a,1) &
                               & + RPH_EFT_Omega(a,1)*( Hubble_temp**2 + RPH_HubbleDot(a) ) &
                               & + a*Hubble_temp**2*RPH_EFT_Omega(a,2) )/(EFT_H0*a*Hubble_temp)
        end select

    end function RPH_EFT_Gamma2

    ! ----------------------------

    function RPH_EFT_Gamma3(a, deriv)
        !Function that returns the EFT function gamma3 for RPH.
        implicit none

        integer , intent(in)  :: deriv
        real(dl), intent (in) :: a
        real(dl) :: RPH_EFT_Gamma3
        real(dl) :: EFT_H0

        select case (deriv)
            case (0)
                RPH_EFT_Gamma3 = -RPH_Tensor(a,0)*(1._dl +RPH_PlanckMass(a,0))
            case (1)
                RPH_EFT_Gamma3 = -RPH_PlanckMass(a,1)*RPH_Tensor(a,0) &
                               & -(1._dl +RPH_PlanckMass(a,0))*RPH_Tensor(a,1)
            case (2)
                RPH_EFT_Gamma3 = -(1._dl +RPH_PlanckMass(a,0))*RPH_Tensor(a,2) &
                               & -RPH_PlanckMass(a,2)*RPH_Tensor(a,0) &
                               & -2._dl*RPH_PlanckMass(a,1)*RPH_Tensor(a,1)
        end select

    end function RPH_EFT_Gamma3

end module RPH_to_EFT

! -------------------------------------------------------------------------------------------------
