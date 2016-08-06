! -------------------------------------------------------------------------------------------------
!
!   EFTCAMB
!
!   Developed and implemented by:
!       Bin Hu (hu@lorentz.leidenuniv.nl), Marco Raveri (mraveri@sissa.it)
!       Noemi Frusciante (fruscian@iap.fr), Alessandra Silvestri (silvestri@lorentz.leidenuniv.nl)
!
!
!   EFTCAMB module that contains the code that looks for the return to GR of a model,
!   checks the stability of the model and initializes the EFTCAMB code.
!
!   For more informations about the methods contained in this file see
!   the documentation: arXiv:1405.3590
!
! -------------------------------------------------------------------------------------------------

! Code used to detect the return to GR (RGR) of the considered model.
module EFTreturntoGR
    use EFTDef
    use EFTfunctions
    use MassiveNu
    implicit none

    !We need to define this to interface the EFT functions with the RGR code.
    abstract interface
        function EFTfunction(a,ind) result(y)
            use precision
            implicit none
            real(dl), intent(in) :: a
            integer , intent(in) :: ind
            real(dl)  y
        end function EFTfunction
    end interface

    ! This two variables stores the pointer to the specific EFT function which is considered
    ! by the RGR code.
    integer :: deriv
    procedure(EFTfunction), pointer :: EFTfuncTemp

contains

    ! ----------------------------

    subroutine EFTCheckReturnToGR
        ! RGR test main subroutine.
        implicit none

        real(dl) :: OmegaTime=0._dl, OmegaDotTime=0._dl
        real(dl) :: EFTcTime=0._dl, EFTcDotTime=0._dl, EFTLambdaTime=0._dl, EFTLambdaDotTime=0._dl
        real(dl) :: EFTGamma1Time=0._dl, EFTGamma1DotTime=0._dl
        real(dl) :: EFTGamma2Time=0._dl, EFTGamma2DotTime=0._dl
        real(dl) :: EFTGamma3Time=0._dl, EFTGamma3DotTime=0._dl
        real(dl) :: EFTGamma4Time=0._dl, EFTGamma4DotTime=0._dl
        real(dl) :: EFTGamma5Time=0._dl, EFTGamma5DotTime=0._dl
        real(dl) :: EFTGamma6Time=0._dl, EFTGamma6DotTime=0._dl

        real(dl) :: temp

        if (Feedbacklevel>1) write(*,*) 'EFTCAMB: checking return to GR'

        ! 1) Test Omega:
        deriv = 0
        EFTfuncTemp => EFTOmega
        call EFTdetectionRGR(EFTfunctionRGR, EFTturnonpiInitial, 1._dl, OmegaTime)
        ! 2) Test Omega Dot:
        call EFTdetectionRGR(RGR_EFTOmegaDot, EFTturnonpiInitial, 1._dl, OmegaDotTime)
        ! 3) Test c:
        call EFTdetectionRGR(RGR_EFTc, EFTturnonpiInitial, 1._dl, EFTcTime)
        ! 4) Test c Dot:
        call EFTdetectionRGR(RGR_EFTcDot, EFTturnonpiInitial, 1._dl, EFTcDotTime)
        ! 5) Test Lambda:
        call EFTdetectionRGR(RGR_EFTLambda, EFTturnonpiInitial, 1._dl, EFTLambdaTime)
        ! 6) Test Lambda Dot:
        call EFTdetectionRGR(RGR_EFTLambdaDot, EFTturnonpiInitial, 1._dl, EFTLambdaDotTime)
        ! 7) Test Gamma1:
        deriv = 0
        EFTfuncTemp => EFTGamma1
        call EFTdetectionRGR(EFTfunctionRGR, EFTturnonpiInitial, 1._dl, EFTGamma1Time)
        ! 8) Test Gamma2:
        deriv = 0
        EFTfuncTemp => EFTGamma2
        call EFTdetectionRGR(EFTfunctionRGR, EFTturnonpiInitial, 1._dl, EFTGamma2Time)
        ! 9) Test Gamma3:
        deriv = 0
        EFTfuncTemp => EFTGamma3
        call EFTdetectionRGR(EFTfunctionRGR, EFTturnonpiInitial, 1._dl, EFTGamma3Time)
        ! 10) Test Gamma4:
        deriv = 0
        EFTfuncTemp => EFTGamma4
        call EFTdetectionRGR(EFTfunctionRGR, EFTturnonpiInitial, 1._dl, EFTGamma4Time)
        ! 11) Test Gamma5:
        deriv = 0
        EFTfuncTemp => EFTGamma5
        call EFTdetectionRGR(EFTfunctionRGR, EFTturnonpiInitial, 1._dl, EFTGamma5Time)
        ! 12) Test Gamma6:
        deriv = 0
        EFTfuncTemp => EFTGamma6
        call EFTdetectionRGR(EFTfunctionRGR, EFTturnonpiInitial, 1._dl, EFTGamma6Time)
        ! 13) Test Gamma1 prime:
        deriv = 1
        EFTfuncTemp => EFTGamma1
        call EFTdetectionRGR(EFTfunctionRGR, EFTturnonpiInitial, 1._dl, EFTGamma1DotTime)
        ! 14) Test Gamma2 prime:
        deriv = 1
        EFTfuncTemp => EFTGamma2
        call EFTdetectionRGR(EFTfunctionRGR, EFTturnonpiInitial, 1._dl, EFTGamma2DotTime)
        ! 15) Test Gamma3 prime:
        deriv = 1
        EFTfuncTemp => EFTGamma3
        call EFTdetectionRGR(EFTfunctionRGR, EFTturnonpiInitial, 1._dl, EFTGamma3DotTime)
        ! 16) Test Gamma4 prime:
        deriv = 1
        EFTfuncTemp => EFTGamma4
        call EFTdetectionRGR(EFTfunctionRGR, EFTturnonpiInitial, 1._dl, EFTGamma4DotTime)
        ! 17) Test Gamma5 prime:
        deriv = 1
        EFTfuncTemp => EFTGamma5
        call EFTdetectionRGR(EFTfunctionRGR, EFTturnonpiInitial, 1._dl, EFTGamma5DotTime)
        ! 18) Test Gamma6 prime:
        deriv = 1
        EFTfuncTemp => EFTGamma6
        call EFTdetectionRGR(EFTfunctionRGR, EFTturnonpiInitial, 1._dl, EFTGamma6DotTime)

        ! 19) Now determine EFTturnonpi

        temp = min(OmegaTime,OmegaDotTime, EFTcTime, EFTcDotTime, EFTLambdaTime, EFTLambdaDotTime,&
            &EFTGamma1Time, EFTGamma1DotTime, EFTGamma2Time, EFTGamma2DotTime,&
            &EFTGamma3Time, EFTGamma3DotTime, EFTGamma4Time, EFTGamma4DotTime,&
            &EFTGamma5Time, EFTGamma5DotTime, EFTGamma6Time, EFTGamma6DotTime)

        if (Feedbacklevel>2) then
            write(*,*) 'EFTCAMB: RGR times'
            write(*,*) OmegaTime,OmegaDotTime, EFTcTime, EFTcDotTime, EFTLambdaTime, EFTLambdaDotTime,&
                &EFTGamma1Time, EFTGamma1DotTime, EFTGamma2Time, EFTGamma2DotTime,&
                &EFTGamma3Time, EFTGamma3DotTime, EFTGamma4Time, EFTGamma4DotTime,&
                &EFTGamma5Time, EFTGamma5DotTime, EFTGamma6Time, EFTGamma6DotTime
        end if

        if (temp.le.EFTturnonpiInitial) then
            EFTturnonpi = EFTturnonpiInitial
        else if (temp>EFTturnonpiInitial) then
            EFTturnonpi = temp
        end if

        if (Feedbacklevel>1) write(*,*) 'RGR: a=', EFTturnonpi

    end subroutine EFTCheckReturnToGR

    ! ----------------------------

    function EFTfunctionRGR(a)
        ! This function is needed to make the two input variable EFT functions into a single input functions.
        implicit none
        real(dl), intent(in) :: a
        real(dl) :: EFTfunctionRGR
        EFTfunctionRGR = ABS(EFTfuncTemp(a,deriv))
    end function EFTfunctionRGR

    ! ----------------------------

    subroutine EFTdetectionRGR(EFTfunction,TimeStart, TimeEnd, RGRtime)
        ! Subroutine that test the return to GR of the EFT functions.
        implicit none
        external EFTfunction
        real(dl) :: EFTfunction
        real(dl), intent(in)  :: TimeStart, TimeEnd
        real(dl), intent(out) :: RGRtime
        real(dl) :: zbrent
        external zbrent

        real(dl) :: temp1, temp2, temp3
        integer  :: j
        logical :: succes

        temp1 = EFTfunction(TimeStart) - EFTtoGR
        temp2 = EFTfunction(TimeEnd) - EFTtoGR

        if (temp1*temp2>0.and.temp2>0) then         ! always not GR
            RGRtime = TimeStart
            return
        else if (temp1*temp2>0.and.temp2<0) then    ! Omega is always GR
            RGRtime = 1.1_dl*TimeEnd
            return
        else if (temp1*temp2<0 .and. temp2<0) then  ! Omega is not GR at early times
            RGRtime = TimeStart
            return
        else if (temp1*temp2<0.and. temp2>0) then   ! There is a transition
            RGRtime = zbrent(EFTfunction,TimeStart,TimeEnd,1.d-10,EFTtoGR,succes)
            if (.not.succes) RGRtime = TimeStart
            return
        end if

        RGRtime = TimeStart
        return

    end subroutine EFTdetectionRGR

    ! ----------------------------

    function RGR_EFTOmegaDot(a)
        !Dummy function for the absolute value of EFTOmegaDot
        implicit none
        real(dl) a,RGR_EFTOmegaDot, EFT_H0
        real(dl) a2, grhob_t,grhoc_t,grhor_t,grhog_t,grhov_t,grho,adotoa
        real(dl) EFT_grhonu, EFT_gpinu, grhormass_t, EFT_grhonu_tot
        integer  nu_i
        a2 = a*a
        grhob_t=grhob/a
        grhoc_t=grhoc/a
        grhor_t=grhornomass/a2
        grhog_t=grhog/a2
        grhov_t=grhov*EFTw(a,3)
        grho=grhob_t+grhoc_t+grhor_t+grhog_t+grhov_t
        EFT_grhonu_tot = 0._dl
        if (CP%Num_Nu_Massive /= 0) then
            do nu_i = 1, CP%Nu_mass_eigenstates
                EFT_grhonu    = 0._dl
                EFT_gpinu     = 0._dl
                grhormass_t=grhormass(nu_i)/a**2
                call Nu_background(a*nu_masses(nu_i),EFT_grhonu,EFT_gpinu)
                EFT_grhonu_tot = EFT_grhonu_tot + grhormass_t*EFT_grhonu
            end do
        end if
        grho  = grho  + EFT_grhonu_tot
        adotoa=sqrt(grho/3)
        EFT_H0 = (CP%h0/c_EFT)*1000._dl
        RGR_EFTOmegaDot = ABS(a*adotoa*EFTOmega(a,1))
        return
    end function RGR_EFTOmegaDot

    ! ----------------------------

    function RGR_EFTc(a)
        !Dummy function for the absolute value of EFTc
        implicit none
        real(dl) a,RGR_EFTc
        real(dl) a2, grhob_t,grhoc_t,grhor_t,grhog_t,grhov_t,grho, gpres, EFT_H0
        real(dl) adotoa, adotdota, Hdot
        real(dl) EFTOmegaV, EFTOmegaP, EFTOmegaPP
        real(dl) EFTc
        real(dl) EFT_grhonu, EFT_gpinu, EFT_grhonu_tot, EFT_gpinu_tot, grhormass_t
        integer  nu_i
        EFT_H0 = (CP%h0/c_EFT)*1000._dl
        a2 = a*a
        grhob_t=grhob/a
        grhoc_t=grhoc/a
        grhor_t=grhornomass/a2
        grhog_t=grhog/a2
        grhov_t=grhov*EFTw(a,3)
        if ( CP%EFTflag==4 ) then ! Full mapping EFT overwrites some of these quantities:
            if ( CP%FullMappingEFTmodel==1 ) then ! Horava gravity
                grhov_t = ( grhov &
                    &   +3._dl*EFT_H0**2*((CP%Horava_eta +3._dl*CP%Horava_lambda-2._dl*CP%Horava_xi)/(2._dl*CP%Horava_xi+2._dl-CP%Horava_eta)))*EFTw(a,3)
            end if
        end if
        grho=grhob_t+grhoc_t+grhor_t+grhog_t+grhov_t
        gpres= (grhog_t+grhor_t)/3._dl +EFTw(a,0)*grhov_t
        EFT_gpinu_tot = 0._dl
        EFT_grhonu_tot = 0._dl
        if (CP%Num_Nu_Massive /= 0) then
            do nu_i = 1, CP%Nu_mass_eigenstates
                EFT_grhonu    = 0._dl
                EFT_gpinu     = 0._dl
                grhormass_t=grhormass(nu_i)/a**2
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
        if ( CP%EFTflag==4 ) then ! Full mapping EFT overwrites some of these quantities:
            if ( CP%FullMappingEFTmodel==1 ) then ! Horava gravity
                adotoa   = adotoa*sqrt(( 2._dl +2._dl*CP%Horava_xi -CP%Horava_eta )/(3._dl*CP%Horava_lambda +2._dl))
                adotdota = adotdota*( 2._dl +2._dl*CP%Horava_xi -CP%Horava_eta )/(3._dl*CP%Horava_lambda +2._dl)
                Hdot     = Hdot*( 2._dl +2._dl*CP%Horava_xi -CP%Horava_eta )/(3._dl*CP%Horava_lambda +2._dl)
            end if
        end if
        ! 4) EFT functions:
        EFTOmegaV  = EFTOmega(a,0)
        EFTOmegaP  = EFTOmega(a,1)
        EFTOmegaPP = EFTOmega(a,2)

        select case (CP%EFTflag) ! 8*pi*G*c*a^2
            case (1) ! Pure EFT
                EFTc = (adotoa*adotoa - Hdot)*(EFTOmegaV + a*EFTOmegaP*0.5_dl) &
                    & - 0.5_dl*a2*adotoa*adotoa*EFTOmegaPP&
                    & + 0.5_dl*grhov_t*(1._dl+EFTw(a,0))
            case (2) ! Designer mapping
                select case (CP%DesignerEFTModel)
                    case (2) ! Designer mc quintessence
                        EFTc = (adotoa*adotoa - Hdot)*(EFTOmegaV + a*EFTOmegaP*0.5_dl) &
                            & - 0.5_dl*a2*adotoa*adotoa*EFTOmegaPP&
                            & + 0.5_dl*grhov_t*(1._dl+EFTw(a,0))
                    case default ! normally we can just use the designed one.
                        EFTc = EFTcTemp(a,0)
                end select
            case (3) ! EFT alternative parametrizations
                select case (CP%AltParEFTmodel)
                    case (1) ! RPH
                        EFTc = (adotoa*adotoa - Hdot)*(EFTOmegaV + a*EFTOmegaP*0.5_dl) &
                            & - 0.5_dl*a2*adotoa*adotoa*EFTOmegaPP&
                            & + 0.5_dl*grhov_t*(1._dl+EFTw(a,0))
                    case default
                        stop 'Wrong selection of model'
                end select
            case (4) ! Full mapping EFT
                select case (CP%FullMappingEFTmodel)
                    case (1) ! Horava gravity
                        EFTc = -( 2._dl*CP%Horava_xi - 3._dl*CP%Horava_lambda )*( Hdot - adotoa**2 )/( 2._dl +2._dl*CP%Horava_xi -CP%Horava_eta )
                end select
        end select

        RGR_EFTc = ABS(EFTc/a2)
        return
    end function RGR_EFTc

    ! ----------------------------

    function RGR_EFTcDot(a)
        ! Dummy function for the absolute value of EFTcdot
        implicit none
        real(dl) a,RGR_EFTcDot
        real(dl) a2, grhob_t,grhoc_t,grhor_t,grhog_t,grhov_t,grho, gpres, EFT_H0
        real(dl) adotoa, adotdota, Hdot,Hdotdot
        real(dl) EFTOmegaV, EFTOmegaP, EFTOmegaPP,EFTOmegaPPP
        real(dl) EFTc, EFTcdot, EFTLambda, EFTLambdadot
        real(dl) EFT_grhonu, EFT_gpinu, EFT_grhonu_tot, EFT_gpinu_tot, EFT_gpinudot_tot, grhormass_t
        integer  nu_i
        EFT_H0 = (CP%h0/c_EFT)*1000._dl
        a2 = a*a
        grhob_t=grhob/a
        grhoc_t=grhoc/a
        grhor_t=grhornomass/a2
        grhog_t=grhog/a2
        grhov_t=grhov*EFTw(a,3)
        if ( CP%EFTflag==4 ) then ! Full mapping EFT overwrites some of these quantities:
            if ( CP%FullMappingEFTmodel==1 ) then ! Horava gravity
                grhov_t = ( grhov &
                    &   +3._dl*EFT_H0**2*((CP%Horava_eta +3._dl*CP%Horava_lambda-2._dl*CP%Horava_xi)/(2._dl*CP%Horava_xi+2._dl-CP%Horava_eta)))*EFTw(a,3)
            end if
        end if
        grho=grhob_t+grhoc_t+grhor_t+grhog_t+grhov_t
        gpres= (grhog_t+grhor_t)/3._dl +EFTw(a,0)*grhov_t

        EFT_grhonu_tot = 0._dl
        EFT_gpinu_tot = 0._dl
        if (CP%Num_Nu_Massive /= 0) then
            do nu_i = 1, CP%Nu_mass_eigenstates
                EFT_grhonu    = 0._dl
                EFT_gpinu     = 0._dl
                grhormass_t=grhormass(nu_i)/a**2
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

        if ( CP%EFTflag==4 ) then ! Full mapping EFT overwrites some of these quantities:
            if ( CP%FullMappingEFTmodel==1 ) then ! Horava gravity
                adotoa   = adotoa*sqrt(( 2._dl +2._dl*CP%Horava_xi -CP%Horava_eta )/(3._dl*CP%Horava_lambda +2._dl))
                adotdota = adotdota*( 2._dl +2._dl*CP%Horava_xi -CP%Horava_eta )/(3._dl*CP%Horava_lambda +2._dl)
                Hdot     = Hdot*( 2._dl +2._dl*CP%Horava_xi -CP%Horava_eta )/(3._dl*CP%Horava_lambda +2._dl)
            end if
        end if

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

        if ( CP%EFTflag==4 ) then ! Full mapping EFT overwrites some of these quantities:
            if ( CP%FullMappingEFTmodel==1 ) then ! Horava gravity
                Hdotdot = ( 2._dl +2._dl*CP%Horava_xi -CP%Horava_eta )/(3._dl*CP%Horava_lambda +2._dl)*&
                    &( +adotoa*grhov_t*2._dl/3._dl                                                 &
                    &  +adotoa*( -grhob_t -grhoc_t -4._dl*(grhog_t+grhor_t) )/3._dl                &
                    &  +adotoa/6._dl*EFT_grhonu_tot -0.5_dl*adotoa*EFT_gpinu_tot -0.5_dl*EFT_gpinudot_tot )
            end if
        end if

        ! 4) EFT functions:
        EFTOmegaV = EFTOmega(a,0)
        EFTOmegaP = EFTOmega(a,1)
        EFTOmegaPP = EFTOmega(a,2)
        EFTOmegaPPP = EFTOmega(a,3)

        select case (CP%EFTflag) ! 8*pi*G*cdot*a^2
            case (1) ! Pure EFT
                EFTcdot = -EFTOmegaV*(Hdotdot-4._dl*adotoa*Hdot+2._dl*adotoa*adotoa*adotoa) &
                    & + 0.5_dl*a*EFTOmegaP*(-Hdotdot+adotoa*Hdot+adotoa*adotoa*adotoa)&
                    & +0.5_dl*a2*adotoa*EFTOmegaPP*(adotoa*adotoa-3._dl*Hdot)&
                    & -0.5_dl*a*a2*adotoa*adotoa*adotoa*EFTOmegaPPP&
                    & +0.5_dl*adotoa*grhov_t*(-3._dl*(1._dl+EFTw(a,0))**2 + a*EFTw(a,1))
            case (2) ! Designer mapping
                select case (CP%DesignerEFTModel)
                    case (2) ! Designer mc quintessence:
                        EFTcdot = -EFTOmegaV*(Hdotdot-4._dl*adotoa*Hdot+2._dl*adotoa*adotoa*adotoa) &
                            & + 0.5_dl*a*EFTOmegaP*(-Hdotdot+adotoa*Hdot+adotoa*adotoa*adotoa)&
                            & +0.5_dl*a2*adotoa*EFTOmegaPP*(adotoa*adotoa-3._dl*Hdot)&
                            & -0.5_dl*a*a2*adotoa*adotoa*adotoa*EFTOmegaPPP&
                            & +0.5_dl*adotoa*grhov_t*(-3._dl*(1._dl+EFTw(a,0))**2 + a*EFTw(a,1))
                    case default
                        EFTcdot = EFTcTemp(a,1)
                end select
            case (3) ! EFT alternative parametrizations
                select case (CP%AltParEFTmodel)
                    case (1) ! RPH:
                        EFTcdot = -EFTOmegaV*(Hdotdot-4._dl*adotoa*Hdot+2._dl*adotoa*adotoa*adotoa) &
                            & + 0.5_dl*a*EFTOmegaP*(-Hdotdot+adotoa*Hdot+adotoa*adotoa*adotoa)&
                            & +0.5_dl*a2*adotoa*EFTOmegaPP*(adotoa*adotoa-3._dl*Hdot)&
                            & -0.5_dl*a*a2*adotoa*adotoa*adotoa*EFTOmegaPPP&
                            & +0.5_dl*adotoa*grhov_t*(-3._dl*(1._dl+EFTw(a,0))**2 + a*EFTw(a,1))
                    case default
                        stop 'Wrong selection of model'
                end select
            case (4) ! Full mapping EFT
                select case (CP%FullMappingEFTmodel)
                    case (1) ! Horava gravity
                        EFTcdot = ( 3._dl*CP%Horava_lambda - 2._dl*CP%Horava_xi )*&
                            &( Hdotdot -4._dl*adotoa*Hdot + 2._dl*adotoa**3 )/( 2._dl +2._dl*CP%Horava_xi -CP%Horava_eta )
                end select
        end select

        RGR_EFTcDot = ABS(EFTcdot/a2)
        return
    end function RGR_EFTcDOT

    ! ----------------------------

    function RGR_EFTLambda(a)
        ! Dummy function for the absolute value of EFTLambda
        implicit none
        real(dl) a,RGR_EFTLambda
        real(dl) a2, grhob_t,grhoc_t,grhor_t,grhog_t,grhov_t,grho, gpres, EFT_H0
        real(dl) adotoa, adotdota, Hdot,Hdotdot
        real(dl) EFTOmegaV, EFTOmegaP, EFTOmegaPP,EFTOmegaPPP
        real(dl)  EFTc, EFTcdot, EFTLambda, EFTLambdadot
        real(dl) EFT_grhonu, EFT_gpinu, grhormass_t
        real(dl) EFT_grhonu_tot, EFT_gpinu_tot, EFT_gpinudot_tot
        integer  nu_i
        EFT_H0 = (CP%h0/c_EFT)*1000._dl
        a2 = a*a
        grhob_t=grhob/a
        grhoc_t=grhoc/a
        grhor_t=grhornomass/a2
        grhog_t=grhog/a2
        grhov_t=grhov*EFTw(a,3)
        if ( CP%EFTflag==4 ) then ! Full mapping EFT overwrites some of these quantities:
            if ( CP%FullMappingEFTmodel==1 ) then ! Horava gravity
                grhov_t = ( grhov &
                    &   +3._dl*EFT_H0**2*((CP%Horava_eta +3._dl*CP%Horava_lambda-2._dl*CP%Horava_xi)/(2._dl*CP%Horava_xi+2._dl-CP%Horava_eta)))*EFTw(a,3)
            end if
        end if
        grho=grhob_t+grhoc_t+grhor_t+grhog_t+grhov_t
        gpres= (grhog_t+grhor_t)/3._dl +EFTw(a,0)*grhov_t

        EFT_gpinu_tot = 0._dl
        EFT_grhonu_tot = 0._dl

        if (CP%Num_Nu_Massive /= 0) then
            do nu_i = 1, CP%Nu_mass_eigenstates
                EFT_grhonu    = 0._dl
                EFT_gpinu     = 0._dl
                grhormass_t=grhormass(nu_i)/a**2
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

        if ( CP%EFTflag==4 ) then ! Full mapping EFT overwrites some of these quantities:
            if ( CP%FullMappingEFTmodel==1 ) then ! Horava gravity
                adotoa   = adotoa*sqrt(( 2._dl +2._dl*CP%Horava_xi -CP%Horava_eta )/(3._dl*CP%Horava_lambda +2._dl))
                adotdota = adotdota*( 2._dl +2._dl*CP%Horava_xi -CP%Horava_eta )/(3._dl*CP%Horava_lambda +2._dl)
                Hdot     = Hdot*( 2._dl +2._dl*CP%Horava_xi -CP%Horava_eta )/(3._dl*CP%Horava_lambda +2._dl)
            end if
        end if

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

        if ( CP%EFTflag==4 ) then ! Full mapping EFT overwrites some of these quantities:
            if ( CP%FullMappingEFTmodel==1 ) then ! Horava gravity
                Hdotdot = ( 2._dl +2._dl*CP%Horava_xi -CP%Horava_eta )/(3._dl*CP%Horava_lambda +2._dl)*&
                    &( +adotoa*grhov_t*2._dl/3._dl                                                 &
                    &  +adotoa*( -grhob_t -grhoc_t -4._dl*(grhog_t+grhor_t) )/3._dl                &
                    &  +adotoa/6._dl*EFT_grhonu_tot -0.5_dl*adotoa*EFT_gpinu_tot -0.5_dl*EFT_gpinudot_tot )
            end if
        end if

        ! 4) EFT functions:
        EFTOmegaV = EFTOmega(a,0)
        EFTOmegaP = EFTOmega(a,1)
        EFTOmegaPP = EFTOmega(a,2)
        EFTOmegaPPP = EFTOmega(a,3)

        select case (CP%EFTflag) !8*pi*G*Lambda*a^2
            case (1) ! Pure EFT
                EFTLambda = +EFTw(a,0)*grhov_t &
                    &-EFTOmegaV*(2._dl*Hdot+adotoa**2._dl) &
                    &-a*EFTOmegaP*(2._dl*adotoa**2._dl + Hdot) &
                    &-a2*adotoa**2._dl*EFTOmegaPP
            case (2) ! Designer mapping
                select case (CP%DesignerEFTModel)
                    case (2) ! Designer mc quintessence
                        EFTLambda = +EFTw(a,0)*grhov_t &
                            &-EFTOmegaV*(2._dl*Hdot+adotoa**2._dl) &
                            &-a*EFTOmegaP*(2._dl*adotoa**2._dl + Hdot) &
                            &-a2*adotoa**2._dl*EFTOmegaPP
                    case default
                        EFTLambda = EFTLambdaTemp(a,0)
                end select
            case (3) ! EFT alternative parametrizations
                select case (CP%AltParEFTmodel)
                    case (1) ! RPH:
                        EFTLambda = +EFTw(a,0)*grhov_t &
                            &-EFTOmegaV*(2._dl*Hdot+adotoa**2._dl) &
                            &-a*EFTOmegaP*(2._dl*adotoa**2._dl + Hdot) &
                            &-a2*adotoa**2._dl*EFTOmegaPP
                    case default
                        stop 'Wrong selection of model'
                end select
            case (4) ! Full mapping EFT
                select case (CP%FullMappingEFTmodel)
                    case (1) ! Horava gravity
                        EFTLambda = +EFTw(a,0)*grhov_t + 2._dl*( 3._dl*CP%Horava_lambda - 2._dl*CP%Horava_xi)*&
                            &( 0.5_dl*adotoa**2 + Hdot)/( 2._dl +2._dl*CP%Horava_xi -CP%Horava_eta )
                end select
        end select

        RGR_EFTLambda = ABS(EFTLambda/a2 + grhov)
        return
    end function RGR_EFTLambda

    ! ----------------------------

    function RGR_EFTLambdaDot(a)
        ! Dummy function for the absolute value of EFTLambdaDot
        implicit none
        real(dl) a,RGR_EFTLambdaDot
        real(dl) a2, grhob_t,grhoc_t,grhor_t,grhog_t,grhov_t,grho, gpres, EFT_H0
        real(dl) adotoa, adotdota, Hdot,Hdotdot
        real(dl) EFTOmegaV, EFTOmegaP, EFTOmegaPP,EFTOmegaPPP
        real(dl) EFTc, EFTcdot, EFTLambda, EFTLambdadot
        real(dl) EFT_grhonu, EFT_gpinu, grhormass_t
        real(dl) EFT_grhonu_tot, EFT_gpinu_tot, EFT_gpinudot_tot
        integer  nu_i
        EFT_H0 = (CP%h0/c_EFT)*1000._dl
        a2 = a*a
        grhob_t=grhob/a
        grhoc_t=grhoc/a
        grhor_t=grhornomass/a2
        grhog_t=grhog/a2
        grhov_t=grhov*EFTw(a,3)
        if ( CP%EFTflag==4 ) then ! Full mapping EFT overwrites some of these quantities:
            if ( CP%FullMappingEFTmodel==1 ) then ! Horava gravity
                grhov_t = ( grhov &
                    &   +3._dl*EFT_H0**2*((CP%Horava_eta +3._dl*CP%Horava_lambda-2._dl*CP%Horava_xi)/(2._dl*CP%Horava_xi+2._dl-CP%Horava_eta)))*EFTw(a,3)
            end if
        end if
        grho=grhob_t+grhoc_t+grhor_t+grhog_t+grhov_t
        gpres= (grhog_t+grhor_t)/3._dl +EFTw(a,0)*grhov_t

        EFT_gpinu_tot = 0._dl
        EFT_grhonu_tot = 0._dl

        if (CP%Num_Nu_Massive /= 0) then
            do nu_i = 1, CP%Nu_mass_eigenstates
                EFT_grhonu    = 0._dl
                EFT_gpinu     = 0._dl
                grhormass_t=grhormass(nu_i)/a**2
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

        if ( CP%EFTflag==4 ) then ! Full mapping EFT overwrites some of these quantities:
            if ( CP%FullMappingEFTmodel==1 ) then ! Horava gravity
                adotoa   = adotoa*sqrt(( 2._dl +2._dl*CP%Horava_xi -CP%Horava_eta )/(3._dl*CP%Horava_lambda +2._dl))
                adotdota = adotdota*( 2._dl +2._dl*CP%Horava_xi -CP%Horava_eta )/(3._dl*CP%Horava_lambda +2._dl)
                Hdot     = Hdot*( 2._dl +2._dl*CP%Horava_xi -CP%Horava_eta )/(3._dl*CP%Horava_lambda +2._dl)
            end if
        end if

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

        if ( CP%EFTflag==4 ) then ! Full mapping EFT overwrites some of these quantities:
            if ( CP%FullMappingEFTmodel==1 ) then ! Horava gravity
                Hdotdot = ( 2._dl +2._dl*CP%Horava_xi -CP%Horava_eta )/(3._dl*CP%Horava_lambda +2._dl)*&
                    &( +adotoa*grhov_t*2._dl/3._dl                                                 &
                    &  +adotoa*( -grhob_t -grhoc_t -4._dl*(grhog_t+grhor_t) )/3._dl                &
                    &  +adotoa/6._dl*EFT_grhonu_tot -0.5_dl*adotoa*EFT_gpinu_tot -0.5_dl*EFT_gpinudot_tot )
            end if
        end if

        ! 4) EFT functions:
        EFTOmegaV = EFTOmega(a,0)
        EFTOmegaP = EFTOmega(a,1)
        EFTOmegaPP = EFTOmega(a,2)
        EFTOmegaPPP = EFTOmega(a,3)

        select case (CP%EFTflag) ! 8*pi*G*Ldot*a^2
            case (1) ! Pure EFT
                EFTLambdadot = -2._dl*EFTOmegaV*(Hdotdot-adotoa*Hdot-adotoa*adotoa*adotoa)&
                    & - a*EFTOmegaP*(4._dl*adotoa*Hdot+Hdotdot)&
                    & -a2*EFTOmegaPP*adotoa*(3._dl*Hdot+2._dl*adotoa*adotoa)&
                    & -a*a2*EFTOmegaPPP*adotoa*adotoa*adotoa&
                    & +grhov_t*adotoa*(a*EFTw(a,1)-3._dl*EFTw(a,0)*(1._dl+EFTw(a,0)))
            case (2) ! Designer mapping
                select case (CP%DesignerEFTModel)
                    case (2) ! Designer mc quintessence
                        EFTLambdadot = -2._dl*EFTOmegaV*(Hdotdot-adotoa*Hdot-adotoa*adotoa*adotoa)&
                            & - a*EFTOmegaP*(4._dl*adotoa*Hdot+Hdotdot)&
                            & -a2*EFTOmegaPP*adotoa*(3._dl*Hdot+2._dl*adotoa*adotoa)&
                            & -a*a2*EFTOmegaPPP*adotoa*adotoa*adotoa&
                            & +grhov_t*adotoa*(a*EFTw(a,1)-3._dl*EFTw(a,0)*(1._dl+EFTw(a,0)))
                    case default
                        EFTLambdadot = EFTLambdaTemp(a,1)
                end select
            case (3) ! EFT alternative parametrizations
                select case (CP%AltParEFTmodel)
                    case (1) ! RPH:
                        EFTLambdadot = -2._dl*EFTOmegaV*(Hdotdot-adotoa*Hdot-adotoa*adotoa*adotoa)&
                            & - a*EFTOmegaP*(4._dl*adotoa*Hdot+Hdotdot)&
                            & -a2*EFTOmegaPP*adotoa*(3._dl*Hdot+2._dl*adotoa*adotoa)&
                            & -a*a2*EFTOmegaPPP*adotoa*adotoa*adotoa&
                            & +grhov_t*adotoa*(a*EFTw(a,1)-3._dl*EFTw(a,0)*(1._dl+EFTw(a,0)))
                    case default
                        stop 'Wrong selection of model'
                end select
            case (4) ! Full mapping EFT
                select case (CP%FullMappingEFTmodel)
                    case (1) ! Horava gravity
                        EFTLambdadot = + 2._dl*( 3._dl*CP%Horava_lambda - 2._dl*CP%Horava_xi)*&
                            &( Hdotdot -adotoa*Hdot -adotoa**3 )/( 2._dl +2._dl*CP%Horava_xi -CP%Horava_eta )
                end select
        end select

        RGR_EFTLambdaDot = ABS(EFTLambdadot/a2)
        return
    end function RGR_EFTLambdaDot

end module EFTreturntoGR

! -------------------------------------------------------------------------------------------------

!Performs the stability check for the considered model.
module EFTstability
    use EFTDef
    use EFTfunctions
    use MassiveNu
    implicit none

    ! storage for some utility values that needs to be stored for the stability check.
    real(dl), save :: PastA1 = 0._dl
    real(dl), save :: PastAT = 0._dl

contains

    ! ----------------------------

    subroutine EFTCheck_Stability(success, astart, aend)
        ! CheckStability subroutine.
        ! Three scanning strategies: linear sampling from astart to aend and two log samplig on the borders.
        implicit none

        logical  :: success
        real(dl), intent(in) :: astart, aend

        ! parameters of the stability sampler:
        integer , parameter :: indMax            = 10000   ! Number of points sampled.
        real(dl), parameter :: LogSamplingScale  = -10._dl ! Where to start with the log sampling

        real(dl) :: Atest, y
        integer  :: ind

        if (Feedbacklevel>1) write(*,*) 'EFTCAMB: checking stability of the theory'

        ! 1) Parallel stability code:
        success = .true.
        !$omp parallel do private(ind, Atest, y) shared(success) schedule(static)
        do ind=1, indMax

            if (success) then
                ! linear sampling:
                Atest = astart + REAL(ind-1)*(aend-astart)/REAL(indMax-1)
                if (.not.EFTStabilityComputation(Atest)) then
                    success = .false.
                    if (Feedbacklevel>2) write(*,*) 'EFTCAMB: instability detected at a=', Atest
                end if
                ! log sampling close to astart
                y = LogSamplingScale + REAL(ind-1)*(0._dl-LogSamplingScale)/REAL(indMax-1)

                Atest = astart +(aend-astart)*10._dl**y
                if (.not.EFTStabilityComputation(Atest)) then
                    success = .false.
                    if (Feedbacklevel>2) write(*,*) 'EFTCAMB: instability detected at a=', Atest
                end if
                ! log sampling close to aend
                Atest = aend +(astart-aend)*10._dl**y
                if (.not.EFTStabilityComputation(Atest)) then
                    success = .false.
                    if (Feedbacklevel>2) write(*,*) 'EFTCAMB: instability detected at a=', Atest
                end if
            end if

        end do
        !$omp end parallel do

        ! 2) If the model is stable do stability cleanup
        call EFTStability_cleanup()

        return
    end subroutine EFTCheck_Stability

    ! ----------------------------

    subroutine EFTStability_cleanup()
        ! EFTStability_cleanup subroutine.
        ! Restores the values stored in the module to the default. Needed for successive calls.
        implicit none

        PastA1  = 0._dl
        PastAT  = 0._dl

        return
    end subroutine EFTStability_cleanup

    ! ----------------------------

    function EFTStabilityComputation(a)
        !This function computes if the stability requirements are fullfilled in a given time.
        implicit none
        ! 1) Definitions of variables:
        logical :: EFTStabilityComputation, EFT_HaveNan
        real(dl), intent(in) :: a
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

        ! 1) Stability check initialization
        EFTStabilityComputation = .true.

        ! 2) Computation of various quantities:
        EFTtemp_H0     = CP%H0
        EFTtemp_Omegac = CP%omegac
        EFTtemp_Omegab = CP%omegab
        EFTtemp_Omegav = CP%omegav
        EFTtemp_TCMB   = CP%tcmb
        EFTtemp_nu_massless_degeneracy = CP%Num_Nu_massless
        EFTtemp_grhom = 3*EFTtemp_H0**2/c_EFT**2*1000**2
        EFTtemp_grhog = kappa_EFT/c_EFT**2*4*sigma_boltz_EFT/c_EFT**3*EFTtemp_TCMB**4*Mpc_EFT**2 !8*pi*G/c^2*4*sigma_B/c^3 T^4
        EFTtemp_grhor = 7._dl/8*(4._dl/11)**(4._dl/3)*EFTtemp_grhog
        EFT_H0 = (CP%h0/c_EFT)*1000._dl
        ! 2) Computation of various quantities:
        a2=a*a
        grhob_t=EFTtemp_grhom*EFTtemp_Omegab/a
        grhoc_t=EFTtemp_grhom*EFTtemp_Omegac/a
        grhor_t=grhornomass/a2
        grhog_t=EFTtemp_grhog/a2
        grhov_t=EFTtemp_grhom*EFTtemp_Omegav*EFTw(a,3)

        if ( CP%EFTflag==4 ) then ! Full mapping EFT overwrites some of these quantities:
            if ( CP%FullMappingEFTmodel==1 ) then ! Horava gravity
                grhov_t = ( EFTtemp_grhom*EFTtemp_Omegav &
                    &   +3._dl*EFT_H0**2*((CP%Horava_eta +3._dl*CP%Horava_lambda-2._dl*CP%Horava_xi)/(2._dl*CP%Horava_xi+2._dl-CP%Horava_eta)))*EFTw(a,3)
            end if
        end if

        grho_matter  = grhob_t +grhoc_t +grhor_t +grhog_t
        gpres_matter = (grhog_t+grhor_t)/3._dl

        EFT_grhonu_tot = 0._dl
        EFT_gpinu_tot  = 0._dl
        if (CP%Num_Nu_Massive /= 0) then
            do nu_i = 1, CP%Nu_mass_eigenstates
                EFT_grhonu    = 0._dl
                EFT_gpinu     = 0._dl
                grhormass_t=grhormass(nu_i)/a**2
                call Nu_background(a*nu_masses(nu_i),EFT_grhonu,EFT_gpinu)
                EFT_grhonu_tot = EFT_grhonu_tot + grhormass_t*EFT_grhonu
                EFT_gpinu_tot  = EFT_gpinu_tot + grhormass_t*EFT_gpinu
            end do
        end if

        grho_matter  = grho_matter  +EFT_grhonu_tot
        gpres_matter = gpres_matter +EFT_gpinu_tot
        grho         = grho_matter  +grhov_t
        gpres        = gpres_matter +EFTw(a,0)*grhov_t

        ! 3) Hubble and derivatives:
        adotoa=sqrt(grho/3)
        adotdota=(adotoa*adotoa-gpres)/2._dl
        Hdot =adotdota-adotoa**2._dl

        if ( CP%EFTflag==4 ) then ! Full mapping EFT overwrites some of these quantities:
            if ( CP%FullMappingEFTmodel==1 ) then ! Horava gravity
                adotoa   = adotoa*sqrt(( 2._dl +2._dl*CP%Horava_xi -CP%Horava_eta )/(3._dl*CP%Horava_lambda +2._dl))
                adotdota = adotdota*( 2._dl +2._dl*CP%Horava_xi -CP%Horava_eta )/(3._dl*CP%Horava_lambda +2._dl)
                Hdot     = Hdot*( 2._dl +2._dl*CP%Horava_xi -CP%Horava_eta )/(3._dl*CP%Horava_lambda +2._dl)
            end if
        end if

        EFT_gpinu_tot    = 0._dl
        EFT_grhonu_tot   = 0._dl
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

        if ( CP%EFTflag==4 ) then ! Full mapping EFT overwrites some of these quantities:
            if ( CP%FullMappingEFTmodel==1 ) then ! Horava gravity
                Hdotdot = ( 2._dl +2._dl*CP%Horava_xi -CP%Horava_eta )/(3._dl*CP%Horava_lambda +2._dl)*&
                    &( +adotoa*grhov_t*2._dl/3._dl                                                 &
                    &  +adotoa*( -grhob_t -grhoc_t -4._dl*(grhog_t+grhor_t) )/3._dl                &
                    &  +adotoa/6._dl*EFT_grhonu_tot -0.5_dl*adotoa*EFT_gpinu_tot -0.5_dl*EFT_gpinudot_tot )
            end if
        end if

        ! 4) EFT functions:
        EFTOmegaV   = EFTOmega(a,0)
        EFTOmegaP   = EFTOmega(a,1)
        EFTOmegaPP  = EFTOmega(a,2)
        EFTOmegaPPP = EFTOmega(a,3)
        EFTGamma1V  = EFTGamma1(a,0)
        EFTGamma1P  = EFTGamma1(a,1)
        EFTGamma2V  = EFTGamma2(a,0)
        EFTGamma2P  = EFTGamma2(a,1)
        EFTGamma3V  = EFTGamma3(a,0)
        EFTGamma3P  = EFTGamma3(a,1)
        EFTGamma4V  = EFTGamma4(a,0)
        EFTGamma4P  = EFTGamma4(a,1)
        EFTGamma5V  = EFTGamma5(a,0)
        EFTGamma5P  = EFTGamma5(a,1)
        EFTGamma6V  = EFTGamma6(a,0)
        EFTGamma6P  = EFTGamma6(a,1)

        !8*pi*G*c*a^2
        EFTc = (adotoa*adotoa - Hdot)*(EFTOmegaV + a*EFTOmegaP*0.5_dl) &
            & - 0.5_dl*a2*adotoa*adotoa*EFTOmegaPP&
            & + 0.5_dl*grhov_t*(1._dl+EFTw(a,0))
        !8*pi*G*Lambda*a^2
        EFTLambda = +EFTw(a,0)*grhov_t &
            &-EFTOmegaV*(2._dl*Hdot+adotoa**2._dl) &
            &-a*EFTOmegaP*(2._dl*adotoa**2._dl + Hdot) &
            &-a2*adotoa**2._dl*EFTOmegaPP
        !EFT C DOT: 8*pi*G*cdot*a^2
        EFTcdot = -EFTOmegaV*(Hdotdot-4._dl*adotoa*Hdot+2._dl*adotoa*adotoa*adotoa) &
            & + 0.5_dl*a*EFTOmegaP*(-Hdotdot+adotoa*Hdot+adotoa*adotoa*adotoa)&
            & +0.5_dl*a2*adotoa*EFTOmegaPP*(adotoa*adotoa-3._dl*Hdot)&
            & -0.5_dl*a*a2*adotoa*adotoa*adotoa*EFTOmegaPPP&
            & +0.5_dl*adotoa*grhov_t*(-3._dl*(1._dl+EFTw(a,0))**2 + a*EFTw(a,1))
        !EFT LAMBDA DOT: 8*pi*G*Ldot*a^2
        EFTLambdadot = -2._dl*EFTOmegaV*(Hdotdot-adotoa*Hdot-adotoa*adotoa*adotoa)&
            & - a*EFTOmegaP*(4._dl*adotoa*Hdot+Hdotdot)&
            & -a2*EFTOmegaPP*adotoa*(3._dl*Hdot+2._dl*adotoa*adotoa)&
            & -a*a2*EFTOmegaPPP*adotoa*adotoa*adotoa&
            & +grhov_t*adotoa*(a*EFTw(a,1)-3._dl*EFTw(a,0)*(1._dl+EFTw(a,0)))

        if (CP%EFTflag==2) then
            if (CP%DesignerEFTmodel==1) then
                EFTc = EFTcTemp(a,0)
                EFTcdot = EFTcTemp(a,1)
                EFTLambda = EFTLambdaTemp(a,0)
                EFTLambdadot = EFTLambdaTemp(a,1)
            else if (CP%DesignerEFTmodel==2) then
                ! nothing to do for mc5e
            end if
        else if (CP%EFTflag==4) then
            if (CP%FullMappingEFTmodel==1) then
                EFTc         = -( 2._dl*CP%Horava_xi - 3._dl*CP%Horava_lambda )*( Hdot - adotoa**2 )/( 2._dl +2._dl*CP%Horava_xi -CP%Horava_eta )
                EFTcdot      = ( 3._dl*CP%Horava_lambda - 2._dl*CP%Horava_xi )*&
                    &( Hdotdot -4._dl*adotoa*Hdot + 2._dl*adotoa**3 )/( 2._dl +2._dl*CP%Horava_xi -CP%Horava_eta )
                EFTLambda    = +EFTw(a,0)*grhov_t + 2._dl*( 3._dl*CP%Horava_lambda - 2._dl*CP%Horava_xi)*&
                    &( 0.5_dl*adotoa**2 + Hdot)/( 2._dl +2._dl*CP%Horava_xi -CP%Horava_eta )
                EFTLambdadot = + 2._dl*( 3._dl*CP%Horava_lambda - 2._dl*CP%Horava_xi)*&
                    &( Hdotdot -adotoa*Hdot -adotoa**3 )/( 2._dl +2._dl*CP%Horava_xi -CP%Horava_eta )
                EFTGamma1V   = +0.5_dl*( 2._dl*CP%Horava_xi -3._dl*CP%Horava_lambda )*(Hdot-adotoa**2)/(a**2*EFT_H0**2)/( 2._dl +2._dl*CP%Horava_xi -CP%Horava_eta )
                EFTGamma1P   = +0.5_dl*( 2._dl*CP%Horava_xi -3._dl*CP%Horava_lambda )*(+2._dl*adotoa**2 -4._dl*Hdot +Hdotdot/adotoa)/(a**3*EFT_H0**2)/( 2._dl +2._dl*CP%Horava_xi -CP%Horava_eta )
            end if
        end if

        ! 5) Effective EFT quantities:
        EFTgrhoq = 2._dl*EFTc -EFTLambda -3._dl*a*adotoa*adotoa*EFTOmegaP
        EFTgpresq = EFTLambda + a2*adotoa*adotoa*EFTOmegaPP&
            & +a*EFTOmegaP*(Hdot + 2._dl*adotoa*adotoa)
        EFTgrhodotq = -3._dl*adotoa*(EFTgrhoq+EFTgpresq) + 3._dl*a*adotoa**3._dl*EFTOmegaP
        EFTgpresdotq = EFTLambdadot +a2*a*adotoa**3*EFTOmegaPPP + 3._dl*a2*adotoa*Hdot*EFTOmegaPP&
            & +a*EFTOmegaP*Hdotdot +3._dl*a*adotoa*Hdot*EFTOmegaP +2._dl*a2*adotoa**3*EFTOmegaPP&
            & -2._dl*a*adotoa**3*EFTOmegaP

        ! 6) Pi field equation coefficients:
        k = 0._dl
        k2 = 0._dl
        ! Coefficient A of tensors:
        EFTAT = 1._dl + EFTOmegaV - EFTGamma4V
        ! First part of the A coefficient:
        EFTpiA1 = EFTc +2._dl*a2*EFT_H0**2*EFTGamma1V +1.5_dl*a2*(adotoa*EFTOmegaP+EFT_H0*EFTGamma2V)**2&
            &/(2._dl*(1+EFTOmegaV)+EFTGamma3V+EFTGamma4V)
        ! Second part of the A coefficient, the k^2 term:
        EFTpiA2 = +4._dl*EFTGamma6V
        ! First part of the B coefficient:
        EFTpiB1 = EFTcdot +4._dl*adotoa*EFTc +8._dl*a2*adotoa*EFT_H0**2*(EFTGamma1V+ 0.25_dl*a*EFTGamma1P)&
            & -a*(adotoa*EFTOmegaP+EFT_H0*EFTGamma2V)/(4._dl*(1._dl+EFTOmegaV)+6._dl*EFTGamma3V +2._dl*EFTGamma4V)*&
            &(-3._dl*(EFTgrhoQ + EFTgpresQ) -3._dl*a*adotoa**2*EFTOmegaP*(4._dl +Hdot/(adotoa**2)) -3._dl*a2*adotoa**2*EFTOmegaPP&
            & -3._dl*a*adotoa*EFT_H0*(4._dl*EFTGamma2V + a*EFTGamma2P) -(9._dl*EFTGamma3V -3._dl*EFTGamma4V)*&
            &(Hdot-adotoa**2))&
            & +1._dl/(1._dl+EFTOmegaV+2._dl*EFTGamma5V)*(a*adotoa*EFTOmegaP+2._dl*adotoa*(EFTGamma5V + EFTGamma5P)&
            & -(1._dl+EFTOmegaV)*(a*adotoa*EFTOmegaP+a*EFT_H0*EFTGamma2V)/(2._dl*(1._dl+EFTOmegaV) +3._dl*EFTGamma3V +EFTGamma4V))*&
            &(-EFTc +1.5_dl*a*adotoa**2*EFTOmegaP -2._dl*a2*EFT_H0*EFTGamma1V +1.5_dl*a*adotoa*EFT_H0*EFTGamma2V)
        ! Second part of the B coefficient, the k^2 term:
        k2      = 1._dl
        EFTpiB2 = +4._dl*k2*adotoa*(2._dl*EFTGamma6V+ a*EFTGamma6P) +a*k2*(EFTGamma4V&
            & +2._dl*EFTGamma5V)/(2._dl*(1._dl+EFTOmegaV)-2._dl*EFTGamma4V)*(adotoa*EFTOmegaP+EFT_H0*EFTGamma2V)&
            & -a*(adotoa*EFTOmegaP+EFT_H0*EFTGamma2V)/(4._dl*(1._dl+EFTOmegaV)+6._dl*EFTGamma3V +2._dl*EFTGamma4V)*&
            &( +k2*(3._dl*EFTGamma3V -EFTGamma4V +4._dl*EFTGamma5V))&
            & +1._dl/(1._dl+EFTOmegaV+2._dl*EFTGamma5V)*(a*adotoa*EFTOmegaP+2._dl*adotoa*(EFTGamma5V + EFTGamma5P)&
            & -(1._dl+EFTOmegaV)*(a*adotoa*EFTOmegaP+a*EFT_H0*EFTGamma2V)/(2._dl*(1._dl+EFTOmegaV) +3._dl*EFTGamma3V +EFTGamma4V))*&
            &(-4._dl*EFTGamma6V*k2)
        ! The C coefficient:
        EFTpiC = +adotoa*EFTcdot +(6._dl*adotoa**2-2._dl*Hdot)*EFTc +1.5_dl*a*adotoa*EFTOmegaP*(Hdotdot-2._dl*adotoa**3) &
            & +6._dl*a2*adotoa**2*EFT_H0**2*EFTGamma1V +2._dl*a2*Hdot*EFT_H0**2*EFTGamma1V &
            & +2._dl*a2*a*adotoa**2*EFT_H0**2*EFTGamma1P +1.5_dl*(Hdot-adotoa**2)**2*(EFTGamma4V +3._dl*EFTGamma3V )&
            & +4.5_dl*adotoa*EFT_H0*a*(Hdot-adotoa**2)*(EFTGamma2V + a*EFTGamma2P/3._dl )&
            & +0.5_dl*a*EFT_H0*EFTGamma2V*(3._dl*Hdotdot -12._dl*Hdot*adotoa +6._dl*adotoa**3) &
            & -a*(adotoa*EFTOmegaP+EFT_H0*EFTGamma2V)/(4._dl*(1._dl+EFTOmegaV)+6._dl*EFTGamma3V +2._dl*EFTGamma4V)*&
            &(-3._dl*EFTgpresdotQ -3._dl*adotoa*(EFTgrhoQ +EFTgpresQ)-3._dl*a*adotoa**3*(a*EFTOmegaPP +6._dl*EFTOmegaP) &
            & -6._dl*a*adotoa*Hdot*EFTOmegaP +3._dl*(Hdotdot-2._dl*adotoa*Hdot)*(EFTGamma4V +3._dl*EFTGamma3V)&
            & +6._dl*adotoa*(Hdot-adotoa**2)*(3._dl*EFTGamma3V +1.5_dl*a*EFTGamma3P +EFTGamma4V + 0.5_dl*a*EFTGamma4P)&
            & -3._dl*a*EFT_H0*(3._dl*adotoa**2*EFTGamma2V +Hdot*EFTGamma2V +a*adotoa**2*EFTGamma2P))&
            & +1._dl/(1._dl+EFTOmegaV+2._dl*EFTGamma5V)*(a*adotoa*EFTOmegaP+2._dl*adotoa*(EFTGamma5V +a*EFTGamma5P)&
            & -(1._dl+EFTOmegaV)*(a*adotoa*EFTOmegaP+a*EFT_H0*EFTGamma2V)/(2._dl*(1._dl+EFTOmegaV)+3._dl*EFTGamma3V +EFTGamma4V))*&
            &(-0.5*EFTgrhodotQ -adotoa*EFTc +1.5_dl*a*adotoa*EFTOmegaP*(3._dl*adotoa**2-Hdot) -2._dl*a2*adotoa*EFT_H0**2*EFTGamma1V&
            & -1.5_dl*a*EFT_H0*EFTGamma2V*(Hdot-2._dl*adotoa**2) -3._dl*adotoa*(Hdot-adotoa**2)*(1.5_dl*EFTGamma3V +0.5_dl*EFTGamma4V))
        ! First part of the D coefficient, the k^2 term:
        EFTpiD1 = EFTc -0.5_dl*a*adotoa*EFT_H0*(EFTGamma2V +a*EFTGamma2P) +(adotoa**2-Hdot)*(3._dl*EFTGamma3V +EFTGamma4V)&
            & +4._dl*(Hdot*EFTGamma6V + adotoa**2*EFTGamma6V + a*adotoa**2*EFTGamma6P)&
            & +2._dl*(Hdot*EFTGamma5V +a*adotoa**2*EFTGamma5P)&
            & -a*(adotoa*EFTOmegaP+EFT_H0*EFTGamma2V)/(4._dl*(1._dl+EFTOmegaV)+6._dl*EFTGamma3V +2._dl*EFTGamma4V)*&
            &(-2._dl*a*adotoa*EFTOmegaP +4._dl*adotoa*EFTGamma5V -2._dl*adotoa*(3._dl*EFTGamma3V +1.5_dl*a*EFTGamma3P &
            & +EFTGamma4V +0.5_dl*a*EFTGamma4P))&
            & +1._dl/(1._dl+EFTOmegaV+2._dl*EFTGamma5V)*(a*adotoa*EFTOmegaP+2._dl*adotoa*(EFTGamma5V +a*EFTGamma5P)&
            & -(1._dl+EFTOmegaV)*(a*adotoa*EFTOmegaP +a*EFT_H0*EFTGamma2V)/(2._dl*(1._dl+EFTOmegaV)+3._dl*EFTGamma3V +EFTGamma4V))*&
            &(+0.5_dl*a*adotoa*EFTOmegaP -2._dl*adotoa*EFTGamma5V +0.5_dl*a*EFT_H0*EFTGamma2V +1.5_dl*adotoa*EFTGamma3V&
            & +0.5_dl*adotoa*EFTGamma4V -4._dl*adotoa*EFTGamma6V)&
            & +(EFTGamma4V +2._dl*EFTGamma5V)/(2._dl*(1._dl+EFTOmegaV) -2._dl*EFTGamma4V)*(EFTgrhoQ +EFTgpresQ +a*adotoa**2*EFTOmegaP&
            & -EFTGamma4V*(Hdot-adotoa**2) +a*adotoa*EFT_H0*EFTGamma2V +3._dl*EFTGamma3V*(adotoa**2-Hdot))
        ! Second part of the D coefficient, the k^4 term:
        EFTpiD2 = +(+0.5_dl*EFTGamma3V +0.5_dl*EFTGamma4V &
            & +(EFTGamma4V +2._dl*EFTGamma5V)/(2._dl*(1._dl+EFTOmegaV) -2._dl*EFTGamma4V)*(EFTGamma3V +EFTGamma4V))

        if ( CP%EFTflag==4 ) then ! Full mapping EFT overwrites some of these quantities:
            if ( CP%FullMappingEFTmodel==1 ) then ! Horava gravity
                EFTpiA1 = 0._dl
                EFTpiA2 = 0.5_dl*CP%Horava_eta
                EFTpiB1 = 0._dl
                EFTpiB2 = adotoa*CP%Horava_eta
                EFTpiC  = 0._dl
                EFTpiD1 = 0.5_dl*(3._dl*CP%Horava_lambda-2._dl*CP%Horava_xi)*(adotoa**2-Hdot) + 0.5_dl*CP%Horava_eta*(adotoa**2+Hdot)
                EFTpiD2 = 0.5_dl*CP%Horava_lambda*(1._dl+CP%Horava_xi)
            end if
        end if

        EFT_kinetic  = 9._dl*( 1._dl +EFTOmegaV -EFTGamma4V )*( 4._dl*EFTc*( 1._dl +EFTOmegaV -EFTGamma4V ) &
            & +3._dl*adotoa**2*EFTOmegaP**2*a2 + &
            & a2*EFT_H0*( EFT_H0*( 3._dl*EFTGamma2V**2 +8._dl*EFTGamma1V*( 1._dl +EFTOmegaV -EFTGamma4V ) +6._dl*adotoa*EFTGamma2V*EFTOmegaP ) ) )
        ! Gradient term:
        EFT_gradient = 9._dl*(8._dl*a*adotoa**2*EFTGamma5P - 16._dl*adotoa**2*EFTGamma5V**2 + 16._dl*EFTc*EFTGamma5V**2 - 2._dl*a**2*adotoa**2*EFTGamma4P*EFTOmegaP + 4._dl*a**2*adotoa**2*EFTGamma5P*EFTOmegaP -&
            &4._dl*a*adotoa**2*EFTGamma5V*EFTOmegaP - 4._dl*a**2*adotoa**2*EFTGamma4P*EFTGamma5V*EFTOmegaP - 8._dl*a*adotoa**2*EFTGamma5V**2*EFTOmegaP + 3._dl*a**2*adotoa**2*EFTOmegaP**2 +&
            &4._dl*a**2*adotoa**2*EFTGamma5V*EFTOmegaP**2 + 16._dl*a*adotoa**2*EFTGamma5P*EFTOmegaV + 16._dl*EFTc*EFTGamma5V*EFTOmegaV - 16._dl*adotoa**2*EFTGamma5V**2*EFTOmegaV -&
            &2._dl*a**2*adotoa**2*EFTGamma4P*EFTOmegaP*EFTOmegaV + 4._dl*a**2*adotoa**2*EFTGamma5P*EFTOmegaP*EFTOmegaV - 4._dl*a*adotoa**2*EFTGamma5V*EFTOmegaP*EFTOmegaV +&
            &3._dl*a**2*adotoa**2*EFTOmegaP**2*EFTOmegaV + 8._dl*a*adotoa**2*EFTGamma5P*EFTOmegaV**2 - a**2*EFTGamma2V**2*EFT_H0**2*(1 + EFTOmegaV) +&
            &4._dl*a**2*adotoa**2*EFTGamma5V*EFTOmegaPP*(1._dl + 2._dl*EFTGamma5V + EFTOmegaV) + 4._dl*EFTc*(4._dl*EFTGamma5V + (1._dl + EFTOmegaV)**2) -&
            &2._dl*a*adotoa*EFT_H0*(a*EFTGamma2P*(1._dl - EFTGamma4V + EFTOmegaV)*(1._dl + 2._dl*EFTGamma5V + EFTOmegaV) +&
            &EFTGamma2V*(EFTGamma4V*(-1._dl + 2._dl*a*EFTGamma5P + 2._dl*EFTGamma5V + a*EFTOmegaP - EFTOmegaV) +&
            &(1._dl + EFTOmegaV)*(1._dl + a*(EFTGamma4P - 2._dl*EFTGamma5P - EFTOmegaP) + EFTOmegaV) - 2._dl*EFTGamma5V*(1._dl - a*EFTGamma4P + a*EFTOmegaP + EFTOmegaV))) +&
            &8._dl*EFTGamma5V*Hdot + 16._dl*EFTGamma5V**2*Hdot + 4._dl*a*EFTGamma5V*EFTOmegaP*Hdot + 8._dl*a*EFTGamma5V**2*EFTOmegaP*Hdot + 16._dl*EFTGamma5V*EFTOmegaV*Hdot +&
            &16._dl*EFTGamma5V**2*EFTOmegaV*Hdot + 4._dl*a*EFTGamma5V*EFTOmegaP*EFTOmegaV*Hdot + 8._dl*EFTGamma5V*EFTOmegaV**2*Hdot +&
            &4._dl*EFTGamma4V**2*(adotoa**2*(1._dl + 2._dl*a*EFTGamma5P + 4._dl*EFTGamma5V + a*EFTOmegaP + EFTOmegaV) - (1._dl + 2._dl*EFTGamma5V + EFTOmegaV)*Hdot) +&
            &2._dl*EFTGamma4V*(adotoa**2*(-(a**2*EFTOmegaP**2) + a**2*EFTOmegaPP*(1._dl + 2._dl*EFTGamma5V + EFTOmegaV) -&
            &4._dl*(1._dl + EFTOmegaV)*(1._dl + 2._dl*a*EFTGamma5P + 4._dl*EFTGamma5V + EFTOmegaV) - a*EFTOmegaP*(3._dl + 2._dl*a*EFTGamma5P + 2._dl*EFTGamma5V + 3._dl*EFTOmegaV)) +&
            &(1._dl + 2._dl*EFTGamma5V + EFTOmegaV)*(4._dl + a*EFTOmegaP + 4._dl*EFTOmegaV)*Hdot))

        ! 7) Compute k_max_CAMB. This is slightly overestimated to be safe...

        !    This is just to be safe if doing only background stuff. We still require stability of linear scales.
        kmax = 0.1_dl
        !    This is if we want scalar Cls.
        if (CP%WantCls) kmax = CP%Max_eta_k/CP%tau0
        !    This is if we want also (or only) transfer functions.
        if (CP%WantTransfer) then
            kmax = CP%Max_eta_k/CP%tau0*exp((int(log(CP%Transfer%kmax/(CP%Max_eta_k/CP%tau0))*(3*AccuracyBoost))+2)/(3*AccuracyBoost))
        end if

        ! All the coefficients should not be Nan. This can happen for strange values of the parameters
        ! for which a division by zero may occur.

        EFT_HaveNan = IsNaN(EFTAT).or.IsNaN(EFTpiA1).or.IsNaN(EFTpiA2).or.IsNaN(EFTpiB1).or.IsNaN(EFTpiB2)&
            & .or.IsNaN(EFTpiC).or.IsNaN(EFTpiD1).or.IsNaN(EFTpiD2).or.IsNaN(EFT_kinetic).or.IsNaN(EFT_gradient)

        if (EFT_HaveNan) then
            EFTStabilityComputation = .false.
            if (Feedbacklevel > 0) write(*,*) 'EFTCAMB: model has Nan.'
            return
        end if

        ! Mathematical stability:
        if (CP%EFT_mathematical_stability) then

            ! 1- the A coefficient should not change sign in time and in k, i.e. it shall not be zero.
            !    This is the stronger stability constraint since violating it would violate the mathematical
            !    consistency of the pi field equation.

            !    The first condition is A1/=0. Implemented by detecting sign changes in A1.
            if ( EFTpiA1*PastA1 < 0._dl ) then
                EFTStabilityComputation = .false.
                if (Feedbacklevel > 0) write(*,*) 'EFTCAMB: mathematical instability, A is zero in time.'
            end if
            PastA1 = EFTpiA1
            !    The second one is the condition on k.
            if ( (EFTpiA1 > 0 .and. EFTpiA1 + kmax**2*EFTpiA2 < 0) .or. &
                &(EFTpiA1 < 0 .and. EFTpiA1 + kmax**2*EFTpiA2 > 0) ) then
                EFTStabilityComputation = .false.
                if (Feedbacklevel > 0) write(*,*) 'EFTCAMB: mathematical instability, A is zero in k.'
            end if

            ! 2- the AT coefficient should not change sign in time, i.e. it shall not be zero.
            !    This is the second stronger stability constraint since violating it would
            !    violate the mathematical consistency of the tensor perturbation equation.
            !    Implemented by detecting sign changes in AT.
            if ( EFTAT*PastAT < 0._dl ) then
                EFTStabilityComputation = .false.
                if (Feedbacklevel > 0) write(*,*) 'EFTCAMB: mathematical instability, AT is zero in time.'
            end if
            PastAT = EFTAT

            ! 3- we do not want (fast) growing exponential modes.
            !    This condition prevents the pi field from growing exponentially and destroying everything.
            !    Even though this condition is not completely related to physics not mathematics, violating it
            !    would completely screw up cosmological observables.

            !    This is the maximum allowed rate of instability. Units shall be Mpc^-1.
            EFT_instability_rate = 0._dl

            !    This condition needs to be tested in k. Sample in k.

            ind_max = 10

            do ind = 1, ind_max
                ! kmode to test. Linear sampling. Should suffice... (??)
                tempk = 0._dl + REAL(ind-1)*(kmax)/REAL(ind_max-1)
                ! vaule that discriminates between different cases:
                temp1 = (EFTpiB1 +EFTpiB2*tempk**2)
                temp2 = (EFTpiA1 +EFTpiA2*tempk**2)
                temp3 = temp1**2 -4._dl*temp2*(EFTpiC +EFTpiD1*tempk**2 + EFTpiD2*tempk**4)

                !    case 1:
                if (temp3 > 0._dl .and. temp2 /= 0._dl) then
                    temp4 = +0.5_dl*(-temp1 +sqrt(temp3))/temp2
                    temp5 = +0.5_dl*(-temp1 -sqrt(temp3))/temp2
                    if (temp4>EFT_instability_rate .or. temp5>EFT_instability_rate) then
                        EFTStabilityComputation = .false.
                        if (Feedbacklevel > 0) write(*,*) 'EFTCAMB: mathematical instability. Growing exponential at k', tempk, temp4, temp5
                        exit
                    end if
                !    case 2:
                else if ( temp2 /= 0._dl ) then
                    temp4 = -0.5_dl*temp1/temp2
                    if (temp4>EFT_instability_rate) then
                        EFTStabilityComputation = .false.
                        if (Feedbacklevel > 0) write(*,*) 'EFTCAMB: mathematical instability. Growing exponential at k', tempk, temp4
                        exit
                    end if
                end if

            end do

        end if

        ! Additional priors:
        if (CP%EFTAdditionalPriors) then
            if (EFTw(a,0)>-1._dl/3._dl) EFTStabilityComputation = .false.
        end if

        ! Minkowsky prior: some theories have known stability properties on Minkowsky background:
        if (CP%MinkowskyPriors) then
            if (CP%EFTflag==4) then
                if (CP%FullMappingEFTmodel==1) then ! Horava gravity

                    if ( CP%Horava_lambda > -2._dl/3._dl .and. CP%Horava_lambda < 0._dl ) then
                        EFTStabilityComputation = .false.
                        if (Feedbacklevel > 0) write(*,*) 'EFTCAMB: Instability on Minkowsky backgound'
                    end if

                    if ( CP%Horava_eta < 0._dl .or. CP%Horava_eta > 2._dl*CP%Horava_xi +2._dl ) then
                        EFTStabilityComputation = .false.
                        if (Feedbacklevel > 0) write(*,*) 'EFTCAMB: Instability on Minkowsky backgound'
                    end if

                end if
            end if
        end if



        ! Physical viability:
        if (CP%EFT_physical_stability) then

            if ( .not. EFT_old_stability .and. &
                & CP%EFTFlag /= 4 .and. ( &
                & (EFTGamma6V /= 0._dl) .or.      &
                & ((EFTGamma3V + EFTGamma4V) /= 0._dl) ) ) then

                write(*,*) 'EFTCAMB WARNING: stability for model beyond GLPV has not been worked out.'
                write(*,*) 'It will be added in a future release.'
                write(*,*) 'If you want to run this model disable EFT_physical_stability.'

                EFTStabilityComputation = .false.
                return
            end if

            ! 1- Positive gravitational constant:
            if (1._dl +EFTOmegaV <= 0) then
                EFTStabilityComputation = .false.
                if (Feedbacklevel > 0) write(*,*) 'EFTCAMB: negative gravitational constant', 1._dl +EFTOmegaV
            end if

            ! 2- Old ghost and gradient conditions:
            if ( EFT_old_stability .or. &
                & (EFTGamma6V /= 0._dl) .or. &
                & ((EFTGamma3V + EFTGamma4V) /= 0._dl) ) then

                ! Ghost instability:
                if (EFTpiA1 < 0 .or. ( EFTpiA1 + kmax**2*EFTpiA2 < 0)) then
                    EFTStabilityComputation = .false.
                    if (Feedbacklevel > 0) write(*,*) 'EFTCAMB: ghost instability', EFTpiA1, EFTpiA1 + kmax**2*EFTpiA2
                end if
                ! Gradient instability 1:
                if (EFTpiD1 < 0) then
                    EFTStabilityComputation = .false.
                    if (Feedbacklevel > 0) write(*,*) 'EFTCAMB: gradient instability k^2', EFTpiD1
                end if
                ! Gradient instability 2:
                if (EFTpiD2 < 0) then
                    EFTStabilityComputation = .false.
                    if (Feedbacklevel > 0) write(*,*) 'EFTCAMB: gradient instability k^4', EFTpiD2
                end if

            else
                ! New ghost and gradient conditions:
                ! ghost condition:
                if ( EFT_kinetic < 0._dl ) then
                    EFTStabilityComputation = .false.
                    if (Feedbacklevel > 0) write(*,*) 'EFTCAMB new stability: ghost instability. Kinetic: ', EFT_kinetic
                end if
                ! gradient instability:
                if ( EFT_gradient < 0._dl ) then
                    EFTStabilityComputation = .false.
                    if (Feedbacklevel > 0) write(*,*) 'EFTCAMB new stability: gradient instability. Gradient: ', EFT_gradient
                end if

            end if

            !5- Positive effective mass of pi:
            if (EFTpiC < 0.and.EFTpiMassPrior) then
                EFTStabilityComputation = .false.
                if (Feedbacklevel > 0) write(*,*) 'EFTCAMB: negative mass'
            end if
            ! 6- No tensor ghosts:
            if (EFTAT < 0) then
                EFTStabilityComputation = .false.
                if (Feedbacklevel > 0) write(*,*) 'EFTCAMB: tensor ghost instability'
            end if
            ! 7- Sub-luminal propagation:
            if (EFTlightspeedPrior) then
                if (EFTpiA2==0.and.EFTpiD2==0.and.(EFTpiD1/EFTpiA1)>1.001_dl) then
                    EFTStabilityComputation = .false.
                    if (Feedbacklevel > 0) write(*,*) 'EFTCAMB: tachion perturbations'
                else if (EFTpiA2/=0.and.EFTpiD2/=0.and.(EFTpiD2/EFTpiA2)>1.001_dl) then
                    EFTStabilityComputation = .false.
                    if (Feedbacklevel > 0) write(*,*) 'EFTCAMB: tachion perturbations'
                end if
            end if
            ! 8- Every theory has it's own peculiarities...
            ! 1) F(R): for this model it is easy to show that the positive mass condition requires that OmegaPrime
            !    should be positive. We add this test for this models as it is numerically easier to check.
            if (CP%EFTflag==2.and.CP%DesignerEFTmodel==1) then
                if (EFTOmega(0.11_dl*EFTturnonpiInitial,1)*EFTOmegaP<0) EFTStabilityComputation = .false.
            end if
        end if

        return
    end function EFTStabilityComputation

end module EFTstability

! -------------------------------------------------------------------------------------------------

!Performs EFT initialization just after CAMB is called.
module EFTinitialization
    use EFTDef
    use EFTdeEOS
    use EFTfunctions
    use EFTreturntoGR
    use EFTstability

    implicit none

contains

    ! ----------------------------

    subroutine EFTCAMB_initialization(success)
        implicit none
        logical  :: success
        real(dl) :: astart, aend

        integer :: ind

        success = .true.

        ! 1) Feedback and check of the flags:
        if (Feedbacklevel>0) write(*,*) 'EFTCAMB initialization'
        call EFTCAMB_check_flags_consistency
        if (Feedbacklevel>0) call EFTCAMB_print_parameters

        ! 2) EFT hard priors
        call EFTCAMB_HardPriors(success)
        if (.not.success) return

        ! 3) Call the designer
        if (CP%EFTflag==2) then
            call EFTCAMB_Designer(success)
            if (.not.success) return
        end if

        ! 4) Check stability
        if (EarlyTimeStability) then
            astart = 10._dl**(-7.8_dl)
        else
            astart = EFTturnonpiInitial
        end if

        aend   = 1._dl
        call EFTCheck_Stability(success, astart, aend)
        if (.not.success) return

        ! 5) Calculate Return to GR
        call EFTCheckReturnToGR

    end subroutine EFTCAMB_initialization

    ! ----------------------------

    subroutine EFTCAMB_HardPriors(success)
        ! This subroutine enforces hard priors on cosmological parameters to avoid EFTCAMB
        ! being called for weird values of cosmological parameters.
        implicit none
        logical  :: success
        real(dl) :: a, amin, amax
        integer  :: i

        ! 1) There should be dark energy:
        if (CP%Omegav<0._dl.and.CP%EFTAdditionalPriors) then
            success=.false.
            if (Feedbacklevel > 0) write(*,*) 'EFTCAMB: Hard priors rejecting the model.'
            return
        end if

        ! 2) Dark energy should source cosmic acceleration:
        if (CP%EFTAdditionalPriors) then
            do i=1, 100
                amin = 0._dl
                amax = 1._dl
                a = amin +REAL(i-1)*(amax-amin)/REAL(100-1)
                if (EFTw(a,0)>-1._dl/3._dl) then
                    success=.false.
                    if (Feedbacklevel > 0) write(*,*) 'EFTCAMB: Hard priors rejecting the model.'
                    return
                end if
            end do
        end if

        ! 3) Protection against weird input from CosmoMC:
        if ( ( CP%Omegav<0._dl .or. CP%Omegav>1.0_dl).and.CP%EFTAdditionalPriors) then
            success=.false.
            if (Feedbacklevel > 0) write(*,*) 'EFTCAMB: Hard priors rejecting the model.'
            return
        end if
        if ( ( CP%omegac<0._dl .or. CP%omegac>1.0_dl).and.CP%EFTAdditionalPriors) then
            success=.false.
            if (Feedbacklevel > 0) write(*,*) 'EFTCAMB: Hard priors rejecting the model.'
            return
        end if
        if ( ( CP%omegab<0._dl .or. CP%omegab>1.0_dl).and.CP%EFTAdditionalPriors) then
            success=.false.
            if (Feedbacklevel > 0) write(*,*) 'EFTCAMB: Hard priors rejecting the model.'
            return
        end if
        if ( ( CP%H0<1._dl ).and.CP%EFTAdditionalPriors ) then
            success=.false.
            if (Feedbacklevel > 0) write(*,*) 'EFTCAMB: Hard priors rejecting the model.'
            return
        end if

        ! 4) One Horava condition:
        if ( CP%EFTflag==4 ) then ! Full mapping EFT overwrites some of these quantities:
            if ( CP%FullMappingEFTmodel==1 ) then ! Horava gravity
                if ( ( CP%Omegav+(CP%Horava_eta +3._dl*CP%Horava_lambda-2._dl*CP%Horava_xi)/(2._dl*CP%Horava_xi+2._dl-CP%Horava_eta) < 0._dl ).and.CP%EFTAdditionalPriors ) then
                    success=.false.
                    if (Feedbacklevel > 0) write(*,*) 'EFTCAMB: Hard priors rejecting the model.'
                    return
                end if
            !                if ( ( CP%Omegav+(CP%Horava_eta +3._dl*CP%Horava_lambda-2._dl*CP%Horava_xi)/(2._dl*CP%Horava_xi+2._dl-CP%Horava_eta) > 1._dl ).and.CP%EFTAdditionalPriors ) then
            !                    success=.false.
            !                    if (Feedbacklevel > 0) write(*,*) 'EFTCAMB: Hard priors rejecting the model.'
            !                    return
            !                end if
            end if
        end if

    end subroutine EFTCAMB_HardPriors

    ! ----------------------------

    subroutine EFTCAMB_print_parameters
        ! This is used to print EFT parameters on the screen.

        use compile_time_eft

        implicit none

        logical, parameter :: print_all = .false.

        if ( compile_time_eftcamb ) then
            write(*,*)
            write(*,*) 'WARNING: EFTCAMB running with compile time model selection.'
        end if

        write(*,*)
        write(*,*) 'EFTCAMB stability flags:'

        write(*,*) ' Mathematical stability = ', CP%EFT_mathematical_stability
        write(*,*) ' Physical stability     = ', CP%EFT_physical_stability
        write(*,*) ' Additional priors      = ', CP%EFTAdditionalPriors
        write(*,*) ' Minkowsky priors       = ', CP%MinkowskyPriors
        write(*,*)

        write(*,*) 'EFTCAMB model flags:'

        write(*,"(A24,I3)") '   EFTflag             =', CP%EFTflag

        if (CP%EFTflag==1.or.CP%EFTflag==2.or.print_all) &
            write(*,"(A24,I3)") '   w_DE                =', CP%EFTwDE

        if (CP%EFTflag==2.or.print_all) &
            write(*,"(A24,I3)") '   DesignerEFTmodel    =', CP%DesignerEFTmodel

        if (CP%EFTflag==3.or.print_all) &
            write(*,"(A24,I3)") '   AltParEFTmodel      =', CP%DesignerEFTmodel

        if (CP%EFTflag==4.or.print_all) &
            write(*,"(A24,I3)") '   FullMappingEFTmodel =', CP%DesignerEFTmodel

        if ( (CP%EFTflag==1 .and. CP%PureEFTHorndeski) .or. print_all ) &
            write(*,"(a)") '   Pure EFT Horndeski'

        if (CP%EFTflag==1.or.print_all) then
            if (CP%PureEFTmodelOmega/=0.or.print_all) &
                write(*,"(A24,I3)") '   PureEFTmodelOmega   =', CP%PureEFTmodelOmega
            if (CP%PureEFTmodelGamma1/=0.or.print_all) &
                write(*,"(A24,I3)") '   PureEFTmodelGamma1  =', CP%PureEFTmodelGamma1
            if (CP%PureEFTmodelGamma2/=0.or.print_all) &
                write(*,"(A24,I3)") '   PureEFTmodelGamma2  =', CP%PureEFTmodelGamma2
            if (CP%PureEFTmodelGamma3/=0.or.print_all) &
                write(*,"(A24,I3)") '   PureEFTmodelGamma3  =', CP%PureEFTmodelGamma3
            if ((CP%PureEFTmodelGamma4/=0.and..not.CP%PureEFTHorndeski).or.print_all) &
                write(*,"(A24,I3)") '   PureEFTmodelGamma4  =', CP%PureEFTmodelGamma4
            if ((CP%PureEFTmodelGamma5/=0.and..not.CP%PureEFTHorndeski).or.print_all) &
                write(*,"(A24,I3)") '   PureEFTmodelGamma5  =', CP%PureEFTmodelGamma5
            if ((CP%PureEFTmodelGamma6/=0.and..not.CP%PureEFTHorndeski).or.print_all) &
                write(*,"(A24,I3)") '   PureEFTmodelGamma6  =', CP%PureEFTmodelGamma6
        end if

        if (CP%EFTflag==3.and.CP%AltParEFTmodel==1.or.print_all) then
            if (CP%RPHmassPmodel/=0.or.print_all) &
                write(*,"(A24,I3)") '   RPHmassPmodel       =', CP%RPHmassPmodel
            if (CP%RPHkineticitymodel/=0.or.print_all) &
                write(*,"(A24,I3)") '   RPHkineticitymodel  =', CP%RPHkineticitymodel
            if (CP%RPHbraidingmodel/=0.or.print_all) &
                write(*,"(A24,I3)") '   RPHbraidingmodel    =', CP%RPHbraidingmodel
            if (CP%RPHtensormodel/=0.or.print_all) &
                write(*,"(A24,I3)") '   RPHtensormodel      =', CP%RPHtensormodel
        end if

        write(*,*) 'EFTCAMB parameters:'

        if ((CP%EFTflag==1.or.CP%EFTflag==2.or.CP%EFTflag==3).and.CP%EFTwDE/=0.or.print_all) then
            write(*,"(A18,F9.6)") '   EFTw0        = ', CP%EFTw0
            if (CP%EFTwDE>=2.or.print_all) &
                write(*,"(A18,F9.6)") '   EFTwa        = ', CP%EFTwa
            if (CP%EFTwDE==3.or.print_all) &
                write(*,"(A18,F9.6)") '   EFTwn        = ', CP%EFTwn
            if (CP%EFTwDE==4.or.print_all) &
                write(*,"(A18,F9.6)") '   EFTwat       = ', CP%EFTwat
            if (CP%EFTwDE>=5.or.print_all) &
                write(*,"(A18,F9.6)") '   EFtw2        = ', CP%EFtw2
            if (CP%EFTwDE>=5.or.print_all) &
                write(*,"(A18,F9.6)") '   EFTw3        = ', CP%EFTw3
        end if

        if (CP%EFTflag==1.or.print_all) then
            if (CP%PureEFTmodelOmega/=0.or.print_all)  &
                write(*,"(A18,F12.6)") '   EFTOmega0    = ', CP%EFTOmega0
            if (CP%PureEFTmodelOmega >2.or.print_all)  &
                write(*,"(A18,F12.6)") '   EFTOmegaExp  = ', CP%EFTOmegaExp
            if (CP%PureEFTmodelGamma1/=0.or.print_all) &
                write(*,"(A18,F12.6)") '   EFTGamma10   = ', CP%EFTGamma10
            if (CP%PureEFTmodelGamma1 >2.or.print_all) &
                write(*,"(A18,F12.6)") '   EFTGamma1Exp = ', CP%EFTGamma1Exp
            if (CP%PureEFTmodelGamma2/=0.or.print_all) &
                write(*,"(A18,F12.6)") '   EFTGamma20   = ', CP%EFTGamma20
            if (CP%PureEFTmodelGamma2 >2.or.print_all) &
                write(*,"(A18,F12.6)") '   EFTGamma2Exp = ', CP%EFTGamma2Exp
            if (CP%PureEFTmodelGamma3/=0.or.print_all) &
                write(*,"(A18,F12.6)") '   EFTGamma30   = ', CP%EFTGamma30
            if (CP%PureEFTmodelGamma3 >2.or.print_all) &
                write(*,"(A18,F12.6)") '   EFTGamma3Exp = ', CP%EFTGamma3Exp
            if ((CP%PureEFTmodelGamma4/=0.and..not.CP%PureEFTHorndeski).or.print_all) &
                write(*,"(A18,F12.6)") '   EFTGamma40   = ', CP%EFTGamma40
            if ((CP%PureEFTmodelGamma4 >2.and..not.CP%PureEFTHorndeski).or.print_all) &
                write(*,"(A18,F12.6)") '   EFTGamma4Exp = ', CP%EFTGamma4Exp
            if ((CP%PureEFTmodelGamma5/=0.and..not.CP%PureEFTHorndeski).or.print_all) &
                write(*,"(A18,F12.6)") '   EFTGamma50   = ', CP%EFTGamma50
            if ((CP%PureEFTmodelGamma5 >2.and..not.CP%PureEFTHorndeski).or.print_all) &
                write(*,"(A18,F12.6)") '   EFTGamma5Exp = ', CP%EFTGamma5Exp
            if ((CP%PureEFTmodelGamma6/=0.and..not.CP%PureEFTHorndeski).or.print_all) &
                write(*,"(A18,F12.6)") '   EFTGamma60   = ', CP%EFTGamma60
            if ((CP%PureEFTmodelGamma6 >2.and..not.CP%PureEFTHorndeski).or.print_all) &
                write(*,"(A18,F12.6)") '   EFTGamma6Exp = ', CP%EFTGamma6Exp
        end if

        if (CP%EFTflag==2.or.print_all) then
            if (CP%DesignerEFTmodel==1.or.print_all) write(*,"(A18,F12.6)")  '   F(R) B0      = ', CP%EFTB0
        end if

        if (CP%EFTflag==3.or.print_all) then
            if (CP%AltParEFTmodel==1.or.print_all) then
                if (CP%RPHmassPmodel/=0.or.print_all) write(*,"(A22,F12.6)") &
                    & '   RPHmassP0        = ', CP%RPHmassP0
                if (CP%RPHmassPmodel >1.or.print_all) write(*,"(A22,F12.6)") &
                    & '   RPHmassPexp      = ', CP%RPHmassPexp
                if (CP%RPHkineticitymodel/=0.or.print_all) write(*,"(A22,F12.6)") &
                    & '   RPHkineticity0   = ', CP%RPHkineticity0
                if (CP%RPHkineticitymodel >1.or.print_all) write(*,"(A22,F12.6)") &
                    & '   RPHkineticityexp = ', CP%RPHkineticityexp
                if (CP%RPHbraidingmodel/=0.or.print_all) write(*,"(A22,F12.6)")  &
                    & '   RPHbraiding0     = ', CP%RPHbraiding0
                if (CP%RPHbraidingmodel >1.or.print_all) write(*,"(A22,F12.6)")  &
                    & '   RPHbraidingexp   = ', CP%RPHbraidingexp
                if (CP%RPHtensormodel/=0.or.print_all) write(*,"(A22,F12.6)") &
                    & '   RPHtensor0       = ', CP%RPHtensor0
                if (CP%RPHtensormodel >1.or.print_all) write(*,"(A22,F12.6)") &
                    & '   RPHtensorexp     = ', CP%RPHtensorexp
            end if
        end if

        if (CP%EFTflag==4.or.print_all) then
            if (CP%FullMappingEFTmodel==1.or.print_all) then
                if ( CP%HoravaSolarSystem .or. print_all )  write(*,"(a)") '   Solar System Horava'
                write(*,"(A22,F12.6)") '   Horava_xi        = ', CP%Horava_xi
                write(*,"(A22,F12.6)") '   Horava_lambda    = ', CP%Horava_lambda
                write(*,"(A22,F12.6)") '   Horava_eta       = ', CP%Horava_eta
            end if
        end if

    end subroutine EFTCAMB_print_parameters

    ! ----------------------------

    subroutine EFTCAMB_check_flags_consistency
        ! This subroutine enforces the consistency of the EFTCAMB selection flags.
        implicit none

        ! check Horava gravity. If Horava is selected the equation of state for dark energy should be -1
        if ( CP%EFTflag==4 .and. CP%FullMappingEFTmodel==1 ) CP%EFTwDE=0

    end subroutine EFTCAMB_check_flags_consistency

end module EFTinitialization

! -------------------------------------------------------------------------------------------------
