! -------------------------------------------------------------------------------------------------
!
!   EFTCAMB
!
!   Developed and implemented by:
!       Bin Hu (hu@lorentz.leidenuniv.nl), Marco Raveri (mraveri@sissa.it)
!       Noemi Frusciante (fruscian@iap.fr), Alessandra Silvestri (silvestri@lorentz.leidenuniv.nl)
!
!
!   EFTCAMB module for designer models. Includes the designer f(R) code.
!
!   For more informations about the methods contained in this file see
!   the documentation: arXiv:1405.3590
!
! -------------------------------------------------------------------------------------------------

! Solve the designer differential equations for the chosen model and computes,
! from the solution, the values of the EFT functions
module EFTdesigner

    use EFTDef
    use EFTdeEOS
    use MassiveNu

    implicit none

    ! 1) Definitions of the variables common to all the designer code.
    integer , parameter :: des_nstep = 1000        ! Number of points sampled by the designer code.
    integer , parameter :: des_ninterpol = 2       ! Number of points used by the interpolator. Only even.
    real(dl), parameter :: xInitial = log(10._dl**(-8._dl))    ! log(a start)
    real(dl), parameter :: xFinal   = 0._dl                    ! Log(a final)
    integer , parameter :: DesignerInFuture = 0    ! How many points in the future will the code sample
    real(dl), save, dimension(des_nstep+DesignerInFuture+1) :: xp_des   ! Array with the time sampling

    ! 2) Tables containing the sampled values of the EFT functions.
    real(dl), save, dimension(des_nstep+DesignerInFuture+1) :: OmegaDes, OmegapDes, OmegappDes,Omega3pDes
    real(dl), save, dimension(des_nstep+DesignerInFuture+1) :: EFT_Lambda_Des, EFT_Lambdadot_Des, EFT_c_Des, EFT_cdot_Des
    real(dl), save, dimension(des_nstep+DesignerInFuture+1) :: Des_Gamma1, Des_Gamma1p, Des_Gamma2, Des_Gamma2p
    real(dl), save, dimension(des_nstep+DesignerInFuture+1) :: Des_Gamma3, Des_Gamma3p, Des_Gamma4, Des_Gamma4p
    real(dl), save, dimension(des_nstep+DesignerInFuture+1) :: Des_Gamma5, Des_Gamma5p, Des_Gamma6, Des_Gamma6p
    real(dl), save, dimension(des_nstep+DesignerInFuture+1) :: Des_Gamma4pp

    ! 3) Definitions used by a specific model.
    !    1- designer f(R) models:
    integer , parameter  :: nvar_fR = 2
    real(dl), save, dimension(des_nstep+DesignerInFuture+1) :: BFunc, BPrime, fRFunc
    real(dl), save, dimension(des_nstep+DesignerInFuture+1) :: zp_fR ! zp_fR = Ricci/H0^2
    real(dl), save, dimension(nvar_fR, des_nstep+DesignerInFuture+1) :: yp_fR, dyp_fR ! yp_fR = f(R)/H0^2

    logical :: DesignerPrintFlag = .False.

contains

    ! ----------------------------

    subroutine EFTCAMB_Designer(success)
        ! This subroutine calls the correct designer subroutine depending on the selected model.
        implicit none

        logical, intent(inout) :: success

        if (Feedbacklevel>1) write(*,*) 'EFTCAMB: calling designer code'

        select case (CP%DesignerEFTmodel)
            case (1) ! designer f(R) models:
                call EFTCAMB_Designer_FR(success)
            case (2) ! designer mc quintessence:
                if (CP%EFTwDE==0) write(*,*) 'EFTCAMB: WARNING designer quintessence makes sense only if wDE/=-1'
                success = .true.
            case (3) ! designer RPH models:
                success = .true.
            case default
                write(*,*) 'EFTCAMB: wrong choice of designer model'
                stop
        end select

        return
    end subroutine EFTCAMB_Designer

    ! ----------------------------

    subroutine EFTCAMB_Designer_FR(success)
        ! This subroutine finds the appropriate initial conditions and solves
        ! the designer equation for f(R) models.
        implicit none

        logical :: success

        real(dl) zbrent

        integer  :: Debug_n=1, Debug_MaxNum
        real(dl) :: B0_Temp, Temp_A, Temp_h, TempMin, TempMax
        integer  :: ind=1 , i=1
        real(dl) :: HorizAsyntB, VertAsyntA, BLeftAsynt, BRightAsynt
        real(dl) :: realAp
        real(dl) :: BTemp1, BTemp2, ATemp1 = 1._dl, Atemp2 = 1._dl, ATemp3, ATemp4, BTemp3, BTemp4

        !  Find initial conditions for designer F(R).
        !  Initial conditions are found solving B0(A)=B0wanted.
        !  The next part of code is complicated by the fact that B0(A) is not continuous and we have to adopt ad-hoc strategies.

        ! 1) Find the horizontal asymptote of B0(A)
        do ind=10, 100,1
            ATemp1 = ATemp2
            ATemp2 = 10._dl**REAL(ind)
            BTemp1 = DesFR_BfuncA(ATemp1)
            BTemp2 = DesFR_BfuncA(ATemp2)
            if (Abs(BTemp1-BTemp2)/Abs(ATemp1-ATemp2).lt.1.d-100) exit
        end do
        HorizAsyntB = BTemp2

        if (FeedbackLevel>2) write(*,*) 'f(R) designer: horizontal asymptote =',  HorizAsyntB

        ! 2) Check that the value of B0 given is not the forbidden one: the one corresponding to the horizontal asynt.
        if (ABS(HorizAsyntB-CP%EFTB0)<1.d-15) then
            success=.false.
            return
        end if

        ! 3) Bracket the vertical asyntote
        ATemp1=-10._dl
        ATemp2=10._dl
        call zbrac(DesFR_BfuncA,ATemp1,ATemp2,success,HorizAsyntB)
        if (.not.success) then
            write(*,*) 'f(R) designer: failure of vert asynt bracketing'
            return
        end if

        ! 4) Find the vertical asyntote by tricking the root finding algorithm.
        VertAsyntA = zbrent(DesFR_BfuncA,ATemp1,ATemp2,1.d-100,HorizAsyntB,success)
        if (.not.success) then
            write(*,*) 'f(R) designer: vertical asyntote not found.'
            return
        end if

        if (FeedbackLevel>2) write(*,*) 'f(R) designer: vertical asymptote =',  VertAsyntA

        ! 5) Find values for A that are on the left and on the right of the asyntote.
        do ind=-10, -1, 1
            ATemp1 = VertAsyntA+10._dl**ind
            ATemp2 = VertAsyntA-10._dl**ind
            BTemp1 = DesFR_BfuncA(ATemp1)
            BTemp2 = DesFR_BfuncA(ATemp2)
            if (BTemp1*BTemp2<0._dl) exit
        end do

        ! 6) Extablish on which side of the asyntote to look for the solution.
        ATemp3 = ATemp1
        ATemp4 = ATemp2
        do ind=1, 10, 1
            ATemp3 = ATemp3 + 10._dl**ind
            ATemp4 = ATemp4 - 10._dl**ind
            BTemp3 = DesFR_BfuncA(ATemp3)
            BTemp4 = DesFR_BfuncA(ATemp4)
            if ((BTemp1-CP%EFTB0)*(BTemp3-CP%EFTB0)<0.or.(BTemp2-CP%EFTB0)*(BTemp4-CP%EFTB0)<0) exit
        end do
        if ((BTemp1-CP%EFTB0)*(BTemp3-CP%EFTB0)<0.and.(BTemp2-CP%EFTB0)*(BTemp4-CP%EFTB0)<0) then
            write(*,*) 'f(R) designer: the root is on both side ???'
            stop
        end if

        ! 7) Solve the equation B0(A)=B0_wanted.
        if ((BTemp1-CP%EFTB0)*(BTemp3-CP%EFTB0)<0) then
            realAp = zbrent(DesFR_BfuncA,ATemp1,ATemp3,1.d-50,CP%EFTB0,success)
            if (.not.success) then
                write(*,*) 'f(R) designer: right side solution not found.'
                stop
            end if
        else if ((BTemp2-CP%EFTB0)*(BTemp4-CP%EFTB0)<0) then
            realAp = zbrent(DesFR_BfuncA,ATemp2,ATemp4,1.d-50,CP%EFTB0,success)
            if (.not.success) then
                write(*,*) 'f(R) designer: right side solution not found'
                stop
            end if
        else
            write(*,*) 'f(R) designer: the root was not on the right side nor the left one...'
            stop
        end if

        if (FeedbackLevel>2) write(*,*) 'f(R) designer: initial condition A =', realAp

        ! 8) Check if the result found is compatible with the requested one. This is required only for debug.
        BTemp1 = DesFR_BfuncA(realAp)
        if (Feedbacklevel > 1) write(*,*) 'f(R) designer: B0 found=', BTemp1, 'B0 given= ',CP%EFTB0

        if (ABS(BTemp1-CP%EFTB0)/ABS(CP%EFTB0)>0.1_dl.and.ABS(BTemp1-CP%EFTB0)>1.d-8.and..false.) then
            write(*,*)  'EFTCAMB: designer code unable to find appropriate initial conditions'
            success = .false.
            return
        end if

        ! 9) Finally: solve the F(R) background equation.
        call DesFR_SolveDesignerEq(realAp)

        ! Code to Debug the f(R) designer code.
        if (DebugDesigner) then
            ! 1) Prints the function B0(A). This is used to debug the initial conditions part.
            call CreateTxtFile('./DebugDesignerB.dat',34)
            print*, 'f(R) designer: Printing B(A) results'
            TempMin     = -10._dl
            TempMax     = +10._dl
            Debug_MaxNum = 1000
            do Debug_n = 1, Debug_MaxNum
                Temp_h = TempMin +REAL(Debug_n-1)*(TempMax-TempMin)/REAL(Debug_MaxNum-1)
                write(34,*) Temp_h, DesFR_BFuncA(Temp_h)
            end do
            close(34)
            ! 2) Prints F(R) quantities.
            DesignerPrintFlag = .True.
            print*, 'f(R) designer: Printing F(R) results'
            call CreateTxtFile('./DebugDesigner.dat',33)
            call DesFR_SolveDesignerEq(realAp)
            close(33)
        end if

        return
    end subroutine EFTCAMB_Designer_FR

    ! ----------------------------

    function DesFR_BfuncA(Ap)
        ! This function gives the present day value of B given the initial conditions.
        implicit none

        real(dl), intent(in) :: Ap
        real(dl) :: DesFR_BfuncA
        call DesFR_SolveDesignerEq(Ap)
        DesFR_BfuncA = BFunc(des_nstep)

        return
    end function DesFR_BfuncA

    ! ----------------------------

    subroutine DesFR_SolveDesignerEq(coefA)
        ! This subroutine solves the designer equation for f(R) models for a given value
        ! of the initial conditions (coefA).
        implicit none

        ! Variable definitions:
        real(dl), intent(in) :: coefA
        real(dl) :: H0_EFT, Omegam_EFT, Omegavac_EFT, OmegaGamma_EFT, OmegaNu_EFT, Omegarad_EFT
        real(dl) :: EquivalenceScale_fR, Ratio_fR, Initial_B_fR, Initial_C_fR, PPlus
        real(dl) :: yStar, CoeffA_Part, yPlus
        real(dl) :: yv_fR(nvar_fR), dyv_fR(nvar_fR)
        real(dl) :: EFT_E_g, EFT_E_gp, EFT_E_gpp, EFT_E_gppp, EFun, EP, EPP, E3P, Ede, Edep, Edepp, fsubRR
        ! Massive neutrinos variables
        real(dl) :: rhonu, presnu, rhonu_tot, presnu_tot, grhormass_t, EFT_E_nu, EFT_EP_nu, EFT_EPP_nu, EFT_E3P_nu
        real(dl) :: presnudot_tot, adotoa, Hdot, presnudotdot_tot, presnudot, presnudotdot
        integer  :: nu_i
        real(dl) :: a, OmegaMassiveNu_EFT

        real(dl) :: x   ! x=log(a) is the time variable.
        real(dl) :: h   ! h is the timestep.
        integer  :: k,i

        real(dl) :: temp1, temp2, temp3, temp4  ! Debug variables

        ! 1) Cosmological parameters:
        H0_EFT = (CP%h0/c_EFT)*1000._dl
        Omegam_EFT = CP%omegab + CP%omegac
        Omegavac_EFT = CP%omegav
        OmegaMassiveNu_EFT = CP%omegan
        OmegaGamma_EFT = kappa_EFT/c_EFT**2*4*sigma_boltz_EFT/c_EFT**3*CP%tcmb**4*Mpc_EFT**2/(3._dl*H0_EFT**2) !8*pi*G/c^2*4*sigma_B/c^3 T^4
        OmegaNu_EFT = 7._dl/8*(4._dl/11)**(4._dl/3)*OmegaGamma_EFT*CP%Num_Nu_massless !7/8*(4/11)^(4/3)*grhog (per neutrino species)
        OmegaNu_EFT = grhornomass/(3*CP%H0**2/c_EFT**2*1000**2)
        Omegarad_EFT = OmegaGamma_EFT + OmegaNu_EFT
        EquivalenceScale_fR = (Omegarad_EFT+ OmegaMassiveNu_EFT)/Omegam_EFT
        Ratio_fR = EquivalenceScale_fR/Exp(xInitial)

        ! 2) Growing mode solution:
        Initial_B_fR = (7._dl + 8._dl*Ratio_fR)/(2._dl*(1._dl + Ratio_fR))
        Initial_C_fR = -3._dl/(2._dl*(1._dl + Ratio_fR))
        PPlus = 0.5_dl*(-Initial_B_fR + Sqrt(Initial_B_fR**2 - 4._dl*Initial_C_fR))
        yPlus = exp(PPlus*xInitial)
        !    Construction of the particolar solution:
        CoeffA_Part = (-6._dl*Initial_C_fR)/(-3._dl*Exp(xInitial)*EFTw(Exp(xInitial),1)&
            & +9._dl*EFTw(Exp(xInitial),0)**2&
            & +(18._dl-3._dl*Initial_B_fR)*EFTw(Exp(xInitial),0)&
            & +9._dl -3._dl*Initial_B_fR +Initial_C_fR)
        yStar = CoeffA_Part*Omegavac_EFT*Exp(-2._dl*xInitial)*EFTw(Exp(xInitial),3)

        !    Initial conditions:
        yv_fR(1)=coefA*yPlus + yStar
        yv_fR(2)=PPlus*coefA*yPlus - 3._dl*(1._dl+EFTw(Exp(xInitial),0))*yStar
        xp_des(1)=xInitial
        yp_fR(1,1)=yv_fR(1)
        yp_fR(2,1)=yv_fR(2)
        x=xInitial
        !    Function g(x) at initial time:
        EFT_E_g    = -(Log(EFTw(Exp(x),3)) -2._dl*x)/3._dl
        EFT_E_gp   = 1._dl +EFTw(Exp(x),0)
        EFT_E_gpp  = Exp(x)*EFTw(Exp(x),1)
        EFT_E_gppp = Exp(x)*EFTw(Exp(x),1) +Exp(2._dl*x)*EFTw(Exp(x),2)
        !   Massive neutrinos at initial conditions:
        rhonu_tot  = 0._dl
        presnu_tot = 0._dl
        EFT_E_nu   = 0._dl
        EFT_EP_nu  = 0._dl
        if (CP%Num_Nu_Massive /= 0) then
            do nu_i = 1, CP%Nu_mass_eigenstates
                a= Exp(x)
                rhonu  = 0._dl
                presnu = 0._dl
                grhormass_t=grhormass(nu_i)/a**2
                call Nu_background(a*nu_masses(nu_i),rhonu,presnu)
                rhonu_tot  = rhonu_tot + grhormass_t*rhonu
                presnu_tot = presnu_tot + grhormass_t*presnu

                EFT_E_nu   = EFT_E_nu + grhormass(nu_i)/3._dl/a**4/H0_EFT**2*rhonu
                EFT_EP_nu  = EFT_EP_nu - grhormass(nu_i)/H0_EFT**2/a**4*(rhonu +presnu)
            end do
        end if
        !    Initial values of energy, energy derivative and Ricci.
        EFun = +OmegaRad_EFT*exp(-4._dl*x)&
            & +Omegam_EFT*exp(-3._dl*x)&
            & +Omegavac_EFT*exp(-3._dl*EFT_E_g) + EFT_E_nu
        EP  = -4._dl*OmegaRad_EFT*exp(-4._dl*x)&
            & -3._dl*Omegam_EFT*exp(-3._dl*x)&
            & -3._dl*Omegavac_EFT*EFT_E_gp*exp(-3._dl*EFT_E_g) +EFT_EP_nu
        zp_fR(1) = 3._dl*(4._dl*EFun+EP)

        ! 3) Do the timestep loop:
        h=(xFinal-xInitial)/real(des_nstep-1)
        do k=1, des_nstep+DesignerInFuture

            call  DesFR_derivs(nvar_fR, x, yv_fR, dyv_fR)
            call  EFT_rk4(nvar_fR, yv_fR, dyv_fR, x, h, yv_fR, DesFR_derivs)

            x=x+h
            if (k.le.des_nstep+DesignerInFuture) then
                xp_des(k+1)=x
                yp_fR(1,k+1) = yv_fR(1)
                yp_fR(2,k+1) = yv_fR(2)
                dyp_fR(1,k+1) = dyv_fR(1)
                dyp_fR(2,k+1) = dyv_fR(2)
                ! Compute the g(x) functions
                EFT_E_g    = -(Log(EFTw(Exp(x),3)) -2._dl*x)/3._dl
                EFT_E_gp   = 1._dl +EFTw(Exp(x),0)
                EFT_E_gpp  = Exp(x)*EFTw(Exp(x),1)
                EFT_E_gppp = Exp(x)*EFTw(Exp(x),1) +Exp(2._dl*x)*EFTw(Exp(x),2)
                ! Compute Energy and its derivatives wrt to ln(a):

                ! First compute massive neutrinos contribution:
                rhonu_tot  = 0._dl
                presnu_tot = 0._dl
                EFT_E_nu   = 0._dl
                EFT_EP_nu  = 0._dl
                if (CP%Num_Nu_Massive /= 0) then
                    do nu_i = 1, CP%Nu_mass_eigenstates
                        a= Exp(x)
                        rhonu  = 0._dl
                        presnu = 0._dl
                        grhormass_t=grhormass(nu_i)/a**2
                        call Nu_background(a*nu_masses(nu_i),rhonu,presnu)
                        rhonu_tot  = rhonu_tot + grhormass_t*rhonu
                        presnu_tot = presnu_tot + grhormass_t*presnu

                        EFT_E_nu   = EFT_E_nu + grhormass(nu_i)/3._dl/a**4/H0_EFT**2*rhonu
                        EFT_EP_nu  = EFT_EP_nu - grhormass(nu_i)/H0_EFT**2/a**4*(rhonu +presnu)

                    end do
                end if
                ! Add its contribution to E and E':
                EFun = +OmegaRad_EFT*exp(-4._dl*x)&
                    & +Omegam_EFT*exp(-3._dl*x)&
                    & +Omegavac_EFT*exp(-3._dl*EFT_E_g) + EFT_E_nu
                EP  = -4._dl*OmegaRad_EFT*exp(-4._dl*x)&
                    & -3._dl*Omegam_EFT*exp(-3._dl*x)&
                    & -3._dl*Omegavac_EFT*EFT_E_gp*exp(-3._dl*EFT_E_g) +EFT_EP_nu
                ! Compute everything of massive nu again to get the time derivatives:
                rhonu_tot  = 0._dl
                presnu_tot = 0._dl
                presnudot_tot = 0._dl
                presnudotdot_tot = 0._dl
                EFT_E_nu   = 0._dl
                EFT_EP_nu  = 0._dl
                EFT_EPP_nu = 0._dl
                EFT_E3P_nu = 0._dl
                if (CP%Num_Nu_Massive /= 0) then
                    do nu_i = 1, CP%Nu_mass_eigenstates
                        a= Exp(x)
                        adotoa = +a*H0_EFT*sqrt(EFun)
                        Hdot   = +0.5_dl*H0_EFT**2*a**2*EP +adotoa**2

                        rhonu  = 0._dl
                        presnu = 0._dl
                        presnudot = 0._dl
                        presnudotdot = 0._dl

                        grhormass_t=grhormass(nu_i)/a**2

                        call Nu_background(a*nu_masses(nu_i),rhonu,presnu)
                        presnudot = Nu_pidot(a*nu_masses(nu_i),adotoa,presnu)
                        presnudotdot = Nu_pidotdot(a*nu_masses(nu_i),adotoa,Hdot,presnu,presnudot)

                        rhonu_tot  = rhonu_tot + grhormass_t*rhonu
                        presnu_tot = presnu_tot + grhormass_t*presnu
                        presnudot_tot  = presnudot_tot + grhormass_t*(presnudot -4._dl*adotoa*presnu)
                        presnudotdot_tot = presnudotdot_tot + grhormass_t*(presnudotdot &
                            & -8._dl*adotoa*presnudot +4._dl*presnu*(+4._dl*adotoa**2-Hdot))

                        EFT_E_nu   = EFT_E_nu + grhormass(nu_i)/3._dl/a**4/H0_EFT**2*rhonu
                        EFT_EP_nu  = EFT_EP_nu - grhormass(nu_i)/H0_EFT**2/a**4*(rhonu +presnu)
                        EFT_EPP_nu = EFT_EPP_nu + 3._dl/H0_EFT**2*grhormass(nu_i)/a**4*(rhonu +presnu)&
                            & -grhormass_t*(presnudot -4._dl*adotoa*presnu)/H0_EFT**3/sqrt(EFun)/a**3
                        EFT_E3P_nu = EFT_E3P_nu -9._dl/H0_EFT**2*grhormass(nu_i)/a**4*(rhonu +presnu)&
                            & +(3._dl/adotoa/H0_EFT**2/a**2+Hdot/adotoa**3/H0_EFT**2/a**2)&
                            &*grhormass_t*(presnudot -4._dl*adotoa*presnu)&
                            & -grhormass_t*(presnudotdot &
                            & -8._dl*adotoa*presnudot +4._dl*presnu*(+4._dl*adotoa**2-Hdot))/adotoa**2/H0_EFT**2/a**2
                    end do
                end if

                EPP = 16._dl*OmegaRad_EFT*exp(-4._dl*x)&
                    & +9._dl*Omegam_EFT*exp(-3._dl*x)&
                    & -3._dl*Omegavac_EFT*exp(-3._dl*EFT_E_g)*(EFT_E_gpp-3._dl*EFT_E_gp**2) +EFT_EPP_nu
                E3P = -64._dl*OmegaRad_EFT*exp(-4._dl*x)&
                    & -27._dl*Omegam_EFT*exp(-3._dl*x)&
                    & -3._dl*Omegavac_EFT*exp(-3._dl*EFT_E_g)*(EFT_E_gppp-9._dl*EFT_E_gp*EFT_E_gpp+9._dl*EFT_E_gp**3)&
                    & + EFT_E3P_nu
                Ede= +Omegavac_EFT*exp(-3._dl*EFT_E_g)
                Edep = -3._dl*Omegavac_EFT*EFT_E_gp*exp(-3._dl*EFT_E_g)
                Edepp = -3._dl*Omegavac_EFT*exp(-3._dl*EFT_E_g)*(EFT_E_gpp -3._dl*EFT_E_gp**2)
                ! Compute the Ricci:
                zp_fR(k+1) = 3._dl*(4._dl*EFun+EP)
                ! Compute fsubR:
                fRFunc(k+1) = dyp_fR(1,k+1)/3._dl/(4._dl*EP+EPP)
                ! Compute function B:
                BFunc(k+1) = 2._dl/3._dl/(1._dl+fRFunc(k+1))*1._dl/(4._dl*EP+EPP)*EFun/EP*&
                    &(dyp_fR(2,k+1)-dyp_fR(1,k+1)*(4._dl*EPP+E3P)/(4._dl*EP+EPP))
                ! Compute the EFT functions
                OmegaDes(k+1)   = fRFunc(k+1)
                OmegapDes(k+1)  = Exp(-x)/(6._dl*Efun*(EPP+4._dl*EP))*(-EPP*(6._dl*Ede+yp_fR(1,k+1))&
                    &+EP*(dyp_fR(1,k+1)-4._dl*(6._dl*Ede+yp_fR(1,k+1)))+2._dl*Efun*dyp_fR(1,k+1))
                OmegappDes(k+1) = Exp(-2._dl*x)/(12._dl*Efun**2*(EPP+4._dl*EP))*&
                    &(EP*(EP*(24._dl*Ede-dyp_fR(1,k+1)+4._dl*yp_fR(1,k+1))-6._dl*EFun*(8._dl*Edep+dyp_fR(1,k+1)))&
                    &+EPP*(EP*(6._dl*Ede+yp_fR(1,k+1))-12._dl*Efun*Edep))
                Omega3pDes(k+1) = Exp(-3._dl*x)/(24._dl*Efun**3*(EPP+4._dl*EP))*&
                    &(2._dl*Efun*(EPP**2*(6._dl*Ede+yp_fR(1,k+1))&
                    & +4._dl*EP**2*(18._dl*Edep+6._dl*Ede+2._dl*dyp_fR(1,k+1)+yp_fR(1,k+1))&
                    & +EP*EPP*(18._dl*Edep+30._dl*Ede-dyp_fR(1,k+1)+5._dl*yp_fR(1,k+1)))&
                    & +12._dl*Efun**2*(EPP*(-2._dl*Edepp+4._dl*Edep-dyp_fR(1,k+1))&
                    & +EP*(-8._dl*Edepp+16._dl*Edep+dyp_fR(1,k+1)))&
                    & +3._dl*EP**2*(EP*(dyp_fR(1,k+1)-4._dl*(6._dl*Ede+yp_fR(1,k+1)))-EPP*(6._dl*Ede+yp_fR(1,k+1))))
                EFT_Lambda_Des(k+1)    = 0.5_dl*H0_EFT**2*(yp_fR(1,k+1)-fRFunc(k+1)*zp_fR(k+1))*Exp(x)**2
                EFT_Lambdadot_Des(k+1) = -1.5_dl*H0_EFT**3*Sqrt(Efun)*(Exp(x)**4*OmegapDes(k+1)*(4._dl*Efun+EP))
                ! Compute other f(R) quantities:
                fsubRR = 0.5_dl*EP/Efun/3._dl*BFunc(k+1)*(1._dl+fRFunc(k+1))/(4._dl*EP+EPP)
                ! Debug output:
                if (DebugDesigner.and.DesignerPrintFlag) then
                !     write (33,'(20E15.5)') x, Exp(x), zp_fR(k+1), yp_fR(1,k+1), &
                !         &OmegaDes(k+1), OmegapDes(k+1), OmegappDes(k+1), Omega3pDes(k+1), EFT_Lambda_Des(k+1), EFT_Lambdadot_Des(k+1),&
                !         &BFunc(k+1), fsubRR
                end if
            end if
        end do

        return
    end subroutine DesFR_SolveDesignerEq

    ! ----------------------------

    subroutine DesFR_derivs(n,x,y_fR,dydx_fR)
        ! This subroutine computes the designer f(R) equation of motion.
        implicit none

        integer n
        real(dl) :: x
        real(dl), dimension(n) :: y_fR,dydx_fR
        real(dl) :: EFunction, EFunPrime, EFunPrime2, EFunPrime3
        real(dl) :: EFTwde0, EFTwde1, EFTwde2, EFTwde3
        real(dl) :: EFT_E_gfun, EFT_E_gfunp, EFT_E_gfunpp, EFT_E_gfunppp
        real(dl) :: H0_EFT, Omegam_EFT, Omegavac_EFT, OmegaGamma_EFT, OmegaNu_EFT, Omegarad_EFT
        ! Massive neutrinos variables
        real(dl) :: rhonu, presnu, rhonu_tot, presnu_tot, grhormass_t, EFT_E_nu, EFT_EP_nu, EFT_EPP_nu, EFT_E3P_nu
        real(dl) :: presnudot_tot, adotoa, Hdot, presnudotdot_tot, presnudot, presnudotdot
        integer  :: nu_i
        real(dl) :: a
        ! 1) Cosmological parameters:
        H0_EFT = (CP%h0/c_EFT)*1000._dl
        Omegam_EFT = CP%omegab + CP%omegac
        Omegavac_EFT = CP%omegav
        OmegaGamma_EFT = kappa_EFT/c_EFT**2*4*sigma_boltz_EFT/c_EFT**3*CP%tcmb**4*Mpc_EFT**2/(3._dl*H0_EFT**2) !8*pi*G/c^2*4*sigma_B/c^3 T^4
        OmegaNu_EFT = 7._dl/8*(4._dl/11)**(4._dl/3)*OmegaGamma_EFT !7/8*(4/11)^(4/3)*grhog (per neutrino species)
        OmegaNu_EFT = grhornomass/(3*CP%H0**2/c_EFT**2*1000**2)
        Omegarad_EFT = OmegaGamma_EFT + OmegaNu_EFT
        ! 2) Compute the function g(x) and its derivatives:
        EFT_E_gfun    = -(Log(EFTw(Exp(x),3)) -2._dl*x)/3._dl
        EFT_E_gfunp   = 1._dl +EFTw(Exp(x),0)
        EFT_E_gfunpp  = Exp(x)*EFTw(Exp(x),1)
        EFT_E_gfunppp = Exp(x)*EFTw(Exp(x),1) +Exp(2._dl*x)*EFTw(Exp(x),2)
        ! 3) Compute energy and its derivatives:

        ! First compute massive neutrinos contribution:
        rhonu_tot  = 0._dl
        presnu_tot = 0._dl
        EFT_E_nu   = 0._dl
        EFT_EP_nu  = 0._dl
        if (CP%Num_Nu_Massive /= 0) then
            do nu_i = 1, CP%Nu_mass_eigenstates
                a= Exp(x)
                rhonu  = 0._dl
                presnu = 0._dl
                grhormass_t=grhormass(nu_i)/a**2
                call Nu_background(a*nu_masses(nu_i),rhonu,presnu)
                rhonu_tot  = rhonu_tot + grhormass_t*rhonu
                presnu_tot = presnu_tot + grhormass_t*presnu

                EFT_E_nu   = EFT_E_nu + grhormass(nu_i)/3._dl/a**4/H0_EFT**2*rhonu
                EFT_EP_nu  = EFT_EP_nu - grhormass(nu_i)/H0_EFT**2/a**4*(rhonu +presnu)

            end do
        end if

        ! Add its contribution to E and E':

        EFunction = +OmegaRad_EFT*exp(-4._dl*x)&
            & +Omegam_EFT*exp(-3._dl*x)&
            & +Omegavac_EFT*exp(-3._dl*EFT_E_gfun) + EFT_E_nu
        EFunPrime = -4._dl*OmegaRad_EFT*exp(-4._dl*x)&
            & -3._dl*Omegam_EFT*exp(-3._dl*x)&
            & -3._dl*Omegavac_EFT*EFT_E_gfunp*exp(-3._dl*EFT_E_gfun) +EFT_EP_nu

        ! Compute everything of massive nu again to get the time derivatives:
        rhonu_tot  = 0._dl
        presnu_tot = 0._dl
        presnudot_tot = 0._dl
        presnudotdot_tot = 0._dl
        EFT_E_nu   = 0._dl
        EFT_EP_nu  = 0._dl
        EFT_EPP_nu = 0._dl
        EFT_E3P_nu = 0._dl
        if (CP%Num_Nu_Massive /= 0) then
            do nu_i = 1, CP%Nu_mass_eigenstates
                a= Exp(x)
                adotoa = +a*H0_EFT*sqrt(EFunction)
                Hdot   = +0.5_dl*H0_EFT**2*a**2*EFunPrime +adotoa**2

                rhonu  = 0._dl
                presnu = 0._dl
                presnudot = 0._dl
                presnudotdot = 0._dl

                grhormass_t=grhormass(nu_i)/a**2

                call Nu_background(a*nu_masses(nu_i),rhonu,presnu)
                presnudot = Nu_pidot(a*nu_masses(nu_i),adotoa,presnu)
                presnudotdot = Nu_pidotdot(a*nu_masses(nu_i),adotoa,Hdot,presnu,presnudot)

                rhonu_tot  = rhonu_tot + grhormass_t*rhonu
                presnu_tot = presnu_tot + grhormass_t*presnu
                presnudot_tot  = presnudot_tot + grhormass_t*(presnudot -4._dl*adotoa*presnu)
                presnudotdot_tot = presnudotdot_tot + grhormass_t*(presnudotdot &
                    & -8._dl*adotoa*presnudot +4._dl*presnu*(+4._dl*adotoa**2-Hdot))

                EFT_E_nu   = EFT_E_nu + grhormass(nu_i)/3._dl/a**4/H0_EFT**2*rhonu
                EFT_EP_nu  = EFT_EP_nu - grhormass(nu_i)/H0_EFT**2/a**4*(rhonu +presnu)
                EFT_EPP_nu = EFT_EPP_nu + 3._dl/H0_EFT**2*grhormass(nu_i)/a**4*(rhonu +presnu)&
                    & -grhormass_t*(presnudot -4._dl*adotoa*presnu)/H0_EFT**3/sqrt(EFunction)/a**3
                EFT_E3P_nu = EFT_E3P_nu -9._dl/H0_EFT**2*grhormass(nu_i)/a**4*(rhonu +presnu)&
                    & +(3._dl/adotoa/H0_EFT**2/a**2+Hdot/adotoa**3/H0_EFT**2/a**2)&
                    &*grhormass_t*(presnudot -4._dl*adotoa*presnu)&
                    & -grhormass_t*(presnudotdot &
                    & -8._dl*adotoa*presnudot +4._dl*presnu*(+4._dl*adotoa**2-Hdot))/adotoa**2/H0_EFT**2/a**2
            end do
        end if

        EFunPrime2 = 16._dl*OmegaRad_EFT*exp(-4._dl*x)&
            & +9._dl*Omegam_EFT*exp(-3._dl*x)&
            & -3._dl*Omegavac_EFT*exp(-3._dl*EFT_E_gfun)*(EFT_E_gfunpp -3._dl*EFT_E_gfunp**2) + EFT_EPP_nu
        EFunPrime3 = -64._dl*OmegaRad_EFT*exp(-4._dl*x)&
            & -27._dl*Omegam_EFT*exp(-3._dl*x)&
            & -3._dl*Omegavac_EFT*exp(-3._dl*EFT_E_gfun)*&
            &(EFT_E_gfunppp-9._dl*EFT_E_gfunp*EFT_E_gfunpp+9._dl*EFT_E_gfunp**3) + EFT_E3P_nu

        ! 4) Get the equation of motion:
        dydx_fR(1) = y_fR(2)
        dydx_fR(2) = (1._dl+0.5_dl*EFunPrime/EFunction+(4._dl*EFunPrime2+EFunPrime3)/(4._dl*EFunPrime+EFunPrime2))*y_fR(2) &
            & -0.5_dl*(4._dl*EFunPrime+EFunPrime2)/EFunction*y_fR(1) &
            & -3._dl*Omegavac_EFT*exp(-3._dl*EFT_E_gfun)*(4._dl*EFunPrime+EFunPrime2)/EFunction
        return
    end subroutine DesFR_derivs

    ! ----------------------------

    function Designed_EFT_Function(EFTFunctionTable, a, GR_value)
        ! The designer code will provide a table of sampled values for the EFT functions.
        ! This function is called in the EFTfunctions module to interpolate those tables.
        implicit none

        real(dl), intent(in) :: a, GR_value
        real(dl) :: EFTFunctionTable(des_nstep+DesignerInFuture+1)
        real(dl) :: Designed_EFT_Function

        real(dl) :: x, temp, dtemp
        real(dl) :: xb(des_ninterpol),yb(des_ninterpol)
        integer  :: i, jlo, stint

        x = log(a)

        if (x.lt.xp_des(1)) then
            temp = GR_value

        else if (x.ge.xp_des(1) .and. x.le.xp_des(des_nstep+DesignerInFuture)) then
            ! 1) find the point corresponing to the requested input.
            call hunt(xp_des, des_nstep+DesignerInFuture, x, jlo)
            ! 2) construct the table that will be interpolated
            stint = des_ninterpol/2
            if((des_nstep+DesignerInFuture-jlo-stint).lt.0) then
                ! Last des_ninterpol points
                do i=1, des_ninterpol
                    xb(i)=xp_des(des_nstep+DesignerInFuture-des_ninterpol+i)
                    yb(i)=EFTFunctionTable(des_nstep+DesignerInFuture-des_ninterpol+i)
                end do
            else if (jlo.eq.1) then
                ! First des_ninterpol points
                do i=1, des_ninterpol
                    xb(i)=xp_des(i)
                    yb(i)=EFTFunctionTable(i)
                end do
            else
                ! Surrounding des_ninterpol points
                do i=1, des_ninterpol
                    xb(i)=xp_des(jlo-stint+i)
                    yb(i)=EFTFunctionTable(jlo-stint+i)
                end do
            endif
            ! 3) Call the Neville interpolator.
            call Polint(des_ninterpol, xb, yb, x, temp, dtemp)

        else
            temp = EFTFunctionTable(des_nstep+DesignerInFuture)

        end if

        Designed_EFT_Function = temp

    end function Designed_EFT_Function

end module EFTdesigner

! -------------------------------------------------------------------------------------------------
