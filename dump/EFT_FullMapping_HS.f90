! -------------------------------------------------------------------------------------------------
!
!   EFT_FullMapping_HS module
!
!   Developed and implemented by:
!       Matteo Rizzato (matteo.rizzato@studenti.unipd.it) (Hu-Sawicki f(R) model, 0705.1158)
!       Bin HU (hu@lorentz.leidenuniv.nl)
!
!   EFT_FullMapping module provides EFTCAMB with consistent background (EFTw(a)) and
!   perturbation (EFT functions, e.g. EFTOmega, EFTLambda, ...) quantities.
!
!   New developers should check the consistency of between background and perturbation quantities.
!
! -------------------------------------------------------------------------------------------------

module EFT_FullMapping_HS
    use EFT_FullMapping_common
    implicit none

    !-----------definitions-and-setting-parameters-values---------------------------------------!
    real(dl) :: fm_hs_m2,fm_hs_c1,fm_hs_c2,fm_hs_n,fm_hs_fro           !HuSawicki param!
    real(dl) :: fm_hs_Omegacdmh2,fm_hs_Omegagammah2                    !Comos param!
    real(dl) :: fm_hs_Omeganuh2,fm_hs_Omegadeh2,fm_hs_h                !Comos param!
    real(dl), save, dimension(fm_nstep+1) :: fm_hs_ourguess ,fm_hs_ourguessprime
    real(dl), save, dimension(fm_nstep+1) :: fm_hs_eqstate,fm_hs_v,fm_hs_yH,fm_hs_vtrue
    real(dl), save, dimension(fm_nstep+1) :: fm_hs_yHtrue,fm_hs_c,fm_hs_Lambda,fm_hs_Omega,fm_hs_rhoeff
    real(dl), save, dimension(fm_nstep+1) :: fm_hs_efold
    real(dl) :: fm_hs_llim,fm_hs_ulim,fm_hs_Zetastar,fm_hs_Totmass
    real(dl) :: fm_hs_Neff,fm_hs_mgamma2,fm_hs_rnge,fm_hs_step
    real(dl) :: fm_hs_z_ast_threshold
    integer :: fm_hs_nnosteps

    real(dl), parameter :: solver_switch = -5._dl

contains

    !---------------------------------------------------------------------------!
    subroutine FM_HS_main(success)

        implicit none
        logical, intent(inout) :: success

        call FM_HS_Initialization(success)
        call FM_HS_solver(success) !... solved yH,rho_eff, EoS !bh
        call FM_HS_output(success) !... calculate background and EFT functions !bh

        return
    end subroutine FM_HS_main

    !---------------------------------------------------------------------------!
    subroutine FM_HS_Initialization(success)

        implicit none

        logical, intent(inout) :: success

        !-----------parameters-fixed-by-the-user-in-full-options-version-------!
        !bh: this part need to update with params.ini
        fm_hs_n = CP%EFT_FM_HS_n
        fm_hs_fro = CP%EFT_FM_HS_fR_0
        !fm_hs_Zetastar = 0.0 !-2.6
        fm_hs_Totmass = 94.*CP%omegan*(CP%H0/100)**2
        fm_hs_h = CP%h0/100.
        fm_hs_Omegacdmh2 = CP%omegab*(CP%H0/100)**2 + CP%omegac*(CP%H0/100)**2
        fm_hs_Omegagammah2 = 2.469D-5
        fm_hs_Omeganuh2 = fm_hs_Totmass/94.
        fm_hs_Omegadeh2 = (1-fm_hs_Omegacdmh2/fm_hs_h**2-fm_hs_Omegagammah2/fm_hs_h**2-fm_hs_Omeganuh2/fm_hs_h**2)*fm_hs_h**2
        fm_hs_llim = fm_xInitial
        fm_hs_ulim = fm_xFinal
        fm_hs_nnosteps = fm_nstep+1
        fm_hs_Neff = CP%Num_Nu_massless + CP%Num_Nu_Massive
        fm_hs_z_ast_threshold = 5.d-5

        if (Feedbacklevel>0) then
            print*,'n=',fm_hs_n
            print*,'fR_0=',fm_hs_fro
            print*,'fm_hs_mnu',fm_hs_Totmass
            print*,'h0=',fm_hs_h
            print*,'fm_hs_Omegacdmh2',fm_hs_Omegacdmh2
            print*,'Neff=',fm_hs_Neff
        end if


        !------------rho,rhop,rhopp,rhoppp--------!
        !allocate (fm_hs_ourguess(fm_hs_nnosteps))
        !allocate (fm_hs_ourguessprime(fm_hs_nnosteps))
        !allocate (fm_hs_eqstate(fm_hs_nnosteps))
        !allocate (fm_hs_rhoeff(fm_hs_nnosteps))
        !allocate (fm_hs_v(fm_hs_nnosteps))
        !allocate (fm_hs_yH(fm_hs_nnosteps))
        !allocate (fm_hs_yHtrue(fm_hs_nnosteps))
        !allocate (fm_hs_vtrue(fm_hs_nnosteps))
        !allocate (fm_hs_c(fm_hs_nnosteps))
        !allocate (fm_hs_Omega(fm_hs_nnosteps))
        !allocate (fm_hs_Lambda(fm_hs_nnosteps))
        !allocate (fm_hs_efold(fm_hs_nnosteps))

        fm_hs_c2 = -fm_hs_n*(6.*fm_hs_Omegadeh2/fm_hs_Omegacdmh2)*( 3.+12.*fm_hs_Omegadeh2/fm_hs_Omegacdmh2)**(-fm_hs_n-1.)*(1./fm_hs_fro)
        fm_hs_c1 = 6.*(fm_hs_Omegadeh2/fm_hs_Omegacdmh2)*fm_hs_c2
        fm_hs_m2 = ((8320.5)**(-2.))*(fm_hs_Omegacdmh2/0.13)
        fm_hs_mgamma2 = fm_hs_Omegagammah2*(1./9.)*10**(-6.)
        fm_hs_rnge = ABS(fm_hs_llim-fm_hs_ulim)
        fm_hs_step = fm_hs_rnge/real(fm_nstep)

        success = .true.

        return
    end subroutine FM_HS_Initialization


    subroutine FM_HS_derivs( num_eq, x, y, ydot )

        implicit none

        integer  :: num_eq
        real(dl) :: x
        real(dl), dimension(num_eq) :: y
        real(dl), dimension(num_eq) :: ydot

        real(dl) :: dummy

        dummy = fm_hs_derivativeForV( x, y(1), y(2), fm_hs_m2, fm_hs_c1, fm_hs_c2, fm_hs_n)

        ydot(1) = y(2)
        ydot(2) = dummy

    end subroutine FM_HS_derivs

    subroutine FM_HS_Jac( num_eq, x, y, ml, mu, pd, nrowpd )

        implicit none

        integer  :: num_eq
        integer  :: ml ! ignored
        integer  :: mu ! ignored
        integer  :: nrowpd
        real(dl) :: x
        real(dl), dimension(num_eq) :: y
        real(dl), dimension(nrowpd,num_eq) :: pd

    end subroutine FM_HS_Jac

    !---------------------------------------------------------------------------!
    !... solve the background defferential eq. in HS model
    subroutine FM_HS_solver(success)

        implicit none

        integer :: fm_hs_ni
        real(dl) :: fm_hs_i
        real(dl) :: fm_hs_h1,fm_hs_h2,fm_hs_h3,fm_hs_h4,fm_hs_g1,fm_hs_g2,fm_hs_g3,fm_hs_g4
        real(dl) :: fm_hs_t
        logical, intent(inout) :: success

        integer  :: neq, itol, itask, istate, iopt, LRN, LRS, LRW, LIS, LIN, LIW, JacobianMode
        real(dl) :: t0, tfin, rtol, atol, t1, t2
        real(dl), allocatable :: rwork(:), y(:)
        integer,  allocatable :: iwork(:)

        ! Initialization of LSODA:
        neq  = 2
        t0   = fm_hs_llim
        tfin = fm_hs_ulim
        ! set-up the relative and absolute tollerances:
        itol = 1
        rtol = 1.d-10
        atol = 1.d-14
        ! initialize task to do:
        itask  = 1
        istate = 1
        iopt   = 1
        ! initialize the work space:
        LRN = 20 + 16*NEQ
        LRS = 22 + 9*NEQ + NEQ**2
        LRW = max(LRN,LRS)
        LIS = 20 + NEQ
        LIN = 20
        LIW = max(LIS,LIN)
        ! allocate the arrays:
        allocate(rwork(LRW))
        allocate(iwork(LIW))
        allocate(y(neq))
        ! optional lsoda input:
        RWORK(5) = 0._dl  ! the step size to be attempted on the first step. The default value is determined by the solver.
        RWORK(6) = 0._dl  ! the maximum absolute step size allowed. The default value is infinite.
        RWORK(7) = 0._dl  ! the minimum absolute step size allowed. The default value is 0.
        IWORK(5) = 0  !  flag to generate extra printing at method switches. IXPR = 0 means no extra printing (the default). IXPR = 1 means print data on each switch.
        IWORK(6) = 100  !  maximum number of (internally defined) steps allowed during one call to the solver. The default value is 500.
        IWORK(7) = 0  !  maximum number of messages printed (per problem) warning that T + H = T on a step (H = step size). This must be positive to result in a non-default value.  The default value is 10.
        IWORK(8) = 0  !  the maximum order to be allowed for the nonstiff (Adams) method.  the default value is 12.
        IWORK(9) = 0  !  the maximum order to be allowed for the stiff (BDF) method.  The default value is 5.
        ! additional lsoda stuff:
        CALL XSETF(0) ! suppress odepack printing
        ! Jacobian mode: 1=fullJacobian, 2=not provided
        JacobianMode = 2

        !-----------runge-kutta-(cdm)-------------------------------------!
        fm_hs_i = 1.

        ! set initial conditions:
        ! initial time:
        fm_hs_t = fm_hs_llim
        ! guess and guess first derivative:
        fm_hs_yH(1) = fm_hs_sourceFunction(fm_hs_m2,fm_hs_c1,fm_hs_c2,fm_hs_n,fm_hs_t)/fm_hs_massFunction(fm_hs_m2,fm_hs_c1,fm_hs_c2,fm_hs_n,fm_hs_t)
        fm_hs_v(1) = fm_hs_firstderive(fm_hs_m2,fm_hs_c1,fm_hs_c2,fm_hs_n,fm_hs_t)
        ! initial conditions:
        fm_hs_yHtrue(1) = fm_hs_yH(1)*EXP(-3.*fm_hs_t)+fm_hs_c1/(fm_hs_c2*6.)
        fm_hs_vtrue(1) = -3.*EXP(-3.*fm_hs_t)*fm_hs_yH(1)+EXP(-3.*fm_hs_t)*fm_hs_v(1)
        ! derived quantities:
        fm_hs_eqstate(1) = -1-(fm_hs_vtrue(1)/(3*fm_hs_yHtrue(1)))
        fm_hs_rhoeff(1) = 3*fm_hs_m2*fm_hs_yHtrue(1)
        fm_hs_i = fm_hs_i+1.
        ! initial conditions for LSODA:
        y(1) = fm_hs_yHtrue(1)
        y(2) = fm_hs_vtrue(1)

        ! solve the ODE:
        DO fm_hs_ni = 2,fm_hs_nnosteps,1

            t1 = fm_hs_t
            t2 = fm_hs_llim+fm_hs_step*(fm_hs_i-1.)
            fm_hs_t = t2

            if ( .false. ) then

                !... evolving full non-linear eq. using rk4 method
                fm_hs_h1 = fm_hs_step*fm_hs_vtrue(fm_hs_ni-1)
                fm_hs_g1 = fm_hs_step*fm_hs_derivativeForV(fm_hs_llim+(fm_hs_i-2.)*fm_hs_step,fm_hs_yHtrue(fm_hs_ni-1),fm_hs_vtrue(fm_hs_ni-1),fm_hs_m2,fm_hs_c1,fm_hs_c2,fm_hs_n)
                fm_hs_h2 = fm_hs_step*(fm_hs_vtrue(fm_hs_ni-1)+fm_hs_h1/2)
                fm_hs_g2 = fm_hs_step*fm_hs_derivativeForV(fm_hs_llim+(fm_hs_i-2.)*fm_hs_step+fm_hs_step/2.,fm_hs_yHtrue(fm_hs_ni-1)+fm_hs_h1/2,&
                    &fm_hs_vtrue(fm_hs_ni-1)+fm_hs_g1/2.,fm_hs_m2,fm_hs_c1,fm_hs_c2,fm_hs_n)
                fm_hs_h3 = fm_hs_step*(fm_hs_vtrue(fm_hs_ni-1)+fm_hs_h2/2)
                fm_hs_g3 = fm_hs_step*fm_hs_derivativeForV(fm_hs_llim+fm_hs_step*(fm_hs_i-2.)+fm_hs_step/2.,fm_hs_yHtrue(fm_hs_ni-1)+fm_hs_h2/2,fm_hs_vtrue(fm_hs_ni-1)+&
                    &fm_hs_g2/2.,fm_hs_m2,fm_hs_c1,fm_hs_c2,fm_hs_n)
                fm_hs_h4 = fm_hs_step*(fm_hs_vtrue(fm_hs_ni-1)+fm_hs_h3)
                fm_hs_g4 = fm_hs_step*fm_hs_derivativeForV(fm_hs_llim+(fm_hs_i-2.)*fm_hs_step+fm_hs_step,fm_hs_yHtrue(fm_hs_ni-1)+fm_hs_h3,fm_hs_vtrue(fm_hs_ni-1)+fm_hs_g3,&
                    &fm_hs_m2,fm_hs_c1,fm_hs_c2,fm_hs_n)

                fm_hs_v(fm_hs_ni) = fm_hs_firstderive(fm_hs_m2,fm_hs_c1,fm_hs_c2,fm_hs_n,fm_hs_t)
                fm_hs_yH(fm_hs_ni) = fm_hs_sourceFunction(fm_hs_m2,fm_hs_c1,fm_hs_c2,fm_hs_n,fm_hs_t)/fm_hs_massFunction(fm_hs_m2,fm_hs_c1,fm_hs_c2,fm_hs_n,fm_hs_t)
                fm_hs_vtrue(fm_hs_ni) = fm_hs_vtrue(fm_hs_ni-1) + (fm_hs_g1+fm_hs_g2*2.+fm_hs_g3*2.+fm_hs_g4)/6.
                fm_hs_yHtrue(fm_hs_ni) = fm_hs_yHtrue(fm_hs_ni-1) + (fm_hs_h1+fm_hs_h2*2.+fm_hs_h3*2.+fm_hs_h4)/6.

            else if ( fm_hs_t < solver_switch ) then

                !... using the particular guess
                fm_hs_yH(fm_hs_ni) = fm_hs_sourceFunction(fm_hs_m2,fm_hs_c1,fm_hs_c2,fm_hs_n,fm_hs_t)/fm_hs_massFunction(fm_hs_m2,fm_hs_c1,fm_hs_c2,fm_hs_n,fm_hs_t)
                fm_hs_v(fm_hs_ni) = fm_hs_firstderive(fm_hs_m2,fm_hs_c1,fm_hs_c2,fm_hs_n,fm_hs_t)

                fm_hs_yHtrue(fm_hs_ni) = fm_hs_yH(fm_hs_ni)*EXP(-3.*fm_hs_t)+fm_hs_c1/(fm_hs_c2*6.)
                fm_hs_vtrue(fm_hs_ni) = -3.*EXP(-3.*fm_hs_t)*fm_hs_yH(fm_hs_ni)+EXP(-3.*fm_hs_t)*fm_hs_v(fm_hs_ni)

                y(1) = fm_hs_yHtrue(fm_hs_ni)
                y(2) = fm_hs_vtrue(fm_hs_ni)

            else if ( fm_hs_t >= solver_switch ) then

                call DLSODA ( FM_HS_derivs, neq, y, t1, t2, itol, rtol, atol, itask, istate, iopt, RWORK, LRW, IWORK, LIW, FM_HS_Jac, JacobianMode)

                ! check the quality of the output:
                if ( istate<0 ) then

                    fm_hs_yH(fm_hs_ni) = fm_hs_sourceFunction(fm_hs_m2,fm_hs_c1,fm_hs_c2,fm_hs_n,fm_hs_t)/fm_hs_massFunction(fm_hs_m2,fm_hs_c1,fm_hs_c2,fm_hs_n,fm_hs_t)
                    fm_hs_v(fm_hs_ni) = fm_hs_firstderive(fm_hs_m2,fm_hs_c1,fm_hs_c2,fm_hs_n,fm_hs_t)

                    fm_hs_yHtrue(fm_hs_ni) = fm_hs_yH(fm_hs_ni)*EXP(-3.*fm_hs_t)+fm_hs_c1/(fm_hs_c2*6.)
                    fm_hs_vtrue(fm_hs_ni) = -3.*EXP(-3.*fm_hs_t)*fm_hs_yH(fm_hs_ni)+EXP(-3.*fm_hs_t)*fm_hs_v(fm_hs_ni)

                    y(1) = fm_hs_yHtrue(fm_hs_ni)
                    y(2) = fm_hs_vtrue(fm_hs_ni)

                    istate = 1

                end if

                fm_hs_yHtrue(fm_hs_ni) = y(1)
                fm_hs_vtrue(fm_hs_ni)  = y(2)

            end if


            fm_hs_eqstate(fm_hs_ni) = -1-(fm_hs_vtrue(fm_hs_ni)/(3*fm_hs_yHtrue(fm_hs_ni)))
            fm_hs_rhoeff(fm_hs_ni) = 3*fm_hs_m2*fm_hs_yHtrue(fm_hs_ni)
            fm_hs_i = fm_hs_i+1.

        END DO

        success = .true.

        return
    end subroutine FM_HS_solver

    !---------------------------------------------------------------------------!
    !... calculate background and EFT functions
    subroutine FM_HS_output(success)

        implicit none

        logical, intent(inout) :: success
        integer :: fm_hs_ni
        real(dl) :: fm_hs_i
        real(dl) :: fm_hs_rhoeffprime
        real(dl) :: fm_hs_t,fm_hs_rho,fm_hs_rhop,fm_hs_Ricci,fm_hs_Hubble,fm_hs_Hprime

        !--------------------------------Hybride-Part---------------------------------!
        fm_hs_i=1.
        DO fm_hs_ni = 1,fm_hs_nnosteps,1
            fm_hs_t = fm_hs_llim+fm_hs_step*(fm_hs_i-1.)
            fm_hs_efold(fm_hs_ni) = fm_hs_t
            fm_hs_rhoeffprime = 3*fm_hs_m2*fm_hs_vtrue(fm_hs_ni)
            fm_hs_rho = fm_hs_rhofunction(fm_hs_t,fm_hs_mgamma2,fm_hs_Neff,fm_hs_Omeganuh2)
            fm_hs_rhop = fm_hs_rhoprime(fm_hs_t,fm_hs_mgamma2,fm_hs_Neff,fm_hs_Omeganuh2)
            fm_hs_Hubble = SQRT(fm_hs_m2*EXP(-3.*fm_hs_t) + fm_hs_mgamma2*EXP(-4.*fm_hs_t) + fm_hs_rho/3. +&
                fm_hs_rhoeff(fm_hs_ni)/3.)
            fm_hs_Hprime = ((-3*fm_hs_m2)/EXP(3*fm_hs_t) - (4*fm_hs_mgamma2)/EXP(4*fm_hs_t) + fm_hs_rhop/3. + fm_hs_rhoeffprime/3.)/(2*fm_hs_Hubble)
            fm_hs_Ricci = 6*fm_hs_Hprime*fm_hs_Hubble + 12*(fm_hs_Hubble**2)
            fm_hs_Lambda(fm_hs_ni) = (1./2.)*(-(fm_hs_c1/fm_hs_c2)*fm_hs_m2 + (fm_hs_c1/(fm_hs_c2**2.))*((fm_hs_m2/fm_hs_Ricci)**(fm_hs_n+1.))*fm_hs_Ricci -&
                fm_hs_Ricci*(-(fm_hs_c1/(fm_hs_c2**2.))*fm_hs_n*(fm_hs_m2/fm_hs_Ricci)**(1.+fm_hs_n)))
            fm_hs_Omega(fm_hs_ni) = -(fm_hs_c1/(fm_hs_c2**2.))*fm_hs_n*((fm_hs_m2/fm_hs_Ricci)**(1.+fm_hs_n))
            fm_hs_c(fm_hs_ni) = 0.
            fm_hs_i = fm_hs_i+1.
        END DO
        success = .true.

        return

    end subroutine FM_HS_output

    !---------------------------------------------------------------------------!
    function fm_hs_rhofunction(t,mgamma2,Neff,Omeganuh2)

        implicit none
        real(dl)::t,mgamma2,Neff,Omeganuh2
        real(dl) :: fm_hs_rhofunction

        IF (Omeganuh2==0.) then
            fm_hs_rhofunction = 0.
        else
            fm_hs_rhofunction = (0.6813*mgamma2*Neff*(1 + 5.434509745635465e8*(EXP(t)*Omeganuh2)**1.83)**0.5464480874316939)/EXP(4*t)
        end if

    end function fm_hs_rhofunction

    !---------------------------------------------------------------------------!
    function fm_hs_rhoprime(t,mgamma2,Neff,Omeganuh2)

        implicit none
        real(dl)::t,mgamma2,Neff,Omeganuh2
        real(dl) :: fm_hs_rhoprime

        if (Omeganuh2==0.) then
            fm_hs_rhoprime = 0.
        else
            fm_hs_rhoprime = (3.7025315D8*mgamma2*Neff*Omeganuh2*(EXP(t)*Omeganuh2)**83D-2)/&
                (EXP(3.*t)*(1.+ 5.43451D8*(EXP(t)*Omeganuh2)**1.83)**453552D-6)&
                - (2.7252*mgamma2*Neff*EXP(-4.*t))*((1. + 5.43451D8*(EXP(t)*Omeganuh2)**1.83)**(1./1.83))
        end if

    end function fm_hs_rhoprime

    !---------------------------------------------------------------------------!
    function fm_hs_derivativeForV(t,yh,v,m2,c1,c2,n)

        implicit none
        real(dl)::t,yh,v,yhp,m2,c1,c2,n
        real(dl)::firstpart,secondpart,thirdpart,fourthpart
        real(dl) :: fm_hs_derivativeForV

        firstpart = (c1*n*(1/(3/EXP(3*t) + 3*v + 12*yH))**(1 + n)*(-1/(2.*EXP(3*t)) - yH + (3*v + 12*yH)/6.))/c2**2
        secondpart = (-((c1*m2)/c2) + (c1*m2*(1/(3/EXP(3*t) + 3*v + 12*yH))**n)/c2**2)/(6.*m2)
        thirdpart = (2*c1*n*(1/(3/EXP(3*t) + 3*v + 12*yH))**(2 + n))/(c2**2*m2)
        fourthpart = (c1*(-1 + n)*n*(1/(3/EXP(3*t) + 3*v + 12*yH))**(2 + n))/(c2**2*m2)
        fm_hs_derivativeForV =  -4*v+(1./3.)*(9*EXP(-3*t) - (yH+firstpart + secondpart)/(m2*(EXP(-3*t)+yH)*(thirdpart+fourthpart)))

    end function fm_hs_derivativeForV

    !-------------------------------------------------------------------------------------------!
    function fm_hs_sourceFunction(m2,c1,c2,n,t)

        implicit none
        real(dl),intent(in)::m2,c1,c2,n,t
        real(dl) :: fm_hs_sourceFunction

        fm_hs_sourceFunction =3 + (3*c2)/((6*c2 + c1*EXP(3*t))*(1 + n)) +&
            (c1*EXP(3*t))/((6*c2 + c1*EXP(3*t))*(1 + n)) - &
            (2*c1**2*EXP(6*t))/(3.*c2*(6*c2 + c1*EXP(3*t))*(1 + n)) - &
            (3*c2)/((6*c2 + c1*EXP(3*t))*n*(1 + n)) - &
            (4*c1*EXP(3*t))/((6*c2 + c1*EXP(3*t))*n*(1 + n)) - &
            (4*c1**2*EXP(6*t))/(3*c2*(6*c2 + c1*EXP(3*t))*n*(1 + n))

    end function fm_hs_sourceFunction

    !-------------------------------------------------------------------------------------------!
    function fm_hs_massFunction(m2,c1,c2,n,t)

        implicit none
        real(dl),intent(in) :: m2,c1,c2,n,t
        real(dl) :: fm_hs_massFunction

        fm_hs_massFunction =-3 - (18*c2**2)/((6*c2 + c1*EXP(3*t))**2*(1 + n)) -&
            (6*c1*c2*EXP(3*t))/((6*c2 + c1*EXP(3*t))**2*(1 + n)) - &
            (5*c1**2*EXP(6*t))/((6*c2 + c1*EXP(3*t))**2*(1 + n)) + &
            (18*c2**2)/((6*c2 + c1*EXP(3*t))**2*n*(1 + n)) + &
            (6*c1*c2*EXP(3*t))/((6*c2 + c1*EXP(3*t))**2*n*(1 + n)) - &
            (4*c1**2*EXP(6*t))/((6*c2 + c1*EXP(3*t))**2*n*(1 + n)) + &
            (162*c2**3)/&
            ((1/((2*c1)/c2 + 3/EXP(3*t)))**n*(6*c2 + c1*EXP(3*t))**2*n*&
            (1 + n)) + (108*c2**4)/&
            (c1*EXP(3*t)*(1/((2*c1)/c2 + 3/EXP(3*t)))**n*&
            (6*c2 + c1*EXP(3*t))**2*n*(1 + n)) + &
            (72*c1*c2**2*EXP(3*t))/&
            ((1/((2*c1)/c2 + 3/EXP(3*t)))**n*(6*c2 + c1*EXP(3*t))**2*n*&
            (1 + n)) + (8*c1**2*c2*EXP(6*t))/&
            ((1/((2*c1)/c2 + 3/EXP(3*t)))**n*(6*c2 + c1*EXP(3*t))**2*n*&
            (1 + n))

    end function fm_hs_massFunction

    !-------------------------------------------------------------------------------------------!
    function fm_hs_fFunction(m2,c1,c2,n,t)

        implicit none
        real(dl),intent(in) :: m2,c1,c2,n,t
        real(dl) :: fm_hs_fFunction

        fm_hs_fFunction = (EXP(-3.*t)*(3.*c1*n*EXP(3.*t)+4.*c1*EXP(3.*t)+6.*c2))/(n*(n+1.)*(c1*EXP(3.*t)+6.*c2))
        fm_hs_fFunction = 2.*EXP(-3.*t) - fm_hs_fFunction
        fm_hs_fFunction = -EXP(-3.*t)*fm_hs_fFunction

    end function fm_hs_fFunction

    !-------------------------------------------------------------------------------------------!
    function fm_hs_est_z_ast(success)

        implicit none
        real(dl) :: fm_hs_est_z_ast
        integer fm_hs_ni
        real(dl) fm_hs_i,fm_hs_t,fm_hs_ratio
        logical, intent(inout) :: success

        if (ABS(fm_hs_fro) .gt. 1.d-2) then
            fm_hs_i = 1.
            DO fm_hs_ni = 1,fm_hs_nnosteps,1
                fm_hs_t = fm_hs_llim+fm_hs_step*(fm_hs_i-1.)
                fm_hs_yH(fm_hs_ni) = fm_hs_sourceFunction(fm_hs_m2,fm_hs_c1,fm_hs_c2,fm_hs_n,fm_hs_t)/fm_hs_massFunction(fm_hs_m2,fm_hs_c1,fm_hs_c2,fm_hs_n,fm_hs_t)
                fm_hs_v(fm_hs_ni) = fm_hs_firstderive(fm_hs_m2,fm_hs_c1,fm_hs_c2,fm_hs_n,fm_hs_t)
                if ((fm_hs_yH(fm_hs_ni) .ne. 0.0) .and. (fm_hs_massFunction(fm_hs_m2,fm_hs_c1,fm_hs_c2,fm_hs_n,fm_hs_t) .ne. 0.0)) then
                    !fm_hs_ratio = fm_hs_fFunction(fm_hs_m2,fm_hs_c1,fm_hs_c2,fm_hs_n,fm_hs_t)*fm_hs_v(fm_hs_ni)
                    fm_hs_ratio = fm_hs_v(fm_hs_ni)
                    fm_hs_ratio = fm_hs_ratio/fm_hs_massFunction(fm_hs_m2,fm_hs_c1,fm_hs_c2,fm_hs_n,fm_hs_t)/fm_hs_yH(fm_hs_ni)
                    fm_hs_ratio = ABS(fm_hs_ratio)
                endif
                if ( fm_hs_ratio .gt. fm_hs_z_ast_threshold) then
                    fm_hs_est_z_ast = fm_hs_t
                    return
                endif
                fm_hs_i = fm_hs_i+1.
            END DO
            fm_hs_est_z_ast = fm_xFinal !bh
        else
            fm_hs_est_z_ast = fm_xFinal !bh
        endif
        success = .true.

    end function fm_hs_est_z_ast

    !-------------------------------------------------------------------------------------------!
    function fm_hs_firstderive(m2,c1,c2,n,t)

        implicit none
        real(dl),intent(in)::m2,c1,c2,n,t
        real(dl) :: fm_hs_firstderive

        fm_hs_firstderive = -(-27/EXP(3*t) - (6*c1*(3*c2 + 2*c1*EXP(3*t) - 3*c2*n + c1*EXP(3*t)*n))/(c2*(6*c2 + c1*EXP(3*t))*n*(1 + n)) +&
            (3*c1*(3*c2 + 2*c1*EXP(3*t))*(3*c2 + 2*c1*EXP(3*t) - 3*c2*n + c1*EXP(3*t)*n))/&
            (c2*(6*c2 + c1*EXP(3*t))**2*n*(1 + n)) +&
            (3*(3*c2 + 2*c1*EXP(3*t))*(3*c2 + 2*c1*EXP(3*t) - 3*c2*n + c1*EXP(3*t)*n))/&
            (c2*EXP(3*t)*(6*c2 + c1*EXP(3*t))*n*(1 + n)) -&
            ((3*c2 + 2*c1*EXP(3*t))*(6*c1*EXP(3*t) + 3*c1*EXP(3*t)*n))/(c2*EXP(3*t)*(6*c2 + c1*EXP(3*t))*n*(1 + n)))/&
            (3*(3/EXP(3*t) + (-108*c2**4 - 162*c1*c2**3*EXP(3*t) - 72*c1**2*c2**2*EXP(6*t) - 8*c1**3*c2*EXP(9*t) -&
            18*c1*c2**2*EXP(3*t)*(1/((2*c1)/c2 + 3/EXP(3*t)))**n - 6*c1**2*c2*EXP(6*t)*(1/((2*c1)/c2 + 3/EXP(3*t)))**n +&
            4*c1**3*EXP(9*t)*(1/((2*c1)/c2 + 3/EXP(3*t)))**n + 18*c1*c2**2*EXP(3*t)*(1/((2*c1)/c2 + 3/EXP(3*t)))**n*n +&
            6*c1**2*c2*EXP(6*t)*(1/((2*c1)/c2 + 3/EXP(3*t)))**n*n + 5*c1**3*EXP(9*t)*(1/((2*c1)/c2 + 3/EXP(3*t)))**n*n)/&
            (c1*EXP(6*t)*(1/((2*c1)/c2 + 3/EXP(3*t)))**n*(6*c2 + c1*EXP(3*t))**2*n*(1 + n)))) +&
            ((9/EXP(3*t) - ((3*c2 + 2*c1*EXP(3*t))*(3*c2 + 2*c1*EXP(3*t) - 3*c2*n + c1*EXP(3*t)*n))/&
            (c2*EXP(3*t)*(6*c2 + c1*EXP(3*t))*n*(1 + n)))*&
            (-9/EXP(3*t) - (9*(1/((2*c1)/c2 + 3/EXP(3*t)))**(1 - n)*&
            (-108*c2**4 - 162*c1*c2**3*EXP(3*t) - 72*c1**2*c2**2*EXP(6*t) - 8*c1**3*c2*EXP(9*t) -&
            18*c1*c2**2*EXP(3*t)*(1/((2*c1)/c2 + 3/EXP(3*t)))**n -&
            6*c1**2*c2*EXP(6*t)*(1/((2*c1)/c2 + 3/EXP(3*t)))**n + 4*c1**3*EXP(9*t)*(1/((2*c1)/c2 + 3/EXP(3*t)))**n +&
            18*c1*c2**2*EXP(3*t)*(1/((2*c1)/c2 + 3/EXP(3*t)))**n*n +&
            6*c1**2*c2*EXP(6*t)*(1/((2*c1)/c2 + 3/EXP(3*t)))**n*n + 5*c1**3*EXP(9*t)*(1/((2*c1)/c2 + 3/EXP(3*t)))**n*n))&
            /(c1*EXP(9*t)*(6*c2 + c1*EXP(3*t))**2*(1 + n)) -&
            (6*(-108*c2**4 - 162*c1*c2**3*EXP(3*t) - 72*c1**2*c2**2*EXP(6*t) - 8*c1**3*c2*EXP(9*t) -&
            18*c1*c2**2*EXP(3*t)*(1/((2*c1)/c2 + 3/EXP(3*t)))**n -&
            6*c1**2*c2*EXP(6*t)*(1/((2*c1)/c2 + 3/EXP(3*t)))**n + 4*c1**3*EXP(9*t)*(1/((2*c1)/c2 + 3/EXP(3*t)))**n +&
            18*c1*c2**2*EXP(3*t)*(1/((2*c1)/c2 + 3/EXP(3*t)))**n*n +&
            6*c1**2*c2*EXP(6*t)*(1/((2*c1)/c2 + 3/EXP(3*t)))**n*n + 5*c1**3*EXP(9*t)*(1/((2*c1)/c2 + 3/EXP(3*t)))**n*n))&
            /(EXP(3*t)*(1/((2*c1)/c2 + 3/EXP(3*t)))**n*(6*c2 + c1*EXP(3*t))**3*n*(1 + n)) -&
            (6*(-108*c2**4 - 162*c1*c2**3*EXP(3*t) - 72*c1**2*c2**2*EXP(6*t) - 8*c1**3*c2*EXP(9*t) -&
            18*c1*c2**2*EXP(3*t)*(1/((2*c1)/c2 + 3/EXP(3*t)))**n -&
            6*c1**2*c2*EXP(6*t)*(1/((2*c1)/c2 + 3/EXP(3*t)))**n + 4*c1**3*EXP(9*t)*(1/((2*c1)/c2 + 3/EXP(3*t)))**n +&
            18*c1*c2**2*EXP(3*t)*(1/((2*c1)/c2 + 3/EXP(3*t)))**n*n +&
            6*c1**2*c2*EXP(6*t)*(1/((2*c1)/c2 + 3/EXP(3*t)))**n*n + 5*c1**3*EXP(9*t)*(1/((2*c1)/c2 + 3/EXP(3*t)))**n*n))&
            /(c1*EXP(6*t)*(1/((2*c1)/c2 + 3/EXP(3*t)))**n*(6*c2 + c1*EXP(3*t))**2*n*(1 + n)) +&
            (-486*c1*c2**3*EXP(3*t) - 432*c1**2*c2**2*EXP(6*t) - 72*c1**3*c2*EXP(9*t) -&
            54*c1*c2**2*EXP(3*t)*(1/((2*c1)/c2 + 3/EXP(3*t)))**n - 36*c1**2*c2*EXP(6*t)*(1/((2*c1)/c2 + 3/EXP(3*t)))**n +&
            36*c1**3*EXP(9*t)*(1/((2*c1)/c2 + 3/EXP(3*t)))**n + 54*c1*c2**2*EXP(3*t)*(1/((2*c1)/c2 + 3/EXP(3*t)))**n*n +&
            36*c1**2*c2*EXP(6*t)*(1/((2*c1)/c2 + 3/EXP(3*t)))**n*n +&
            45*c1**3*EXP(9*t)*(1/((2*c1)/c2 + 3/EXP(3*t)))**n*n - 162*c1*c2**2*(1/((2*c1)/c2 + 3/EXP(3*t)))**(1 + n)*n -&
            54*c1**2*c2*EXP(3*t)*(1/((2*c1)/c2 + 3/EXP(3*t)))**(1 + n)*n +&
            36*c1**3*EXP(6*t)*(1/((2*c1)/c2 + 3/EXP(3*t)))**(1 + n)*n +&
            162*c1*c2**2*(1/((2*c1)/c2 + 3/EXP(3*t)))**(1 + n)*n**2 +&
            54*c1**2*c2*EXP(3*t)*(1/((2*c1)/c2 + 3/EXP(3*t)))**(1 + n)*n**2 +&
            45*c1**3*EXP(6*t)*(1/((2*c1)/c2 + 3/EXP(3*t)))**(1 + n)*n**2)/&
            (c1*EXP(6*t)*(1/((2*c1)/c2 + 3/EXP(3*t)))**n*(6*c2 + c1*EXP(3*t))**2*n*(1 + n))))/&
            (3*(3/EXP(3*t) + (-108*c2**4 - 162*c1*c2**3*EXP(3*t) - 72*c1**2*c2**2*EXP(6*t) - 8*c1**3*c2*EXP(9*t) -&
            18*c1*c2**2*EXP(3*t)*(1/((2*c1)/c2 + 3/EXP(3*t)))**n - 6*c1**2*c2*EXP(6*t)*(1/((2*c1)/c2 + 3/EXP(3*t)))**n +&
            4*c1**3*EXP(9*t)*(1/((2*c1)/c2 + 3/EXP(3*t)))**n + 18*c1*c2**2*EXP(3*t)*(1/((2*c1)/c2 + 3/EXP(3*t)))**n*n +&
            6*c1**2*c2*EXP(6*t)*(1/((2*c1)/c2 + 3/EXP(3*t)))**n*n + 5*c1**3*EXP(9*t)*(1/((2*c1)/c2 + 3/EXP(3*t)))**n*n)/&
            (c1*EXP(6*t)*(1/((2*c1)/c2 + 3/EXP(3*t)))**n*(6*c2 + c1*EXP(3*t))**2*n*(1 + n)))**2)

    end function fm_hs_firstderive

end module EFT_FullMapping_HS
!--------------------------------------------------------------------------------------------!
