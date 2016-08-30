module EFT_HS_fR

    use EFTDef
    use EFTdeEOS
    use MassiveNu

    ! 1) Definitions of the variables common to all the designer code.
    integer , parameter :: MG_nstep = 1000      ! Number of points sampled by the designer code.
    integer , parameter :: MG_ninterpol = 2       ! Number of points used by the interpolator. Only even.
    real(dl), parameter :: xInitial = log(10._dl**(-8._dl))    ! log(a start)
    real(dl), parameter :: xFinal   = 0._dl                    ! Log(a final)

    ! 2) Tables containing the sampled values of the EFT functions.
    real(dl), save, dimension( MG_nstep ) :: xp_MG                      ! Array with the time sampling
    real(dl), save, dimension( MG_nstep ) :: EFT_Omega_MG, EFT_OmegaP_MG, EFT_OmegaPP_MG, EFT_OmegaPPP_MG      ! Array with the functions
    real(dl), save, dimension( MG_nstep ) :: EFT_c_MG, EFT_cdot_MG      ! Array with the functions
    real(dl), save, dimension( MG_nstep ) :: EFT_Lambda_MG, EFT_Lambdadot_MG      ! Array with the functions

    real(dl), save, dimension( MG_nstep ) :: MG_w_DE      ! Array with derived functions

    logical :: MGPrintFlag = .False.

contains

    !***************************************************
    ! Function f(R)

    function EFT_foR( R, deriv )
        ! f(R) function
        implicit none

        real(dl) :: EFT_foR
        real(dl) , intent(in) :: R      ! the Ricci scalar
        integer  , intent(in) :: deriv  ! order of derivative wrt the Ricci Scalar

        if (CP%EFTflag /= 5) stop 'Calling f(R) when EFTflag is not 5'

        select case (CP%FullMappingfoRmodel)
            case (0) ! 0) Operator neglected f(R)=0
                EFT_foR = 0._dl
            case (1) ! 1) HS model:
                if (deriv==0) then
                    EFT_foR = CP%EFT_foR_mass**2*( R/CP%EFT_foR_mass**2 -CP%EFT_foR_c1*(R/CP%EFT_foR_mass**2)**CP%EFT_foR_n/(1._dl +CP%EFT_foR_c2*(R/CP%EFT_foR_mass**2)**CP%EFT_foR_n) )
                else if (deriv==1) then
                    EFT_foR = (-(CP%EFT_foR_c1*CP%EFT_foR_mass**2*CP%EFT_foR_n*&
                            &(R/CP%EFT_foR_mass**2)**CP%EFT_foR_n) + &
                            &R*(1._dl + CP%EFT_foR_c2*(R/CP%EFT_foR_mass**2)**CP%EFT_foR_n)**2)/&
                            &(R*(1._dl + CP%EFT_foR_c2*(R/CP%EFT_foR_mass**2)**CP%EFT_foR_n)**2)
                else if (deriv==2) then
                    EFT_foR = (CP%EFT_foR_c1*CP%EFT_foR_mass**2*CP%EFT_foR_n*&
                            &(R/CP%EFT_foR_mass**2)**CP%EFT_foR_n*&
                            &(1._dl - CP%EFT_foR_n + CP%EFT_foR_c2*(1._dl + CP%EFT_foR_n)*&
                            &(R/CP%EFT_foR_mass**2)**CP%EFT_foR_n))/&
                            &(R**2*(1._dl + CP%EFT_foR_c2*(R/CP%EFT_foR_mass**2)**CP%EFT_foR_n)**3)
                end if
            case default
                stop 'Wrong model for f(R)'
        end select

    end function EFT_foR

    !***************************************************
    ! Derivs of the ODE system for the f(R) model

    subroutine EFT_foR_derivs( num_eq, x, y, ydot )

        implicit none

        integer  :: num_eq
        real(dl) :: x
        real(dl), dimension(num_eq) :: y
        real(dl), dimension(num_eq) :: ydot

        real(dl) :: coeff_1, coeff_2, coeff_3
        real(dl) :: a, a2, gpres, EFT_H0
        real(dl) :: EFT_grhonu_tot, EFT_gpinu_tot, EFT_grhonu, EFT_gpinu, grhormass_t
        real(dl) :: grhor_t, grhog_t

        real(dl) :: grhob_t, grhoc_t, grhov_t, grho, adotoa, adotdota, Hdot !DB
        real(dl) :: EFT_c, EFT_Lambda, adotoa_rec, grho_matter

        integer  :: nu_i

        ! initial computations:
        a  = exp(x)
        a2 = a*a

        EFT_H0 = (CP%h0/c_EFT)*1000._dl

        ! compute total pressure
        ! massive neutrinos:
        EFT_gpinu_tot     = 0._dl
        EFT_grhonu_tot    = 0._dl
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
        ! other species:
        grhor_t = grhornomass/a2
        grhog_t = grhog/a2
        ! collect all pressure terms:
        gpres= EFT_gpinu_tot + (grhog_t+grhor_t)/3._dl

        coeff_1 = 1._dl + EFT_OL_Omega(a,0) +0.5_dl*a*EFT_OL_Omega(a,1)
        coeff_2 = 1._dl + EFT_OL_Omega(a,0) +2._dl*a*EFT_OL_Omega(a,1) +a2*EFT_OL_Omega(a,2)
        coeff_3 = gpres -3._dl*EFT_H0**2*CP%omegav*( 1._dl+EFT_OL_Lambda(a,0) )*a2

        ydot(1) = 1._dl/(coeff_1)*( -coeff_2*y(1) -coeff_3 )

    end subroutine EFT_foR_derivs

    !***************************************************
    ! Takes the solution of the ODE and computes the other quantities that we need:

    subroutine EFT_foR_output( num_eq, x, y,  EFT_c, EFT_cdot, w_DE )

        implicit none

        integer  :: num_eq
        real(dl) :: x
        real(dl), dimension(num_eq) :: y
        real(dl), intent(out) :: EFT_c, EFT_cdot, w_DE

        real(dl), dimension(num_eq) :: ydot
        real(dl) :: coeff_1, coeff_2, coeff_3
        real(dl) :: a, a2, gpres, EFT_H0
        real(dl) :: EFT_grhonu_tot, EFT_gpinu_tot, EFT_grhonu, EFT_gpinu, grhormass_t
        real(dl) :: grhor_t, grhog_t

        real(dl) :: grhob_t, grhoc_t, grhov_t, grho, adotoa, adotdota, Hdot !DB
        real(dl) :: grho_matter

        integer  :: nu_i

        call Omega_Lambda_derivs( num_eq, x, y, ydot )

        ! initial computations:
        a  = exp(x)
        a2 = a*a

        EFT_H0 = (CP%h0/c_EFT)*1000._dl

        ! compute total pressure and density

        grhob_t=grhob/a     !DB
        grhoc_t=grhoc/a     !DB
        grhov_t=grhov*a2    !DB

        ! massive neutrinos:
        EFT_gpinu_tot     = 0._dl
        EFT_grhonu_tot    = 0._dl
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
        ! other species:
        grhor_t = grhornomass/a2
        grhog_t = grhog/a2
        ! collect all pressure terms:

        gpres= EFT_gpinu_tot + (grhog_t+grhor_t)/3._dl
        grho_matter = grhob_t +grhoc_t +grhor_t +grhog_t +EFT_grhonu_tot

        ! compute the output:

        adotoa = sqrt(y(1))
        Hdot   = 0.5_dl*ydot(1)

        EFT_c      = +1.5_dl*(1._dl +EFT_OL_Omega(a,0) +a*EFT_OL_Omega(a,1) )*adotoa**2 -0.5_dl*(grho_matter) -1.5_dl*EFT_H0**2*CP%omegav*( 1._dl+EFT_OL_Lambda(a,0) )*a2
        EFT_cdot   = +3._dl*adotoa*( &
            & -( 1._dl +EFT_OL_Omega(a,0) -0.5_dl*a2*EFT_OL_Omega(a,2) )*adotoa**2 &
            & +( 1._dl +EFT_OL_Omega(a,0) +a*EFT_OL_Omega(a,1) )*Hdot &
            & +0.5_dl*( grho_matter +gpres ) &
            & -0.5_dl*a**3*EFT_H0**2*CP%omegav**EFT_OL_Lambda(a,1) )

        w_DE       = ( -2._dl*Hdot -adotoa**2 -gpres )/( 3._dl*adotoa**2 -grho_matter )

    end subroutine EFT_foR_output

    !***************************************************
    ! Jacobian for the ODE system of the f(R) model.

    subroutine EFT_foR_Jac( num_eq, x, y, ml, mu, pd, nrowpd )

        implicit none

        integer  :: num_eq
        integer  :: ml ! ignored
        integer  :: mu ! ignored
        integer  :: nrowpd
        real(dl) :: x
        real(dl), dimension(num_eq) :: y
        real(dl), dimension(nrowpd,num_eq) :: pd

        real(dl) :: coeff_1, coeff_2, coeff_3
        real(dl) :: a, a2, gpres, EFT_H0
        real(dl) :: EFT_grhonu_tot, EFT_gpinu_tot, EFT_grhonu, EFT_gpinu, grhormass_t
        real(dl) :: grhor_t, grhog_t

        real(dl) :: grhob_t, grhoc_t, grhov_t, grho, adotoa, adotdota, Hdot !DB
        real(dl) :: EFT_c, EFT_Lambda, adotoa_rec, grho_matter

        integer  :: nu_i

        ! initial computations:
        a  = exp(x)
        a2 = a*a

        coeff_1 = 1._dl + EFT_OL_Omega(a,0) +0.5_dl*a*EFT_OL_Omega(a,1)
        coeff_2 = 1._dl + EFT_OL_Omega(a,0) +2._dl*a*EFT_OL_Omega(a,1) +a2*EFT_OL_Omega(a,2)

        pd(1,1) = -1._dl/(coeff_1)*( coeff_2 )

    end subroutine EFT_foR_Jac

    !*******************************************
    ! This function takes a table of values and returns the interpolation of them

    function Interpolate_EFT_Function(EFTFunctionTable, a, GR_value)
        ! The designer code will provide a table of sampled values for the EFT functions.
        ! This function is called in the EFTfunctions module to interpolate those tables.
        implicit none

        real(dl), intent(in) :: a, GR_value
        real(dl) :: EFTFunctionTable(MG_nstep)
        real(dl) :: Interpolate_EFT_Function

        real(dl) :: x, temp, dtemp, width
        real(dl) :: xb(MG_ninterpol),yb(MG_ninterpol)
        integer  :: i, jlo, stint

        x = log(a)

        if (x.lt.xp_MG(1)) then
            temp = GR_value
        else if (x.ge.xp_MG(1) .and. x.le.xp_MG(MG_nstep)) then
            ! 1) find the point corresponing to the requested input.
            !call hunt(xp_MG, MG_nstep, x, jlo)
            ! we use the equispaced hash function:
            width = (xFinal-xInitial)/REAL(MG_nstep-1)
            jlo = int( ( x-xp_MG(1) - Mod(x-xp_MG(1),width))/(width) ) +1
            ! 2) construct the table that will be interpolated
            stint = MG_ninterpol/2

            if((MG_nstep-jlo-stint).lt.0) then
                ! Last MG_ninterpol points
                do i=1, MG_ninterpol
                    xb(i)=xp_MG(MG_nstep-MG_ninterpol+i)
                    yb(i)=EFTFunctionTable(MG_nstep-MG_ninterpol+i)
                end do
            else if (jlo.eq.1) then
                ! First MG_ninterpol points
                do i=1, MG_ninterpol
                    xb(i)=xp_MG(i)
                    yb(i)=EFTFunctionTable(i)
                end do
            else
                ! Surrounding MG_ninterpol points
                do i=1, MG_ninterpol
                    xb(i)=xp_MG(jlo-stint+i)
                    yb(i)=EFTFunctionTable(jlo-stint+i)
                end do
            endif
            ! 3) Call the Neville interpolator.
            call Polint(MG_ninterpol, xb, yb, x, temp, dtemp)

        else
            temp = EFTFunctionTable(MG_nstep)
        end if

        Interpolate_EFT_Function = temp

    end function Interpolate_EFT_Function

    !****************************************
    ! Does initialization of the model

    subroutine Initialize_FullMapping_foR(success)
        ! This subroutine calls the correct designer subroutine depending on the selected model.
        implicit none

        logical, intent(inout) :: success

        real(dl) :: y(1)
        integer  :: i, num_points
        real(dl) :: t0, tfin, t1, t2
        real(dl), allocatable :: RWORK(:)
        integer , allocatable :: IWORK(:)

        real(dl) :: rtol, atol, EFT_H0
        integer  :: LRN, LRS, LRW, LIS, LIN, LIW, NEQ
        integer  :: itol, itask, istate, iopt, JacobianMode

        real(dl) :: EFT_c, EFT_cdot, w_DE

        success = .True.

        if (Feedbacklevel>1) write(*,*) 'EFTCAMB: initializing full mapping f(R)'

        ! initial computations:

        EFT_H0 = (CP%h0/c_EFT)*1000._dl

        EFT_c_MG    = 0._dl
        EFT_cdot_MG = 0._dl
        OL_w_DE     = 0._dl

        ! dimension of the system:
        neq = 1
        ! set the initial conditions:
        y = EFT_H0**2
        ! set time boundary:
        t0   = xFinal
        tfin = xInitial
        ! set number of points:
        num_points = MG_nstep
        ! set-up the relative and absolute tollerances:
        itol = 1
        rtol = 1.d-10
        atol = 1.d-16
        ! initialize task to do:
        itask  = 1
        istate = 1
        iopt   = 0
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
        ! Jacobian mode: 1=fullJacobian, 2=not provided
        JacobianMode = 1

        ! solve the system and store the solution:

        ! compute the output for the first point (notice that the system is solved in the past):
        call Omega_Lambda_output( neq, t0, y,  EFT_c, EFT_cdot, w_DE )
        !
        xp_MG(MG_nstep)       = t0
        EFT_c_MG(MG_nstep)    = EFT_c
        EFT_cdot_MG(MG_nstep) = EFT_cdot
        OL_w_DE(MG_nstep)     = w_DE

        do i=1, num_points-1

            t1 = t0 + REAL(i-1)/REAL(num_points-1)*(tfin-t0)
            t2 = t0 + REAL(i)/REAL(num_points-1)*(tfin-t0)

            ! solve the system:
            call DLSODA ( Omega_Lambda_derivs, neq, y, t1, t2, itol, rtol, atol, itask, istate, iopt, RWORK, LRW, IWORK, LIW, Omega_Lambda_Jac, JacobianMode)
            ! check the quality of the output:
            if ( istate<0 ) then
                success = .False.
                if (Feedbacklevel>=1) write(*,*) 'Problems with initializing Lambda Omega'
                return
            end if

            ! compute the desired output:
            call Omega_Lambda_output( neq, t2, y,  EFT_c, EFT_cdot, w_DE )
            ! store the output:
            xp_MG(MG_nstep-i)       = t2
            EFT_c_MG(MG_nstep-i)    = EFT_c
            EFT_cdot_MG(MG_nstep-i) = EFT_cdot
            OL_w_DE(MG_nstep-i)     = w_DE

        end do

        return
    end subroutine Initialize_Omega_Lambda

    !*****************************
    ! computes c

    function Omega_Lambda_c( a, deriv )

        implicit none

        real(dl), intent(in) :: a
        integer , intent(in) :: deriv

        real(dl) :: Omega_Lambda_c

        if ( deriv==0 ) then
            Omega_Lambda_c = Interpolate_EFT_Function(EFT_c_MG, a, 0._dl)
        else if (deriv==1 ) then
            Omega_Lambda_c = Interpolate_EFT_Function(EFT_cdot_MG, a, 0._dl)
        else
            stop 'Wrong deriv in call to Omega_Lambda_c'
        end if

    end function Omega_Lambda_c

    !*****************************
    ! computes w DE

    function Omega_Lambda_wDE( a )

        implicit none

        real(dl), intent(in) :: a
        real(dl) :: Omega_Lambda_wDE

        Omega_Lambda_wDE = Interpolate_EFT_Function(OL_w_DE, a, -1._dl)

    end function Omega_Lambda_wDE


end module EFT_HS_fR
