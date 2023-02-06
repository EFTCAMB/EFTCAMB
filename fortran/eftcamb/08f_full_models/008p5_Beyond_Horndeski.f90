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

!> @file 10p5_Beyond_Horndeski.f90
!! This file contains the definition of the Beyond Horndeski model consistent with GW170817.
!! Please refer to the numerical notes for details.


!----------------------------------------------------------------------------------------
!> This module contains the definition of the Beyond Horndeski model consistent with GW170817.
!! Please refer to the numerical notes for details.

!> @author Simone Peirone, Bin Hu, Marco Raveri

module EFTCAMB_full_Beyond_Horndeski

    use precision
    use IniObjects
    use MpiUtils
    use FileUtils
    use EFTCAMB_cache
    use EFT_def
    use EFTCAMB_abstract_model_full
    use EFTCAMB_abstract_parametrizations_1D
    use EFTCAMB_constant_parametrization_1D
    use equispaced_linear_interpolation_1D

    implicit none

    private

    public EFTCAMB_Beyond_Horndeski

    !----------------------------------------------------------------------------------------
    !> This is the type that contains the definition of the Beyond Horndeski model.
    type, extends ( EFTCAMB_full_model ) :: EFTCAMB_Beyond_Horndeski

        ! the model parameters:
        real(dl)  :: bHor_x10
        real(dl)  :: bHor_x20
        real(dl)  :: bHor_x30
        real(dl)  :: bHor_x40
        real(dl)  :: bHor_x1today
        real(dl)  :: bHor_x2today
        real(dl)  :: bHor_x3today
        real(dl)  :: bHor_x4today

        ! the interpolated EFT functions that come out of the background sover:
        type(equispaced_linear_interpolate_function_1D) :: fun_x1       !< The interpolated function x1
        type(equispaced_linear_interpolate_function_1D) :: fun_x2       !< The interpolated function x2
        type(equispaced_linear_interpolate_function_1D) :: fun_x3       !< The interpolated function x3
        type(equispaced_linear_interpolate_function_1D) :: fun_x4       !< The interpolated function x4

        ! some parameters for the ODE solver:
        integer  :: designer_num_points = 10000                            !< Number of points sampled by the designer code.
        real(dl) :: x_initial           = log(1._dl/(1.5*10._dl**(5._dl)+1._dl))                   !< log(a start)
        real(dl) :: x_final             = 0._dl         !< log(a final)

        ! some additional parameters
        logical :: debug_flag = .true.
        logical :: LCDM_background = .false.

    contains

        ! initialization of the model:
        procedure :: read_model_selection            => EFTCAMBBeyondHorndeskiReadModelSelectionFromFile  !< subroutine that reads the parameters of the model from file
        procedure :: allocate_model_selection        => EFTCAMBBeyondHorndeskiAllocateModelSelection      !< subroutine that allocates the model selection. For Horava this is a dummy procedure.
        procedure :: init_model_parameters           => EFTCAMBBeyondHorndeskiInitModelParameters         !< subroutine that initializes the model parameters based on the values found in an input array.
        procedure :: init_model_parameters_from_file => EFTCAMBBeyondHorndeskiInitModelParametersFromFile !< subroutine that reads the parameters of the model from file.

        ! background solver:
        procedure :: initialize_background           => EFTCAMBBeyondHorndeskiInitBackground              !< subroutine that initializes the background of Quintic Galileon.
        procedure :: solve_background                => EFTCAMBBeyondHorndeskiSolveBackground             !< subroutine that solves the background equations.

        ! utility functions:
        procedure :: compute_param_number  => EFTCAMBBeyondHorndeskiComputeParametersNumber    !< subroutine that computes the number of parameters of the model.
        procedure :: feedback              => EFTCAMBBeyondHorndeskiFeedback                   !< subroutine that prints on the screen feedback information about the model.
        procedure :: parameter_names       => EFTCAMBBeyondHorndeskiParameterNames             !< subroutine that returns the i-th parameter name of the model.
        procedure :: parameter_names_latex => EFTCAMBBeyondHorndeskiParameterNamesLatex        !< subroutine that returns the i-th parameter name of the model.
        procedure :: parameter_values      => EFTCAMBBeyondHorndeskiParameterValues            !< subroutine that returns the i-th parameter value.

        ! CAMB related procedures:
        procedure :: compute_background_EFT_functions  => EFTCAMBBeyondHorndeskiBackgroundEFTFunctions   !< subroutine that computes the value of the background EFT functions at a given time.
        procedure :: compute_secondorder_EFT_functions => EFTCAMBBeyondHorndeskiSecondOrderEFTFunctions  !< subroutine that computes the value of the second order EFT functions at a given time.

        ! stability procedures:
        procedure :: additional_model_stability        => EFTCAMBBeyondHorndeskiAdditionalModelStability !< function that computes model specific stability requirements.

    end type EFTCAMB_Beyond_Horndeski

    ! ---------------------------------------------------------------------------------------------

contains

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that reads the parameters of the model from file.
    subroutine EFTCAMBBeyondHorndeskiReadModelSelectionFromFile( self, Ini, eft_error )

        implicit none

        class(EFTCAMB_Beyond_Horndeski)         :: self   !< the base class
        type(TIniFile)                          :: Ini    !< Input ini file
        integer                                 :: eft_error     !< error code: 0 all fine, 1 initialization

        ! Nothing needs to be done but procedure present because its definition is deferred.

    end subroutine EFTCAMBBeyondHorndeskiReadModelSelectionFromFile

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that allocates the model selection.
    subroutine EFTCAMBBeyondHorndeskiAllocateModelSelection( self, Ini, eft_error)

        implicit none

        class(EFTCAMB_Beyond_Horndeski)  :: self         !< the base class
        type(TIniFile)                   :: Ini          !< Input ini file
        integer                          :: eft_error    !< error code: 0 all fine, 1 initialization failed

        ! Nothing needs to be done but procedure present because its definition is deferred.

    end subroutine EFTCAMBBeyondHorndeskiAllocateModelSelection

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that initializes the model parameters based on the values found in an input array.

    subroutine EFTCAMBBeyondHorndeskiInitModelParameters( self, array )

        implicit none

        class(EFTCAMB_Beyond_Horndeski)                            :: self   !< the base class
        real(dl), dimension(self%parameter_number), intent(in)     :: array  !< input array with the values of the parameters of the model.

        self%bHor_x10                = array(1)
        self%bHor_x30                = array(2)
        self%bHor_x40                = array(3)


    end subroutine EFTCAMBBeyondHorndeskiInitModelParameters

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that reads the parameters of the model from file.
    subroutine EFTCAMBBeyondHorndeskiInitModelParametersFromFile( self, Ini, eft_error )

        implicit none

        class(EFTCAMB_Beyond_Horndeski)  :: self         !< the base class
        type(TIniFile)                   :: Ini          !< Input ini file
        integer                          :: eft_error    !< error code: 0 all fine, 1 initialization failed

        self%bHor_x10                = Ini%Read_Double('Beyond_Horndeski_x10', 0._dl )
        self%bHor_x30                = Ini%Read_Double('Beyond_Horndeski_x30', 0._dl )
        self%bHor_x40                = Ini%Read_Double('Beyond_Horndeski_x40', 0._dl )
        self%debug_flag              = Ini%Read_logical( 'want_debug', .false. )
        self%LCDM_background         = Ini%Read_logical( 'want_LCDM',.false. )

    end subroutine EFTCAMBBeyondHorndeskiInitModelParametersFromFile

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that initializes the background of Quintic Galileon.
    subroutine EFTCAMBBeyondHorndeskiInitBackground( self, params_cache, feedback_level, success, outroot )

        implicit none

        class(EFTCAMB_Beyond_Horndeski)               :: self            !< the base class
        type(TEFTCAMB_parameter_cache), intent(inout) :: params_cache    !< a EFTCAMB parameter cache containing cosmological parameters
        integer                       , intent(in)    :: feedback_level  !< level of feedback from the background code. 0=none; 1=some; 2=chatty.
        logical                       , intent(out)   :: success         !< wether the background initialization succeded or not
        character(LEN=*), optional    , intent(in)    :: outroot         !< the output root for the debug files

        real(dl)      :: OmegaDE_0,OmegaDE, OmegaM_0, OmegaRad_0,x2_i, x2_il, x2_ir, tol,tol_r,tol_l
        real(dl)      :: aj, xj, mu, x1,x2,x3,x4, OmegaRad
        real(dl)      :: H2, H_0, grhom, grhorad, epsilon, h, w_DE, H2_LCDM
        integer       :: j, ind, i

        grhom      = params_cache%grhob + params_cache%grhoc
        grhorad    = params_cache%grhornomass +  params_cache%grhog
        OmegaDE_0  = 1._dl-params_cache%omegac -params_cache%omegab -params_cache%omegag- params_cache%omegar
        OmegaRad_0 = params_cache%omegag+params_cache%omegar
        OmegaM_0   = params_cache%omegac +params_cache%omegab

       ! some feedback:
        if ( feedback_level>0 ) then
            write(*,'(a)') "***************************************************************"
            write(*,'(a)') ' EFTCAMB beyond Horndeski background solver'
            write(*,'(a)')
        end if

        success= .true.

        !call self%feedback()

        ! initialize interpolating functions:
        call self%fun_x1%initialize      ( self%designer_num_points, self%x_initial, self%x_final )
        call self%fun_x2%initialize      ( self%designer_num_points, self%x_initial, self%x_final )
        call self%fun_x3%initialize      ( self%designer_num_points, self%x_initial, self%x_final )
        call self%fun_x4%initialize      ( self%designer_num_points, self%x_initial, self%x_final )

	    !SP: check omegav
	    if (params_cache%omegav < 0.3_dl) then
		    success = .false.
		    return
	    end if

        ! Initial values of x2 used in the bisection algorithm

        if (self%LCDM_background .or. self%bHor_x10*1.d-16 == -2._dl* OmegaDE_0 /(OmegaRad_0/exp(self%x_initial)**4 + OmegaM_0/exp(self%x_initial)**3 + OmegaDE_0)) then
        ! In this case x2=-1.5 x1 and then x1_i= -2 * Omega_DE_i
            self%bHor_x10= -2._dl* OmegaDE_0 /(OmegaRad_0/exp(self%x_initial)**4 + OmegaM_0/exp(self%x_initial)**3 + OmegaDE_0)
            self%bHor_x20= -1.5*self%bHor_x10
            if (feedback_level>0) then
                write(*,*)  '   x1 initial   =   -2 *Omega_DE initial   =   ', self%bHor_x10
                write(*,*)  '   x2 initial   =    3 *Omega_DE initial   =   ', self%bHor_x20
                write(*,*)  '   x3 initial   =    0'
                write(*,*)  '   x4 initial   =    0'
            end if
            call self%solve_background( params_cache, tol, success=success )
        else
            self%bHor_x20= -1.5*self%bHor_x10*1.d-16
            x2_i= self%bHor_x20
            x2_il= 0._dl
            x2_ir= 10._dl*x2_i*(abs(self%bHor_x10*1.d-16/(-2._dl* OmegaDE_0 /(OmegaRad_0/exp(self%x_initial)**4 + OmegaM_0/exp(self%x_initial)**3 + OmegaDE_0))))
            call self%solve_background( params_cache, tol, success=success )
            i=0
            do while (abs(tol)>1.d-6)
                self%bHor_x20= x2_ir
                call self%solve_background( params_cache, tol_r, success=success )
                self%bHor_x20= x2_il
                call self%solve_background( params_cache, tol_l, success=success )
                if (sign(1._dl,tol_r) == sign(1._dl,tol_l)) then
                    if (sign(1._dl,tol_r)<0._dl) then
                        write(*,*) 'Error in the bisection algorithm '
                        write(*,*) 'Negative right tolerance:', tol_r
                        write(*,*) 'Use an higher value for x2_ir'
                        write(*,*) 'Omegav = ',params_cache%omegav
                    else
                        write(*,*) 'Error in the bisection algorithm '
                        write(*,*) 'Positive left tolerance:', tol_l
                        write(*,*) 'No positive values of x2_i allow to recover Omega_{DE}. Use a lower initial value for x3 and/or x4 '
                        write(*,*) 'Omegav = ',params_cache%omegav
                    end if
                    call MpiStop('EFTCAMB error')
                else
                    if (sign(1._dl,tol)==sign(1._dl,tol_r)) then
                        x2_ir=x2_i
                        x2_i= (x2_il+x2_i)/2._dl
                    end if
                    if (sign(1._dl,tol)==sign(1._dl,tol_l)) then
                        x2_il= x2_i
                        x2_i= (x2_ir+x2_i)/2._dl
                    end if
                end if
                self%bHor_x20= x2_i
                call self%solve_background( params_cache, tol, success=success )
                i=i+1
            end do
            if (feedback_level>0) then
                write(*,*) '---------------------------------------------------'
                write(*,'(a, E13.3,a,I3)') 'Bisection algorithm ended successfully, with tollerance',  tol , '  and steps  ', i
                write(*,'(a)')
                write(*,'(a, E13.3)') 'x2 initial    =     ',  x2_i
                write(*,'(a)')
                write(*,*) '---------------------------------------------------'
            end if
        end if
        if (feedback_level>0) then
            xj = 0._dl
            call self%fun_x1%precompute(xj, ind,mu)
            write(*,'(a)') 'today values'
            write(*,'(a,F11.8)') 'x1 = ', self%fun_x1%value(xj, index=ind, coeff=mu)
            write(*,'(a,F10.8)') 'x2 = ', self%fun_x2%value(xj, index=ind, coeff=mu)
            write(*,'(a,F10.8)') 'x3 = ', self%fun_x3%value(xj, index=ind, coeff=mu)
            write(*,'(a,F10.8)') 'x4 = ', self%fun_x4%value(xj, index=ind, coeff=mu)
        end if
        if (self%debug_flag) then
          xj = 0._dl
          call self%fun_x1%precompute(xj, ind,mu)
          OmegaDE_0 = self%fun_x1%value(xj, index=ind, coeff=mu)+self%fun_x2%value(xj, index=ind, coeff=mu)+&
              &self%fun_x3%value(xj, index=ind, coeff=mu) +self%fun_x4%value(xj, index=ind, coeff=mu)
          H_0= sqrt((grhorad +grhom )/((1._dl-OmegaDE_0)*3._dl))
          OmegaRad_0 =  grhorad/(3._dl*H_0**2)
          open(unit=111, action = 'write', file='evolution_x1.dat')
          open(unit=222, action = 'write', file='evolution_x2.dat')
          open(unit=333, action = 'write', file='evolution_x3.dat')
          open(unit=444, action = 'write', file='evolution_x4.dat')
          open(unit=555, action = 'write', file='evolution_Omega.dat')
          open(unit=666, action = 'write', file='specific_background_quantities.dat')
          do j =1, self%fun_x1%num_points
            aj = exp(self%fun_x1%x(j))
            xj = self%fun_x1%x(j)
            call self%fun_x1%precompute(xj, ind,mu)
            write(111, *) aj, self%fun_x1%value(xj, index=ind, coeff=mu)
            write(222, *) aj, self%fun_x2%value(xj, index=ind, coeff=mu)
            write(333, *) aj, self%fun_x3%value(xj, index=ind, coeff=mu)
            write(444, *) aj, self%fun_x4%value(xj, index=ind, coeff=mu)

            x1=self%fun_x1%value(xj, index=ind, coeff=mu)
            x2=self%fun_x2%value(xj, index=ind, coeff=mu)
            x3=self%fun_x3%value(xj, index=ind, coeff=mu)
            x4=self%fun_x4%value(xj, index=ind, coeff=mu)
            OmegaDE = self%fun_x1%value(xj, index=ind, coeff=mu)+self%fun_x2%value(xj, index=ind, coeff=mu)+&
                &self%fun_x3%value(xj, index=ind, coeff=mu) +self%fun_x4%value(xj, index=ind, coeff=mu)
            H2= (grhorad/aj**4+grhom/aj**3)/((1._dl-OmegaDE)*3._dl)
            H2_LCDM= (grhorad/aj**4+grhom/aj**3)/((1._dl-OmegaDE_0)*3._dl)
            OmegaRad= (grhorad/aj**4)/(3._dl*(H2))
            epsilon = -((20._dl*(3._dl*x1 + 2._dl*x2) - 5._dl*x3*(-3._dl + 3._dl*x1 + x2 + OmegaRad) -x4*(36._dl*x1 + 16._dl*x2 + 3._dl*x3 + 8._dl*OmegaRad))/(5._dl*x3**2._dl + 20._dl*(x1 + 2._dl*x2 + x3) + 4._dl*(6._dl &
                      - x1 - 2._dl*x2 + 3._dl*x3)*x4 + 8._dl*x4**2._dl))
            h = -((15._dl*x3**2._dl + 12._dl*x4**2._dl + 10._dl*(x1 + 2._dl*x2)*(3._dl + 3._dl*x1 + x2 + OmegaRad) + 10._dl*x3*(3._dl + 6._dl*x1+ 3._dl*x2 + OmegaRad) + x4*(36._dl + 78._dl*x1 + 32._dl*x2 + 30._dl*x3     &
                + 12._dl*OmegaRad))/(5._dl*x3**2._dl + 20._dl*(x1 + 2._dl*x2 + x3) + 4._dl*(6._dl - x1 - 2._dl*x2 + 3._dl*x3)*x4 + 8._dl*x4**2._dl))

            w_DE= (5._dl*(3._dl*x1+x2-epsilon*x3)-(3._dl+8._dl*epsilon+2._dl*h)*x4)/(15._dl* (x1+x2+x3+x4))

            write(555, '(20E15.5)')1._dl/aj-1._dl, sqrt(H2), (grhorad/aj**4)/(3._dl*H2), (grhom/aj**3)/(3._dl*H2), OmegaDE
            write(666, '(20E15.5)')1._dl/aj-1._dl, h, epsilon, w_DE, (sqrt(H2)-sqrt(H2_LCDM))/sqrt(H2_LCDM)
          end do
          close(111)
          close(222)
          close(333)
          close(444)
          close(555)
          close(666)

        end if

        !store energy densities today
        xj = 0._dl

        call self%fun_x1%precompute(xj, ind,mu)

        OmegaDE_0  = self%fun_x1%value(xj, index=ind, coeff=mu)+self%fun_x2%value(xj, index=ind, coeff=mu)+ self%fun_x3%value(xj, index=ind, coeff=mu) + self%fun_x4%value(xj, index=ind, coeff=mu)
        H_0        = sqrt((grhorad +grhom )/((1._dl-OmegaDE_0)*3._dl))
        OmegaM_0   = grhom/(3._dl*H_0**2)
        OmegaRAD_0 =  grhorad/(3._dl*H_0**2)

        if (feedback_level>0) then
            write(*,*) '---------------------------------------------------'
            write(*,'(a,E13.4)') 'relative error on H_0     ',  (H_0 - params_cache%h0_Mpc)/H_0
            write(*,'(a,E13.4)') 'relative error on OmegaM_0',  (OmegaM_0-(params_cache%omegac +params_cache%omegab))/OmegaM_0
            write(*,'(a,F8.6)') 'OmegaRad_0 = ',  OmegaRad_0
            write(*,'(a,F8.6)') 'OmegaDE_0  = ',  OmegaDE_0
            write(*,'(a,F8.6)') 'OmegaM_0   = ',  OmegaM_0
            write(*,*) '---------------------------------------------------'
        end if

        !store x values today
        self%bHor_x1today = self%fun_x1%value(xj, index=ind, coeff=mu)
        self%bHor_x2today = self%fun_x2%value(xj, index=ind, coeff=mu)
        self%bHor_x3today = self%fun_x3%value(xj, index=ind, coeff=mu)
        self%bHor_x4today = self%fun_x4%value(xj, index=ind, coeff=mu)
        return

    end subroutine EFTCAMBBeyondHorndeskiInitBackground

    ! ---------------------------------------------------------------------------------------------

    subroutine EFTCAMBBeyondHorndeskiSolveBackground( self, params_cache, tol, success)

        implicit none

        class(EFTCAMB_Beyond_Horndeski)              :: self          !< the base class.
        type(TEFTCAMB_parameter_cache), intent(in)    :: params_cache  !< a EFTCAMB parameter cache containing cosmological parameters.
        real(dl), intent(out)                        :: tol            !< present day relative difference between H_0 input and H_0 computed.
        logical , intent(out)                        :: success       !< whether the calculation ended correctly or not

        real(dl) :: PPlus, yPlus, CoeffA_Part, yStar, x, OmegaDE_0 ,OmegaNu_EFT,Omegarad_EFT, grhom, grhorad

        integer, parameter :: num_eq = 4  !<  Number of equations of the dynamical system

        real(dl) :: y(num_eq), ydot(num_eq)

        integer  :: itol, itask, istate, iopt, LRN, LRS, LRW, LIS, LIN, LIW, JacobianMode, i
        real(dl) :: rtol, atol, t1, t2, B
        real(dl), allocatable :: rwork(:)
        integer,  allocatable :: iwork(:)

        ! 1) Cosmological densities:
        grhom   = params_cache%grhob + params_cache%grhoc
        grhorad = params_cache%grhornomass +  params_cache%grhog
        OmegaDE_0        = 1._dl-params_cache%omegac -params_cache%omegab -params_cache%omegag- params_cache%omegar

        if (self%LCDM_background) then
            ! 2) Set initial conditions to reproduce LCDM bacground
            x    = self%x_initial
            y(1) = self%bHor_x10
            y(2) = self%bHor_x20
            y(3) = 0.
            y(4) = 0.
        else
            ! 2) Set initial conditions:
            x    = self%x_initial
            y(1) = self%bHor_x10*1.d-16
            y(2) = self%bHor_x20
            y(3) = self%bHor_x30*1.d-9
            y(4) = self%bHor_x40*1.d-6
        end if

        ydot = 0._dl
        ! 3) Initialize DLSODA:
        ! set-up the relative and absolute tollerances:
        itol = 1
        rtol = 1.d-13
        atol = 1.d-18
        ! initialize task to do:
        itask  = 1
        istate = 1
        iopt   = 1
        ! initialize the work space:
        LRN = 20 + 16*num_eq
        LRS = 22 + 9*num_eq + num_eq**2
        LRW = max(LRN,LRS)
        LIS = 20 + num_eq
        LIN = 20
        LIW = max(LIS,LIN)
        ! allocate the arrays:
        allocate(rwork(LRW))
        allocate(iwork(LIW))
        ! optional lsoda input:
        RWORK(5) = 0._dl  ! the step size to be attempted on the first step. The default value is determined by the solver.
        RWORK(6) = 0._dl  ! the maximum absolute step size allowed. The default value is infinite.
        RWORK(7) = 0._dl  ! the minimum absolute step size allowed. The default value is 0.
        IWORK(5) = 0      ! flag to generate extra printing at method switches. IXPR = 0 means no extra printing (the default). IXPR = 1 means print data on each switch.
        IWORK(6) = 100    ! maximum number of (internally defined) steps allowed during one call to the solver. The default value is 500.
        IWORK(7) = 0      ! maximum number of messages printed (per problem) warning that T + H = T on a step (H = step size). This must be positive to result in a non-default value.  The default value is 10.
        IWORK(8) = 0      ! the maximum order to be allowed for the nonstiff (Adams) method.  the default value is 12.
        IWORK(9) = 0      ! the maximum order to be allowed for the stiff (BDF) method.  The default value is 5.
        ! additional lsoda stuff:
        CALL XSETF(0) ! suppress odepack printing
        ! Jacobian mode: 1=fullJacobian, 2=not provided
        JacobianMode = 2

        ! store the first value of the EFT functions:
        t1  = self%fun_x1%x(1)
        ! compute output EFT functions if needed:
        call output( num_eq, 1, t1, y )

        ! 4) solve the equations:
        do i=1, self%fun_x1%num_points-1

            ! set the time step:
            t1 = self%fun_x1%x(i)
            t2 = self%fun_x1%x(i+1)
            ! solve the system:
            call DLSODA ( derivs, num_eq, y, t1, t2, itol, rtol, atol, itask, istate, iopt, RWORK, LRW, IWORK, LIW, jacobian, JacobianMode)
            ! check istate for LSODA good completion:
            if ( istate < 0 ) then
                success = .False.
                return
            end if
            call output( num_eq, i+1, t2, y )
            if ( i== self%fun_x1%num_points-1 ) then
                tol = (OmegaDE_0-(y(1)+y(2)+y(3)+y(4)))/OmegaDE_0
            end if

        end do

    contains

      ! ---------------------------------------------------------------------------------------------
      !> Subroutine that computes y' given y for the background of designer f(R)
      subroutine derivs( num_eq, x, y, ydot )

          implicit none

          integer , intent(in)                     :: num_eq !< number of equations in the ODE system
          real(dl), intent(in)                     :: x      !< time at which the derivatives of the system are computed
          real(dl), intent(in) , dimension(num_eq) :: y      !< input status of the system
          real(dl), intent(out), dimension(num_eq) :: ydot   !< value of the derivative at time x

          real(dl) :: qs, epsilon_phi, h, H2, OmegaRad

          !1) compute the functionals entering the ODEs

          H2= (grhorad/exp(x)**4+grhom/exp(x)**3)/((1._dl-y(1)-y(2)-y(3)-y(4))*3._dl)

          OmegaRad= (grhorad/exp(x)**4)/(3._dl*(H2))

          qs = 5._dl*y(3)**2._dl + 20._dl*(y(1) + 2._dl*y(2) + y(3)) + 4._dl*(6._dl - y(1) - 2._dl*y(2) + 3._dl*y(3))*y(4) &
              &+ 8._dl*y(4)**2._dl

          h = -((15._dl*y(3)**2._dl + 12._dl*y(4)**2._dl + 10._dl*(y(1) + 2._dl*y(2))*(3._dl + 3._dl*y(1) + y(2) + OmegaRad) &
              &+ 10._dl*y(3)*(3._dl + 6._dl*y(1)+ 3._dl*y(2) + OmegaRad) + y(4)*(36._dl + 78._dl*y(1) + 32._dl*y(2) + 30._dl*y(3)&
              & + 12._dl*OmegaRad))/(5._dl*y(3)**2._dl + 20._dl*(y(1) + 2._dl*y(2) + y(3)) + 4._dl*(6._dl - y(1) - 2._dl*y(2) &
              &+ 3._dl*y(3))*y(4) + 8._dl*y(4)**2._dl))

          epsilon_phi = -((20._dl*(3._dl*y(1) + 2._dl*y(2)) - 5._dl*y(3)*(-3._dl + 3._dl*y(1) + y(2) + OmegaRad) - y(4)*(36._dl*y(1) &
              &+ 16._dl*y(2) + 3._dl*y(3) + 8._dl*OmegaRad))/(5._dl*y(3)**2._dl + 20._dl*(y(1) + 2._dl*y(2) + y(3)) + 4._dl*(6._dl &
              &- y(1) - 2._dl*y(2) + 3._dl*y(3))*y(4) + 8._dl*y(4)**2._dl))

          ! 2) Get the equations of motion:
          ydot(1) = 2._dl*y(1)*( epsilon_phi -h)
          ydot(2) = 2._dl*y(2)*( 2._dl*epsilon_phi -h)
          ydot(3) = y(3)*( 3._dl*epsilon_phi -h)
          ydot(4) = 4._dl*y(4)*epsilon_phi

      end subroutine

      ! ---------------------------------------------------------------------------------------------
      !> Subroutine that computes the Jacobian of the system. Now a dummy function.
      !! Implementing it might increase performances.
      subroutine jacobian( num_eq, x, y, ml, mu, pd, nrowpd )

          implicit none

          integer                            :: num_eq !< number of components of the Jacobian
          integer                            :: ml     !< ignored
          integer                            :: mu     !< ignored
          integer                            :: nrowpd !< ignored
          real(dl)                           :: x      !< time at which the Jacobian is computed
          real(dl), dimension(num_eq)        :: y      !< input status of the system
          real(dl), dimension(nrowpd,num_eq) :: pd     !< output Jacobian

      end subroutine jacobian

      ! ---------------------------------------------------------------------------------------------
      !> Subroutine that takes the solution of the background f(R) equations and computes the value of
      !! B and stores the values of the EFT functions.
      subroutine output( num_eq, ind, x, y )

          implicit none

          integer , intent(in)                     :: num_eq !< number of equations in the ODE system.
          integer , intent(in)                     :: ind    !< index of the EFT functions interpolation tables to fill.
          real(dl), intent(in)                     :: x      !< time at which the derivatives of the system are computed.
          real(dl), intent(in) , dimension(num_eq) :: y      !< input status of the system.

          real(dl) :: qs, epsilon, h, OmegaDE, epsilon_p, epsilon_pp, hprime
          real(dl) :: ydot(num_eq),alpha_K, alpha_B, Q_s, Q_t, Sigma, cs2

          real(dl) :: x1,x2,x3,x4,x1p,x2p,x3p,x4p, x1pp,x2pp, x3pp,x4pp, x4ppp,OmegaRad, OmegaRadprime,OmegaRadprimeprime

          real(dl) :: a , H2,adotoa,Hdot,temp

          a= exp(x)

          ! 1) call derivs to compute ydot at a given time:
          call derivs( num_eq, x, y, ydot )

          ! 2) store some quantities
          x1 = y(1)
          x2 = y(2)
          x3 = y(3)
          x4 = y(4)

          x1p = ydot(1)
          x2p = ydot(2)
          x3p = ydot(3)
          x4p = ydot(4)

          if (x <self%x_initial) then
              x1 = self%bHor_x10*1.d-16*(a/exp(self%x_initial))**(14._dl/3._dl)
              x2 = self%bHor_x20*(a/exp(self%x_initial))**(16._dl/3._dl)
              x3 = self%bHor_x30*1.d-9*(a/exp(self%x_initial))**(3._dl)
              x4 = self%bHor_x40*1.d-6*(a/exp(self%x_initial))**(4._dl/3._dl)
          end if


          H2= (grhorad/exp(x)**4+grhom/exp(x)**3)/((1._dl-x1-x2-x3-x4)*3._dl)

          OmegaRad= (grhorad/exp(x)**4)/(3._dl*(H2))

          ! 3) compute the functionals entering the ODEs

          epsilon = -((20._dl*(3._dl*x1 + 2._dl*x2) - 5._dl*x3*(-3._dl + 3._dl*x1 + x2 + OmegaRad) -x4*(36._dl*x1 + 16._dl*x2 + 3._dl*x3 + 8._dl*OmegaRad))/(5._dl*x3**2._dl + 20._dl*(x1 + 2._dl*x2 + x3) + 4._dl*(6._dl     &
                    - x1 - 2._dl*x2 + 3._dl*x3)*x4 + 8._dl*x4**2._dl))

          h = -((15._dl*x3**2._dl + 12._dl*x4**2._dl + 10._dl*(x1 + 2._dl*x2)*(3._dl + 3._dl*x1 + x2 + OmegaRad) + 10._dl*x3*(3._dl + 6._dl*x1+ 3._dl*x2 + OmegaRad) + x4*(36._dl + 78._dl*x1 + 32._dl*x2 + 30._dl*x3         &
              + 12._dl*OmegaRad))/(5._dl*x3**2._dl + 20._dl*(x1 + 2._dl*x2 + x3) + 4._dl*(6._dl - x1 - 2._dl*x2 + 3._dl*x3)*x4 + 8._dl*x4**2._dl))

          if (x <self%x_initial) then
              x1p = 2._dl*x1*( epsilon -h)
              x2p  = 2._dl*x2*( 2._dl*epsilon -h)
              x3p  = x3*( 3._dl*epsilon -h)
              x4p  = 4._dl*x4*epsilon
          end if

          OmegaRadprime = -2._dl*OmegaRad*(2._dl+h)

          hprime = ((15._dl*x3**2._dl + 12._dl*x4**2._dl + 10._dl*(x1 + 2._dl*x2)*(3._dl + 3._dl*x1 + x2 + OmegaRad) &
              &+ 10._dl*x3*(3._dl + 6._dl*x1 + 3._dl*x2 + OmegaRad) + x4*(36._dl + 78._dl*x1 + 32._dl*x2 + 30._dl*x3 &
              &+ 12._dl*OmegaRad))*(10._dl*x3*x3p + 20._dl*(x1p + 2._dl*x2p + x3p) + 4._dl*x4*(-x1p - 2._dl*x2p + 3._dl*x3p) &
              &+4._dl*(6._dl - x1 - 2._dl*x2 + 3._dl*x3)*x4p + 16._dl*x4*x4p))/(5._dl*x3**2._dl + 20._dl*(x1 + 2._dl*x2 + x3)&
              & + 4._dl*(6._dl - x1 - 2._dl*x2 + 3._dl*x3)*x4 + 8._dl*x4**2._dl)**2._dl -(10._dl*(3._dl + 3._dl*x1 + x2 &
              &+ OmegaRad)*(x1p + 2._dl*x2p) + 30._dl*x3*x3p + 10._dl*(3._dl + 6._dl*x1 + 3._dl*x2 + OmegaRad)*x3p &
              &+24._dl*x4*x4p + (36._dl + 78._dl*x1 + 32._dl*x2 + 30._dl*x3 + 12._dl*OmegaRad)*x4p +10._dl*(x1 + 2._dl*x2)*(&
              &3._dl*x1p + x2p + OmegaRadprime) + 10._dl*x3*(6._dl*x1p + 3._dl*x2p + OmegaRadprime) +x4*(78._dl*x1p &
              &+ 32._dl*x2p + 30._dl*x3p + 12._dl*OmegaRadprime))/(5._dl*x3**2. + 20._dl*(x1 + 2._dl*x2 + x3) + 4._dl*(6._dl &
              &- x1 - 2._dl*x2 + 3._dl*x3)*x4 + 8._dl*x4**2._dl)

          epsilon_p = ((20._dl*(3._dl*x1 + 2._dl*x2) - 5._dl*x3*(-3._dl +3._dl*x1 +x2 +OmegaRad) -x4*(36._dl*x1 +16._dl*x2 &
              &+3._dl*x3 +8._dl*OmegaRad))*(10._dl*x3*x3p + 20._dl*(x1p + 2._dl*x2p + x3p) + 4._dl*x4*(-x1p - 2._dl*x2p + 3._dl*x3p) &
              &+4._dl*(6._dl - x1 - 2._dl*x2 + 3._dl*x3)*x4p + 16._dl*x4*x4p))/(5._dl*x3**2._dl + 20._dl*(x1 + 2._dl*x2 + x3) + 4._dl*(&
              &6._dl - x1 - 2._dl*x2 + 3._dl*x3)*x4 + 8._dl*x4**2._dl)**2._dl -(20._dl*(3._dl*x1p + 2._dl*x2p) - 5._dl*(-3._dl + 3._dl&
              &*x1 + x2 + OmegaRad)*x3p - (36._dl*x1 + 16._dl*x2 + 3._dl*x3 + 8._dl*OmegaRad)*x4p -5._dl*x3*(3._dl*x1p + x2p &
              &+ OmegaRadprime) - x4*(36._dl*x1p + 16._dl*x2p + 3._dl*x3p + 8._dl*OmegaRadprime))/(5._dl*x3**2._dl + 20._dl*(x1 &
              &+ 2._dl*x2 + x3) + 4._dl*(6._dl - x1 - 2._dl*x2 + 3._dl*x3)*x4 + 8._dl*x4**2._dl)

          x1pp = 2._dl*( x1p*( epsilon- h ) +x1*( epsilon_p -hprime ) )
          x2pp = 2._dl*( x2p*( 2._dl*epsilon- h ) +x2*( 2._dl*epsilon_p -hprime ) )
          x3pp = x3p*( 3._dl*epsilon- h ) +x3*( 3._dl*epsilon_p -hprime )
          x4pp = 4._dl*( x4p*epsilon +x4*epsilon_p )
          OmegaRadprimeprime = -2._dl*OmegaRadprime*(2._dl +h) -2._dl*OmegaRad*hprime

          epsilon_pp = (-2._dl*(60._dl*x1 + 40._dl*x2 - 5._dl*x3*(-3._dl + 3._dl*x1 + x2 + OmegaRad) - x4*(36._dl*x1 + 16._dl*x2 &
              &+ 3._dl*x3 + 8._dl*OmegaRad))*(-4._dl*x4*(x1p + 2._dl*x2p - 3._dl*x3p) + 10._dl*x3*x3p + 20._dl*(x1p &
              &+ 2._dl*x2p + x3p) -4._dl*(x1 + 2._dl*x2 - 3._dl*(2._dl + x3))*x4p + 16._dl*x4*x4p)**2._dl +2._dl*(5._dl*x3**2._dl &
              &+ 20._dl*(x1 + 2._dl*x2 + x3) - 4._dl*(x1 + 2._dl*x2 - 3._dl*(2._dl + x3))*x4 + 8._dl*x4**2._dl)*(-4._dl*x4*(x1p &
              &+ 2._dl*x2p - 3._dl*x3p) + 10._dl*x3*x3p + 20._dl*(x1p + 2._dl*x2p + x3p) -4._dl*(x1 + 2._dl*x2 - 3._dl*(2._dl + x3))*x4p &
              &+ 16._dl*x4*x4p)*(60._dl*x1p + 40._dl*x2p - 5._dl*(-3._dl + 3._dl*x1 + x2 + OmegaRad)*x3p - (36._dl*x1 + 16._dl*x2 &
              &+ 3._dl*x3 + 8._dl*OmegaRad)*x4p -5._dl*x3*(3._dl*x1p + x2p + OmegaRadprime) - x4*(36._dl*x1p + 16._dl*x2p + 3._dl*x3p &
              &+ 8._dl*OmegaRadprime)) +2._dl*(5._dl*x3**2._dl + 20._dl*(x1 + 2._dl*x2 + x3) - 4._dl*(x1 + 2._dl*x2 - 3._dl*(2._dl &
              &+ x3))*x4 + 8._dl*x4**2._dl)*(60._dl*x1 + 40._dl*x2 - 5._dl*x3*(-3._dl + 3._dl*x1 + x2 + OmegaRad) - x4*(36._dl*x1 &
              &+ 16._dl*x2 + 3._dl*x3 + 8._dl*OmegaRad))*(5._dl*x3p**2._dl - 4._dl*(x1p + 2._dl*x2p - 3._dl*x3p)*x4p + 8._dl*x4p**2._dl &
              &-2._dl*x4*(x1pp + 2._dl*x2pp - 3._dl*x3pp) + 5._dl*x3*x3pp + 10._dl*(x1pp + 2._dl*x2pp + x3pp) -2._dl*(x1 + 2._dl*x2 &
              &- 3._dl*(2._dl + x3))*x4pp + 8._dl*x4*x4pp) -(5._dl*x3**2._dl + 20._dl*(x1 + 2._dl*x2 + x3) - 4._dl*(x1 + 2._dl*x2 &
              &- 3._dl*(2._dl + x3))*x4 + 8._dl*x4**2._dl)**2._dl*(-10._dl*x3p*(3._dl*x1p + x2p + OmegaRadprime) -2._dl*x4p*(36._dl*x1p &
              &+ 16._dl*x2p + 3._dl*x3p + 8._dl*OmegaRadprime) + 60._dl*x1pp + 40._dl*x2pp -5._dl*(-3._dl + 3._dl*x1 + x2 &
              &+ OmegaRad)*x3pp - (36._dl*x1 + 16._dl*x2 + 3._dl*x3 + 8._dl*OmegaRad)*x4pp -5._dl*x3*(3._dl*x1pp + x2pp &
              &+ OmegaRadprimeprime) - x4*(36._dl*x1pp + 16._dl*x2pp + 3._dl*x3pp + 8._dl*OmegaRadprimeprime)))/(5._dl*x3**2._dl &
              &+ 20._dl*(x1 + 2._dl*x2 + x3) - 4._dl*(x1 + 2._dl*x2 - 3._dl*(2._dl + x3))*x4 + 8._dl*x4**2._dl)**3._dl

          x4ppp = 4._dl*( x4pp*epsilon +2._dl*X4p*epsilon_p +x4*epsilon_pp )

          ! 3) compute the x_i functions: more derivatives needed
          self%fun_x1%y(ind) = x1
          self%fun_x1%yp(ind) = x1p
          self%fun_x1%ypp(ind) = x1pp
          self%fun_x2%y(ind) = x2
          self%fun_x2%yp(ind) = x2p
          self%fun_x2%ypp(ind) = x2pp
          self%fun_x3%y(ind) = x3
          self%fun_x3%yp(ind) = x3p
          self%fun_x3%ypp(ind) = x3pp
          self%fun_x4%y(ind) = x4
          self%fun_x4%yp(ind) = x4p
          self%fun_x4%ypp(ind) = x4pp
          self%fun_x4%yppp(ind) = x4ppp

      end subroutine

    end subroutine EFTCAMBBeyondHorndeskiSolveBackground

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes the number of parameters of the model.
    subroutine EFTCAMBBeyondHorndeskiComputeParametersNumber( self )

        implicit none

        class(EFTCAMB_Beyond_Horndeski)  :: self   !< the base class

        self%parameter_number = 3

    end subroutine EFTCAMBBeyondHorndeskiComputeParametersNumber

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that prints on the screen feedback information about the model.
    subroutine EFTCAMBBeyondHorndeskiFeedback( self, print_params )

        implicit none

        class(EFTCAMB_Beyond_Horndeski)    :: self         !< the base class
        logical, optional                  :: print_params !< optional flag that decised whether to print numerical values
                                                       !! of the parameters.

        logical                            :: print_params_temp

        ! print general model informations:
        if (self%LCDM_background) then

            write(*,*) 'LambdaCDM background flag selected'
            write(*,'(a,a)')    '   Model               =  ', self%name
            write(*,'(a,I3)')   '   Number of params    ='  , self%parameter_number
        else
            write(*,*)
            write(*,'(a,a)')    '   Model               =  ', self%name
            write(*,'(a,I3)')   '   Number of params    ='  , self%parameter_number
            write(*,'(a,E13.3)')  '   x1 initial          ='  , self%bHor_x10*10.**(-16)
            write(*,'(a,E13.3)')  '   x3 initial          ='  , self%bHor_x30*10.**(-9)
            write(*,'(a,E13.3)')  '   x4 initial          ='  , self%bHor_x40*10.**(-6)
        end if

        ! print the values of the parameters:
        if ( present(print_params) ) then
            print_params_temp = print_params
        else
            print_params_temp = .True.
        end if

    end subroutine EFTCAMBBeyondHorndeskiFeedback

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the i-th parameter name of the model
    subroutine EFTCAMBBeyondHorndeskiParameterNames( self, i, name )

        implicit none

        class(EFTCAMB_Beyond_Horndeski)   :: self   !< the base class
        integer     , intent(in)          :: i      !< the index of the parameter
        character(*), intent(out)         :: name   !< the output name of the i-th parameter

        ! check the input index:
        if ( i>self%parameter_number ) then
            write(*,*) 'Illegal index for parameter_names.'
            write(*,*) 'Maximum value is:', self%parameter_number
            call MpiStop('EFTCAMB error')
        end if
        ! return the appropriate name:
        if ( i==1 ) then
            name = 'Beyond_Horndeski_x10'
            return
        end if
        if ( i==2 ) then
            name = 'Beyond_Horndeski_x30'
            return
        end if
        if ( i==3 ) then
            name = 'Beyond_Horndeski_x40'
            return
        end if
        if ( i==0 ) then
            name = 'no_name'
            return
        end if

    end subroutine EFTCAMBBeyondHorndeskiParameterNames

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the i-th parameter name of the model
    subroutine EFTCAMBBeyondHorndeskiParameterNamesLatex( self, i, latexname )

        implicit none

        class(EFTCAMB_Beyond_Horndeski)   :: self       !< the base class
        integer     , intent(in)          :: i          !< The index of the parameter
        character(*), intent(out)         :: latexname  !< the output latex name of the i-th parameter

        ! check the input index:
        if ( i>self%parameter_number ) then
            write(*,*) 'Illegal index for parameter_names_latex.'
            write(*,*) 'Maximum value is:', self%parameter_number
            call MpiStop('EFTCAMB error')
        end if
        ! return the appropriate name:
        if ( i==1 ) then
            latexname = 'x_1^0'
            return
        end if
        if ( i==2 ) then
            latexname = 'x_3^0'
            return
        end if
        if ( i==3 ) then
            latexname = 'x_4^0'
            return
        end if
        if ( i==0 ) then
            latexname = 'noname'
            return
        end if

    end subroutine EFTCAMBBeyondHorndeskiParameterNamesLatex

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the i-th parameter value of the model
    subroutine EFTCAMBBeyondHorndeskiParameterValues( self, i, value )

        implicit none

        class(EFTCAMB_Beyond_Horndeski)   :: self   !< the base class
        integer , intent(in)              :: i      !< The index of the parameter
        real(dl), intent(out)             :: value  !< the output value of the i-th parameter

        ! check the input index:
        if ( i>self%parameter_number ) then
            write(*,*) 'Illegal index for parameter_value.'
            write(*,*) 'Maximum value is:', self%parameter_number
            call MpiStop('EFTCAMB error')
        end if
        ! return the appropriate name:
        if ( i==0 ) then
            value = 0._dl
            return
        end if
        if ( i==1 ) then
            value = self%bHor_x10
            return
        end if
        if ( i==2 ) then
            value = self%bHor_x30
            return
        end if
        if ( i==3 ) then
            value = self%bHor_x40
            return
        end if

        if ( i==-1 ) then
            value = self%bHor_x1today
            return
        end if
        if ( i==-2 ) then
            value = self%bHor_x2today
            return
        end if
        if ( i==-3 ) then
            value = self%bHor_x3today
            return
        end if
        if ( i==-4 ) then
            value = self%bHor_x4today
            return
        end if

    end subroutine EFTCAMBBeyondHorndeskiParameterValues

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes the value of the background EFT functions at a given time.
    subroutine EFTCAMBBeyondHorndeskiBackgroundEFTFunctions( self, a, eft_par_cache, eft_cache )

        implicit none

        class(EFTCAMB_Beyond_Horndeski)              :: self          !< the base class
        real(dl), intent(in)                         :: a             !< the input scale factor
        type(TEFTCAMB_parameter_cache), intent(inout) :: eft_par_cache !< the EFTCAMB parameter cache that contains all the physical parameters.
        type(TEFTCAMB_timestep_cache ), intent(inout) :: eft_cache     !< the EFTCAMB timestep cache that contains all the physical values.
        real(dl) :: x, mu
        integer  :: ind

        real(dl) :: x1,x2,x3,x4, x1p,x2p,x3p,x4p, x1pp,x2pp,x4pp,x3pp,x3ppp, x4ppp,qs,qsprime
        real(dl) :: OmegaRad,OmegaRadprime,h,hprime,hpprime, grhom, grhorad, H2, qspprime, adotoa, Hdot,Hdotdot,Hdotdotdot
        real(dl) :: epsilon,epsilon_p,epsilon_pp,OmegaRadprimeprime, EFTcdot, EFTLAMBDAdot, EFTc, EFTLAMBDA

        if(a==0._dl) return

        x   = log(a)

        grhom   = eft_par_cache%grhob + eft_par_cache%grhoc
        grhorad = eft_par_cache%grhornomass +  eft_par_cache%grhog

        call self%fun_x4%precompute(x, ind, mu )
        x1    = self%fun_x1%value(x, index=ind, coeff=mu)
        x2    = self%fun_x2%value(x, index=ind, coeff=mu)
        x3    = self%fun_x3%value(x, index=ind, coeff=mu)
        x4    = self%fun_x4%value(x, index=ind, coeff=mu)

        if (x <self%x_initial) then
            x1 = self%bHor_x10*1.d-16*(a/exp(self%x_initial))**(14._dl/3._dl)
            x2 = self%bHor_x20*(a/exp(self%x_initial))**(16._dl/3._dl)
            x3 = self%bHor_x30*1.d-9*(a/exp(self%x_initial))**(3._dl)
            x4 = self%bHor_x40*1.d-6*(a/exp(self%x_initial))**(4._dl/3._dl)
        end if

        H2= (grhorad/exp(x)**4+grhom/exp(x)**3)/((1._dl-x1-x2-x3-x4)*3._dl)
        OmegaRad= (grhorad/exp(x)**4)/(3._dl*(H2))

          epsilon = -((20._dl*(3._dl*x1 + 2._dl*x2) - 5._dl*x3*(-3._dl + 3._dl*x1 + x2 + OmegaRad) -x4*(36._dl*x1 &
              &+ 16._dl*x2 + 3._dl*x3 + 8._dl*OmegaRad))/(5._dl*x3**2._dl + 20._dl*(x1 + 2._dl*x2 + x3) + 4._dl*(6._dl &
              &- x1 - 2._dl*x2 + 3._dl*x3)*x4 + 8._dl*x4**2._dl))

          h = -((15._dl*x3**2._dl + 12._dl*x4**2._dl + 10._dl*(x1 + 2._dl*x2)*(3._dl + 3._dl*x1 + x2 + OmegaRad) &
              &+ 10._dl*x3*(3._dl + 6._dl*x1+ 3._dl*x2 + OmegaRad) + x4*(36._dl + 78._dl*x1 + 32._dl*x2 + 30._dl*x3&
              & + 12._dl*OmegaRad))/(5._dl*x3**2._dl + 20._dl*(x1 + 2._dl*x2 + x3) + 4._dl*(6._dl - x1 - 2._dl*x2 &
              &+ 3._dl*x3)*x4 + 8._dl*x4**2._dl))

          x1p = 2._dl*x1*( epsilon -h)
          x2p = 2._dl*x2*( 2._dl*epsilon -h)
          x3p = x3*( 3._dl*epsilon -h)
          x4p = 4._dl*x4*epsilon

          OmegaRadprime = -2._dl*OmegaRad*(2._dl+h)

          hprime = ((15._dl*x3**2._dl + 12._dl*x4**2._dl + 10._dl*(x1 + 2._dl*x2)*(3._dl + 3._dl*x1 + x2 + OmegaRad) &
              &+ 10._dl*x3*(3._dl + 6._dl*x1 + 3._dl*x2 + OmegaRad) + x4*(36._dl + 78._dl*x1 + 32._dl*x2 + 30._dl*x3 &
              &+ 12._dl*OmegaRad))*(10._dl*x3*x3p + 20._dl*(x1p + 2._dl*x2p + x3p) + 4._dl*x4*(-x1p - 2._dl*x2p + 3._dl*x3p) &
              &+4._dl*(6._dl - x1 - 2._dl*x2 + 3._dl*x3)*x4p + 16._dl*x4*x4p))/(5._dl*x3**2._dl + 20._dl*(x1 + 2._dl*x2 + x3)&
              & + 4._dl*(6._dl - x1 - 2._dl*x2 + 3._dl*x3)*x4 + 8._dl*x4**2._dl)**2._dl -(10._dl*(3._dl + 3._dl*x1 + x2 &
              &+ OmegaRad)*(x1p + 2._dl*x2p) + 30._dl*x3*x3p + 10._dl*(3._dl + 6._dl*x1 + 3._dl*x2 + OmegaRad)*x3p &
              &+24._dl*x4*x4p + (36._dl + 78._dl*x1 + 32._dl*x2 + 30._dl*x3 + 12._dl*OmegaRad)*x4p +10._dl*(x1 + 2._dl*x2)*(&
              &3._dl*x1p + x2p + OmegaRadprime) + 10._dl*x3*(6._dl*x1p + 3._dl*x2p + OmegaRadprime) +x4*(78._dl*x1p &
              &+ 32._dl*x2p + 30._dl*x3p + 12._dl*OmegaRadprime))/(5._dl*x3**2. + 20._dl*(x1 + 2._dl*x2 + x3) + 4._dl*(6._dl &
              &- x1 - 2._dl*x2 + 3._dl*x3)*x4 + 8._dl*x4**2._dl)

        epsilon_p = ((20._dl*(3._dl*x1 + 2._dl*x2) - 5._dl*x3*(-3._dl +3._dl*x1 +x2 +OmegaRad) -x4*(36._dl*x1 +16._dl*x2 &
              &+3._dl*x3 +8._dl*OmegaRad))*(10._dl*x3*x3p + 20._dl*(x1p + 2._dl*x2p + x3p) + 4._dl*x4*(-x1p - 2._dl*x2p + 3._dl*x3p) &
              &+4._dl*(6._dl - x1 - 2._dl*x2 + 3._dl*x3)*x4p + 16._dl*x4*x4p))/(5._dl*x3**2._dl + 20._dl*(x1 + 2._dl*x2 + x3) + 4._dl*(&
              &6._dl - x1 - 2._dl*x2 + 3._dl*x3)*x4 + 8._dl*x4**2._dl)**2._dl -(20._dl*(3._dl*x1p + 2._dl*x2p) - 5._dl*(-3._dl + 3._dl&
              &*x1 + x2 + OmegaRad)*x3p - (36._dl*x1 + 16._dl*x2 + 3._dl*x3 + 8._dl*OmegaRad)*x4p -5._dl*x3*(3._dl*x1p + x2p &
              &+ OmegaRadprime) - x4*(36._dl*x1p + 16._dl*x2p + 3._dl*x3p + 8._dl*OmegaRadprime))/(5._dl*x3**2._dl + 20._dl*(x1 &
              &+ 2._dl*x2 + x3) + 4._dl*(6._dl - x1 - 2._dl*x2 + 3._dl*x3)*x4 + 8._dl*x4**2._dl)

        x1pp = 2._dl*( x1p*( epsilon- h ) +x1*( epsilon_p -hprime ) )
        x2pp = 2._dl*( x2p*( 2._dl*epsilon- h ) +x2*( 2._dl*epsilon_p -hprime ) )
        x3pp = x3p*( 3._dl*epsilon- h ) +x3*( 3._dl*epsilon_p -hprime )
        x4pp = 4._dl*( x4p*epsilon +x4*epsilon_p )
        OmegaRadprimeprime = -2._dl*OmegaRadprime*(2._dl +h) -2._dl*OmegaRad*hprime

        epsilon_pp = (-2._dl*(60._dl*x1 + 40._dl*x2 - 5._dl*x3*(-3._dl + 3._dl*x1 + x2 + OmegaRad) - x4*(36._dl*x1 + 16._dl*x2 &
              &+ 3._dl*x3 + 8._dl*OmegaRad))*(-4._dl*x4*(x1p + 2._dl*x2p - 3._dl*x3p) + 10._dl*x3*x3p + 20._dl*(x1p &
              &+ 2._dl*x2p + x3p) -4._dl*(x1 + 2._dl*x2 - 3._dl*(2._dl + x3))*x4p + 16._dl*x4*x4p)**2._dl +2._dl*(5._dl*x3**2._dl &
              &+ 20._dl*(x1 + 2._dl*x2 + x3) - 4._dl*(x1 + 2._dl*x2 - 3._dl*(2._dl + x3))*x4 + 8._dl*x4**2._dl)*(-4._dl*x4*(x1p &
              &+ 2._dl*x2p - 3._dl*x3p) + 10._dl*x3*x3p + 20._dl*(x1p + 2._dl*x2p + x3p) -4._dl*(x1 + 2._dl*x2 - 3._dl*(2._dl + x3))*x4p &
              &+ 16._dl*x4*x4p)*(60._dl*x1p + 40._dl*x2p - 5._dl*(-3._dl + 3._dl*x1 + x2 + OmegaRad)*x3p - (36._dl*x1 + 16._dl*x2 &
              &+ 3._dl*x3 + 8._dl*OmegaRad)*x4p -5._dl*x3*(3._dl*x1p + x2p + OmegaRadprime) - x4*(36._dl*x1p + 16._dl*x2p + 3._dl*x3p &
              &+ 8._dl*OmegaRadprime)) +2._dl*(5._dl*x3**2._dl + 20._dl*(x1 + 2._dl*x2 + x3) - 4._dl*(x1 + 2._dl*x2 - 3._dl*(2._dl &
              &+ x3))*x4 + 8._dl*x4**2._dl)*(60._dl*x1 + 40._dl*x2 - 5._dl*x3*(-3._dl + 3._dl*x1 + x2 + OmegaRad) - x4*(36._dl*x1 &
              &+ 16._dl*x2 + 3._dl*x3 + 8._dl*OmegaRad))*(5._dl*x3p**2._dl - 4._dl*(x1p + 2._dl*x2p - 3._dl*x3p)*x4p + 8._dl*x4p**2._dl &
              &-2._dl*x4*(x1pp + 2._dl*x2pp - 3._dl*x3pp) + 5._dl*x3*x3pp + 10._dl*(x1pp + 2._dl*x2pp + x3pp) -2._dl*(x1 + 2._dl*x2 &
              &- 3._dl*(2._dl + x3))*x4pp + 8._dl*x4*x4pp) -(5._dl*x3**2._dl + 20._dl*(x1 + 2._dl*x2 + x3) - 4._dl*(x1 + 2._dl*x2 &
              &- 3._dl*(2._dl + x3))*x4 + 8._dl*x4**2._dl)**2._dl*(-10._dl*x3p*(3._dl*x1p + x2p + OmegaRadprime) -2._dl*x4p*(36._dl*x1p &
              &+ 16._dl*x2p + 3._dl*x3p + 8._dl*OmegaRadprime) + 60._dl*x1pp + 40._dl*x2pp -5._dl*(-3._dl + 3._dl*x1 + x2 &
              &+ OmegaRad)*x3pp - (36._dl*x1 + 16._dl*x2 + 3._dl*x3 + 8._dl*OmegaRad)*x4pp -5._dl*x3*(3._dl*x1pp + x2pp &
              &+ OmegaRadprimeprime) - x4*(36._dl*x1pp + 16._dl*x2pp + 3._dl*x3pp + 8._dl*OmegaRadprimeprime)))/(5._dl*x3**2._dl &
              &+ 20._dl*(x1 + 2._dl*x2 + x3) - 4._dl*(x1 + 2._dl*x2 - 3._dl*(2._dl + x3))*x4 + 8._dl*x4**2._dl)**3._dl

        x4ppp = 4._dl*( x4pp*epsilon +2._dl*X4p*epsilon_p +x4*epsilon_pp )
        
        adotoa = sqrt(H2)*a
        Hdot = adotoa**2._dl*( 1._dl+h )
        Hdotdot = 2._dl*adotoa*Hdot*( 1._dl+h )+adotoa**3._dl*hprime
                
        ! compute the background EFT functions:
        eft_cache%EFTOmegaV    = -x4/5._dl
        eft_cache%EFTOmegaP    = -x4p/(5._dl*a)
        eft_cache%EFTOmegaPP   = x4p/(5._dl*a**2._dl) - x4pp/(5._dl*a**2._dl)
        eft_cache%EFTOmegaPPP  = -x4p/(2._dl*a**3._dl)*4._dl/5._dl + (3._dl*x4pp)/(5._dl*a**3._dl) - x4ppp/(5._dl*a**3._dl)

        if (x4 == 0._dl) then
            eft_cache%EFTOmegaPPPP = 0._dl
        end if

        eft_cache%EFTc         = adotoa**2.*(h*x4p/10._dl-x4p/10._dl + x4pp/10._dl + (1._dl/10._dl)*(5._dl*(3._dl*x1+x2-epsilon*x3)-(3._dl+8._dl*epsilon)*x4) + (x1+x2+x3+x4)*3._dl/2._dl)

        eft_cache%EFTcdot      = eft_cache%EFTc*(2._dl*Hdot/adotoa-2._dl*adotoa)  + adotoa**3 * ( h*x4pp/10._dl+hprime*x4p/10._dl-x4pp/10._dl + x4ppp/10._dl + (1._dl/10._dl)*(5._dl*(3._dl*x1p+x2p-epsilon*x3p-epsilon_p*x3)  &
                               - 8._dl*epsilon_p*x4-x4p*(3._dl+8._dl*epsilon)) + (x1p+x2p+x3p+x4p)*3._dl/2._dl)

        eft_cache%EFTLambda    = adotoa**2.*((3._dl+ h)* x4p/5._dl + 6._dl*x4/5._dl -x4p/5._dl + x4pp/5._dl + 3._dl*x1+x2-epsilon*x3 -(3._dl+8._dl*epsilon)*x4/5._dl)

        eft_cache%EFTLambdadot = eft_cache%EFTLambda*(2._dl*Hdot/adotoa-2._dl*adotoa)  + adotoa**3 * ((3._dl+ h)* x4pp/5._dl+hprime* x4p/5._dl+6._dl*x4p/5._dl -x4pp/5._dl + x4ppp/5._dl + 3._dl*x1p + x2p              &
                               - epsilon_p*x3 -epsilon*x3p -(3._dl+8._dl*epsilon)*x4p/5._dl- 8._dl*epsilon_p*x4/5._dl)

        if ( x4 == 0._dl ) then

           qs = 5._dl*x3**2._dl + 20._dl*(x1 + 2._dl*x2 + x3)

           qsprime = 10._dl*(2._dl*x1p + 4._dl*x2p + (x3 + 2._dl)*x3p)

           qspprime = 10._dl*(2._dl*x1pp + 4._dl*x2pp + x3pp*x3 + 2._dl*x3pp + x3p**2._dl)

           hpprime = 1._dl/qs**3._dl*(5._dl*(qs*qspprime*(2._dl*x3*(6._dl*x1 + 3._dl*x2 + OmegaRad + 3._dl) + 2._dl*(x1 + 2._dl*x2)*(3._dl*x1 + x2 + OmegaRad + 3._dl) + 3._dl*x3**2.)                                                &
                     + 4._dl*qs*qsprime*(x3*(6._dl*x1p + 3._dl*x2p + OmegaRadprime) + (x1 + 2._dl*x2)*(3._dl*x1p + x2p + OmegaRadprime) + (3._dl*x1 + x2 + OmegaRad + 3._dl)*(x1p + 2._dl*x2p)                                        &
                     + x3p*(6._dl*x1 + 3._dl*x2 + OmegaRad + 3._dl) + 3._dl*x3*x3p) - 2._dl*qsprime**2._dl*(2._dl*x3*(6._dl*x1 + 3._dl*x2 + OmegaRad + 3._dl) + 2._dl*(x1 + 2._dl*x2)*(3._dl*x1+x2+OmegaRad+3._dl)                    &
                     + 3._dl*x3**2._dl) - 2._dl*qs**2._dl*(x3*(6._dl*x1pp + 3._dl*x2pp + OmegaRadprimeprime) + (x1 + 2._dl*x2)*(3._dl*x1pp + x2pp + OmegaRadprimeprime) + (3._dl*x1 + x2 + OmegaRad + 3._dl)*(x1pp + 2._dl*x2pp)      &
                     + 2._dl*x3p*(6*x1p + 3._dl*x2p + OmegaRadprime) + 2._dl*(x1p + 2._dl*x2p)*(3._dl*x1p + x2p + OmegaRad) + x3pp*(6._dl*x1 + 3._dl*x2 + OmegaRad + 3._dl) + 3._dl*x3*x3p)))

           Hdotdotdot = 2._dl*Hdot**2._dl*(1+h) + 2._dl * adotoa*Hdotdot*(1+h) + 5._dl*adotoa**2._dl*Hdot*hprime + adotoa**4._dl*hpprime
           x3ppp = - x3*(hpprime - 3._dl*epsilon_pp) - 2._dl*x3p*(hprime - 3._dl*epsilon_p) - (h - 3._dl*epsilon)*x3pp

           eft_cache%EFTLambdadotdot = 1._dl/3._dl *(-5._dl * adotoa**2._dl * Hdot *( x3*hprime + x3p*h - 2._dl*h*x3 - 9._dl*x1p - 3._dl*x2p - 2._dl*x3p + x3pp + 18._dl*x1 + 6._dl*x2)                     &
                                       + adotoa**4._dl *(-x3*hpprime + 4._dl*x3*hprime - 2._dl*hprime*x3p + 4._dl*h*x3p - h*x3pp - 4._dl*h*x3 - 36._dl*x1p - 12._dl*x2p - 4._dl*x3p + 9._dl*x2pp            &
                                       + 3._dl*x2pp + 4._dl*x3pp + 36._dl*x1 + 12._dl*x2 - x3ppp) + 2._dl*adotoa*Hdotdot*(-h*x3 - x3p + 9._dl*x1 + 3._dl*x2) + 2._dl*Hdot**2._dl*(                          &
                                       -h*x3 - x3p + 9._dl*x1 + 3._dl*x2))
       end if

    end subroutine EFTCAMBBeyondHorndeskiBackgroundEFTFunctions

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes the value of the second order EFT functions at a given time.
    subroutine EFTCAMBBeyondHorndeskiSecondOrderEFTFunctions( self, a, eft_par_cache, eft_cache )

        implicit none

        class(EFTCAMB_Beyond_Horndeski)              :: self          !< the base class
        real(dl), intent(in)                         :: a             !< the input scale factor
        type(TEFTCAMB_parameter_cache), intent(inout) :: eft_par_cache !< the EFTCAMB parameter cache that contains all the physical parameters.
        type(TEFTCAMB_timestep_cache ), intent(inout) :: eft_cache     !< the EFTCAMB timestep cache that contains all the physical values.

        real(dl) :: x1,x2,x3,x4, x1p,x2p,x3p,x4p, x1pp,x2pp,x4pp,x3pp,x3ppp,x4ppp,OmegaRad,OmegaRadprime
        real(dl) :: grhom, grhorad, H2, epsilon,epsilon_p,epsilon_pp,OmegaRadprimeprime
        real(dl) :: mu, x,h,hprime,hpprime,qs,qsprime,qspprime, adotoa, Hdot,Hdotdot,Hdotdotdot
        integer  :: ind

        x   = log(a)
        grhom   = eft_par_cache%grhob + eft_par_cache%grhoc
        grhorad = eft_par_cache%grhornomass +  eft_par_cache%grhog


        call self%fun_x4%precompute(x, ind, mu )
        x1    = self%fun_x1%value(x, index=ind, coeff=mu)
        x2    = self%fun_x2%value(x, index=ind, coeff=mu)
        x3    = self%fun_x3%value(x, index=ind, coeff=mu)
        x4    = self%fun_x4%value(x, index=ind, coeff=mu)

        if (x <self%x_initial) then
            x1 = self%bHor_x10*1.d-16*(a/exp(self%x_initial))**(14._dl/3._dl)
            x2 = self%bHor_x20*(a/exp(self%x_initial))**(16._dl/3._dl)
            x3 = self%bHor_x30*1.d-9*(a/exp(self%x_initial))**(3._dl)
            x4 = self%bHor_x40*1.d-6*(a/exp(self%x_initial))**(4._dl/3._dl)
        end if

        H2= (grhorad/exp(x)**4+grhom/exp(x)**3)/((1._dl-x1-x2-x3-x4)*3._dl)
        OmegaRad= (grhorad/exp(x)**4)/(3._dl*(H2))

        epsilon = - ((20._dl*(3._dl*x1 + 2._dl*x2) - 5._dl*x3*(-3._dl + 3._dl*x1 + x2 + OmegaRad) -x4*(36._dl*x1 + 16._dl*x2 + 3._dl*x3 + 8._dl*OmegaRad))/ &
                      (5._dl*x3**2._dl + 20._dl*(x1 + 2._dl*x2 + x3) + 4._dl*(6._dl- x1 - 2._dl*x2 + 3._dl*x3)*x4 + 8._dl*x4**2._dl))

        h = -((15._dl*x3**2._dl + 12._dl*x4**2._dl + 10._dl*(x1 + 2._dl*x2)*(3._dl + 3._dl*x1 + x2 + OmegaRad) &
              &+ 10._dl*x3*(3._dl + 6._dl*x1+ 3._dl*x2 + OmegaRad) + x4*(36._dl + 78._dl*x1 + 32._dl*x2 + 30._dl*x3&
              & + 12._dl*OmegaRad))/(5._dl*x3**2._dl + 20._dl*(x1 + 2._dl*x2 + x3) + 4._dl*(6._dl - x1 - 2._dl*x2 &
              &+ 3._dl*x3)*x4 + 8._dl*x4**2._dl))

        x1p = 2._dl*x1*( epsilon -h)
        x2p = 2._dl*x2*( 2._dl*epsilon -h)
        x3p = x3*( 3._dl*epsilon -h)
        x4p = 4._dl*x4*epsilon

        OmegaRadprime = -2._dl*OmegaRad*(2._dl+h)

        hprime = ((15._dl*x3**2._dl + 12._dl*x4**2._dl + 10._dl*(x1 + 2._dl*x2)*(3._dl + 3._dl*x1 + x2 + OmegaRad) + 10._dl*x3*(3._dl + 6._dl*x1 + 3._dl*x2 + OmegaRad) &
                  + x4*(36._dl + 78._dl*x1 + 32._dl*x2 + 30._dl*x3 &
                  + 12._dl*OmegaRad))*(10._dl*x3*x3p + 20._dl*(x1p + 2._dl*x2p + x3p) + 4._dl*x4*(-x1p - 2._dl*x2p + 3._dl*x3p) &
                  + 4._dl*(6._dl - x1 - 2._dl*x2 + 3._dl*x3)*x4p + 16._dl*x4*x4p))/(5._dl*x3**2._dl + 20._dl*(x1 + 2._dl*x2 + x3)&
                  + 4._dl*(6._dl - x1 - 2._dl*x2 + 3._dl*x3)*x4 + 8._dl*x4**2._dl)**2._dl -(10._dl*(3._dl + 3._dl*x1 + x2 &
                  + OmegaRad)*(x1p + 2._dl*x2p) + 30._dl*x3*x3p + 10._dl*(3._dl + 6._dl*x1 + 3._dl*x2 + OmegaRad)*x3p &
                  +24._dl*x4*x4p + (36._dl + 78._dl*x1 + 32._dl*x2 + 30._dl*x3 + 12._dl*OmegaRad)*x4p +10._dl*(x1 + 2._dl*x2)*(&
                  3._dl*x1p + x2p + OmegaRadprime) + 10._dl*x3*(6._dl*x1p + 3._dl*x2p + OmegaRadprime) +x4*(78._dl*x1p &
                  + 32._dl*x2p + 30._dl*x3p + 12._dl*OmegaRadprime))/(5._dl*x3**2. + 20._dl*(x1 + 2._dl*x2 + x3) + 4._dl*(6._dl &
                  - x1 - 2._dl*x2 + 3._dl*x3)*x4 + 8._dl*x4**2._dl)

        epsilon_p = ((20._dl*(3._dl*x1 + 2._dl*x2) - 5._dl*x3*(-3._dl +3._dl*x1 +x2 +OmegaRad) -x4*(36._dl*x1 +16._dl*x2 &
                      &+3._dl*x3 +8._dl*OmegaRad))*(10._dl*x3*x3p + 20._dl*(x1p + 2._dl*x2p + x3p) + 4._dl*x4*(-x1p - 2._dl*x2p + 3._dl*x3p) &
                      &+4._dl*(6._dl - x1 - 2._dl*x2 + 3._dl*x3)*x4p + 16._dl*x4*x4p))/(5._dl*x3**2._dl + 20._dl*(x1 + 2._dl*x2 + x3) + 4._dl*(&
                      &6._dl - x1 - 2._dl*x2 + 3._dl*x3)*x4 + 8._dl*x4**2._dl)**2._dl -(20._dl*(3._dl*x1p + 2._dl*x2p) - 5._dl*(-3._dl + 3._dl&
                      &*x1 + x2 + OmegaRad)*x3p - (36._dl*x1 + 16._dl*x2 + 3._dl*x3 + 8._dl*OmegaRad)*x4p -5._dl*x3*(3._dl*x1p + x2p &
                      &+ OmegaRadprime) - x4*(36._dl*x1p + 16._dl*x2p + 3._dl*x3p + 8._dl*OmegaRadprime))/(5._dl*x3**2._dl + 20._dl*(x1 &
                      &+ 2._dl*x2 + x3) + 4._dl*(6._dl - x1 - 2._dl*x2 + 3._dl*x3)*x4 + 8._dl*x4**2._dl)

        x1pp = 2._dl*( x1p*( epsilon- h ) +x1*( epsilon_p -hprime ) )
        x2pp = 2._dl*( x2p*( 2._dl*epsilon- h ) +x2*( 2._dl*epsilon_p -hprime ) )
        x3pp = x3p*( 3._dl*epsilon- h ) +x3*( 3._dl*epsilon_p -hprime )
        x4pp = 4._dl*( x4p*epsilon +x4*epsilon_p )
        OmegaRadprimeprime = -2._dl*OmegaRadprime*(2._dl +h) -2._dl*OmegaRad*hprime
        epsilon_pp = (-2._dl*(60._dl*x1 + 40._dl*x2 - 5._dl*x3*(-3._dl + 3._dl*x1 + x2 + OmegaRad) - x4*(36._dl*x1 + 16._dl*x2 &
                      &+ 3._dl*x3 + 8._dl*OmegaRad))*(-4._dl*x4*(x1p + 2._dl*x2p - 3._dl*x3p) + 10._dl*x3*x3p + 20._dl*(x1p &
                      &+ 2._dl*x2p + x3p) -4._dl*(x1 + 2._dl*x2 - 3._dl*(2._dl + x3))*x4p + 16._dl*x4*x4p)**2._dl +2._dl*(5._dl*x3**2._dl &
                      &+ 20._dl*(x1 + 2._dl*x2 + x3) - 4._dl*(x1 + 2._dl*x2 - 3._dl*(2._dl + x3))*x4 + 8._dl*x4**2._dl)*(-4._dl*x4*(x1p &
                      &+ 2._dl*x2p - 3._dl*x3p) + 10._dl*x3*x3p + 20._dl*(x1p + 2._dl*x2p + x3p) -4._dl*(x1 + 2._dl*x2 - 3._dl*(2._dl + x3))*x4p &
                      &+ 16._dl*x4*x4p)*(60._dl*x1p + 40._dl*x2p - 5._dl*(-3._dl + 3._dl*x1 + x2 + OmegaRad)*x3p - (36._dl*x1 + 16._dl*x2 &
                      &+ 3._dl*x3 + 8._dl*OmegaRad)*x4p -5._dl*x3*(3._dl*x1p + x2p + OmegaRadprime) - x4*(36._dl*x1p + 16._dl*x2p + 3._dl*x3p &
                      &+ 8._dl*OmegaRadprime)) +2._dl*(5._dl*x3**2._dl + 20._dl*(x1 + 2._dl*x2 + x3) - 4._dl*(x1 + 2._dl*x2 - 3._dl*(2._dl &
                      &+ x3))*x4 + 8._dl*x4**2._dl)*(60._dl*x1 + 40._dl*x2 - 5._dl*x3*(-3._dl + 3._dl*x1 + x2 + OmegaRad) - x4*(36._dl*x1 &
                      &+ 16._dl*x2 + 3._dl*x3 + 8._dl*OmegaRad))*(5._dl*x3p**2._dl - 4._dl*(x1p + 2._dl*x2p - 3._dl*x3p)*x4p + 8._dl*x4p**2._dl &
                      &-2._dl*x4*(x1pp + 2._dl*x2pp - 3._dl*x3pp) + 5._dl*x3*x3pp + 10._dl*(x1pp + 2._dl*x2pp + x3pp) -2._dl*(x1 + 2._dl*x2 &
                      &- 3._dl*(2._dl + x3))*x4pp + 8._dl*x4*x4pp) -(5._dl*x3**2._dl + 20._dl*(x1 + 2._dl*x2 + x3) - 4._dl*(x1 + 2._dl*x2 &
                      &- 3._dl*(2._dl + x3))*x4 + 8._dl*x4**2._dl)**2._dl*(-10._dl*x3p*(3._dl*x1p + x2p + OmegaRadprime) -2._dl*x4p*(36._dl*x1p &
                      &+ 16._dl*x2p + 3._dl*x3p + 8._dl*OmegaRadprime) + 60._dl*x1pp + 40._dl*x2pp -5._dl*(-3._dl + 3._dl*x1 + x2 &
                      &+ OmegaRad)*x3pp - (36._dl*x1 + 16._dl*x2 + 3._dl*x3 + 8._dl*OmegaRad)*x4pp -5._dl*x3*(3._dl*x1pp + x2pp &
                      &+ OmegaRadprimeprime) - x4*(36._dl*x1pp + 16._dl*x2pp + 3._dl*x3pp + 8._dl*OmegaRadprimeprime)))/(5._dl*x3**2._dl &
                      &+ 20._dl*(x1 + 2._dl*x2 + x3) - 4._dl*(x1 + 2._dl*x2 - 3._dl*(2._dl + x3))*x4 + 8._dl*x4**2._dl)**3._dl

        x4ppp = 4._dl*( x4pp*epsilon +2._dl*X4p*epsilon_p +x4*epsilon_pp )
        x3ppp = -x3 *(hpprime - 3._dl*epsilon_pp) - 2._dl*x3p*(hprime - 3._dl*epsilon_p) - (h - 3._dl*epsilon)*x3pp

        adotoa = sqrt(H2)*a
        Hdot = adotoa**2._dl*( 1._dl+h )
        Hdotdot = 2._dl*adotoa*Hdot*( 1._dl+h )+adotoa**3._dl*hprime

        ! compute the second order EFT functions:
        eft_cache%EFTGamma1V  = (2*adotoa**2*x2)/(a**2*eft_par_cache%h0_Mpc**2) +(adotoa**2*((3 + h/3.)*x3 + x3p/3.))/(4.*a**2*eft_par_cache%h0_Mpc**2) + &
                                (adotoa**2*(6*x4 + (3*x4p)/4. - (h*x4p)/4. -x4pp/4.))/(5.*a**2*eft_par_cache%h0_Mpc**2)
        eft_cache%EFTGamma1P  = (2*Hdot*(120*x2 + 5*(9 + h)*x3 + 72*x4 +5*x3p + 9*x4p - 3*h*x4p -3*x4pp) + adotoa**2*(-240*x2 - 144*x4 -5*x3*(18 + 2*h - hprime) &
                                 + 120*x2p +35*x3p + 5*h*x3p + 54*x4p +6*h*x4p - 3*hprime*x4p +5*x3pp + 15*x4pp - 3*h*x4pp -3*x4ppp))/(60.*a**3*eft_par_cache%h0_Mpc**2)
        eft_cache%EFTGamma2V  = -(adotoa*x3)/(a*eft_par_cache%h0_Mpc) + (adotoa*(-8*x4 + x4p))/(5.*a*eft_par_cache%h0_Mpc)
        eft_cache%EFTGamma2P  = -(Hdot*(5*x3 + 8*x4 - x4p) - adotoa**2*(5*x3 + 8*x4 - 5*x3p - 9*x4p +x4pp))/(5.*a**2*adotoa*eft_par_cache%h0_Mpc)

        if (x4 == 0._dl) then
          qs       = 5._dl*x3**2._dl + 20._dl*(x1 + 2._dl*x2 + x3)

          qsprime  = 10._dl*(2._dl*x1p + 4._dl*x2p + (x3 + 2._dl)*x3p)

          qspprime = 10._dl*(2._dl*x1pp + 4._dl*x2pp + x3pp*x3 + 2._dl*x3pp + x3p**2._dl)

          hpprime = 1._dl/qs**3._dl*(5._dl*(qs*qspprime*(2._dl*x3*(6._dl*x1 + 3._dl*x2 + OmegaRad + 3._dl) + 2._dl*(x1 + 2._dl*x2)*(3._dl*x1 + x2 + OmegaRad + 3._dl) + 3._dl*x3**2.)                                                &
                    + 4._dl*qs*qsprime*(x3*(6._dl*x1p + 3._dl*x2p + OmegaRadprime) + (x1 + 2._dl*x2)*(3._dl*x1p + x2p + OmegaRadprime) + (3._dl*x1 + x2 + OmegaRad + 3._dl)*(x1p + 2._dl*x2p)                                        &
                    + x3p*(6._dl*x1 + 3._dl*x2 + OmegaRad + 3._dl) + 3._dl*x3*x3p) - 2._dl*qsprime**2._dl*(2._dl*x3*(6._dl*x1 + 3._dl*x2 + OmegaRad + 3._dl) + 2._dl*(x1 + 2._dl*x2)*(3._dl*x1+x2+OmegaRad+3._dl)                    &
                    + 3._dl*x3**2._dl) - 2._dl*qs**2._dl*(x3*(6._dl*x1pp + 3._dl*x2pp + OmegaRadprimeprime) + (x1 + 2._dl*x2)*(3._dl*x1pp + x2pp + OmegaRadprimeprime) + (3._dl*x1 + x2 + OmegaRad + 3._dl)*(x1pp + 2._dl*x2pp)      &
                    + 2._dl*x3p*(6*x1p + 3._dl*x2p + OmegaRadprime) + 2._dl*(x1p + 2._dl*x2p)*(3._dl*x1p + x2p + OmegaRad) + x3pp*(6._dl*x1 + 3._dl*x2 + OmegaRad + 3._dl) + 3._dl*x3*x3p)))

          Hdotdotdot= 2._dl*Hdot**2._dl*(1+h) + 2._dl * adotoa*Hdotdot*(1+h) + 5._dl*adotoa**2._dl *Hdot*hprime + adotoa**4._dl*hpprime

          x3ppp = - x3*(hpprime - 3._dl*epsilon_pp) - 2._dl*x3p*(hprime - 3._dl*epsilon_p) - (h - 3._dl*epsilon)*x3pp

          eft_cache%EFTGamma2PP = (x3*Hdot**2._dl - x3*adotoa*Hdotdot + adotoa**4._dl*(3._dl*x3p-x3pp-2._dl*x3) + adotoa**2._dl*Hdot*(3._dl*x3 - 2._dl*x3p))/(a**3._dl*eft_par_cache%h0_Mpc*adotoa**3._dl)
                                  
          eft_cache%EFTGamma2PPP= (3._dl*adotoa**2._dl*Hdot**2._dl*(x3p-2._dl*x3) + 4._dl*adotoa*Hdot*Hdotdot*x3 + 3._dl*adotoa**3._dl*Hdotdot*(2._dl*x3 - x3p) - adotoa**2._dl*Hdotdotdot*x3 +                                      &
                                   adotoa**4._dl*Hdot*(-11._dl*x3+12._dl*x3p-3._dl*x3pp) + adotoa**6._dl*(6._dl*x3 - 11._dl*x3p+6._dl*x3pp - x3ppp) - Hdot**3._dl*x3)/(a**4._dl*eft_par_cache%h0_Mpc*adotoa**5._dl)                  

        end if


        eft_cache%EFTGamma3V  = 0._dl
        eft_cache%EFTGamma3P  = 0._dl
        eft_cache%EFTGamma4V  = 0._dl
        eft_cache%EFTGamma4P  = 0._dl
        eft_cache%EFTGamma4PP = 0._dl
        eft_cache%EFTGamma5V  = (2._dl*x4)/5._dl
        eft_cache%EFTGamma5P  = (2._dl*x4p)/(5._dl*a)
        eft_cache%EFTGamma6V  = 0._dl
        eft_cache%EFTGamma6P  = 0._dl



    end subroutine EFTCAMBBeyondHorndeskiSecondOrderEFTFunctions

    ! ---------------------------------------------------------------------------------------------
    !> Function that computes model specific stability requirements.
    function EFTCAMBBeyondHorndeskiAdditionalModelStability( self, a, eft_par_cache, eft_cache )

        implicit none

        class(EFTCAMB_Beyond_Horndeski)              :: self            !< the base class
        real(dl), intent(in)                         :: a               !< the input scale factor.
        type(TEFTCAMB_parameter_cache), intent(inout) :: eft_par_cache  !< the EFTCAMB parameter cache that contains all the physical parameters.
        type(TEFTCAMB_timestep_cache ), intent(inout) :: eft_cache      !< the EFTCAMB timestep cache that contains all the physical values.

        logical :: EFTCAMBBeyondHorndeskiAdditionalModelStability       !< the return value of the stability computation. True if the model specific stability criteria are met

        EFTCAMBBeyondHorndeskiAdditionalModelStability = .True.

    end function EFTCAMBBeyondHorndeskiAdditionalModelStability

end module EFTCAMB_full_Beyond_Horndeski
