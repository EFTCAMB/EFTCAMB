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

!> @file 10p6_Scaling_Cubic.f90
!! This file contains the definition of the Scaling Cubic model consistent with the equations on the mathematica file.
!! Please refer to the numerical notes for details.


!----------------------------------------------------------------------------------------
!> This module contains the definition of the Scaling Cubic model consistent with the mathematica file.
!! Please refer to the numerical notes for details.

!> @author Inês Albuquerque, Noemi Frusciante

module EFTCAMB_full_Scaling_Cubic

    use precision
    use IniObjects
    use MpiUtils
    use FileUtils
    use EFTCAMB_cache
    use EFT_def
    use EFTCAMB_abstract_model_full
    use equispaced_linear_interpolation_1D
    use EFTCAMB_abstract_parametrizations_1D
    use EFTCAMB_constant_parametrization_1D

    implicit none

    private

    public EFTCAMB_Scaling_Cubic

    !----------------------------------------------------------------------------------------
    !> This is the type that contains the definition of the Scaling Cubic model.
    type, extends ( EFTCAMB_full_model ) :: EFTCAMB_Scaling_Cubic

        ! the model parameters:
        real(dl)  :: SCG_A
        real(dl)  :: SCG_beta1
        real(dl)  :: SCG_beta2
        real(dl)  :: SCG_lambda
        real(dl)  :: SCG_x10
        real(dl)  :: SCG_y10
        real(dl)  :: SCG_y20
        real(dl)  :: SCG_x1today
        real(dl)  :: SCG_y1today
        real(dl)  :: SCG_y2today

        ! the interpolated EFT functions that come out of the background sover:
        type(equispaced_linear_interpolate_function_1D) :: fun_x1       !< The interpolated function x1
        type(equispaced_linear_interpolate_function_1D) :: fun_y1       !< The interpolated function y1
        type(equispaced_linear_interpolate_function_1D) :: fun_y2       !< The interpolated function y2

        ! some parameters for the ODE solver:
        integer  :: designer_num_points = 10000                            !< Number of points sampled by the designer code.
        real(dl) :: x_initial           = log(1._dl/(1.5*10._dl**(5._dl)+1._dl))                   !< log(a start)
        real(dl) :: x_final             = 0._dl         !< log(a final)

        ! some additional parameters
        logical :: debug_flag = .true.
        logical :: LCDM_background = .false.

    contains

        ! initialization of the model:
        procedure :: read_model_selection            => EFTCAMBScalingCubicReadModelSelectionFromFile  !< subroutine that reads the parameters of the model from file
        procedure :: allocate_model_selection        => EFTCAMBScalingCubicAllocateModelSelection      !< subroutine that allocates the model selection. For Horava this is a dummy procedure.
        procedure :: init_model_parameters           => EFTCAMBScalingCubicInitModelParameters         !< subroutine that initializes the model parameters based on the values found in an input array.
        procedure :: init_model_parameters_from_file => EFTCAMBScalingCubicInitModelParametersFromFile !< subroutine that reads the parameters of the model from file.

        ! background solver:
        procedure :: initialize_background           => EFTCAMBScalingCubicInitBackground              !< subroutine that initializes the background of the Scaling Cubic Model.
        procedure :: solve_background                => EFTCAMBScalingCubicSolveBackground             !< subroutine that solves the background equations.

        ! utility functions:
        procedure :: compute_param_number  => EFTCAMBScalingCubicComputeParametersNumber    !< subroutine that computes the number of parameters of the model.
        procedure :: feedback              => EFTCAMBScalingCubicFeedback                   !< subroutine that prints on the screen feedback information about the model.
        procedure :: parameter_names       => EFTCAMBScalingCubicParameterNames             !< subroutine that returns the i-th parameter name of the model.
        procedure :: parameter_names_latex => EFTCAMBScalingCubicParameterNamesLatex        !< subroutine that returns the i-th parameter name of the model.
        procedure :: parameter_values      => EFTCAMBScalingCubicParameterValues            !< subroutine that returns the i-th parameter value.

        ! CAMB related procedures:
        procedure :: compute_background_EFT_functions  => EFTCAMBScalingCubicBackgroundEFTFunctions   !< subroutine that computes the value of the background EFT functions at a given time.
        procedure :: compute_secondorder_EFT_functions => EFTCAMBScalingCubicSecondOrderEFTFunctions  !< subroutine that computes the value of the second order EFT functions at a given time.
        procedure :: compute_adotoa                    => EFTCAMBScalingCubicComputeAdotoa        !< subroutine that computes adotoa = H.
        procedure :: compute_H_derivs                  => EFTCAMBScalingCubicComputeHubbleDer     !< subroutine that computes the two derivatives wrt conformal time of H.

        ! stability procedures:
        procedure :: additional_model_stability        => EFTCAMBScalingCubicAdditionalModelStability !< function that computes model specific stability requirements.


    end type EFTCAMB_Scaling_Cubic

    ! ---------------------------------------------------------------------------------------------

contains

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that reads the parameters of the model from file. Nothing needs to be done
    !! but procedure present because it is deferred.
    subroutine EFTCAMBScalingCubicReadModelSelectionFromFile( self, Ini, eft_error )

        implicit none

        class(EFTCAMB_Scaling_Cubic)            :: self      !< the base class
        type(TIniFile)                          :: Ini       !< Input ini file
        integer                                 :: eft_error !< error code: 0 all fine, 1 initialization failed

        ! read precision parameters
        self%designer_num_points = Ini%Read_Int( 'model_background_num_points', 10000 )
        self%x_initial = Log( Ini%Read_Double( 'model_background_a_ini', 1._dl/(1.5*10._dl**(5._dl)+1._dl) ) )
        self%x_final = Log( Ini%Read_Double( 'model_background_a_final', 1._dl ) )
        
    end subroutine EFTCAMBScalingCubicReadModelSelectionFromFile

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that allocates the model selection. Nothing needs to be done
    !! but procedure present because it is deferred.
    subroutine EFTCAMBScalingCubicAllocateModelSelection( self, Ini, eft_error )

        implicit none

        class(EFTCAMB_Scaling_Cubic)            :: self      !< the base class
        type(TIniFile)                          :: Ini       !< Input ini file
        integer                                 :: eft_error !< error code: 0 all fine, 1 initialization failed

    end subroutine EFTCAMBScalingCubicAllocateModelSelection

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that initializes the model parameters based on the values found in an input array.
    !! Nothing needs to be done but procedure present because it is deferred.
    subroutine EFTCAMBScalingCubicInitModelParameters( self, array )

        implicit none

        class(EFTCAMB_Scaling_Cubic)                               :: self   !< the base class
        real(dl), dimension(self%parameter_number), intent(in)     :: array  !< input array with the values of the parameters of the model.

        self%SCG_A                    = array(1)
        self%SCG_beta1                = array(2)
        self%SCG_beta2                = array(3)
        self%SCG_lambda               = array(4)

    end subroutine EFTCAMBScalingCubicInitModelParameters

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that reads the parameters of the model from file. Nothing needs to be done
    !! but procedure present because it is deferred.
    subroutine EFTCAMBScalingCubicInitModelParametersFromFile( self, Ini, eft_error )

        implicit none

        class(EFTCAMB_Scaling_Cubic)     :: self      !< the base class
        type(TIniFile)                   :: Ini       !< Input ini file
        integer                          :: eft_error !< error code: 0 all fine, 1 initialization failed

        self%SCG_A                = Ini%Read_Double( 'Scaling_Cubic_A', -0.3_dl )
        self%SCG_beta1            = Ini%Read_Double( 'Scaling_Cubic_beta1', 100._dl )
        self%SCG_beta2            = Ini%Read_Double( 'Scaling_Cubic_beta2', 0.7_dl )
        self%SCG_lambda           = Ini%Read_Double( 'Scaling_Cubic_lambda', 154._dl )
        self%debug_flag           = Ini%Read_logical( 'want_debug', .false. )
        self%LCDM_background      = Ini%Read_logical( 'want_LCDM',.false. )
	    self%SCG_x10              = (2._dl*sqrt(2._dl/3._dl))/(self%SCG_beta1)
	    self%SCG_y10              = sqrt((4._dl/3._dl) + 2._dl*self%SCG_A*(self%SCG_beta1 - (4._dl/3._dl)*self%SCG_lambda))/(self%SCG_beta1)

    end subroutine EFTCAMBScalingCubicInitModelParametersFromFile

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that initializes the background of Scaling Cubic Galileon.
    subroutine EFTCAMBScalingCubicInitBackground( self, params_cache, feedback_level, success, outroot )

        implicit none

        class(EFTCAMB_Scaling_Cubic)                   :: self           !< the base class
        type(TEFTCAMB_parameter_cache), intent(inout)  :: params_cache   !< a EFTCAMB parameter cache containing cosmological parameters
	    integer                       , intent(in)     :: feedback_level !< level of feedback from the background code. 0=none; 1=some; 2=chatty.
        logical                       , intent(out)    :: success        !< wether the background initialization succeded or not
        character(LEN=*), optional    , intent(in)     :: outroot        !< the output root for the debug files
        real(dl)   :: OmegaDE_0,OmegaDE, OmegaM_0, OmegaRad_0, y2_i, y2_il, y2_ir, tol,tol_r,tol_l ,H2, H_0, grhom, grhorad, f, h, w_DE, Omega_G2
        integer :: j, ind, i
        real(dl) :: aj, xj, mu, x1,y1,y2, OmegaRad

        grhom   = params_cache%grhob + params_cache%grhoc
        grhorad = params_cache%grhornomass +  params_cache%grhog
        OmegaDE_0        = 1._dl-params_cache%omegac -params_cache%omegab -params_cache%omegag- params_cache%omegar
        OmegaRad_0=  params_cache%omegag+params_cache%omegar
        OmegaM_0= params_cache%omegac +params_cache%omegab
       ! some feedback:
        if ( feedback_level>0 ) then
            write(*,'(a)') "***************************************************************"
            write(*,'(a)') ' EFTCAMB Scaling Cubic Galileon v2. background solver'
            write(*,'(a)')
        end if

        success= .true.

        !call self%feedback()
        
        !self%SCG_x10 = (2._dl*sqrt(2._dl/3._dl))/(self%SCG_beta1)
	!self%SCG_y10 = sqrt((4._dl/3._dl) + 2._dl*self%SCG_A*(self%SCG_beta1 - (4._dl/3._dl)*self%SCG_lambda))/(self%SCG_beta1)

        ! initialize interpolating functions:
        call self%fun_x1%initialize      ( self%designer_num_points, self%x_initial, self%x_final )
        call self%fun_y1%initialize      ( self%designer_num_points, self%x_initial, self%x_final )
        call self%fun_y2%initialize      ( self%designer_num_points, self%x_initial, self%x_final )

	!SP: check omegav
	if (params_cache%omegav < 0.3_dl) then
		success = .false.
		return
	end if

        ! Initial values of y2 used in the bisection algorithm

            self%SCG_y20= 4.74*1.d-9
            y2_i= self%SCG_y20
            y2_il= 0._dl
            y2_ir= self%SCG_y20 + (1._dl/2._dl)*self%SCG_y20

            call self%solve_background( params_cache, tol, success=success )
            i=0
            do while (abs(tol)>1.d-6)

                self%SCG_y20= y2_ir
                call self%solve_background( params_cache, tol_r, success=success )

                self%SCG_y20= y2_il
                call self%solve_background( params_cache, tol_l, success=success )

                if (sign(1._dl,tol_r) == sign(1._dl,tol_l)) then
                    if (sign(1._dl,tol_r)<0._dl) then
                        write(*,*) 'Error in the bisection algorithm '
                        write(*,*) 'Negative right tolerance:', tol_r
                        write(*,*) 'Use an higher value for y2_ir'
                        write(*,*) 'Omegav = ',params_cache%omegav
                    else
                        write(*,*) 'Error in the bisection algorithm '
                        write(*,*) 'Positive left tolerance:', tol_l
                        write(*,*) 'No positive values of y2_i allow to recover Omega_{DE}.'
                        write(*,*) 'Omegav = ',params_cache%omegav
                    end if
                    call MpiStop('EFTCAMB error')
                else
                    if (sign(1._dl,tol)==sign(1._dl,tol_r)) then
                        y2_ir=y2_i
                        y2_i= (y2_il+y2_i)/2._dl
                    end if

                    if (sign(1._dl,tol)==sign(1._dl,tol_l)) then
                        y2_il= y2_i
                        y2_i= (y2_ir+y2_i)/2._dl
                    end if
                end if
                self%SCG_y20= y2_i
                call self%solve_background( params_cache, tol, success=success )
                i=i+1
            end do

            if (feedback_level>0) then
                write(*,*) '---------------------------------------------------'
                write(*,'(a, E13.3,a,I3)') 'Bisection algorithm ended successfully, with tollerance',  tol , '  and steps  ', i
                write(*,'(a)')
                write(*,'(a, E13.3)') 'y2 initial    =     ',  y2_i
                write(*,'(a)')
                write(*,*) '---------------------------------------------------'
            end if


        if (feedback_level>0) then
            xj = 0._dl
            call self%fun_x1%precompute(xj, ind,mu)
            write(*,'(a)') 'today values'
            write(*,'(a,F8.2)') 'x1 = ', self%fun_x1%value(xj, index=ind, coeff=mu)
            write(*,'(a,F8.2)') 'y1 = ', self%fun_y1%value(xj, index=ind, coeff=mu)
            write(*,'(a,F8.2)') 'y2 = ', self%fun_y2%value(xj, index=ind, coeff=mu)
        end if

        if (self%debug_flag) then

          xj = 0._dl
          call self%fun_x1%precompute(xj, ind,mu)
          OmegaDE_0 = (self%fun_x1%value(xj, index=ind, coeff=mu))**2._dl + (self%fun_y1%value(xj, index=ind, coeff=mu))**2._dl &
	   &+ (self%fun_y2%value(xj, index=ind, coeff=mu))**2._dl + 2._dl*self%fun_x1%value(xj, index=ind, coeff=mu)*self%SCG_A*(sqrt(6._dl) &
	   &- self%SCG_lambda*self%fun_x1%value(xj, index=ind, coeff=mu)) 
          H_0= sqrt((grhorad +grhom )/((1._dl-OmegaDE_0)*3._dl))
          OmegaRad_0 =  grhorad/(3._dl*H_0**2)

          open(unit=111, action = 'write', file='SCG_evolution_x1.dat')
          open(unit=222, action = 'write', file='SCG_evolution_y1.dat')
          open(unit=333, action = 'write', file='SCG_evolution_y2.dat')
          open(unit=444, action = 'write', file='SCG_evolution_Omega.dat')
          open(unit=555, action = 'write', file='SCG_specific_background_quantities.dat')
          do j =1, 100000
            aj = 1._dl/j
            xj = Log(aj)
            call self%fun_x1%precompute(xj, ind,mu)
            write(111, *)1._dl/aj-1._dl, self%fun_x1%value(xj, index=ind, coeff=mu)
            write(222, *)1._dl/aj-1._dl, self%fun_y1%value(xj, index=ind, coeff=mu)
            write(333, *)1._dl/aj-1._dl, self%fun_y2%value(xj, index=ind, coeff=mu)

            x1=self%fun_x1%value(xj, index=ind, coeff=mu)
            y1=self%fun_y1%value(xj, index=ind, coeff=mu)
            y2=self%fun_y2%value(xj, index=ind, coeff=mu)
            OmegaDE = (self%fun_x1%value(xj, index=ind, coeff=mu))**2._dl + (self%fun_y1%value(xj, index=ind, coeff=mu))**2._dl &
	                + (self%fun_y2%value(xj, index=ind, coeff=mu))**2._dl + 2._dl*self%fun_x1%value(xj, index=ind, coeff=mu)*self%SCG_A*(sqrt(6._dl)&
	                - self%SCG_lambda*self%fun_x1%value(xj, index=ind, coeff=mu))
            H2= (grhorad/aj**4+grhom/aj**3)/((1._dl-OmegaDE)*3._dl)

            OmegaRad= (grhorad/aj**4)/(3._dl*(H2))

            f = (3._dl*(-3._dl*self%SCG_A + sqrt(6._dl)*(-1._dl + 2._dl*self%SCG_A*self%SCG_lambda)*x1 + 3._dl*self%SCG_A*(1._dl &
              &- 2._dl*self%SCG_A*self%SCG_lambda)*x1**2._dl + self%SCG_A*OmegaRad - 3._dl*self%SCG_A*y2**2._dl + self%SCG_beta2*y2**2._dl &
              &- 3._dl*self%SCG_A*y1**2._dl + self%SCG_beta1*y1**2._dl))/(1._dl + 6._dl*self%SCG_A**2._dl - 2._dl*self%SCG_A*self%SCG_lambda)

            h = ((-1._dl + 2._dl*self%SCG_A*self%SCG_lambda)*OmegaRad - 3._dl*(1._dl + 12._dl*self%SCG_A**2._dl &
              &- 2._dl*self%SCG_A*self%SCG_lambda + 2._dl*sqrt(6._dl)*self%SCG_A*(1._dl - 2._dl*self%SCG_A*self%SCG_lambda)*x1 &
              &+ ((1._dl - 2._dl*self%SCG_A*self%SCG_lambda)**2._dl)*x1**2._dl + (-1._dl + 2._dl*self%SCG_A*self%SCG_lambda &
              &- 2._dl*self%SCG_A*self%SCG_beta2)*y2**2._dl - y1**2._dl + 2._dl*self%SCG_A*self%SCG_lambda*y1**2._dl &
              &- 2._dl*self%SCG_A*self%SCG_beta1*y1**2._dl))/(2._dl + 12._dl*self%SCG_A**2._dl - 4._dl*self%SCG_A*self%SCG_lambda)

            w_DE= (x1**2._dl - y1**2._dl - y2**2._dl - 2._dl*self%SCG_A*(f/3._dl + self%SCG_lambda*x1**2._dl))/(x1**2._dl + y1**2._dl + y2**2._dl + 2._dl*self%SCG_A*(sqrt(6._dl)*x1 - self%SCG_lambda*x1**2._dl))

	        Omega_G2 = x1**2._dl + y1**2._dl + y2**2._dl

            write(444, '(20E15.5)')1._dl/aj-1._dl, sqrt(H2), (grhorad/aj**4)/(3._dl*H2)&
                &, (grhom/aj**3)/(3._dl*H2), OmegaDE
            write(555, '(20E15.5)')1._dl/aj-1._dl, h, f, w_DE, Omega_G2
            
	  end do	
          close(111)
          close(222)
          close(333)
          close(444)
          close(555)

        end if


        !store energy densities today
        xj = 0._dl
        call self%fun_x1%precompute(xj, ind,mu)
        OmegaDE_0 = (self%fun_x1%value(xj, index=ind, coeff=mu))**2._dl + (self%fun_y1%value(xj, index=ind, coeff=mu))**2._dl &
	   &+ (self%fun_y2%value(xj, index=ind, coeff=mu))**2._dl + 2._dl*self%fun_x1%value(xj, index=ind, coeff=mu)*self%SCG_A*(sqrt(6._dl) &
	   &- self%SCG_lambda*self%fun_x1%value(xj, index=ind, coeff=mu)) 
        H_0= sqrt((grhorad +grhom )/((1._dl-OmegaDE_0)*3._dl))

        OmegaM_0 = grhom/(3._dl*H_0**2)

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
        self%SCG_x1today = self%fun_x1%value(xj, index=ind, coeff=mu)
        self%SCG_y1today = self%fun_y1%value(xj, index=ind, coeff=mu)
        self%SCG_y2today = self%fun_y2%value(xj, index=ind, coeff=mu)

        return

    end subroutine EFTCAMBScalingCubicInitBackground

    ! ---------------------------------------------------------------------------------------------

    subroutine EFTCAMBScalingCubicSolveBackground( self, params_cache, tol, success)

        implicit none

        class(EFTCAMB_Scaling_Cubic)                  :: self          !< the base class.
        type(TEFTCAMB_parameter_cache), intent(in)    :: params_cache  !< a EFTCAMB parameter cache containing cosmological parameters.
        real(dl), intent(out)                         :: tol           !< present day relative difference between H_0 input and H_0 computed.
        logical , intent(out)                         :: success       !< whether the calculation ended correctly or not

        real(dl) :: PPlus, yPlus, CoeffA_Part, yStar, x, OmegaDE_0 ,OmegaNu_EFT,Omegarad_EFT, grhom, grhorad

        integer, parameter :: num_eq = 3  !<  Number of equations of the dynamical system

        real(dl) :: y(num_eq), ydot(num_eq)

        integer  :: itol, itask, istate, iopt, LRN, LRS, LRW, LIS, LIN, LIW, JacobianMode, i
        real(dl) :: rtol, atol, t1, t2, B
        real(dl), allocatable :: rwork(:)
        integer,  allocatable :: iwork(:)

        ! 1) Cosmological densities:
        grhom   = params_cache%grhob + params_cache%grhoc
        grhorad = params_cache%grhornomass +  params_cache%grhog
        OmegaDE_0        = 1._dl-params_cache%omegac -params_cache%omegab -params_cache%omegag- params_cache%omegar

        ! 2) Set initial conditions:
        x    = self%x_initial
        y(1) = self%SCG_x10
        y(2) = self%SCG_y10
        y(3) = self%SCG_y20

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
                tol = (OmegaDE_0-(y(1)**2._dl + y(2)**2._dl + y(3)**2._dl + 2._dl*y(1)*self%SCG_A*(sqrt(6._dl) &
	            &- self%SCG_lambda*y(1))))/OmegaDE_0
            end if

        end do

    contains

      ! ---------------------------------------------------------------------------------------------
      !> Subroutine that computes y' given y for the background the Scaling Cubic Galileon.
      subroutine derivs( num_eq, x, y, ydot )

          implicit none

          integer , intent(in)                     :: num_eq !< number of equations in the ODE system
          real(dl), intent(in)                     :: x      !< time at which the derivatives of the system are computed
          real(dl), intent(in) , dimension(num_eq) :: y      !< input status of the system
          real(dl), intent(out), dimension(num_eq) :: ydot   !< value of the derivative at time x

          real(dl) :: f, h, H2, OmegaRad

          !1) compute the functionals entering the ODEs

          H2= (grhorad/exp(x)**4+grhom/exp(x)**3)/((1._dl - y(1)**2._dl - y(2)**2._dl - y(3)**2._dl - 2._dl*y(1)*self%SCG_A*(sqrt(6._dl) - self%SCG_lambda*y(1)))*3._dl)

          OmegaRad= (grhorad/exp(x)**4)/(3._dl*(H2))


          h = ((-1._dl + 2._dl*self%SCG_A*self%SCG_lambda)*OmegaRad - 3._dl*(1._dl + 12._dl*self%SCG_A**2._dl &
              &- 2._dl*self%SCG_A*self%SCG_lambda + 2._dl*sqrt(6._dl)*self%SCG_A*(1._dl - 2._dl*self%SCG_A*self%SCG_lambda)*y(1) &
              &+ ((1._dl - 2._dl*self%SCG_A*self%SCG_lambda)**2._dl)*y(1)**2._dl + (-1._dl + 2._dl*self%SCG_A*self%SCG_lambda &
              &- 2._dl*self%SCG_A*self%SCG_beta2)*y(3)**2._dl - y(2)**2._dl + 2._dl*self%SCG_A*self%SCG_lambda*y(2)**2._dl &
              &- 2._dl*self%SCG_A*self%SCG_beta1*y(2)**2._dl))/(2._dl + 12._dl*self%SCG_A**2._dl - 4._dl*self%SCG_A*self%SCG_lambda)


          f = (3._dl*(-3._dl*self%SCG_A + sqrt(6._dl)*(-1._dl + 2._dl*self%SCG_A*self%SCG_lambda)*y(1) + 3._dl*self%SCG_A*(1._dl &
              &- 2._dl*self%SCG_A*self%SCG_lambda)*y(1)**2._dl + self%SCG_A*OmegaRad - 3._dl*self%SCG_A*y(3)**2._dl &
              &+ self%SCG_beta2*y(3)**2._dl - 3._dl*self%SCG_A*y(2)**2._dl + self%SCG_beta1*y(2)**2._dl))/(1._dl + 6._dl*self%SCG_A**2._dl - 2._dl*self%SCG_A*self%SCG_lambda)

          ! 2) Get the equations of motion:
          ydot(1) = (1._dl/sqrt(6._dl))*f - y(1)*h
          ydot(2) = -sqrt(3._dl/2._dl)*self%SCG_beta1*y(1)*y(2) - y(2)*h
          ydot(3) = -sqrt(3._dl/2._dl)*self%SCG_beta2*y(1)*y(3) - y(3)*h

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
      !> Subroutine that takes the solution of the background equations and stores the values of
      !! the EFT functions.
      subroutine output( num_eq, ind, x, y )

          implicit none

          integer , intent(in)                     :: num_eq !< number of equations in the ODE system.
          integer , intent(in)                     :: ind    !< index of the EFT functions interpolation tables to fill.
          real(dl), intent(in)                     :: x      !< time at which the derivatives of the system are computed.
          real(dl), intent(in) , dimension(num_eq) :: y      !< input status of the system.

          real(dl) :: f, h, OmegaDE, fprime, hprime
          real(dl) :: ydot(num_eq),alpha_K, alpha_B, Q_s, Q_t, Sigma, cs2

          real(dl) :: x1,y1,y2,x1p,y1p,y2p,x1pp,OmegaRad,OmegaRadprime

          real(dl) :: a , H2,adotoa,Hdot,temp

          ! 1) call derivs to compute ydot at a given time:
          call derivs( num_eq, x, y, ydot )

          ! 2) store some quantities
          x1 = y(1)
          y1 = y(2)
          y2 = y(3)

          x1p = ydot(1)
          y1p = ydot(2)
          y2p = ydot(3)

          H2= (grhorad/exp(x)**4+grhom/exp(x)**3)/((1._dl - y(1)**2._dl - y(2)**2._dl - y(3)**2._dl - 2._dl*y(1)*self%SCG_A*(sqrt(6._dl) - self%SCG_lambda*y(1)))*3._dl)

          OmegaRad= (grhorad/exp(x)**4)/(3._dl*(H2))

          ! 3) compute the functionals entering the ODEs

          f = (3._dl*(-3._dl*self%SCG_A + sqrt(6._dl)*(-1._dl + 2._dl*self%SCG_A*self%SCG_lambda)*x1 + 3._dl*self%SCG_A*(1._dl &
              &- 2._dl*self%SCG_A*self%SCG_lambda)*x1**2._dl + self%SCG_A*OmegaRad - 3._dl*self%SCG_A*y2**2._dl + self%SCG_beta2*y2**2._dl &
              &- 3._dl*self%SCG_A*y1**2._dl + self%SCG_beta1*y1**2._dl))/(1._dl + 6._dl*self%SCG_A**2._dl - 2._dl*self%SCG_A*self%SCG_lambda)

          h = ((-1._dl + 2._dl*self%SCG_A*self%SCG_lambda)*OmegaRad - 3._dl*(1._dl + 12._dl*self%SCG_A**2._dl &
              &- 2._dl*self%SCG_A*self%SCG_lambda + 2._dl*sqrt(6._dl)*self%SCG_A*(1._dl - 2._dl*self%SCG_A*self%SCG_lambda)*x1 &
              &+ ((1._dl - 2._dl*self%SCG_A*self%SCG_lambda)**2._dl)*x1**2._dl + (-1._dl + 2._dl*self%SCG_A*self%SCG_lambda &
              &- 2._dl*self%SCG_A*self%SCG_beta2)*y2**2._dl - y1**2._dl + 2._dl*self%SCG_A*self%SCG_lambda*y1**2._dl &
              &- 2._dl*self%SCG_A*self%SCG_beta1*y1**2._dl))/(2._dl + 12._dl*self%SCG_A**2._dl - 4._dl*self%SCG_A*self%SCG_lambda)

          OmegaRadprime = -2._dl*OmegaRad*(2._dl+h)

          fprime = (-3._dl*(-self%SCG_A*OmegaRadprime + (-1._dl + 2._dl*self%SCG_A*self%SCG_lambda)*(-sqrt(6._dl) + 6._dl*self%SCG_A*x1)*x1p &
                 &+ 6._dl*self%SCG_A*y2*y2p - 2._dl*self%SCG_beta2*y2*y2p + 6._dl*self%SCG_A*y1*y1p - 2._dl*self%SCG_beta1*y1*y1p))/(1._dl &
                 &+ 6._dl*self%SCG_A**2._dl - 2._dl*self%SCG_A*self%SCG_lambda)

	      hprime = ((1._dl/2._dl)*OmegaRadprime*(-1._dl + 2._dl*self%SCG_A*self%SCG_lambda) + 3._dl*((-1._dl &
                 &+ 2._dl*self%SCG_A*self%SCG_lambda)*(sqrt(6._dl)*self%SCG_A + x1 - 2._dl*self%SCG_A*self%SCG_lambda*x1)*x1p + (1._dl &
                 &- 2._dl*self%SCG_A*self%SCG_lambda + 2._dl*self%SCG_A*self%SCG_beta2)*y2*y2p + (1._dl - 2._dl*self%SCG_A*self%SCG_lambda &
                 &+ 2._dl*self%SCG_A*self%SCG_beta1)*y1*y1p))/(1._dl + 6._dl*self%SCG_A**2._dl - 2._dl*self%SCG_A*self%SCG_lambda)

          x1pp = (1._dl/sqrt(6._dl))*fprime - hprime*x1 - h*x1p         


          ! 3) compute the x_i functions: more derivatives needed
          self%fun_x1%y(ind) = y(1)
          self%fun_x1%yp(ind) = ydot(1)
          self%fun_x1%ypp(ind) = x1pp
          self%fun_y1%y(ind) = y(2)
          self%fun_y1%yp(ind) = ydot(2)
          self%fun_y2%y(ind) = y(3)
          self%fun_y2%yp(ind) = ydot(3)

      end subroutine

    end subroutine EFTCAMBScalingCubicSolveBackground

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes the number of parameters of the model.
    subroutine EFTCAMBScalingCubicComputeParametersNumber( self )

        implicit none

        class(EFTCAMB_Scaling_Cubic)  :: self   !< the base class

        self%parameter_number = 4

    end subroutine EFTCAMBScalingCubicComputeParametersNumber

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that prints on the screen feedback information about the model.
    subroutine EFTCAMBScalingCubicFeedback( self, print_params )

        implicit none

        class(EFTCAMB_Scaling_Cubic)       :: self         !< the base class
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
            write(*,'(a,E13.3)')  '   A              ='  , self%SCG_A
            write(*,'(a,E13.3)')  '   beta1          ='  , self%SCG_beta1
            write(*,'(a,E13.3)')  '   beta2          ='  , self%SCG_beta2
	    write(*,'(a,E13.3)')  '   lambda         ='  , self%SCG_lambda
	    write(*,'(a,E13.3)')  '   x1_0           ='  , self%SCG_x10
	    write(*,'(a,E13.3)')  '   y1_0           ='  , self%SCG_y10
        end if

        ! print the values of the parameters:
        if ( present(print_params) ) then
            print_params_temp = print_params
        else
            print_params_temp = .True.
        end if

    end subroutine EFTCAMBScalingCubicFeedback

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the i-th parameter name of the model
    subroutine EFTCAMBScalingCubicParameterNames( self, i, name )

        implicit none

        class(EFTCAMB_Scaling_Cubic)      :: self   !< the base class
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
            name = 'Scaling_Cubic_A'
            return
        end if
        if ( i==2 ) then
            name = 'Scaling_Cubic_beta1'
            return
        end if
        if ( i==3 ) then
            name = 'Scaling_Cubic_beta2'
            return
        end if
        if ( i==4 ) then
            name = 'Scaling_Cubic_lambda'
            return
        end if
        if ( i==0 ) then
            name = 'no_name'
            return
        end if

    end subroutine EFTCAMBScalingCubicParameterNames

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the i-th parameter name of the model
    subroutine EFTCAMBScalingCubicParameterNamesLatex( self, i, latexname )

        implicit none

        class(EFTCAMB_Scaling_Cubic)      :: self       !< the base class
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
            latexname = 'A'
            return
        end if
        if ( i==2 ) then
            latexname = 'beta_1'
            return
        end if
        if ( i==3 ) then
            latexname = 'beta_2'
            return
        end if
        if ( i==4 ) then
            latexname = 'lambda'
            return
        end if
        if ( i==0 ) then
            latexname = 'noname'
            return
        end if

    end subroutine EFTCAMBScalingCubicParameterNamesLatex

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the i-th parameter value of the model
    subroutine EFTCAMBScalingCubicParameterValues( self, i, value )

        implicit none

        class(EFTCAMB_Scaling_Cubic)      :: self   !< the base class
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
            value = self%SCG_A
            return
        end if
        if ( i==2 ) then
            value = self%SCG_beta1
            return
        end if
        if ( i==3 ) then
            value = self%SCG_beta2
            return
        end if
        if ( i==4 ) then
            value = self%SCG_lambda
            return
        end if

        if ( i==-1 ) then
            value = self%SCG_x10
            return
        end if
        if ( i==-2 ) then
            value = self%SCG_y10
            return
        end if
        if ( i==-3 ) then
            value = self%SCG_y20
            return
        end if
        if ( i==-4 ) then
            value = self%SCG_x1today
            return
        end if
        if ( i==-5 ) then
            value = self%SCG_y1today
            return
        end if 
        if ( i==-6 ) then
            value = self%SCG_y2today
            return
        end if           


    end subroutine EFTCAMBScalingCubicParameterValues

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes the value of the background EFT functions at a given time.
    subroutine EFTCAMBScalingCubicBackgroundEFTFunctions( self, a, eft_par_cache, eft_cache )

        implicit none

        class(EFTCAMB_Scaling_Cubic)                  :: self          !< the base class
        real(dl), intent(in)                          :: a             !< the input scale factor
        type(TEFTCAMB_parameter_cache), intent(inout) :: eft_par_cache !< the EFTCAMB parameter cache that contains all the physical parameters.
        type(TEFTCAMB_timestep_cache ), intent(inout) :: eft_cache     !< the EFTCAMB timestep cache that contains all the physical values.

        integer  :: ind

        real(dl) :: x1, x1p, x1pp, x1ppp, y1, y1p, y1pp, y2, y2p, y2pp, OmegaRad, OmegaRadprime, OmegaRadpprime
        real(dl) :: f, fprime, fpprime, h, hprime, hpprime, grhom, grhorad, H2, EFTcdot, EFTLAMBDAdot, EFTc, EFTLAMBDA, adotoa, Hdot
        real(dl) :: x, mu

        if(a==0._dl) return

        x   = log(a)

        grhom   = eft_par_cache%grhob + eft_par_cache%grhoc
        grhorad = eft_par_cache%grhornomass +  eft_par_cache%grhog

        call self%fun_x1%precompute(x, ind, mu )
        x1    = self%fun_x1%value(x, index=ind, coeff=mu)
        y1    = self%fun_y1%value(x, index=ind, coeff=mu)
        y2    = self%fun_y2%value(x, index=ind, coeff=mu)

        H2= (grhorad/exp(x)**4+grhom/exp(x)**3)/((1._dl - x1**2._dl - y1**2._dl - y2**2._dl - 2._dl*self%SCG_A*x1*(sqrt(6._dl) - self%SCG_lambda*x1))*3._dl)
        OmegaRad= (grhorad/exp(x)**4)/(3._dl*(H2))

        f = (3._dl*(-3._dl*self%SCG_A + sqrt(6._dl)*(-1._dl + 2._dl*self%SCG_A*self%SCG_lambda)*x1 + 3._dl*self%SCG_A*(1._dl &
              &- 2._dl*self%SCG_A*self%SCG_lambda)*x1**2._dl + self%SCG_A*OmegaRad - 3._dl*self%SCG_A*y2**2._dl + self%SCG_beta2*y2**2._dl &
              &- 3._dl*self%SCG_A*y1**2._dl + self%SCG_beta1*y1**2._dl))/(1._dl + 6._dl*self%SCG_A**2._dl - 2._dl*self%SCG_A*self%SCG_lambda)


        h = ((-1._dl + 2._dl*self%SCG_A*self%SCG_lambda)*OmegaRad - 3._dl*(1._dl + 12._dl*self%SCG_A**2._dl &
              &- 2._dl*self%SCG_A*self%SCG_lambda + 2._dl*sqrt(6._dl)*self%SCG_A*(1._dl - 2._dl*self%SCG_A*self%SCG_lambda)*x1 &
              &+ ((1._dl - 2._dl*self%SCG_A*self%SCG_lambda)**2._dl)*x1**2._dl + (-1._dl + 2._dl*self%SCG_A*self%SCG_lambda &
              &- 2._dl*self%SCG_A*self%SCG_beta2)*y2**2._dl - y1**2._dl + 2._dl*self%SCG_A*self%SCG_lambda*y1**2._dl &
              &- 2._dl*self%SCG_A*self%SCG_beta1*y1**2._dl))/(2._dl + 12._dl*self%SCG_A**2._dl - 4._dl*self%SCG_A*self%SCG_lambda)

        ! first derivative 

        x1p = (1._dl/sqrt(6._dl))*f - x1*h
        y1p = -sqrt(3._dl/2._dl)*self%SCG_beta1*x1*y1 - y1*h
        y2p = -sqrt(3._dl/2._dl)*self%SCG_beta2*x1*y2 - y2*h

        OmegaRadprime = -2._dl*OmegaRad*(2._dl+h)

	    fprime = (-3._dl*(-self%SCG_A*OmegaRadprime + (-1._dl + 2._dl*self%SCG_A*self%SCG_lambda)*(-sqrt(6._dl) + 6._dl*self%SCG_A*x1)*x1p &
                 &+ 6._dl*self%SCG_A*y2*y2p - 2._dl*self%SCG_beta2*y2*y2p + 6._dl*self%SCG_A*y1*y1p - 2._dl*self%SCG_beta1*y1*y1p))/(1._dl &
                 &+ 6._dl*self%SCG_A**2._dl - 2._dl*self%SCG_A*self%SCG_lambda)

	    hprime = ((1._dl/2._dl)*OmegaRadprime*(-1._dl + 2._dl*self%SCG_A*self%SCG_lambda) + 3._dl*((-1._dl &
                 &+ 2._dl*self%SCG_A*self%SCG_lambda)*(sqrt(6._dl)*self%SCG_A + x1 - 2._dl*self%SCG_A*self%SCG_lambda*x1)*x1p + (1._dl &
                 &- 2._dl*self%SCG_A*self%SCG_lambda + 2._dl*self%SCG_A*self%SCG_beta2)*y2*y2p + (1._dl - 2._dl*self%SCG_A*self%SCG_lambda &
                 &+ 2._dl*self%SCG_A*self%SCG_beta1)*y1*y1p))/(1._dl + 6._dl*self%SCG_A**2._dl - 2._dl*self%SCG_A*self%SCG_lambda)

        ! second derivatives

        x1pp  = (1._dl/sqrt(6._dl))*fprime - hprime*x1 - h*x1p
        y1pp  = -sqrt(1.5_dl) *self%SCG_beta1*(x1p*y1 + x1*y1p) - hprime*y1 - h*y1p
        y2pp  = -sqrt(1.5_dl) *self%SCG_beta2*(x1p*y2 + x1*y2p) - hprime*y2 - h*y2p

        OmegaRadpprime = -2._dl*OmegaRadprime*(2._dl+h) - 2._dl*OmegaRad*hprime

        fpprime = -3._dl/(1._dl + 2._dl*self%SCG_A * (3._dl *self%SCG_A - self%SCG_lambda) ) * (-self%SCG_A*OmegaRadpprime + 6._dl*self%SCG_A*(2._dl*self%SCG_A*self%SCG_lambda - 1._dl)*x1p**2._dl   &
                  +2._dl *(3._dl*self%SCG_A - self%SCG_beta1) * (y1p**2._dl + y1*y1pp) + 2._dl*(3._dl*self%SCG_A - self%SCG_beta2 )*(y2p**2._dl + y2*y2pp)  )

        hpprime = - 1._dl/(1._dl + 2._dl*self%SCG_A*(3._dl *self%SCG_A - self%SCG_lambda)) * (0.5_dl*OmegaRadpprime*(2._dl*self%SCG_A*self%SCG_lambda - 1._dl) + 3._dl *( - (1 - 2._dl*self%SCG_A*self%SCG_lambda)**2._dl*x1p**2._dl &
                  + (2._dl*self%SCG_A*self%SCG_lambda - 1._dl)*(sqrt(6._dl)*self%SCG_A + x1 - 2._dl*self%SCG_A*self%SCG_lambda*x1)*x1pp + (1._dl-2._dl*self%SCG_A*self%SCG_lambda+2._dl*self%SCG_A*self%SCG_beta1)*(y1p**2._dl+y1*y1pp)  &
                  + (1._dl + 2._dl*self%SCG_A*self%SCG_lambda + 2._dl*self%SCG_A*self%SCG_beta2)*(y2p**2._dl + y2*y2pp) ) )


        ! third derivative

        x1ppp = sqrt(1._dl/6._dl)*fpprime - hpprime*x1 - h*x1pp - 2._dl*hprime*x1p

        adotoa=sqrt(H2)*a
        Hdot=adotoa**2._dl*( 1._dl+h )

        ! compute the background EFT functions:
        eft_cache%EFTOmegaV    = 0._dl
        eft_cache%EFTOmegaP    = 0._dl
        eft_cache%EFTOmegaPP   = 0._dl
        eft_cache%EFTOmegaPPP  = 0._dl
        eft_cache%EFTOmegaPPPP = 0._dl

        eft_cache%EFTc         = 3._dl*(adotoa**2._dl)*x1**2._dl + self%SCG_A*(adotoa**2._dl)*(3._dl*sqrt(6._dl)*x1 - sqrt(6._dl)*(h*x1 + x1p) - 6._dl*self%SCG_lambda*x1**2._dl)

        eft_cache%EFTcdot      = - 1._dl*(adotoa**3._dl)*(2._dl*sqrt(6._dl)*self%SCG_A*x1*h**2._dl+ h*(-6._dl*sqrt(6._dl)*self%SCG_A*x1 + 6._dl*(-1._dl + 2._dl*self%SCG_A*self%SCG_lambda)*x1**2._dl &
                                 + 3._dl*sqrt(6._dl)*self%SCG_A*x1p) + x1*(sqrt(6._dl)*self%SCG_A*hprime + 6._dl*(-1._dl + 2._dl*self%SCG_A*self%SCG_lambda)*x1p) + sqrt(6._dl)*self%SCG_A*(-3._dl*x1p + x1pp))

        eft_cache%EFTLambda    = 3._dl*(adotoa**2._dl)*(x1**2._dl - y1**2._dl - y2**2._dl) - self%SCG_A*(adotoa**2._dl)*(6._dl*self%SCG_lambda*x1**2._dl + 2._dl*sqrt(6._dl)*(x1*h + x1p))

        eft_cache%EFTLambdadot = -2._dl*(adotoa**3._dl)*(2._dl*sqrt(6._dl)*self%SCG_A*x1*h**2._dl + 3._dl*h*((-1._dl + 2._dl*self%SCG_A*self%SCG_lambda)*x1**2._dl + y1**2._dl + y2**2._dl &
                                 + sqrt(6._dl)*self%SCG_A*x1p) + x1*(sqrt(6._dl)*self%SCG_A*hprime + 3._dl*(-1._dl + 2._dl*self%SCG_A*self%SCG_lambda)*x1p) + 3._dl*y1*y1p + 3._dl*y2*y2p + sqrt(6._dl)*self%SCG_A*x1pp)
	    
	    eft_cache%EFTLambdadotdot = - 2._dl*adotoa**4._dl * (6._dl*sqrt(6._dl)*self%SCG_A*x1*h**3._dl + 3._dl*hprime*(y1**2._dl+y2**2._dl + x1**2.*(2._dl*self%SCG_A*self%SCG_lambda - 1._dl) )                   & 
                                    + 4._dl*sqrt(6._dl)*self%SCG_A*x1p*hprime + 3._dl*(2._dl*self%SCG_A*self%SCG_lambda - 1._dl) + self%SCG_A*sqrt(6._dl)*(x1pp + x1ppp)                                          & 
                                    + h**2._dl*(2._dl*sqrt(6._dl)*self%SCG_A*x1 + 9._dl*(y1**2._dl + y2**2._dl) + 9._dl*(2._dl*self%SCG_A*self%SCG_lambda - 1._dl)*x1**2. + 11._dl*sqrt(6._dl)*x1p*self%SCG_A)    &
                                    + x1p*(sqrt(6._dl)*self%SCG_A*hprime + 3._dl*(x1p + x1pp)*(2._dl*self%SCG_A*self%SCG_lambda - 1._dl) + sqrt(6._dl)*self%SCG_A*hpprime)                                        &
                                    + 3._dl*(y1*(y1p + y1pp) + y1p**2._dl + y2*(y2p + y2pp) + y2p**2._dl)                                                                                                         &
                                    + h*(3._dl*x1**2._dl*(2._dl*self%SCG_A*self%SCG_lambda - 1._dl) + x1*(15._dl*(2._dl*self%SCG_A*self%SCG_lambda - 1._dl) + 7._dl*sqrt(6._dl )*self%SCG_A*hprime )              & 
                                    + 3._dl*sqrt(6._dl)*a*(x1p + 2._dl*x1pp) + 3._dl*(y1**2._dl + y2**2._dl + 5._dl*y1*y1p + 5._dl*y2*y2p)))
    end subroutine EFTCAMBScalingCubicBackgroundEFTFunctions

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes the value of the second order EFT functions at a given time.
    subroutine EFTCAMBScalingCubicSecondOrderEFTFunctions( self, a, eft_par_cache, eft_cache )

        implicit none

        class(EFTCAMB_Scaling_Cubic)                  :: self          !< the base class
        real(dl), intent(in)                          :: a             !< the input scale factor
        type(TEFTCAMB_parameter_cache), intent(inout) :: eft_par_cache !< the EFTCAMB parameter cache that contains all the physical parameters.
        type(TEFTCAMB_timestep_cache ), intent(inout) :: eft_cache     !< the EFTCAMB timestep cache that contains all the physical values.

        real(dl) :: x1, x1p, x1pp, x1ppp, y1, y1p, y1pp, y2, y2p, y2pp, OmegaRad, OmegaRadprime, OmegaRadpprime
        real(dl) :: f, fprime, fpprime, h, hprime, hpprime, grhom, grhorad, H2, EFTcdot, EFTLAMBDAdot, EFTc, EFTLAMBDA, adotoa, Hdot
        real(dl):: mu, x
        integer  :: ind

        x   = log(a)
        grhom   = eft_par_cache%grhob + eft_par_cache%grhoc
        grhorad = eft_par_cache%grhornomass +  eft_par_cache%grhog


        call self%fun_x1%precompute(x, ind, mu )
        x1    = self%fun_x1%value(x, index=ind, coeff=mu)
        y1    = self%fun_y1%value(x, index=ind, coeff=mu)
        y2    = self%fun_y2%value(x, index=ind, coeff=mu)

        H2= (grhorad/exp(x)**4+grhom/exp(x)**3)/((1._dl - x1**2._dl - y1**2._dl - y2**2._dl - 2._dl*self%SCG_A*x1*(sqrt(6._dl) - self%SCG_lambda*x1))*3._dl)
        
        OmegaRad= (grhorad/exp(x)**4)/(3._dl*(H2))

        f = (3._dl*(-3._dl*self%SCG_A + sqrt(6._dl)*(-1._dl + 2._dl*self%SCG_A*self%SCG_lambda)*x1 + 3._dl*self%SCG_A*(1._dl &
              &- 2._dl*self%SCG_A*self%SCG_lambda)*x1**2._dl + self%SCG_A*OmegaRad - 3._dl*self%SCG_A*y2**2._dl + self%SCG_beta2*y2**2._dl &
              &- 3._dl*self%SCG_A*y1**2._dl + self%SCG_beta1*y1**2._dl))/(1._dl + 6._dl*self%SCG_A**2._dl - 2._dl*self%SCG_A*self%SCG_lambda)

        h = ((-1._dl + 2._dl*self%SCG_A*self%SCG_lambda)*OmegaRad - 3._dl*(1._dl + 12._dl*self%SCG_A**2._dl &
              &- 2._dl*self%SCG_A*self%SCG_lambda + 2._dl*sqrt(6._dl)*self%SCG_A*(1._dl - 2._dl*self%SCG_A*self%SCG_lambda)*x1 &
              &+ ((1._dl - 2._dl*self%SCG_A*self%SCG_lambda)**2._dl)*x1**2._dl + (-1._dl + 2._dl*self%SCG_A*self%SCG_lambda &
              &- 2._dl*self%SCG_A*self%SCG_beta2)*y2**2._dl - y1**2._dl + 2._dl*self%SCG_A*self%SCG_lambda*y1**2._dl &
              &- 2._dl*self%SCG_A*self%SCG_beta1*y1**2._dl))/(2._dl + 12._dl*self%SCG_A**2._dl - 4._dl*self%SCG_A*self%SCG_lambda)

        x1p = (1._dl/sqrt(6._dl))*f - x1*h
        y1p = -sqrt(3._dl/2._dl)*self%SCG_beta1*x1*y1 - y1*h
        y2p = -sqrt(3._dl/2._dl)*self%SCG_beta2*x1*y2 - y2*h

        OmegaRadprime = -2._dl*OmegaRad*(2._dl+h)

        fprime = (-3._dl*(-self%SCG_A*OmegaRadprime + (-1._dl + 2._dl*self%SCG_A*self%SCG_lambda)*(-sqrt(6._dl) + 6._dl*self%SCG_A*x1)*x1p &
                 &+ 6._dl*self%SCG_A*y2*y2p - 2._dl*self%SCG_beta2*y2*y2p + 6._dl*self%SCG_A*y1*y1p - 2._dl*self%SCG_beta1*y1*y1p))/(1._dl &
                 &+ 6._dl*self%SCG_A**2._dl - 2._dl*self%SCG_A*self%SCG_lambda)

	    hprime = ((1._dl/2._dl)*OmegaRadprime*(-1._dl + 2._dl*self%SCG_A*self%SCG_lambda) + 3._dl*((-1._dl &
                 &+ 2._dl*self%SCG_A*self%SCG_lambda)*(sqrt(6._dl)*self%SCG_A + x1 - 2._dl*self%SCG_A*self%SCG_lambda*x1)*x1p + (1._dl &
                 &- 2._dl*self%SCG_A*self%SCG_lambda + 2._dl*self%SCG_A*self%SCG_beta2)*y2*y2p + (1._dl - 2._dl*self%SCG_A*self%SCG_lambda &
                 &+ 2._dl*self%SCG_A*self%SCG_beta1)*y1*y1p))/(1._dl + 6._dl*self%SCG_A**2._dl - 2._dl*self%SCG_A*self%SCG_lambda)

        ! second derivatives

        x1pp  = (1._dl/sqrt(6._dl))*fprime - hprime*x1 - h*x1p
        y1pp  = -sqrt(1.5_dl) *self%SCG_beta1*(x1p*y1 + x1*y1p) - hprime*y1 - h*y1p
        y2pp  = -sqrt(1.5_dl) *self%SCG_beta2*(x1p*y2 + x1*y2p) - hprime*y2 - h*y2p

        OmegaRadpprime = -2._dl*OmegaRadprime*(2._dl+h) - 2._dl*OmegaRad*hprime

        fpprime = -3._dl/(1._dl + 2._dl*self%SCG_A * (3._dl *self%SCG_A - self%SCG_lambda) ) * (-self%SCG_A*OmegaRadpprime + 6._dl*self%SCG_A*(2._dl*self%SCG_A*self%SCG_lambda - 1._dl)*x1p**2._dl   &
                  +2._dl *(3._dl*self%SCG_A - self%SCG_beta1) * (y1p**2._dl + y1*y1pp) + 2._dl*(3._dl*self%SCG_A - self%SCG_beta2 )*(y2p**2._dl + y2*y2pp)  )

        hpprime = - 1._dl/(1._dl + 2._dl*self%SCG_A*(3._dl *self%SCG_A - self%SCG_lambda)) * (0.5_dl*OmegaRadpprime*(2._dl*self%SCG_A*self%SCG_lambda - 1._dl) + 3._dl *( - (1 - 2._dl*self%SCG_A*self%SCG_lambda)**2._dl*x1p**2._dl &
                  + (2._dl*self%SCG_A*self%SCG_lambda -1._dl)*(sqrt(6._dl)*self%SCG_A + x1 - 2._dl*self%SCG_A*self%SCG_lambda*x1)*x1pp + (1._dl-2._dl*self%SCG_A*self%SCG_lambda+2._dl*self%SCG_A*self%SCG_beta1)*(y1p**2._dl+y1*y1pp) &
                  + (1._dl + 2._dl*self%SCG_A*self%SCG_lambda + 2._dl*self%SCG_A*self%SCG_beta2)*(y2p**2._dl + y2*y2pp) ) )


        ! third derivative

        x1ppp = sqrt(1._dl/6._dl)*fpprime - hpprime*x1 - h*x1pp - 2._dl*hprime*x1p
          
        adotoa=sqrt(H2)*a
        Hdot=adotoa**2._dl*( 1._dl+h )

        ! compute the second order EFT functions:
        eft_cache%EFTGamma1V    = (sqrt(3._dl/2._dl)*self%SCG_A*(adotoa**2._dl)*(-3._dl*x1 + h*x1+ x1p))/((a**2._dl)*(eft_par_cache%h0_Mpc**2._dl))
        eft_cache%EFTGamma1P    = (sqrt(3._dl/2._dl)*self%SCG_A*(adotoa**2._dl)*(2._dl*x1*h**2._dl + x1*hprime - 3._dl*x1p + h*(-6._dl*x1 + 3._dl*x1p) + x1pp))/((a**3._dl)*(eft_par_cache%h0_Mpc**2._dl))
        eft_cache%EFTGamma2V    = (-2._dl*sqrt(6._dl)*self%SCG_A*adotoa*x1)/(a*eft_par_cache%h0_Mpc)
        eft_cache%EFTGamma2P    = (-2._dl*sqrt(6._dl)*self%SCG_A*adotoa*(h*x1 + x1p))/((a**2._dl)*eft_par_cache%h0_Mpc)
        eft_cache%EFTGamma2PP   = - 2._dl*sqrt(6._dl)*self%SCG_A*adotoa/(eft_par_cache%h0_Mpc*a**3._dl)*(-x1p - h*x1 + 2._dl*x1p*h + hprime*x1 + h**2._dl*x1 + x1pp)
        eft_cache%EFTGamma2PPP  = - 2._dl*sqrt(6._dl)*self%SCG_A*adotoa/(eft_par_cache%h0_Mpc*a**4._dl)*(x1 * (2._dl*h - 3._dl*hprime - 3._dl*h**2._dl + 3._dl*h*hprime + h**3._dl + hpprime) &
                                  + x1p*(2._dl - 6._dl*h + 3._dl*h*2._dl + 3._dl*hprime) + x1pp*(3._dl*h - 3._dl) +x1ppp )
        eft_cache%EFTGamma3V    = 0._dl
        eft_cache%EFTGamma3P    = 0._dl
        eft_cache%EFTGamma3PP   = 0._dl
        eft_cache%EFTGamma3PPP  = 0._dl
        eft_cache%EFTGamma3PPPP = 0._dl
        eft_cache%EFTGamma4V    = 0._dl
        eft_cache%EFTGamma4P    = 0._dl
        eft_cache%EFTGamma4PP   = 0._dl
        eft_cache%EFTGamma5V    = 0._dl
        eft_cache%EFTGamma5P    = 0._dl
        eft_cache%EFTGamma6V    = 0._dl
        eft_cache%EFTGamma6P    = 0._dl


    end subroutine EFTCAMBScalingCubicSecondOrderEFTFunctions

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes adotoa = H.
    subroutine EFTCAMBScalingCubicComputeAdotoa( self, a, eft_par_cache, eft_cache )

        implicit none

        class(EFTCAMB_Scaling_Cubic)                  :: self          !< the base class
        real(dl), intent(in)                          :: a             !< the input scale factor
        type(TEFTCAMB_parameter_cache), intent(inout) :: eft_par_cache !< the EFTCAMB parameter cache that contains all the physical parameters.
        type(TEFTCAMB_timestep_cache ), intent(inout) :: eft_cache     !< the EFTCAMB timestep cache that contains all the physical values.

        real(dl)    :: temp

        eft_cache%grhom_t = (eft_par_cache%grhornomass +  eft_par_cache%grhog)/a**2+(eft_par_cache%grhob + eft_par_cache%grhoc)/a
        temp = 1.0_dl/(1.0_dl + eft_cache%EFTOmegaV+ a*eft_cache%EFTOmegaP)*(eft_cache%grhom_t + 2.0_dl*eft_cache%EFTc -eft_cache%EFTLambda )/3.0_dl
        eft_cache%adotoa = sqrt(temp)

    end subroutine EFTCAMBScalingCubicComputeAdotoa

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes the two derivatives wrt conformal time of H.
    subroutine EFTCAMBScalingCubicComputeHubbleDer( self, a, eft_par_cache, eft_cache )

        implicit none

        class(EFTCAMB_Scaling_Cubic)                    :: self          !< the base class
        real(dl), intent(in)                            :: a             !< the input scale factor
        type(TEFTCAMB_parameter_cache), intent(inout)   :: eft_par_cache !< the EFTCAMB parameter cache that contains all the physical parameters.
        type(TEFTCAMB_timestep_cache ), intent(inout)   :: eft_cache     !< the EFTCAMB timestep cache that contains all the physical values.

        real(dl) :: x1, x1p, x1pp, x1ppp, y1, y1p, y1pp, y2, y2p, y2pp, OmegaRad, OmegaRadprime, OmegaRadpprime
        real(dl) :: f, fprime, ft, h, hprime, hpprime, ht, htt, grhom, grhorad, H2, adotoa, Hdot, Hdotdot
        real(dl):: mu, x
        integer  :: ind

        x   = log(a)
        grhom   = eft_par_cache%grhob + eft_par_cache%grhoc
        grhorad = eft_par_cache%grhornomass +  eft_par_cache%grhog


        call self%fun_x1%precompute(x, ind, mu )
        x1    = self%fun_x1%value(x, index=ind, coeff=mu)
        y1    = self%fun_y1%value(x, index=ind, coeff=mu)
        y2    = self%fun_y2%value(x, index=ind, coeff=mu)

        H2= (grhorad/exp(x)**4+grhom/exp(x)**3)/((1._dl - x1**2._dl - y1**2._dl - y2**2._dl - 2._dl*self%SCG_A*x1*(sqrt(6._dl) - self%SCG_lambda*x1))*3._dl)
        
        OmegaRad= (grhorad/exp(x)**4)/(3._dl*(H2))

        f = (3._dl*(-3._dl*self%SCG_A + sqrt(6._dl)*(-1._dl + 2._dl*self%SCG_A*self%SCG_lambda)*x1 + 3._dl*self%SCG_A*(1._dl &
              &- 2._dl*self%SCG_A*self%SCG_lambda)*x1**2._dl + self%SCG_A*OmegaRad - 3._dl*self%SCG_A*y2**2._dl + self%SCG_beta2*y2**2._dl &
              &- 3._dl*self%SCG_A*y1**2._dl + self%SCG_beta1*y1**2._dl))/(1._dl + 6._dl*self%SCG_A**2._dl - 2._dl*self%SCG_A*self%SCG_lambda)

        h = ((-1._dl + 2._dl*self%SCG_A*self%SCG_lambda)*OmegaRad - 3._dl*(1._dl + 12._dl*self%SCG_A**2._dl &
              &- 2._dl*self%SCG_A*self%SCG_lambda + 2._dl*sqrt(6._dl)*self%SCG_A*(1._dl - 2._dl*self%SCG_A*self%SCG_lambda)*x1 &
              &+ ((1._dl - 2._dl*self%SCG_A*self%SCG_lambda)**2._dl)*x1**2._dl + (-1._dl + 2._dl*self%SCG_A*self%SCG_lambda &
              &- 2._dl*self%SCG_A*self%SCG_beta2)*y2**2._dl - y1**2._dl + 2._dl*self%SCG_A*self%SCG_lambda*y1**2._dl &
              &- 2._dl*self%SCG_A*self%SCG_beta1*y1**2._dl))/(2._dl + 12._dl*self%SCG_A**2._dl - 4._dl*self%SCG_A*self%SCG_lambda)

        x1p = (1._dl/sqrt(6._dl))*f - x1*h
        y1p = -sqrt(3._dl/2._dl)*self%SCG_beta1*x1*y1 - y1*h
        y2p = -sqrt(3._dl/2._dl)*self%SCG_beta2*x1*y2 - y2*h

        OmegaRadprime = -2._dl*OmegaRad*(2._dl+h)

        fprime = (-3._dl*(-self%SCG_A*OmegaRadprime + (-1._dl + 2._dl*self%SCG_A*self%SCG_lambda)*(-sqrt(6._dl) + 6._dl*self%SCG_A*x1)*x1p &
                 &+ 6._dl*self%SCG_A*y2*y2p - 2._dl*self%SCG_beta2*y2*y2p + 6._dl*self%SCG_A*y1*y1p - 2._dl*self%SCG_beta1*y1*y1p))/(1._dl &
                 &+ 6._dl*self%SCG_A**2._dl - 2._dl*self%SCG_A*self%SCG_lambda)

	    hprime = ((1._dl/2._dl)*OmegaRadprime*(-1._dl + 2._dl*self%SCG_A*self%SCG_lambda) + 3._dl*((-1._dl &
                 &+ 2._dl*self%SCG_A*self%SCG_lambda)*(sqrt(6._dl)*self%SCG_A + x1 - 2._dl*self%SCG_A*self%SCG_lambda*x1)*x1p + (1._dl &
                 &- 2._dl*self%SCG_A*self%SCG_lambda + 2._dl*self%SCG_A*self%SCG_beta2)*y2*y2p + (1._dl - 2._dl*self%SCG_A*self%SCG_lambda &
                 &+ 2._dl*self%SCG_A*self%SCG_beta1)*y1*y1p))/(1._dl + 6._dl*self%SCG_A**2._dl - 2._dl*self%SCG_A*self%SCG_lambda)

        ! second derivatives

        x1pp  = (1._dl/sqrt(6._dl))*fprime - hprime*x1 - h*x1p
        y1pp  = -sqrt(1.5_dl) *self%SCG_beta1*(x1p*y1 + x1*y1p) - hprime*y1 - h*y1p
        y2pp  = -sqrt(1.5_dl) *self%SCG_beta2*(x1p*y2 + x1*y2p) - hprime*y2 - h*y2p

        OmegaRadpprime = -2._dl*OmegaRadprime*(2._dl+h) - 2._dl*OmegaRad*hprime

        hpprime = - 1._dl/(1._dl + 2._dl*self%SCG_A*(3._dl *self%SCG_A - self%SCG_lambda)) * (0.5_dl*OmegaRadpprime*(2._dl*self%SCG_A*self%SCG_lambda - 1._dl) + 3._dl *( - (1 - 2._dl*self%SCG_A*self%SCG_lambda)**2._dl*x1p**2._dl &
                  + (2._dl*self%SCG_A*self%SCG_lambda -1._dl)*(sqrt(6._dl)*self%SCG_A + x1 - 2._dl*self%SCG_A*self%SCG_lambda*x1)*x1pp + (1._dl-2._dl*self%SCG_A*self%SCG_lambda+2._dl*self%SCG_A*self%SCG_beta1)*(y1p**2._dl+y1*y1pp) &
                  + (1._dl + 2._dl*self%SCG_A*self%SCG_lambda + 2._dl*self%SCG_A*self%SCG_beta2)*(y2p**2._dl + y2*y2pp) ) )
          
        adotoa=sqrt(H2)*a
        ht = adotoa*hprime
        Hdot=adotoa**2._dl*( 1._dl+h )
        htt = adotoa**2*hpprime + Hdot*hprime
        eft_cache%Hdot = Hdot
        Hdotdot = 2*adotoa*Hdot*( 1._dl+h ) + adotoa**2*ht
        eft_cache%Hdotdot = Hdotdot
        eft_cache%Hdotdotdot = adotoa**2*htt + 4._dl*adotoa*Hdot*ht + 2._dl*adotoa*Hdotdot*(1._dl+h) + 2._dl*Hdot**2*(1._dl+h)
    
    end subroutine EFTCAMBScalingCubicComputeHubbleDer

    ! ---------------------------------------------------------------------------------------------
    !> Function that computes model specific stability requirements.
    function EFTCAMBScalingCubicAdditionalModelStability( self, a, eft_par_cache, eft_cache )

        implicit none

        class(EFTCAMB_Scaling_Cubic)                  :: self          !< the base class
        real(dl), intent(in)                          :: a             !< the input scale factor.
        type(TEFTCAMB_parameter_cache), intent(inout) :: eft_par_cache !< the EFTCAMB parameter cache that contains all the physical parameters.
        type(TEFTCAMB_timestep_cache ), intent(inout) :: eft_cache     !< the EFTCAMB timestep cache that contains all the physical values.

        logical :: EFTCAMBScalingCubicAdditionalModelStability       !< the return value of the stability computation. True if the model specific stability criteria are met, false otherwise.

        EFTCAMBScalingCubicAdditionalModelStability = .True.

    end function EFTCAMBScalingCubicAdditionalModelStability

    ! ---------------------------------------------------------------------------------------------


end module EFTCAMB_full_Scaling_Cubic
