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

!> @file 007p7_ShiftSym_alphaB.f90
!! This file contains the relevant code for the shift-symmetric parametrization of alpha_B.


!----------------------------------------------------------------------------------------
!> This module contains the relevant code for the shift-symmetric parametrization of alpha_B.

!> @author Inês Albuquerque, Noemi Frusciante

module EFTCAMB_designer_ShiftSym_alphaB

    use precision
    use IniObjects
    use EFTCAMB_cache
    use MpiUtils
    use EFT_def
    use EFTCAMB_abstract_parametrizations_1D
    use EFTCAMB_parametrizations_1D
    use EFTCAMB_abstract_model_designer
    !use FileUtils
    !use equispaced_linear_interpolation_1D
    !use EFTCAMB_rootfind
    !use MassiveNu

    implicit none

    private

    public EFTCAMB_ShiftSym_alphaB

    !----------------------------------------------------------------------------------------
    !> This is the shift-symmetric parametrization of the alpha_B function. Inherits from the abstract designer model and has the
    !! freedom of defining the expansion history.
    type, extends ( EFTCAMB_designer_model ) :: EFTCAMB_ShiftSym_alphaB

        ! theory parameters:
        real(dl) :: SS_alphaB0                                                    !< The constant pre-factor of alpha_B.
        real(dl) :: SS_m                                                          !< The parameter that enters in the exponent of alpha_B.
        real(dl) :: SS_alphaK0                                                    !< The constant pre-factor of alpha_K.

        ! the pure EFT functions model selection flags:
        integer  :: EFTwDE                                                !< Model selection flag for the Shift-Symmetric choice of w_DE.

        ! the pure EFT functions:
        class( parametrized_function_1D ), allocatable :: ShiftSymwDE      !< The pure EFT function w_DE.

        ! the interpolated EFT functions that come out of the background sover: [I BELIEVE WE DON'T NEED THIS]
        !type(equispaced_linear_interpolate_function_1D) :: EFTOmega       !< The interpolated function Omega (and derivatives).
        !type(equispaced_linear_interpolate_function_1D) :: EFTLambda      !< The interpolated function Lambda (and derivatives).

        ! some designer parameters:
        integer  :: designer_num_points = 1000                            !< Number of points sampled by the designer code.
        real(dl) :: x_initial           = log(10._dl**(-8._dl))           !< log(a start)
        real(dl) :: x_final             = 0.0_dl                          !< log(a final)

    contains

        ! initialization of the model:
        procedure :: read_model_selection            => EFTCAMBShiftSymmetricReadModelSelectionFromFile   !< subroutine that reads the parameters of the model from file
        procedure :: allocate_model_selection        => EFTCAMBShiftSymmetricAllocateModelSelection       !< subroutine that allocates the model selection.
        procedure :: init_model_parameters           => EFTCAMBShiftSymmetricInitModelParameters          !< subroutine taht initializes the model parameters based on the values found in an input array.
        procedure :: init_model_parameters_from_file => EFTCAMBShiftSymmetricInitModelParametersFromFile  !< subroutine that reads the parameters of the model from file.

        ! background solver: [ALSO THINK WE DON'T NEED THESE]
        !procedure :: initialize_background           => EFTCAMBDesignerFRInitBackground               !< subroutine that initializes the background of designer f(R).
        !procedure :: solve_designer_equations        => EFTCAMBDesignerFRSolveDesignerEquations       !< subroutine that solves the designer f(R) background equations.
        !procedure :: find_initial_conditions         => EFTCAMBDesignerFRFindInitialConditions        !< subroutine that solves the background equations several time to determine the values of the initial conditions.

        ! utility functions:
        procedure :: compute_param_number  => EFTCAMBShiftSymmetricComputeParametersNumber     !< subroutine that computes the number of parameters of the model.
        procedure :: feedback              => EFTCAMBShiftSymmetricFeedback                    !< subroutine that prints on the screen feedback information about the model.
        procedure :: parameter_names       => EFTCAMBShiftSymmetricParameterNames              !< subroutine that returns the i-th parameter name of the model.
        procedure :: parameter_names_latex => EFTCAMBShiftSymmetricParameterNamesLatex         !< subroutine that returns the i-th parameter name of the model.
        procedure :: parameter_values      => EFTCAMBShiftSymmetricParameterValues             !< subroutine that returns the i-th parameter value.

        ! CAMB related procedures:
        procedure :: compute_background_EFT_functions  => EFTCAMBShiftSymmetricBackgroundEFTFunctions   !< subroutine that computes the value of the background EFT functions at a given time.
        procedure :: compute_secondorder_EFT_functions => EFTCAMBShiftSymmetricSecondOrderEFTFunctions  !< subroutine that computes the value of the second order EFT functions at a given time.
        procedure :: compute_dtauda                    => EFTCAMBShiftSymmetricComputeDtauda            !< function that computes dtauda = 1/sqrt(a^2H^2). [ARE THESE STILL NEEDED?]
        procedure :: compute_adotoa                    => EFTCAMBShiftSymmetricComputeAdotoa            !< subroutine that computes adotoa = H and its two derivatives wrt conformal time.
        procedure :: compute_H_derivs                  => EFTCAMBShiftSymmetricComputeHubbleDer         !< subroutine that computes the two derivatives wrt conformal time of H.

        ! stability procedures:
        !procedure :: additional_model_stability        => EFTCAMBShiftSymmetricAdditionalModelStability !< function that computes model specific stability requirements.

    end type EFTCAMB_ShiftSym_alphaB

    ! ---------------------------------------------------------------------------------------------

    ! define debug files
    !type(TTextFile) :: file_debug_1, file_debug_2

    ! ---------------------------------------------------------------------------------------------

contains

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that reads the parameters of the model from file.
    subroutine EFTCAMBShiftSymmetricReadModelSelectionFromFile( self, Ini, eft_error )

        implicit none

        class(EFTCAMB_ShiftSym_alphaB)  :: self      !< the base class
        type(TIniFile)                  :: Ini       !< Input ini file
        integer                         :: eft_error !< error code: 0 all fine, 1 initialization failed

        ! read model selection flags:
        self%EFTwDE             = Ini%Read_Int( 'EFTwDE', 0 )

    end subroutine EFTCAMBShiftSymmetricReadModelSelectionFromFile

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that allocates the model selection.
    subroutine EFTCAMBShiftSymmetricAllocateModelSelection( self, Ini, eft_error )

        implicit none

        class(EFTCAMB_ShiftSym_alphaB) :: self      !< the base class
        type(TIniFile)                 :: Ini       !< Input ini file
        integer                        :: eft_error !< error code: 0 all fine, 1 initialization failed

        character, allocatable, dimension(:) :: param_names       !< an array of strings containing the names of the function parameters
        character, allocatable, dimension(:) :: param_names_latex !< an array of strings containing the latex names of the function parameters

        integer                     :: temp_feedback 

        ! get feedback flag:
        temp_feedback = Ini%Read_Int('feedback_level', 0)

        ! allocate wDE:
        if ( allocated(self%ShiftSymwDE) ) deallocate(self%ShiftSymwDE)
        select case ( self%EFTwDE )
            case(0)
                allocate( wDE_LCDM_parametrization_1D::self%ShiftSymwDE )
                call self%ShiftSymwDE%set_name( 'EFTw', 'w' )
            case(1)
                allocate( constant_parametrization_1D::self%ShiftSymwDE )
                call self%ShiftSymwDE%set_name( 'EFTw', 'w' )
            case(2)
                allocate( CPL_parametrization_1D::self%ShiftSymwDE )
                call self%ShiftSymwDE%set_param_names( ['EFTw0', 'EFTwa'], ['w_0', 'w_a'] )
                call self%ShiftSymwDE%set_name( 'EFTw', 'w' )
            case(3)
                allocate( JBP_parametrization_1D::self%ShiftSymwDE )
                call self%ShiftSymwDE%set_param_names( ['EFTw0', 'EFTwa', 'EFTwn'], [ 'w_0', 'w_a', 'n  ' ] )
                call self%ShiftSymwDE%set_name( 'EFTw', 'w' )
            case(4)
                allocate( turning_point_parametrization_1D::self%ShiftSymwDE )
                call self%ShiftSymwDE%set_param_names( ['EFTw0 ', 'EFTwa ', 'EFTwat'], ['w_0', 'w_a', 'a_t'] )
                call self%ShiftSymwDE%set_name( 'EFTw', 'w' )
            case(5)
                allocate( taylor_parametrization_1D::self%ShiftSymwDE )
                call self%ShiftSymwDE%set_param_names( ['EFTw0', 'EFTwa', 'EFTw2', 'EFTw3'], ['w_0', 'w_a', 'w_2', 'w_3'] )
                call self%ShiftSymwDE%set_name( 'EFTw', 'w' )
            case default
                call allocate_parametrized_1D_function( self%ShiftSymwDE, self%EFTwDE, 'EFTw', 'w', eft_error, temp_feedback )
                if ( eft_error == 1 ) return
        end select

        ! additional initialization of the function:
        call self%ShiftSymwDE%init_func_from_file( Ini,eft_error )

    end subroutine EFTCAMBShiftSymmetricAllocateModelSelection

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that initializes the model parameters based on the values found in an input array.
    subroutine EFTCAMBShiftSymmetricInitModelParameters( self, array )

        implicit none

        class(EFTCAMB_ShiftSym_alphaB)                             :: self   !< the base class
        real(dl), dimension(self%parameter_number), intent(in)     :: array  !< input array with the values of the parameters.
        real(dl), dimension(self%parameter_number -1)              :: temp
        integer                                                    :: i

        self%SS_alphaB0 = array(1)
        self%SS_m       = array(2)
        self%SS_alphaK0 = array(3)

        do i = 1, self%parameter_number -1         !< I'M NOT SURE WHAT THIS DOES
            temp(i) = array(i+1)
        end do
        call self%ShiftSymwDE%init_parameters(temp)

    end subroutine EFTCAMBShiftSymmetricInitModelParameters

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that reads the parameters of the model from file.
    subroutine EFTCAMBShiftSymmetricInitModelParametersFromFile( self, Ini, eft_error )

        implicit none

        class(EFTCAMB_ShiftSym_alphaB)  :: self      !< the base class
        type(TIniFile)                  :: Ini       !< Input ini file
        integer                         :: eft_error !< error code: 0 all fine, 1 initialization failed

        ! read the alpha_B parameters:
        self%SS_alphaB0 = Ini%Read_Double( 'Shift_Symmetric_alphaB0', 0._dl )
        self%SS_m       = Ini%Read_Double( 'Shift_Symmetric_m', 0._dl )
        ! read the alpha_B parameter:
        self%SS_alphaK0 = Ini%Read_Double( 'Shift_Symmetric_alphaK0', 0._dl)
        ! read w_DE parameters:
        call self%ShiftSymwDE%init_from_file( Ini,eft_error)

    end subroutine EFTCAMBShiftSymmetricInitModelParametersFromFile

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes the number of parameters of the model.
    subroutine EFTCAMBShiftSymmetricComputeParametersNumber( self )

        implicit none

        class(EFTCAMB_ShiftSym_alphaB)  :: self   !< the base class

        self%parameter_number = 3
        self%parameter_number = self%parameter_number +self%ShiftSymwDE%parameter_number

    end subroutine EFTCAMBShiftSymmetricComputeParametersNumber

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that prints on the screen feedback information about the model.
    subroutine EFTCAMBShiftSymmetricFeedback( self, print_params )

        implicit none

        class(EFTCAMB_ShiftSym_alphaB)  :: self         !< the base class
        logical, optional               :: print_params !< optional flag that decised whether to print numerical values
                                                        !! of the parameters.

        write(*,*)
        write(*,'(a,a)')    '   Model               =  ', self%name
        write(*,'(a,I3)')   '   Number of params    ='  , self%parameter_number
        ! print model functions informations:
        if ( self%EFTwDE /= 0 ) then
            write(*,*)
            write(*,'(a,I3)')  '   EFTwDE              =', self%EFTwDE
        end if
        write(*,*)
        write(*,'(a,E13.3)') '   alpha_B0                  =', self%SS_alphaB0
        write(*,'(a,E13.3)') '   m                         =', self%SS_m
        write(*,'(a,E13.3)') '   alpha_K0                  =', self%SS_alphaK0

        call self%ShiftSymwDE%feedback( print_params )

    end subroutine EFTCAMBShiftSymmetricFeedback

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the i-th parameter name of the model
    subroutine EFTCAMBShiftSymmetricParameterNames( self, i, name )

        implicit none

        class(EFTCAMB_ShiftSym_alphaB)  :: self   !< the base class
        integer     , intent(in)        :: i      !< the index of the parameter
        character(*), intent(out)       :: name   !< the output name of the i-th parameter

        ! check validity of input:
        if ( i<=0 .or. i>self%parameter_number ) then
            write(*,'(a,I3)') 'EFTCAMB error: no parameter corresponding to number ', i
            write(*,'(a,I3)') 'Total number of parameters is ', self%parameter_number
            call MpiStop('EFTCAMB error')
        ! the first parameter is alpha_B0:
        else if ( i==1 ) then
            name = TRIM('alpha_B0')
            return
        ! the second parameter is m:
        else if ( i==2 ) then
            name = TRIM('m')
            return
        ! the third parameter is alpha_K0:
        else if ( i==3 ) then
            name = TRIM('alpha_K0')
            return
        ! the other parameters are the w_DE parameters: [NOT ENTIRELY SURE ABOUT THIS PART]
        else
            call self%ShiftSymwDE%parameter_names( i-1, name )
            return
        end if

    end subroutine EFTCAMBShiftSymmetricParameterNames

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the i-th parameter name of the model
    subroutine EFTCAMBShiftSymmetricParameterNamesLatex( self, i, latexname )

        implicit none

        class(EFTCAMB_ShiftSym_alphaB) :: self       !< the base class
        integer     , intent(in)       :: i         !< The index of the parameter
        character(*), intent(out)      :: latexname !< the output latex name of the i-th parameter

        ! check validity of input:
        if ( i<=0 .or. i>self%parameter_number ) then
            write(*,'(a,I3)') 'EFTCAMB error: no parameter corresponding to number ', i
            write(*,'(a,I3)') 'Total number of parameters is ', self%parameter_number
            call MpiStop('EFTCAMB error')
        ! the first parameter is alpha_B0:
        else if ( i==1 ) then
            latexname = TRIM('alpha_B0')
            return
        ! the second parameter is m:
        else if ( i==2 ) then
            latexname = TRIM('m')
            return
        ! the third parameter is alpha_K0:
        else if ( i==3 ) then
            latexname = TRIM('alpha_K0')
            return
        ! the other parameters are the w_DE parameters:
        else
            call self%ShiftSymwDE%parameter_names_latex( i-1, latexname )
            return
        end if

    end subroutine EFTCAMBShiftSymmetricParameterNamesLatex

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the i-th parameter value of the model
    subroutine EFTCAMBShiftSymmetricParameterValues( self, i, value )

        implicit none

        class(EFTCAMB_ShiftSym_alphaB) :: self   !< the base class
        integer , intent(in)           :: i     !< The index of the parameter
        real(dl), intent(out)          :: value !< the output value of the i-th parameter

        ! check validity of input:
        if ( i<=0 .or. i>self%parameter_number ) then
            write(*,'(a,I3)') 'EFTCAMB error: no parameter corresponding to number ', i
            write(*,'(a,I3)') 'Total number of parameters is ', self%parameter_number
            call MpiStop('EFTCAMB error')
        ! the first parameter is alpha_B0:
        else if ( i==1 ) then
            value = self%SS_alphaB0
            return
        ! the second parameter is m:
        else if ( i==2 ) then
            value = self%SS_m
            return
        ! the third parameter is alpha_K0:
        else if ( i==3 ) then
            value = self%SS_alphaK0
            return
        ! the other parameters are the w_DE parameters:
        else
            call self%ShiftSymwDE%parameter_value( i-1, value )
            return
        end if

    end subroutine EFTCAMBShiftSymmetricParameterValues

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes the value of the background EFT functions at a given time.
    subroutine EFTCAMBShiftSymmetricBackgroundEFTFunctions( self, a, eft_par_cache, eft_cache )

        implicit none

        class(EFTCAMB_ShiftSym_alphaB)                   :: self          !< the base class
        real(dl), intent(in)                             :: a             !< the input scale factor
        type(TEFTCAMB_parameter_cache), intent(inout)    :: eft_par_cache !< the EFTCAMB parameter cache that contains all the physical parameters.
        type(TEFTCAMB_timestep_cache ), intent(inout)    :: eft_cache     !< the EFTCAMB timestep cache that contains all the physical values.


        real(dl)  :: RPH_PM_V, RPH_AT_V, RPH_PM_P, RPH_AT_P, RPH_PM_PP, RPH_AT_PP, RPH_PM_PPP, RPH_AT_PPP, RPH_PM_PPPP, RPH_AT_PPPP


        ! defining the necessary alpha functions:
        RPH_PM_V               = 0._dl
        RPH_AT_V               = 0._dl
        RPH_PM_P               = 0._dl
        RPH_AT_P               = 0._dl
        RPH_PM_PP              = 0._dl
        RPH_AT_PP              = 0._dl
        RPH_PM_PPP             = 0._dl
        RPH_AT_PPP             = 0._dl
        RPH_PM_PPPP            = 0._dl
        RPH_AT_PPPP             = 0._dl

        ! compute the background EFT functions:

        eft_cache%EFTOmegaV    =  RPH_PM_V + RPH_AT_V*(1._dl + RPH_PM_V)
        eft_cache%EFTOmegaP    = (RPH_PM_V + 1.0_dl)*RPH_AT_P + (1._dl + RPH_AT_V)*RPH_PM_P
        eft_cache%EFTOmegaPP   = (RPH_PM_V + 1.0_dl)*RPH_AT_PP + 2._dl*RPH_PM_P*RPH_AT_P + (1._dl + RPH_AT_V)*RPH_PM_PP
        eft_cache%EFTOmegaPPP  = (RPH_PM_V + 1.0_dl)*RPH_AT_PPP + 3._dl*RPH_PM_PP*RPH_AT_P + 3._dl*RPH_PM_P*RPH_AT_PP + (1._dl +RPH_AT_V)*RPH_PM_PPP
        eft_cache%EFTOmegaPPPP = RPH_PM_PPPP*(1._dl + RPH_AT_V) + 4._dl*RPH_PM_PPP*RPH_AT_P + 6._dl*RPH_PM_PP*RPH_AT_PP + 4._dl*RPH_PM_P*RPH_AT_PPP + RPH_AT_PPPP*(1._dl + RPH_PM_V)

        eft_cache%EFTc         = ( eft_cache%adotoa**2._dl - eft_cache%Hdot )*( eft_cache%EFTOmegaV + 0.5_dl*a*eft_cache%EFTOmegaP ) - 0.5_dl*( a*eft_cache%adotoa )**2._dl*eft_cache%EFTOmegaPP &
            & +0.5_dl*eft_cache%grhov_t*( 1._dl + self%ShiftSymwDE%value(a) )

        eft_cache%EFTLambda    = - eft_cache%EFTOmegaV*( 2._dl*eft_cache%Hdot + eft_cache%adotoa**2._dl ) -a*eft_cache%EFTOmegaP*( 2._dl*eft_cache%adotoa**2._dl + eft_cache%Hdot ) &
                               & -( a*eft_cache%adotoa )**2._dl*eft_cache%EFTOmegaPP + self%ShiftSymwDE%value(a)*eft_cache%grhov_t

        eft_cache%EFTcdot      = +0.5_dl*eft_cache%adotoa*eft_cache%grhov_t*( -3._dl*(1._dl + self%ShiftSymwDE%value(a))**2._dl + a*self%ShiftSymwDE%first_derivative(a) ) &
                                & -eft_cache%EFTOmegaV*( eft_cache%Hdotdot - 4._dl*eft_cache%adotoa*eft_cache%Hdot + 2._dl*eft_cache%adotoa**3._dl ) &
                                & +0.5_dl*a*eft_cache%EFTOmegaP*( -eft_cache%Hdotdot +eft_cache%adotoa*eft_cache%Hdot +eft_cache%adotoa**3._dl) &
                                & +0.5_dl*a**2._dl*eft_cache%adotoa*eft_cache%EFTOmegaPP*( eft_cache%adotoa**2 - 3._dl*eft_cache%Hdot ) - 0.5_dl*(a*eft_cache%adotoa)**3._dl*eft_cache%EFTOmegaPPP

        eft_cache%EFTLambdadot = - 2._dl*eft_cache%EFTOmegaV*( eft_cache%Hdotdot - eft_cache%adotoa*eft_cache%Hdot - eft_cache%adotoa**3._dl )                     &
                                 - a*eft_cache%EFTOmegaP*( +eft_cache%Hdotdot + 5._dl*eft_cache%adotoa*eft_cache%Hdot -eft_cache%adotoa**3._dl  )                  &
                                 - a**2._dl*eft_cache%EFTOmegaPP*eft_cache%adotoa*( +2._dl*eft_cache%adotoa**2._dl + 3._dl*eft_cache%Hdot )                        &
                                 - (a*eft_cache%adotoa)**3._dl*eft_cache%EFTOmegaPPP + eft_cache%grhov_t*eft_cache%adotoa*( a*self%ShiftSymwDE%first_derivative(a) - 3._dl*self%ShiftSymwDE%value(a)*(1._dl + self%ShiftSymwDE%value(a)))

        eft_cache%EFTLambdadotdot = - eft_cache%EFTOmegaV*(2._dl*eft_cache%Hdotdotdot-2._dl*eft_cache%Hdot**2._dl-6._dl*eft_cache%adotoa*eft_cache%Hdotdot-2._dl*eft_cache%adotoa**2._dl*eft_cache%Hdot + 4._dl*eft_cache%adotoa**4._dl)&
                                    - a*eft_cache%EFTOmegaP*(6._dl*eft_cache%adotoa*eft_cache%Hdotdot+eft_cache%adotoa**4._dl+5._dl*eft_cache%Hdot**2._dl+eft_cache%Hdotdotdot - 10._dl*eft_cache%adotoa**2._dl*eft_cache%Hdot)      &
                                    - a**2._dl*eft_cache%EFTOmegaPP*(11._dl*eft_cache%adotoa**2._dl*eft_cache%Hdot - eft_cache%adotoa**4._dl + 4._dl*eft_cache%adotoa*eft_cache%Hdotdot + 3._dl*eft_cache%Hdot**2._dl)               &
                                    - a**3._dl*(3._dl*eft_cache%adotoa**4._dl + 6._dl*eft_cache%adotoa**2._dl*eft_cache%Hdot) - a**4._dl*eft_cache%adotoa**4._dl*eft_cache%EFTOmegaPPPP                                              &
                                    + eft_cache%Hdot*eft_cache%grhov_t*(a*self%ShiftSymwDE%first_derivative(a) - 3._dl*self%ShiftSymwDE%value(a)*(1._dl + self%ShiftSymwDE%value(a)) )                                               &
                                    + eft_cache%adotoa**2._dl*eft_cache%grhov_t*(9._dl*self%ShiftSymwDE%value(a)*(1._dl + self%ShiftSymwDE%value(a)) - 5._dl*a*self%ShiftSymwDE%first_derivative(a)                                  &
                                    + 18._dl*self%ShiftSymwDE%value(a)**2._dl + 9._dl*self%ShiftSymwDE%value(a)**3._dl + a**2._dl*self%ShiftSymwDE%second_derivative(a))

    end subroutine EFTCAMBShiftSymmetricBackgroundEFTFunctions

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes the value of the background EFT functions at a given time.
    subroutine EFTCAMBShiftSymmetricSecondOrderEFTFunctions( self, a, eft_par_cache, eft_cache )

        implicit none

        class(EFTCAMB_ShiftSym_alphaB)                   :: self          !< the base class
        real(dl), intent(in)                             :: a             !< the input scale factor
        type(TEFTCAMB_parameter_cache), intent(inout)    :: eft_par_cache !< the EFTCAMB parameter cache that contains all the physical parameters.
        type(TEFTCAMB_timestep_cache ), intent(inout)    :: eft_cache     !< the EFTCAMB timestep cache that contains all the physical values.

        real(dl) :: RPH_PM_V, RPH_AT_V, RPH_AK_V, RPH_AB_V 
        real(dl) :: RPH_PM_P, RPH_AT_P, RPH_AK_P, RPH_AB_P
        real(dl) :: RPH_PM_PP, RPH_AT_PP, RPH_AK_PP, RPH_AB_PP
        real(dl) :: RPH_PM_PPP, RPH_AT_PPP, RPH_AK_PPP, RPH_AB_PPP
        real(dl) :: RPH_AT_PPPP, RPH_PM_PPPP

        ! defining the necessary alpha functions:
        RPH_PM_V               = 0._dl
        RPH_AT_V               = 0._dl
        RPH_AK_V               = self%SS_alphaK0*a
        RPH_AB_V               = self%SS_alphaB0*( (a*eft_par_cache%h0_Mpc)/(eft_cache%adotoa) )**(4._dl/self%SS_m)

        RPH_PM_P               = 0._dl
        RPH_AT_P               = 0._dl
        RPH_AK_P               = self%SS_alphaK0
        RPH_AB_P               = ((4._dl*self%SS_alphaB0)/self%SS_m)*( (a*eft_par_cache%h0_Mpc)/(eft_cache%adotoa) )**(4._dl/self%SS_m - 1._dl)*( eft_par_cache%h0_Mpc/eft_cache%adotoa &
                                 - (eft_par_cache%h0_Mpc*eft_cache%Hdot)/(eft_cache%adotoa**3) )

        RPH_PM_PP              = 0._dl
        RPH_AT_PP              = 0._dl
        RPH_AK_PP              = 0._dl
        RPH_AB_PP              = - 4._dl*self%SS_alphaB0/(a**2._dl*self%SS_m*eft_cache%adotoa**4._dl)*(a*eft_par_cache%h0_Mpc/eft_cache%adotoa)**(4._dl/self%SS_m)*((self%SS_m - 4._dl)*eft_cache%adotoa**4._dl     &
                                 - (self%SS_m - 8._dl)*eft_cache%adotoa**2._dl*eft_cache%Hdot - 2._dl*(2._dl + self%SS_m)*eft_cache%Hdot**2._dl + self%SS_m*eft_cache%adotoa*eft_cache%Hdotdot) 
        
        RPH_PM_PPP             = 0._dl
        RPH_AT_PPP             = 0._dl
        RPH_AK_PPP             = 0._dl
        RPH_AB_PPP             = - 4._dl*self%SS_alphaB0/(a**3._dl*self%SS_m**3._dl*eft_cache%adotoa**6._dl)*(a*eft_par_cache%h0_Mpc/eft_cache%adotoa)**(4._dl/self%SS_m)*(2._dl*(self%SS_m**2._dl - 6._dl*self%SS_m + 8._dl)     & 
                                 * eft_cache%adotoa**6._dl - 2._dl*(self%SS_m**2._dl - 12._dl*self%SS_m +24._dl)*eft_cache%adotoa**4._dl*eft_cache%Hdot - 8._dl*(self%SS_m + 3._dl*self%SS_m + 2._dl)*eft_cache%Hdot**3._dl       &
                                 + 3._dl*(self%SS_m - 4._dl)*self%SS_m*eft_cache%adotoa**3._dl*eft_cache%Hdotdot + self%SS_m*(7._dl*self%SS_m + 12._dl)*eft_cache%adotoa*eft_cache%Hdot*eft_cache%Hdotdot -                       &
                                 - 6._dl*(2._dl*self%SS_m**2._dl - 2._dl*self%SS_m - 8._dl)*eft_cache%adotoa**2._dl*eft_cache%Hdot**2._dl - self%SS_m**2._dl*eft_cache%adotoa**2._dl*eft_cache%Hdotdotdot) 

        RPH_AT_PPPP            = 0._dl 
        RPH_PM_PPPP            = 0._dl
        ! compute the second order EFT functions:
        eft_cache%EFTGamma1V   = 0.25_dl*( RPH_AK_V*(1._dl + RPH_PM_V)*eft_cache%adotoa**2 - 2._dl*eft_cache%EFTc )/(eft_par_cache%h0_Mpc**2*a**2)

        eft_cache%EFTGamma1P   = - 0.5_dl*( RPH_AK_V*(1._dl + RPH_PM_V)*eft_cache%adotoa**2 - 2._dl*eft_cache%EFTc )/(eft_par_cache%h0_Mpc**2*a**3)                                          &
                                + 0.25_dl*( RPH_AK_P*(1._dl + RPH_PM_V)*eft_cache%adotoa**2 + RPH_AK_V*RPH_PM_P*eft_cache%adotoa**2                                                          &
                                + 2._dl*RPH_AK_V*(1._dl + RPH_PM_V)*eft_cache%Hdot/a - 4._dl*eft_cache%EFTc/a - 2._dl*eft_cache%EFTcdot/a/eft_cache%adotoa )/(eft_par_cache%h0_Mpc**2*a**2)



        eft_cache%EFTGamma2V  = (- 1._dl*RPH_AB_V*(1._dl + RPH_PM_V) - a*eft_cache%EFTOmegaP )*eft_cache%adotoa/(eft_par_cache%h0_Mpc*a)

        eft_cache%EFTGamma2P  = - (- RPH_AB_V*(1._dl + RPH_PM_V) - a*eft_cache%EFTOmegaP )*eft_cache%adotoa/(eft_par_cache%h0_Mpc*a**2) - ((1._dl + RPH_PM_V)*( RPH_AB_P*eft_cache%adotoa**2 + RPH_AB_V*eft_cache%Hdot/a)              &
                                + 1._dl*RPH_AB_V*eft_cache%adotoa**2*RPH_PM_P + eft_cache%EFTOmegaP*( eft_cache%adotoa**2 + eft_cache%Hdot ) + a*eft_cache%adotoa**2*eft_cache%EFTOmegaPP )/(eft_par_cache%h0_Mpc*a*eft_cache%adotoa)
        
        eft_cache%EFTGamma2PP = ((eft_cache%Hdot**2._dl - eft_cache%adotoa*eft_cache%Hdotdot)*((1._dl + RPH_PM_V)*RPH_AB_V + a*eft_cache%EFTOmegaP) + eft_cache%adotoa**2._dl*eft_cache%Hdot             & 
                                * (RPH_AB_V*(3._dl*(1._dl + RPH_PM_V) - 2._dl*a*RPH_PM_P) +a*(eft_cache%EFTOmegaP - 2._dl*(1._dl + RPH_PM_V)*RPH_AB_P - 2._dl*a*eft_cache%EFTOmegaPP) )                  &
                                + eft_cache%adotoa**4._dl * (RPH_AB_V * (- 2._dl * (1._dl + RPH_PM_V) + 2._dl*a*RPH_PM_P - a**2._dl*RPH_PM_PP) + a*(1._dl + RPH_PM_V)*(2._dl*RPH_AB_P - a*RPH_AB_PP) )   &
                                - a * (2._dl*a*RPH_PM_P*RPH_AB_P - a**2._dl*eft_cache%EFTOmegaPPP))/(a**3._dl*eft_par_cache%h0_Mpc*eft_cache%adotoa**3._dl)

        eft_cache%EFTGamma2PPP = ((- 3._dl*eft_cache%Hdot**3._dl* - 4._dl*eft_cache%adotoa*eft_cache%Hdot*eft_cache%Hdotdot)*((1._dl + RPH_PM_V)*RPH_AB_V + a*eft_cache%EFTOmegaP)                                   &
                                 + eft_cache%adotoa**3._dl*eft_cache%Hdotdot*(RPH_AB_V*(6._dl*(1._dl + RPH_PM_V) - 3._dl*a*RPH_PM_V )  + a*(3._dl*eft_cache%EFTOmegaP - 3._dl*(1._dl + RPH_PM_V)*RPH_AB_P            & 
                                 - 3._dl*eft_cache%EFTOmegaPP*a)) - (RPH_AB_V*(6._dl*(1._dl + RPH_PM_V) - 3._dl*a*RPH_PM_V )  + a*(3._dl*eft_cache%EFTOmegaP - 3._dl*(1._dl + RPH_PM_V)*RPH_AB_P                     & 
                                 - 3._dl*eft_cache%EFTOmegaPP*a))*eft_cache%adotoa**2._dl*eft_cache%Hdot**2._dl - eft_cache%adotoa**2._dl*eft_cache%Hdotdotdot*((1._dl + RPH_PM_V)*RPH_AB_V + a*eft_cache%EFTOmegaP) &
                                 + eft_cache%adotoa**4._dl*eft_cache%Hdot*(RPH_AB_V*(-11._dl*(1._dl + RPH_PM_V) + 9._dl*a*RPH_PM_P - 3._dl*a**2._dl*RPH_PM_PP) + a*(-3._dl*a*(1._dl + RPH_PM_V)*RPH_AB_PP            & 
                                 + 6._dl*a*RPH_PM_P*RPH_AB_P - 3._dl*a**2._dl*eft_cache%EFTOmegaPPP - 2._dl*eft_cache%EFTOmegaP + 3._dl*a*eft_cache%EFTOmegaPP + 9._dl*(1._dl + RPH_PM_V)*RPH_AB_P))                 &
                                 + eft_cache%adotoa**6._dl*(RPH_AB_V*(6._dl*(1._dl + RPH_PM_V) - 6._dl*a*RPH_PM_P + 3._dl*a**2._dl*RPH_PM_PP - a**3._dl*RPH_PM_PPP) + a*(-6._dl*(1._dl + RPH_PM_V)*RPH_AB_P          &
                                 + 3._dl*a*(1._dl + RPH_PM_V)*RPH_AB_PP + 6._dl*a*RPH_PM_P*RPH_AB_P - 3._dl*a**2._dl*(RPH_PM_PP*RPH_AB_P + RPH_PM_P*RPH_AB_PP) - a**2._dl*(1._dl + RPH_PM_V)*RPH_AB_PPP -            &
                                 - a**3._dl*eft_cache%EFTOmegaPPPP) ))/(a**4._dl*eft_par_cache%h0_Mpc*eft_cache%adotoa**5._dl)


        eft_cache%EFTGamma3V    = -RPH_AT_V*(1._dl + RPH_PM_V)
        eft_cache%EFTGamma3P    = -RPH_PM_P*RPH_AT_V - (1._dl + RPH_PM_V)*RPH_AT_P
        eft_cache%EFTGamma3PP   = -RPH_AT_PP*(1._dl + RPH_PM_V) - 2._dl*RPH_AT_P*RPH_PM_P - RPH_AT_V*RPH_PM_PP
        eft_cache%EFTGamma3PPP  = -RPH_AT_PPP*(1._dl + RPH_PM_V) - 3._dl*RPH_AT_PP*RPH_PM_P - 3._dl*RPH_AT_P*RPH_PM_PP - RPH_AT_V*RPH_PM_PPP
        eft_cache%EFTGamma3PPPP = -RPH_AT_PPPP*(1._dl + RPH_PM_V) - 4._dl*RPH_AT_PPP*RPH_PM_P - 6._dl*RPH_AT_PP*RPH_PM_PP - 4._dl*RPH_AT_P*RPH_PM_PPP - RPH_AT_V*RPH_PM_PPPP
        eft_cache%EFTGamma4V    = -eft_cache%EFTGamma3V
        eft_cache%EFTGamma4P    = -eft_cache%EFTGamma3P
        eft_cache%EFTGamma4PP   = +(1._dl + RPH_PM_V)*RPH_AT_PP + RPH_PM_PP*RPH_AT_V + 2._dl*RPH_PM_P*RPH_AT_P
        eft_cache%EFTGamma5V    = +0.5_dl*eft_cache%EFTGamma3V
        eft_cache%EFTGamma5P    = +0.5_dl*eft_cache%EFTGamma3P
        eft_cache%EFTGamma6V    = 0._dl
        eft_cache%EFTGamma6P    = 0._dl

    end subroutine EFTCAMBShiftSymmetricSecondOrderEFTFunctions

    ! ---------------------------------------------------------------------------------------------
    !> Function that computes dtauda = 1/sqrt(a^2H^2).For pure EFT std this has to be overridden
    !! for performance reasons.
    function EFTCAMBShiftSymmetricComputeDtauda( self, a, eft_par_cache, eft_cache )

        implicit none

        class(EFTCAMB_ShiftSym_alphaB)                   :: self          !< the base class
        real(dl), intent(in)                             :: a             !< the input scale factor
        type(TEFTCAMB_parameter_cache), intent(inout)    :: eft_par_cache !< the EFTCAMB parameter cache that contains all the physical parameters.
        type(TEFTCAMB_timestep_cache ), intent(inout)    :: eft_cache     !< the EFTCAMB timestep cache that contains all the physical values.

        real(dl) :: EFTCAMBShiftSymmetricComputeDtauda               !< the output dtauda

        real(dl) :: temp

        temp = eft_cache%grhoa2 + eft_par_cache%grhov*a*a*self%ShiftSymwDE%integral(a)
        EFTCAMBShiftSymmetricComputeDtauda = sqrt(3/temp)

    end function EFTCAMBShiftSymmetricComputeDtauda

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes adotoa = H and its two derivatives wrt conformal time.
    subroutine EFTCAMBShiftSymmetricComputeAdotoa( self, a, eft_par_cache, eft_cache )

        implicit none

        class(EFTCAMB_ShiftSym_alphaB)                   :: self          !< the base class
        real(dl), intent(in)                             :: a             !< the input scale factor
        type(TEFTCAMB_parameter_cache), intent(inout)    :: eft_par_cache !< the EFTCAMB parameter cache that contains all the physical parameters.
        type(TEFTCAMB_timestep_cache ), intent(inout)    :: eft_cache     !< the EFTCAMB timestep cache that contains all the physical values.

        eft_cache%grhov_t = eft_par_cache%grhov*self%ShiftSymwDE%integral(a)
        eft_cache%adotoa  = sqrt( ( eft_cache%grhom_t + eft_cache%grhov_t )/3._dl )

    end subroutine EFTCAMBShiftSymmetricComputeAdotoa

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes the two derivatives wrt conformal time of H.
    subroutine EFTCAMBShiftSymmetricComputeHubbleDer( self, a, eft_par_cache, eft_cache )

        implicit none

        class(EFTCAMB_ShiftSym_alphaB)                   :: self          !< the base class
        real(dl), intent(in)                             :: a             !< the input scale factor
        type(TEFTCAMB_parameter_cache), intent(inout)    :: eft_par_cache !< the EFTCAMB parameter cache that contains all the physical parameters.
        type(TEFTCAMB_timestep_cache ), intent(inout)    :: eft_cache     !< the EFTCAMB timestep cache that contains all the physical values.

        eft_cache%gpiv_t  = self%ShiftSymwDE%value(a)*eft_cache%grhov_t

        eft_cache%Hdot    = -0.5_dl*( eft_cache%adotoa**2 +eft_cache%gpresm_t +eft_cache%gpiv_t )

        eft_cache%Hdotdot =  eft_cache%adotoa*( ( eft_cache%grhob_t + eft_cache%grhoc_t)/6._dl +2._dl*( eft_cache%grhor_t +eft_cache%grhog_t)/3._dl )                                          &
                           + eft_cache%adotoa*eft_cache%grhov_t*( 1._dl/6._dl +self%ShiftSymwDE%value(a) +1.5_dl*self%ShiftSymwDE%value(a)**2 -0.5_dl*a*self%ShiftSymwDE%first_derivative(a) ) &
                           + eft_cache%adotoa*eft_cache%grhonu_tot/6._dl -0.5_dl*eft_cache%adotoa*eft_cache%gpinu_tot -0.5_dl*eft_cache%gpinudot_tot

        eft_cache%Hdotdotdot = (2._dl*(eft_cache%grhob_t + eft_cache%grhoc_t)*eft_cache%adotoa**2 + eft_cache%Hdot &
                          &* (eft_cache%grhob_t + eft_cache%grhoc_t) - 3._dl*eft_cache%adotoa**2 *(eft_cache%grhob_t + eft_cache%grhoc_t))/(6._dl) &
                          &+ (2._dl/3._dl)*(2._dl*eft_cache%adotoa**2 * (eft_cache%grhor_t + eft_cache%grhog_t) + eft_cache%Hdot &
                          &* (eft_cache%grhor_t + eft_cache%grhog_t) - 4._dl*eft_cache%adotoa**2 * (eft_cache%grhor_t + eft_cache%grhog_t)) &
                          &+ (2._dl*eft_cache%adotoa**2 * eft_cache%grhov_t - 3._dl*eft_cache%adotoa**2 * eft_cache%grhov_t &
                          &* (1._dl + self%ShiftSymwDE%value(a)) + eft_cache%Hdot*eft_cache%grhov_t)*(1._dl/6._dl + self%ShiftSymwDE%value(a) &
                          &+ (3._dl/2._dl)*self%ShiftSymwDE%value(a)**2 - (1._dl/2._dl)*a*self%ShiftSymwDE%first_derivative(a)) &
                          &+ a*eft_cache%adotoa**2 * eft_cache%grhov_t * ((1._dl/2._dl)*self%ShiftSymwDE%first_derivative(a) &
                          &+ 3._dl*self%ShiftSymwDE%value(a)*self%ShiftSymwDE%first_derivative(a) - (a/2._dl)*self%ShiftSymwDE%second_derivative(a)) &
                          &+ (eft_cache%adotoa**2/3._dl)*eft_cache%grhonu_tot - eft_cache%adotoa**2 * eft_cache%gpinu_tot &
                          &- (3._dl/2._dl)*eft_cache%adotoa*eft_cache%gpinudot_tot + (eft_cache%Hdot/6._dl)*eft_cache%grhonu_tot &
                          &+ (eft_cache%adotoa/6._dl)*eft_cache%grhonudot_tot - (eft_cache%Hdot/2._dl)*eft_cache%gpinu_tot &
                          &- (1._dl/2._dl)*eft_cache%gpinudotdot_tot
    end subroutine EFTCAMBShiftSymmetricComputeHubbleDer

    ! ---------------------------------------------------------------------------------------------

end module EFTCAMB_designer_ShiftSym_alphaB

!----------------------------------------------------------------------------------------
