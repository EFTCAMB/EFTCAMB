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

!> @file 008p3_quintessence.f90
!! This file contains the definition of the Quintessence full mapping model.


!----------------------------------------------------------------------------------------
!> This module contains the definition of the Quintessence full mapping model.

!> @author Marco Raveri

module EFTCAMB_FM_quintessence

    use precision
    use IniObjects
    use MpiUtils
    use FileUtils
    use constants, only : c, const_pi
    use equispaced_linear_interpolation_1D
    use EFTCAMB_cache
    use EFT_def
    use EFTCAMB_mixed_algorithms, only : double_NaN
    use EFTCAMB_rootfind
    use EFTCAMB_abstract_model_full
    use EFTCAMB_abstract_parametrizations_1D
    use EFTCAMB_power_law_parametrizations_1D
    use EFTCAMB_constant_parametrization_1D
    use EFTCAMB_taylor_series_1D
    use EFTCAMB_pade_series_1D
    use EFTCAMB_fourier_parametrizations_1D
    use EFTCAMB_exponential_parametrizations_2_1D
    use EFTCAMB_double_exponential_parametrizations_1D
    use EFTCAMB_cosine_parametrizations_1D
    use EFTCAMB_axion_parametrizations_1D
    use EFTCAMB_power_law_sum_parametrizations_1D
    use MassiveNu

    implicit none

    private

    public EFTCAMB_5e

    !----------------------------------------------------------------------------------------
    !> This is the type that contains the definition of the Quintessence
    !! full mapping model.
    type, extends ( EFTCAMB_full_model ) :: EFTCAMB_5e

        ! 5e theory parameters:
        real(dl) :: phidot_ini      !< initial value of the quintessence kinetic term.
        real(dl) :: field_min       !< bounds on the initial value of the field, minimum value. Useful for periodic potentials. Default to +10^10
        real(dl) :: field_max       !< bounds on the initial value of the field, maximum value. Useful for periodic potentials. Default to -10^10

        ! scalar field potential model selection flag:
        integer  :: potential_model         !< Model selection flag for the 5e potential model.
        ! initial conditions flag:
        logical  :: drag_initial_conditions !< Flag to select Hubble dragged initial conditions.

        ! the scalar field potential in case we need it:
        class( parametrized_function_1D ), allocatable  :: potential   !< The JBD scalar field potential.

        ! the interpolated EFT functions that come out of the background solver:
        type(equispaced_linear_interpolate_function_1D) :: EFTLambda      !< The interpolated function Lambda (and derivatives).
        type(equispaced_linear_interpolate_function_1D) :: EFTc           !< The interpolated function c (and derivatives).

        ! some parameters:
        integer  :: interpolation_num_points = 1100                       !< Number of points sampled by the background solver code.
        real(dl) :: x_initial                = log(10._dl**(-8._dl))      !< log(a start)
        real(dl) :: x_final                  = 0.1_dl                     !< log(a final)

    contains

        ! initialization of the model:
        procedure :: read_model_selection            => EFTCAMB5eReadModelSelectionFromFile  !< subroutine that reads the parameters of the model from file
        procedure :: allocate_model_selection        => EFTCAMB5eAllocateModelSelection      !< subroutine that allocates the model selection.
        procedure :: init_model_parameters           => EFTCAMB5eInitModelParameters         !< subroutine that initializes the model parameters based on the values found in an input array.
        procedure :: init_model_parameters_from_file => EFTCAMB5eInitModelParametersFromFile !< subroutine that reads the parameters of the model from file.

        ! utility functions:
        procedure :: compute_param_number  => EFTCAMB5eComputeParametersNumber    !< subroutine that computes the number of parameters of the model.
        procedure :: feedback              => EFTCAMB5eFeedback                   !< subroutine that prints on the screen feedback information about the model.
        procedure :: parameter_names       => EFTCAMB5eParameterNames             !< subroutine that returns the i-th parameter name of the model.
        procedure :: parameter_names_latex => EFTCAMB5eParameterNamesLatex        !< subroutine that returns the i-th parameter name of the model.
        procedure :: parameter_values      => EFTCAMB5eParameterValues            !< subroutine that returns the i-th parameter value.

        ! CAMB related procedures:
        procedure :: compute_background_EFT_functions  => EFTCAMB5eBackgroundEFTFunctions   !< subroutine that computes the value of the background EFT functions at a given time.
        procedure :: compute_secondorder_EFT_functions => EFTCAMB5eSecondOrderEFTFunctions  !< subroutine that computes the value of the second order EFT functions at a given time.

        ! background solver:
        procedure :: initialize_background           => EFTCAMB5eInitBackground               !< subroutine that initializes the background of 5e.
        procedure :: solve_background_equations      => EFTCAMB5eSolveBackgroundEquations     !< subroutine that solves the 5e background equations.
        procedure :: find_initial_conditions         => EFTCAMB5eFindInitialConditions        !< subroutine that solves the background equations several time to determine the values of the initial conditions.

        ! potential:
        procedure :: quintessence_potential          => EFTCAMB5ePotential !< function that specified the quintessence potential.

    end type

    ! ---------------------------------------------------------------------------------------------

    ! define debug files
    type(TTextFile) :: file_debug_1, file_debug_2, file_debug_3, file_debug_4

    ! ---------------------------------------------------------------------------------------------

contains

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that reads the parameters of the model from file.
    subroutine EFTCAMB5eReadModelSelectionFromFile( self, Ini, eft_error )

        implicit none

        class(EFTCAMB_5e) :: self      !< the base class
        type(TIniFile)    :: Ini       !< Input ini file
        integer           :: eft_error !< error code: 0 all fine, 1 initialization failed

        ! read model selection flags:
        self%potential_model = Ini%Read_Int( 'potential_model', 1 )
        ! read bounds on field value:
        self%field_min = Ini%Read_Double( 'field_min', -1.d10 )
        if ( self%field_min > 0._dl ) then
            write(*,*) 'The minimum value of the field has to be negative.'
            eft_error = 1
            return
        end if
        self%field_max = Ini%Read_Double( 'field_max', +1.d10 )
        if ( self%field_max < 0._dl ) then
            write(*,*) 'The maximum value of the field has to be positive.'
            eft_error = 1
            return
        end if
        ! read the initial condition flag:
        self%drag_initial_conditions = Ini%Read_Logical( 'drag_initial_conditions', .False. )

    end subroutine EFTCAMB5eReadModelSelectionFromFile

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that allocates the model selection.
    subroutine EFTCAMB5eAllocateModelSelection( self, Ini, eft_error )

        implicit none

        class(EFTCAMB_5e) :: self      !< the base class
        type(TIniFile)    :: Ini       !< Input ini file
        integer           :: eft_error !< error code: 0 all fine, 1 initialization failed

        integer                     :: temp_feedback 

        ! get feedback flag:
        temp_feedback = Ini%Read_Int('feedback_level', 0)

        ! here we allocate the model for the quintessence potential
        if ( allocated(self%potential) ) deallocate(self%potential)
        select case ( self%potential_model )
            case(1)
                allocate( constant_parametrization_1D::self%potential )
                call self%potential%set_param_names( ['V0  '], ['V_0'] )
            case(2)
                allocate( power_law_parametrization_1D::self%potential )
                call self%potential%set_param_names( ['V0  ', 'p   '], ['V_0', 'p  '] )
            case(3)
                allocate( taylorseries_parametrization_1D::self%potential )
            case(4)
                allocate( padeseries_parametrization_1D::self%potential )
            case(5)
                allocate( fourier_parametrization_1D::self%potential )
            case(6)
                allocate( exponential_parametrization_2_1D::self%potential )
                call self%potential%set_param_names( ['V0         ', 'EFT_lambda '], ['V_0    ', '\lambda'] )
            case(7)
                allocate( double_exponential_parametrization_1D::self%potential )
                call self%potential%set_param_names( ['V1     ', 'lambda1','V2     ', 'lambda2'], ['V_1      ', '\lambda_1','V_2      ', '\lambda_2'] )
            case(8)
                allocate( cosine_parametrization_1D::self%potential )
                call self%potential%set_param_names( ['V0    ', 'm     ','phi0  ', 'offset'], ['V_1   ', 'm     ','phi_0 ', 'offset'] )
            case(9)
                allocate( axion_parametrization_1D::self%potential )
                call self%potential%set_param_names( ['V0', 'm ','n '], ['V_0', 'm  ','n  '] )
            case(10)
                allocate( power_law_sum_parametrization_1D::self%potential )
                call self%potential%set_param_names( ['L0 ', 'V1 ', 'p1 ', 'V2 ', 'p2 ', 'V3 ', 'p3 '], &
                                                   & ['L_0', 'V_1', 'p_1', 'V_2', 'p_2', 'V_3', 'p_3'] ) 
            case default
                if (temp_feedback > 0) then
                    write(*,'(a,I3)') 'No model corresponding to potential_model =', self%potential_model
                    write(*,'(a)')    'Please select an appropriate model.'
                end if
                eft_error = 1
                return
        end select

        ! initialize the names:
        call self%potential%set_name( 'V', 'V(\phi)' )

        ! additional initialization of the function:
        call self%potential%init_func_from_file( Ini, eft_error )

        ! initialize w0wa flag:
        self%effective_w0wa = .true.

    end subroutine EFTCAMB5eAllocateModelSelection

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that initializes the model parameters based on the values found in an input array.
    subroutine EFTCAMB5eInitModelParameters( self, array )

        implicit none

        class(EFTCAMB_5e)                                      :: self   !< the base class
        real(dl), dimension(self%parameter_number), intent(in) :: array  !< input array with the values of the parameters of the model.

        real(dl), allocatable                                  :: temp(:)
        integer                                                :: i

        if ( self%drag_initial_conditions ) then
            ! read the parameters of the potential:
            allocate( temp(self%parameter_number) )
            do i = 1, self%parameter_number
                temp(i) = array(i)
            end do
        else
            ! read the first parameter:
            self%phidot_ini   = array(1)
            ! read the parameters of the potential:
            allocate( temp(self%parameter_number-1) )
            do i = 1, self%parameter_number -1
                temp(i) = array(i+1)
            end do
        end if

        ! pass the potential parameters:
        if ( self%potential_model<=10 ) then
            call self%potential%init_parameters(temp)
        end if

    end subroutine EFTCAMB5eInitModelParameters

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that reads the parameters of the model from file.
    subroutine EFTCAMB5eInitModelParametersFromFile( self, Ini, eft_error )

        implicit none

        class(EFTCAMB_5e)  :: self      !< the base class
        type(TIniFile)     :: Ini       !< Input ini file
        integer            :: eft_error !< error code: 0 all fine, 1 initialization failed

        ! read standard parameters:
        if ( .not. self%drag_initial_conditions ) then
            self%phidot_ini    = Ini%Read_Double( 'phidot_ini', 0._dl )
        end if
        ! read V(\phi) parameters:
        if ( self%potential_model<=10 ) then
            call self%potential%init_from_file( Ini, eft_error )
        end if

    end subroutine EFTCAMB5eInitModelParametersFromFile

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes the number of parameters of the model.
    subroutine EFTCAMB5eComputeParametersNumber( self )

        implicit none

        class(EFTCAMB_5e)  :: self   !< the base class

        if ( self%drag_initial_conditions ) then
            self%parameter_number = 0
        else
            self%parameter_number = 1
        end if

        ! V(\phi) parameters:
        if ( self%potential_model<=10 ) then
            self%parameter_number = self%parameter_number +self%potential%parameter_number
        end if

    end subroutine EFTCAMB5eComputeParametersNumber

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that prints on the screen feedback information about the model.
    subroutine EFTCAMB5eFeedback( self, print_params )

        implicit none

        class(EFTCAMB_5e)     :: self         !< the base class
        logical, optional     :: print_params !< optional flag that decised whether to print numerical values
                                              !! of the parameters.

        write(*,*)
        write(*,'(a,a)')    '   Model               =  ', self%name
        write(*,'(a,I3)')   '   Number of params    ='  , self%parameter_number
        ! print model functions informations:
        if ( self%potential_model /= 0 ) then
            write(*,*)
            write(*,'(a,I3)')  '   potential_model     =', self%potential_model
        end if
        if ( self%field_min /= -1.d10 ) then
            write(*,'(a,F12.6)')  '   field_min           =', self%field_min
        end if
        if ( self%field_max /= +1.d10 ) then
            write(*,'(a,F12.6)')  '   field_max           =', self%field_max
        end if

        write(*,*)
        if ( .not. self%drag_initial_conditions ) then
            write(*,'(a24,F12.6)') '   phidot_ini          =', self%phidot_ini
        end if

        write(*,*)
        if ( self%potential_model<=10 ) then
            call self%potential%feedback( print_params )
        end if

    end subroutine EFTCAMB5eFeedback

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the i-th parameter name of the model
    subroutine EFTCAMB5eParameterNames( self, i, name )

        implicit none

        class(EFTCAMB_5e)           :: self   !< the base class
        integer     , intent(in)    :: i      !< the index of the parameter
        character(*), intent(out)   :: name   !< the output name of the i-th parameter

        ! check validity of input:
        if ( i<=0 .or. i>self%parameter_number ) then
            write(*,'(a,I3)') 'EFTCAMB error: no parameter corresponding to number ', i
            write(*,'(a,I3)') 'Total number of parameters is ', self%parameter_number
            call MpiStop('EFTCAMB error')
        ! other parameters:
        else if ( i==1 .and. .not. self%drag_initial_conditions ) then
            name = TRIM('phidot_ini')
            return
        ! the other parameters are the potential parameters:
        else
            if ( self%potential_model<=10 ) then
                if ( self%drag_initial_conditions ) then
                    call self%potential%parameter_names( i, name )
                else
                    call self%potential%parameter_names( i-1, name )
                end if
            end if
            return
        end if

    end subroutine EFTCAMB5eParameterNames

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the i-th parameter name of the model
    subroutine EFTCAMB5eParameterNamesLatex( self, i, latexname )

        implicit none

        class(EFTCAMB_5e)          :: self       !< the base class
        integer     , intent(in)    :: i          !< The index of the parameter
        character(*), intent(out)   :: latexname  !< the output latex name of the i-th parameter

        ! check validity of input:
        if ( i<=0 .or. i>self%parameter_number ) then
            write(*,'(a,I3)') 'EFTCAMB error: no parameter corresponding to number ', i
            write(*,'(a,I3)') 'Total number of parameters is ', self%parameter_number
            call MpiStop('EFTCAMB error')
        ! other parameters:
        else if ( i==1 .and. .not. self%drag_initial_conditions ) then
            latexname = TRIM('\dot{\phi}_{\rm ini}')
            return
        ! the other parameters are the potential parameters:
        else
            if ( self%potential_model<=10 ) then
                if ( self%drag_initial_conditions ) then
                    call self%potential%parameter_names_latex( i, latexname )
                else
                    call self%potential%parameter_names_latex( i-1, latexname )
                end if
            end if
            return
        end if

    end subroutine EFTCAMB5eParameterNamesLatex

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the i-th parameter name of the model
    subroutine EFTCAMB5eParameterValues( self, i, value )

        implicit none

        class(EFTCAMB_5e)          :: self   !< the base class
        integer , intent(in)        :: i      !< The index of the parameter
        real(dl), intent(out)       :: value  !< the output value of the i-th parameter

        ! check validity of input:
        if ( i<=0 .or. i>self%parameter_number ) then
            write(*,'(a,I3)') 'EFTCAMB error: no parameter corresponding to number ', i
            write(*,'(a,I3)') 'Total number of parameters is ', self%parameter_number
            call MpiStop('EFTCAMB error')
        ! other parameters:
        else if ( i==1 .and. .not. self%drag_initial_conditions ) then
            value = self%phidot_ini
            return
        ! the other parameters are the w_DE parameters:
        else
            if ( self%potential_model<=10 ) then
                if ( self%drag_initial_conditions ) then
                    call self%potential%parameter_value( i, value )
                else
                    call self%potential%parameter_value( i-1, value )
                end if
            end if
            return
        end if

    end subroutine EFTCAMB5eParameterValues

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes the value of the background EFT functions at a given time.
    subroutine EFTCAMB5eBackgroundEFTFunctions( self, a, eft_par_cache, eft_cache )

        implicit none

        class(EFTCAMB_5e)                            :: self          !< the base class
        real(dl), intent(in)                         :: a             !< the input scale factor
        type(TEFTCAMB_parameter_cache), intent(inout) :: eft_par_cache !< the EFTCAMB parameter cache that contains all the physical parameters.
        type(TEFTCAMB_timestep_cache ), intent(inout) :: eft_cache     !< the EFTCAMB timestep cache that contains all the physical values.

        real(dl) :: x, mu
        integer  :: ind

        ! protect against calling zero and convert to log:
        if ( a < 0.1_dl*exp( self%x_initial ) ) then
            x = log( 0.1_dl*exp( self%x_initial ) )
        else
            x = log(a)
        end if

        ! the eft function Omega is doing the precomputations. All the EFT functions are sampled on the same time grid.
        call self%EFTc%precompute(x, ind, mu )

        eft_cache%EFTc         = self%EFTc%value( x, index=ind, coeff=mu )
        eft_cache%EFTLambda    = self%EFTLambda%value( x, index=ind, coeff=mu )
        eft_cache%EFTcdot      = self%EFTc%first_derivative( x, index=ind, coeff=mu )
        eft_cache%EFTLambdadot = self%EFTLambda%first_derivative( x, index=ind, coeff=mu )
        eft_cache%EFTOmegaV    = 0._dl
        eft_cache%EFTOmegaP    = 0._dl
        eft_cache%EFTOmegaPP   = 0._dl
        eft_cache%EFTOmegaPPP  = 0._dl

    end subroutine EFTCAMB5eBackgroundEFTFunctions

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes the value of the background EFT functions at a given time.
    subroutine EFTCAMB5eSecondOrderEFTFunctions( self, a, eft_par_cache, eft_cache )

        implicit none

        class(EFTCAMB_5e)                            :: self          !< the base class
        real(dl), intent(in)                         :: a             !< the input scale factor
        type(TEFTCAMB_parameter_cache), intent(inout) :: eft_par_cache !< the EFTCAMB parameter cache that contains all the physical parameters.
        type(TEFTCAMB_timestep_cache ), intent(inout) :: eft_cache     !< the EFTCAMB timestep cache that contains all the physical values.

        eft_cache%EFTGamma1V  = 0._dl
        eft_cache%EFTGamma1P  = 0._dl
        eft_cache%EFTGamma2V  = 0._dl
        eft_cache%EFTGamma2P  = 0._dl
        eft_cache%EFTGamma3V  = 0._dl
        eft_cache%EFTGamma3P  = 0._dl
        eft_cache%EFTGamma4V  = 0._dl
        eft_cache%EFTGamma4P  = 0._dl
        eft_cache%EFTGamma4PP = 0._dl
        eft_cache%EFTGamma5V  = 0._dl
        eft_cache%EFTGamma5P  = 0._dl
        eft_cache%EFTGamma6V  = 0._dl
        eft_cache%EFTGamma6P  = 0._dl

    end subroutine EFTCAMB5eSecondOrderEFTFunctions

    ! ---------------------------------------------------------------------------------------------
    !> Function that computes the scalar field potential and its two derivatives
    function EFTCAMB5ePotential( self, phi, deriv )

        implicit none

        class(EFTCAMB_5e)    :: self    !< the base class
        real(dl), intent(in) :: phi     !< the scalar field value
        integer , intent(in) :: deriv   !< the order of the derivative of the field

        real(dl)             :: EFTCAMB5ePotential

        EFTCAMB5ePotential = 1._dl
        if ( deriv == 0 ) then
            if ( self%potential_model<=10 ) then
                EFTCAMB5ePotential = self%potential%value(phi)
            end if
        else if (deriv == 1 ) then
            if ( self%potential_model<=10  ) then
                EFTCAMB5ePotential = self%potential%first_derivative(phi)
            end if
        else if (deriv == 2) then
            if ( self%potential_model<=10  ) then
                EFTCAMB5ePotential = self%potential%second_derivative(phi)
            end if
        end if

    end function EFTCAMB5ePotential

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that initializes the background of the quintessence model.
    subroutine EFTCAMB5eInitBackground( self, params_cache, feedback_level, success, outroot )

        implicit none

        class(EFTCAMB_5e)                             :: self           !< the base class
        type(TEFTCAMB_parameter_cache), intent(inout) :: params_cache   !< a EFTCAMB parameter cache containing cosmological parameters
        integer                       , intent(in)    :: feedback_level !< level of feedback from the background code. 0=none; 1=some; 2=chatty.
        logical                       , intent(out)   :: success        !< wether the background initialization succeded or not
        character(LEN=*), optional    , intent(in)    :: outroot        !< the output root for the debug files

        real(dl) :: TempMin, TempMax, debug_phi, phi_ini, H02, alpha
        real(dl) :: temp_c, temp_L, temp_cdot, temp_Ldot
        integer  :: Debug_MaxNum, Debug_n

        ! some feedback:
        if ( feedback_level>1 ) then
            write(*,'(a)') "***************************************************************"
            write(*,'(a)') ' EFTCAMB Quintessence background solver'
            write(*,'(a)')
        end if

        if ( DebugEFTCAMB .or. feedback_level>2 ) then
            call params_cache%print()
        end if

        ! initialize interpolating functions:
        call self%EFTLambda%initialize ( self%interpolation_num_points, self%x_initial, self%x_final )
        call self%EFTc%initialize      ( self%interpolation_num_points, self%x_initial, self%x_final )

        ! if the potential is periodic then set field minimum and maximum from the period:
        if ( self%potential_model==9 .or. self%potential_model==8 ) then
            call self%potential%parameter_value( 2, alpha )
            self%field_min = -const_pi/alpha
            self%field_max = const_pi/alpha
        end if

        ! debug code print the function H0(phi_ini):
        if ( DebugEFTCAMB ) then
            file_debug_1%unit = 32
            call file_debug_1%CreateFile( TRIM(outroot)//'background_5e.dat' )
            print*, 'EFTCAMB DEBUG ( 5e ): Printing H0^2( phi_ini ) results'
            TempMin      = self%field_min
            TempMax      = self%field_max
            Debug_MaxNum = 1000
            do Debug_n = 1, Debug_MaxNum
                debug_phi = TempMin +REAL(Debug_n-1)*(TempMax-TempMin)/REAL(Debug_MaxNum-1)
                call self%solve_background_equations( params_cache, debug_phi, phidot_ini=self%phidot_ini, H02=H02, only_solve=.True., success=success )
                if ( .not. success ) print*, 'Warning Solver failed with debug_phi=', debug_phi
                write(file_debug_1%unit,'(20ES11.2E3)') debug_phi, self%quintessence_potential(debug_phi,0), H02, params_cache%h0_Mpc**2
            end do
            call file_debug_1%close()
        end if

        ! call boundary conditions lookup:
        call self%find_initial_conditions( params_cache, feedback_level, phi_ini, success )
        if ( success ) then
            if ( feedback_level>1 ) write(*,'(a,E13.4)') '   Initial condition phi_ini  = ', phi_ini
        else if ( .not. success ) then
            return
        end if

        ! debug code to print all background quantities:
        if ( DebugEFTCAMB ) then
            print*, 'EFTCAMB DEBUG ( 5e ): Printing background results'
            file_debug_2%unit = 33
            file_debug_3%unit = 34
            file_debug_4%unit = 35
            call file_debug_2%CreateFile( TRIM(outroot)//'background_5e_solution_1.dat' )
            call file_debug_3%CreateFile( TRIM(outroot)//'background_5e_solution_2.dat' )
            call file_debug_4%CreateFile( TRIM(outroot)//'background_5e_solution_3.dat' )
            write (file_debug_2%unit,'(a)')         '# x a z phi phi_prime phi_prime_prime H2 DHoH Hdot grho_matter gpres_matter'
            write (file_debug_3%unit,'(a)')         '# x a z Ek_phi Ep_phi DEk_DN DEp_DN rho_phi pres_phi Drho_phi_DN Dpres_phi_DN energy_conservation w_phi EFT_lambda'
            write (file_debug_4%unit,'(a)')         '# x a z omega_r_t omega_m_t omega_nu_t omega_phi_t omega_tot_t'
            call self%solve_background_equations( params_cache, phi_ini, phidot_ini=self%phidot_ini, H02=H02, only_solve=.False., success=success )
            call file_debug_2%close()
            call file_debug_3%close()
            call file_debug_4%close()
        end if

        ! solve the background equations and store the solution:
        call self%solve_background_equations( params_cache, phi_ini, phidot_ini=self%phidot_ini, H02=H02, only_solve=.False., success=success )

        ! calculate effective w0wa parameters:
        temp_c = self%EFTc%value(0._dl)
        temp_L = self%EFTLambda%value(0._dl)
        temp_cdot = self%EFTc%first_derivative(0._dl)
        temp_Ldot = self%EFTLambda%first_derivative(0._dl)
        self%w0 = temp_L / (2._dl*temp_c - temp_L)
        self%wa = - 2._dl / params_cache%h0_Mpc * (temp_c*temp_Ldot - temp_cdot*temp_L) / (temp_L - 2._dl*temp_c)**2

        ! feedback:
        if ( feedback_level>1 ) then
            write(*,'(a,E13.4)') '   Effective w0 = ', self%w0
            write(*,'(a,E13.4)') '   Effective wa = ', self%wa
        end if

    end subroutine EFTCAMB5eInitBackground

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that solves the 5e background equations.
    subroutine EFTCAMB5eSolveBackgroundEquations( self, params_cache, phi_ini, phidot_ini, H02, only_solve, success )

        implicit none

        class(EFTCAMB_5e)                            :: self          !< the base class.
        type(TEFTCAMB_parameter_cache), intent(in)    :: params_cache  !< a EFTCAMB parameter cache containing cosmological parameters.
        real(dl), intent(in)                         :: phi_ini       !< the initial value of the field. Refer to the numerical notes for the details.
        real(dl), intent(in)                         :: phidot_ini    !< the initial value of the field derivative. Refer to the numerical notes for the details.
        real(dl), intent(out)                        :: H02           !< value of the squared Hubble constant today.
        logical , optional                           :: only_solve    !< logical flag that tells the code wether to compute only the solution or also the EFT functions.
        logical , intent(out)                        :: success       !< whether the calculation ended correctly or not

        integer, parameter :: num_eq = 2   !<  Number of equations

        real(dl) :: y(num_eq), ydot(num_eq)

        ! odepack quantities:
        integer  :: itol, itask, istate, iopt, LRN, LRS, LRW, LIS, LIN, LIW, JacobianMode, i
        real(dl) :: rtol, atol, t1, t2, t2_temp, B
        real(dl), allocatable :: rwork(:)
        integer,  allocatable :: iwork(:)

        ! background quantities that are shared among the derivs and output routines:
        real(dl) :: a, a2, grhob_t, grhoc_t, grhor_t, grhog_t, grhonu_tot, gpinu_tot
        real(dl) :: grhonu, gpinu, grhormass_t, grho_matter, gpres_matter, grho_m_temp
        real(dl) :: phi, phi_prime, phi_prime_prime, coeff, H2, Hdot, omega_m
        real(dl) :: Ek_phi, Ep_phi, DEk_DN, DEp_DN, rho_phi, pres_phi, w_phi, Drho_phi_DN, Dpres_phi_DN, energy_conservation
        real(dl) :: omega_r_t, omega_m_t, omega_phi_t, omega_nu_t, omega_tot_t, EFT_lambda
        integer  :: nu_i

        ! digest the input:
        if ( .not. present(only_solve) ) only_solve = .False.

        ! set omega_m, notice that we do not include radiation for compatibility with CAMB
        !omega_m = params_cache%omegac +params_cache%omegab +params_cache%omegan +params_cache%omegag +params_cache%omegar
        omega_m = params_cache%omegac +params_cache%omegab +params_cache%omegan

        ! store the first value of the EFT functions:
        t1  = self%EFTc%x(1)

        ! Set initial conditions:
        y(1) = phi_ini
        ! compute drag initial conditions:
        if ( self%drag_initial_conditions ) then
            ! call derivs to compute everything:
            call derivs( num_eq, t1, y, ydot )
            y(2) = +( 6._dl*(1._dl-omega_m)*a2*params_cache%h0_Mpc**2*self%quintessence_potential(phi_ini,deriv=1) )/( +gpres_matter -grho_matter -6._dl*(1._dl-omega_m)*a2*params_cache%h0_Mpc**2*self%quintessence_potential(phi_ini,deriv=0) )
        else
            y(2) = phidot_ini
        end if

        ! now set output vaue:
        H02 = 0._dl

        ! Initialize DLSODA:
        ! set-up the relative and absolute tollerances:
        itol = 1
        rtol = 1.d-10
        atol = 1.d-14
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
        IWORK(6) = 1000   ! maximum number of (internally defined) steps allowed during one call to the solver. The default value is 500.
        IWORK(7) = 0      ! maximum number of messages printed (per problem) warning that T + H = T on a step (H = step size). This must be positive to result in a non-default value.  The default value is 10.
        IWORK(8) = 0      ! the maximum order to be allowed for the nonstiff (Adams) method.  the default value is 12.
        IWORK(9) = 0      ! the maximum order to be allowed for the stiff (BDF) method.  The default value is 5.
        ! additional lsoda stuff:
        CALL XSETF(0) ! suppress odepack printing
        ! Jacobian mode: 1=fullJacobian, 2=not provided
        JacobianMode = 1

        ! compute output EFT functions if needed:
        if ( .not. only_solve ) then
            call output( num_eq, 1, t1, y )
        end if

        ! solve the equations:
        do i=1, self%EFTc%num_points-1

            ! set the time step:
            t1 = self%EFTc%x(i)
            t2 = self%EFTc%x(i+1)
            ! check if the boundary has been reached:
            if ( t1<0._dl .and. t2>=0._dl ) then
                ! solve until zero:
                t2_temp = 0._dl
                call DLSODA ( derivs, num_eq, y, t1, t2_temp, itol, rtol, atol, itask, istate, iopt, RWORK, LRW, IWORK, LIW, jacobian, JacobianMode)
                ! check istate for LSODA good completion:
                if ( istate < 0 ) then
                    if ( istate == -1 ) then
                        istate = 1
                    else
                        success = .False.
                        return
                    end if
                end if
                ! save H0:
                call derivs( num_eq, t2_temp, y, ydot )
                H02 = H2
                ! solve until t2:
                call DLSODA ( derivs, num_eq, y, t2_temp, t2, itol, rtol, atol, itask, istate, iopt, RWORK, LRW, IWORK, LIW, jacobian, JacobianMode)
            else
                call DLSODA ( derivs, num_eq, y, t1, t2, itol, rtol, atol, itask, istate, iopt, RWORK, LRW, IWORK, LIW, jacobian, JacobianMode)
            end if
            ! check istate for LSODA good completion:
            if ( istate < 0 ) then
                if ( istate == -1 ) then
                    istate = 1
                else
                    success = .False.
                    return
                end if
            end if

            ! compute output EFT functions if needed:
            if ( .not. only_solve ) then
                call output( num_eq, i+1, t2, y )
                ! check the results:
                if ( IsNan( self%EFTLambda%yp(i+1) ) .or. IsNan( self%EFTc%yp(i+1) ) ) then
                    success = .False.
                    return
                end if
            end if

        end do

        ! check and return:
        if ( H02 /= 0._dl ) then
            success = .True.
        else
            success = .False.
        end if

        return

    contains

        ! ---------------------------------------------------------------------------------------------
        !> Subroutine that computes y' given y for the background of 5e
        subroutine derivs( num_eq, x, y, ydot )

            implicit none

            integer , intent(in)                     :: num_eq !< number of equations in the ODE system
            real(dl), intent(in)                     :: x      !< time at which the derivatives of the system are computed
            real(dl), intent(in) , dimension(num_eq) :: y      !< input status of the system
            real(dl), intent(out), dimension(num_eq) :: ydot   !< value of the derivative at time x

            ! 0) store variables:
            phi       = y(1)
            phi_prime = y(2)

            ! 1) convert x in a:
            a  = Exp(x)
            a2 = a*a

            ! 2) compute matter densities:
            grhob_t = params_cache%grhob/a         ! 8\pi G_N \rho_b a^2: baryons background density
            grhoc_t = params_cache%grhoc/a         ! 8\pi G_N \rho_{cdm} a^2: cold dark matter background density
            grhor_t = params_cache%grhornomass/a2  ! 8\pi G_N \rho_{\nu} a^2: massless neutrinos background density
            grhog_t = params_cache%grhog/a2        ! 8\pi G_N \rho_{\gamma} a^2: radiation background density

            grhonu_tot = 0._dl
            gpinu_tot  = 0._dl
            if ( params_cache%Num_Nu_Massive /= 0 ) then
                do nu_i = 1, params_cache%Nu_mass_eigenstates
                    grhonu      = 0._dl
                    gpinu       = 0._dl
                    grhormass_t = params_cache%grhormass(nu_i)/a2
                    call ThermalNuBack%rho_P( a*params_cache%nu_masses(nu_i), grhonu, gpinu)
                    grhonu_tot = grhonu_tot + grhormass_t*grhonu ! 8\pi G_N \rho_mnu a^2: massive neutrino background density
                    gpinu_tot  = gpinu_tot  + grhormass_t*gpinu  ! 8\pi G_N P_mnu a^2   : massive neutrino background pressure
                end do
            end if

            grho_matter  = grhob_t +grhoc_t +grhog_t +grhor_t +grhonu_tot
            gpres_matter = gpinu_tot + ( +grhog_t +grhor_t )/3._dl

            ! 3) compute Hubble and derivatives:
            H2   = ( grho_matter +3._dl*(1._dl-omega_m)*params_cache%h0_Mpc**2*a2*self%quintessence_potential(phi,deriv=0) )/3._dl/(1._dl-phi_prime**2/6._dl)
            Hdot = 0.5_dl*( -grho_matter/3._dl -gpres_matter -2._dl/3._dl*H2*phi_prime**2 +2._dl*(1._dl-omega_m)*params_cache%h0_Mpc**2*a2*self%quintessence_potential(phi,deriv=0) )

            ! 4) Get the equations of motion:
            phi_prime_prime = -( Hdot +2._dl*H2 )/H2*phi_prime -3._dl*(1._dl-omega_m)*a2*params_cache%h0_Mpc**2/H2*self%quintessence_potential(phi,deriv=1)

            ydot(1) = phi_prime
            ydot(2) = phi_prime_prime

        end subroutine

        ! ---------------------------------------------------------------------------------------------
        !> Subroutine that computes the Jacobian of the system.
        subroutine jacobian( num_eq, x, y, ml, mu, pd, nrowpd )

            implicit none

            integer                            :: num_eq !< number of components of the Jacobian
            integer                            :: ml     !< ignored
            integer                            :: mu     !< ignored
            integer                            :: nrowpd !< ignored
            real(dl)                           :: x      !< time at which the Jacobian is computed
            real(dl), dimension(num_eq)        :: y      !< input status of the system
            real(dl), dimension(nrowpd,num_eq) :: pd     !< output Jacobian

            real(dl) :: u, v, V0, Vphi, PV, PPV

            ! call derivs to compute everything at the appropriate time step:
            call derivs( num_eq, x, y, ydot )
            ! auxiliary definitions:
            a    = Exp(x)
            a2   = a*a
            u    = phi
            v    = phi_prime
            V0   = 3._dl*(1._dl-omega_m)*a2*params_cache%h0_Mpc**2
            Vphi = self%quintessence_potential(phi,deriv=0)
            PV   = self%quintessence_potential(phi,deriv=1)
            PPV  = self%quintessence_potential(phi,deriv=2)

            ! compute the Jacobian:
            pd(1,1) = 0._dl
            pd(2,1) = 1._dl
            pd(1,2) = (v**2-6._dl)*V0/4._dl/( grho_matter +V0*Vphi )**2*( &
                & +v*PV*( grho_matter +gpres_matter ) -2._dl*V0*PV**2 &
                & +2._dl*PPV*( grho_matter +V0*Vphi ) )
            pd(2,2) = 1._dl/4._dl/( grho_matter +V0*Vphi )*( &
                & -6._dl*( v**2 -2._dl )*V0*Vphi -4._dl*v*V0*PV &
                & +3._dl*(v**2-2._dl)*( gpres_matter -grho_matter) )

        end subroutine jacobian

        ! ---------------------------------------------------------------------------------------------
        !> Subroutine that takes the solution of the 5e equations and computes additional quantities
        !! and the values of the EFT functions.
        subroutine output( num_eq, ind, x, y )

            implicit none

            integer , intent(in)                     :: num_eq !< number of equations in the ODE system.
            integer , intent(in)                     :: ind    !< index of the EFT functions interpolation tables to fill.
            real(dl), intent(in)                     :: x      !< time at which the derivatives of the system are computed.
            real(dl), intent(in) , dimension(num_eq) :: y      !< input status of the system.

            logical :: is_open

            ! 1) call derivs to make sure everything is initialized at the correct time step:
            call derivs( num_eq, x, y, ydot )
            phi             = y(1)
            phi_prime       = y(2)
            phi_prime_prime = ydot(2)

            ! 2) compute the EFT functions:
            self%EFTLambda%y(ind)  = 0.5_dl*H2*phi_prime**2 -3._dl*(1._dl-omega_m)*params_cache%h0_Mpc**2*a2*self%quintessence_potential(phi,deriv=0)
            self%EFTc%y(ind)       = 0.5_dl*H2*phi_prime**2
            ! we have to safeguard against H2 < 0
            if ( H2 > 0._dl ) then
                self%EFTLambda%yp(ind) = sqrt(H2)*( (Hdot-H2)*phi_prime**2 +H2*phi_prime*phi_prime_prime -3._dl*(1._dl-omega_m)*params_cache%h0_Mpc**2*a2*self%quintessence_potential(phi,deriv=1)*phi_prime )
                self%EFTc%yp(ind)      = sqrt(H2)*( (Hdot-H2)*phi_prime**2 +H2*phi_prime*phi_prime_prime )
            else
                self%EFTLambda%yp(ind) = double_NaN
                self%EFTc%yp(ind)      = double_NaN
            end if

            ! 3) compute auxiliary quantities:
            Ek_phi   = 0.5_dl*H2*phi_prime**2/a2
            Ep_phi   = 3._dl*(1._dl-omega_m)*params_cache%h0_Mpc**2*self%quintessence_potential(phi,deriv=0)

            DEk_DN   = +Hdot*phi_prime**2/a2 +H2*phi_prime*phi_prime_prime/a2 -H2*phi_prime**2/a2
            DEp_DN   = 3._dl*(1._dl-omega_m)*params_cache%h0_Mpc**2*self%quintessence_potential(phi,deriv=1)*phi_prime

            rho_phi  = Ek_phi +Ep_phi
            pres_phi = Ek_phi -Ep_phi

            Drho_phi_DN  = DEk_DN +DEp_DN
            Dpres_phi_DN = DEk_DN -DEp_DN

            energy_conservation = ( Drho_phi_DN +3._dl*( rho_phi+pres_phi ) )/rho_phi
            w_phi               = pres_phi/rho_phi
            EFT_lambda          = abs(self%quintessence_potential(phi,deriv=1))/abs(self%quintessence_potential(phi,deriv=0))

            omega_r_t   = (grhog_t +grhor_t)/3._dl/H2
            omega_m_t   = (grhob_t +grhoc_t)/3._dl/H2
            omega_nu_t  = (grhonu_tot)/3._dl/H2
            omega_phi_t = rho_phi*a2/3._dl/H2
            omega_tot_t = omega_r_t +omega_m_t +omega_nu_t +omega_phi_t

            ! 4) debug code:
            if ( DebugEFTCAMB ) then
                inquire( unit=file_debug_2%unit, opened=is_open )
                if ( is_open ) then
                    write (file_debug_2%unit,'(200ES15.4E3)') x, a, 1._dl/a-1._dl, phi, phi_prime, phi_prime_prime, H2, ( H2-( grho_matter +params_cache%grhov*a2 )/3._dl)/H2, Hdot, grho_matter, gpres_matter
                end if
                inquire( unit=file_debug_3%unit, opened=is_open )
                if ( is_open ) then
                    write (file_debug_3%unit,'(200ES15.4E3)') x, a, 1._dl/a-1._dl, Ek_phi, Ep_phi, DEk_DN, DEp_DN, rho_phi, pres_phi, Drho_phi_DN, Dpres_phi_DN, energy_conservation, w_phi, EFT_lambda
                end if
                inquire( unit=file_debug_4%unit, opened=is_open )
                if ( is_open ) then
                    write (file_debug_4%unit,'(200ES15.4E3)') x, a, 1._dl/a-1._dl, omega_r_t, omega_m_t, omega_nu_t, omega_phi_t, omega_tot_t
                end if
            end if

        end subroutine

    end subroutine EFTCAMB5eSolveBackgroundEquations

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that finds the initial conditions for the quintessence field.
    subroutine EFTCAMB5eFindInitialConditions( self, params_cache, feedback_level, phi_ini, success )

        implicit none

        class(EFTCAMB_5e)                            :: self           !< the base class
        type(TEFTCAMB_parameter_cache), intent(in)    :: params_cache   !< a EFTCAMB parameter cache containing cosmological parameters.
        integer                      , intent(in)    :: feedback_level !< level of feedback from the background code. 0=none; 1=some; 2=chatty.
        real(dl), intent(out)                        :: phi_ini        !< the initial value of the field. Refer to the numerical notes for the details.
        logical , intent(out)                        :: success        !< whether the calculation ended correctly or not

        real(dl) :: H0_wanted, H0_temp_1, H0_temp_2, phi_temp, tol, min_field, min_search, phi_exp
        integer  :: ind, num_search

        success = .False.

        ! absolute tollerance in the H0 solution:
        tol       = 1.d-12
        ! other algorithm variables:
        min_field  = 1.d-16
        min_search = 1.d-10
        num_search = 20
        ! copy the desired value, notice that we are solving for the squared Hubble constant:
        H0_wanted = params_cache%h0_Mpc**2

        ! 2) look for the root on the positive side:
        phi_temp  = min_field
        H0_temp_1 = helperH0( phi_temp )

        ! do we already have the solution:
        if ( abs( H0_temp_1-H0_wanted ) <= tol ) then
            phi_ini = phi_temp
            success = .True.
            return
        end if
        ! bracket the root:
        do ind=1, num_search
            phi_exp  = log10(min_search) +REAL(ind-1)/REAL(num_search-1)*( log10(self%field_max) -log10(min_search) )
            phi_temp = 10._dl**phi_exp
            H0_temp_2 = helperH0( phi_temp )
            if ( (H0_temp_1-H0_wanted)*(H0_temp_2-H0_wanted)<=0._dl ) exit
        end do
        ! call the root finder:
        phi_ini = zbrent(helperH0, min_field, phi_temp, tol, H0_wanted, success)
        if ( success ) then
            if ( feedback_level>2 ) then
                write(*,*) '  H0 given:', sqrt(H0_wanted)*(c/1000._dl)
                write(*,*) '  H0 found:', sqrt(helperH0(phi_ini))*(c/1000._dl)
            end if
            return
        else if (.not.success) then
            if ( feedback_level>1 ) write(*,'(a)') '   right side solution not found'
        end if

        ! 3) look for the root on the negative side:
        phi_temp  = -min_field
        H0_temp_1 = helperH0( phi_temp )
        ! do we already have the solution:
        if ( abs( H0_temp_1-H0_wanted ) <= tol ) then
            phi_ini = phi_temp
            success = .True.
            return
        end if
        ! bracket the root:
        do ind=1, num_search
            phi_exp  = log10(min_search) +REAL(ind-1)/REAL(num_search-1)*( log10(abs(self%field_min)) -log10(min_search) )
            phi_temp = -10._dl**phi_exp
            H0_temp_2 = helperH0( phi_temp )
            if ( (H0_temp_1-H0_wanted)*(H0_temp_2-H0_wanted)<=0._dl ) exit
        end do
        ! call the root finder:
        phi_ini = zbrent(helperH0, -min_field, phi_temp, tol, H0_wanted, success)
        if ( success ) then
            if ( feedback_level>2 ) then
                write(*,*) 'H0 given:', sqrt(H0_wanted)*(c/1000._dl)
                write(*,*) 'H0 found:', sqrt(helperH0(phi_ini))*(c/1000._dl)
            end if
            return
        else if (.not.success) then
            if ( feedback_level>1 ) write(*,'(a)') '   left side solution not found'
        end if

        ! 4) look in between:
        phi_ini = zbrent(helperH0, -min_field, +min_field, tol, H0_wanted, success)
        if ( success ) then
            if ( feedback_level>2 ) then
                write(*,*) 'H0 given:', sqrt(H0_wanted)*(c/1000._dl)
                write(*,*) 'H0 found:', sqrt(helperH0(phi_ini))*(c/1000._dl)
            end if
            return
        else if (.not.success) then
            if ( feedback_level>1 ) write(*,'(a)') '   center solution not found'
        end if

        ! 5) in case of failure:
        if (.not.success) then
            if ( feedback_level>1 ) write(*,'(a)') '   initial condition not found'
            return
        end if

    contains

        ! small interface function to feed the root finder:
        function helperH0(phi_ini)
            implicit none
            real(dl) :: phi_ini, helperH0
            call self%solve_background_equations( params_cache, phi_ini, phidot_ini=self%phidot_ini, H02=helperH0, only_solve=.True., success=success )
        end function helperH0

    end subroutine EFTCAMB5eFindInitialConditions

    ! ---------------------------------------------------------------------------------------------

end module EFTCAMB_FM_quintessence
