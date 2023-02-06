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

!> @file 007p5_Omega_Lambda_gamma.f90
!! This file contains the definition of the Omega-Lambda parametrization of Horndeski models.
!! The perturbations are described wuth the gamma functions.
!! Please refer to the numerical notes for details.

!----------------------------------------------------------------------------------------
!! This module contains the definition of the Omega-Lambda parametrization of Horndeski models.
!! The perturbations are described wuth the gamma functions.
!! Please refer to the numerical notes for details.

!> @author Marco Raveri, Simone Peirone

module EFTCAMB_Omega_Lambda_gamma

    use precision
    use IniObjects
    use MpiUtils
    use equispaced_linear_interpolation_1D
    use EFT_def
    use EFTCAMB_mixed_algorithms, only : double_NaN
    use EFTCAMB_cache
    use EFTCAMB_rootfind
    use EFTCAMB_abstract_parametrizations_1D
    use EFTCAMB_parametrizations_1D
    use EFTCAMB_abstract_model_full
    use MassiveNu
    use FileUtils

    implicit none

    private

    public EFTCAMB_OmegaLambda_gamma

    type(TTextFile) :: file_1
    type(TTextFile) :: file_2
    type(TTextFile) :: file_3

    !----------------------------------------------------------------------------------------
    !> This is the type that contains the definition of OL.
    type, extends ( EFTCAMB_full_model ) :: EFTCAMB_OmegaLambda_gamma

        ! the OL functions model selection flags:
        integer  :: OLLambdamodel      !< Model selection flag for OL Lambda.
        integer  :: OLOmegamodel       !< Model selection flag for OL Omega.
        integer  :: OLGamma1model      !< Model selection flag for OL gamma_1.
        integer  :: OLGamma2model      !< Model selection flag for OL gamma_2.
        integer  :: OLGamma3model      !< Model selection flag for OL gamma_3.

        ! the OL functions:
        class( parametrized_function_1D ), allocatable :: OL_Lambda      !< The OL function Lambda.
        class( parametrized_function_1D ), allocatable :: OL_Omega       !< The OL function Omega.
        class( parametrized_function_1D ), allocatable :: OL_Gamma1      !< The OL function Gamma1.
        class( parametrized_function_1D ), allocatable :: OL_Gamma2      !< The OL function Gamma2.
        class( parametrized_function_1D ), allocatable :: OL_Gamma3      !< The OL function Gamma3.

        ! the interpolated EFT functions that come out of the background solver:
        type(equispaced_linear_interpolate_function_1D) :: EFTc          !< The interpolated function c (and derivatives).
        type(equispaced_linear_interpolate_function_1D) :: EFTLambda     !< The interpolated function Lambda (and derivatives).

        ! some parameters:
        integer  :: interpolation_num_points = 2500                      !< Number of points sampled by the background solver code.
        real(dl) :: x_initial                = log(10._dl**(-9._dl))     !< log(a start)
        real(dl) :: x_final                  = log(1.1_dl)               !< log(a final)

    contains

        ! initialization of the model:
        procedure :: read_model_selection            => EFTCAMBOLReadModelSelectionFromFile  !< subroutine that reads the parameters of the model from file
        procedure :: allocate_model_selection        => EFTCAMBOLAllocateModelSelection      !< subroutine that allocates the model selection.
        procedure :: init_model_parameters           => EFTCAMBOLInitModelParameters         !< subroutine that initializes the model parameters based on the values found in an input array.
        procedure :: init_model_parameters_from_file => EFTCAMBOLInitModelParametersFromFile !< subroutine that reads the parameters of the model from file.

        ! utility functions:
        procedure :: compute_param_number  => EFTCAMBOLComputeParametersNumber    !< subroutine that computes the number of parameters of the model.
        procedure :: feedback              => EFTCAMBOLFeedback                   !< subroutine that prints on the screen feedback information about the model.
        procedure :: parameter_names       => EFTCAMBOLParameterNames             !< subroutine that returns the i-th parameter name of the model.
        procedure :: parameter_names_latex => EFTCAMBOLParameterNamesLatex        !< subroutine that returns the i-th parameter name of the model.
        procedure :: parameter_values      => EFTCAMBOLParameterValues            !< subroutine that returns the i-th parameter value.

        ! CAMB related procedures:
        procedure :: compute_background_EFT_functions  => EFTCAMBOLBackgroundEFTFunctions   !< subroutine that computes the value of the background EFT functions at a given time.
        procedure :: compute_secondorder_EFT_functions => EFTCAMBOLSecondOrderEFTFunctions  !< subroutine that computes the value of the second order EFT functions at a given time.

        ! background solver:
        procedure :: initialize_background             => EFTCAMBOLInitBackground              !< subroutine that initializes the background of the Omega Lambda model.
        procedure :: solve_background_equations        => EFTCAMBOLSolveBackgroundEquations    !< subroutine that solves the OmegaLambda background equations.

    end type EFTCAMB_OmegaLambda_gamma

    ! ---------------------------------------------------------------------------------------------

contains

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that reads the parameters of the model from file.
    subroutine EFTCAMBOLReadModelSelectionFromFile( self, Ini, eft_error )

        implicit none

        class(EFTCAMB_OmegaLambda_gamma) :: self      !< the base class
        type(TIniFile)                   :: Ini       !< Input ini file
        integer                          :: eft_error !< error code: 0 all fine, 1 initialization failed

        ! read model selection flags:
        self%OLLambdamodel     = Ini%Read_Int( 'OLLambdamodel', 0 )
        self%OLOmegamodel      = Ini%Read_Int( 'OLOmegamodel' , 0 )
        self%OLGamma1model     = Ini%Read_Int( 'OLGamma1model', 0 )
        self%OLGamma2model     = Ini%Read_Int( 'OLGamma2model', 0 )
        self%OLGamma3model     = Ini%Read_Int( 'OLGamma3model', 0 )

    end subroutine EFTCAMBOLReadModelSelectionFromFile

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that allocates the model selection.
    subroutine EFTCAMBOLAllocateModelSelection( self, Ini, eft_error )

        implicit none

        class(EFTCAMB_OmegaLambda_gamma) :: self      !< the base class
        type(TIniFile)                   :: Ini       !< Input ini file
        integer                          :: eft_error !< error code: 0 all fine, 1 initialization failed

        integer                     :: temp_feedback 

        ! get feedback flag:
        temp_feedback = Ini%Read_Int('feedback_level', 0)

        ! allocate OL_Lambda:
        call allocate_parametrized_1D_function( self%OL_Lambda, self%OLLambdamodel, 'OLLambda', '\tilde{\lambda}', eft_error, temp_feedback )
        if ( eft_error == 1 ) return
        ! allocate OL_Omega:
        call allocate_parametrized_1D_function( self%OL_Omega, self%OLOmegamodel, 'OLOmega', '\Omega', eft_error, temp_feedback )
        if ( eft_error == 1 ) return
        ! allocate OL_Gamma1:
        call allocate_parametrized_1D_function( self%OL_Gamma1, self%OLGamma1model, 'OLGamma1', '{\gamma_1}', eft_error, temp_feedback )
        if ( eft_error == 1 ) return
        ! allocate OL_Gamma2:
        call allocate_parametrized_1D_function( self%OL_Gamma2, self%OLGamma2model, 'OLGamma2', '{\gamma_2}', eft_error, temp_feedback )
        if ( eft_error == 1 ) return
        ! allocate OL_Gamma3:
        call allocate_parametrized_1D_function( self%OL_Gamma3, self%OLGamma3model, 'OLGamma3', '{\gamma_3}', eft_error, temp_feedback )
        if ( eft_error == 1 ) return

        ! additional initialization of the function:
        call self%OL_Lambda%init_func_from_file ( Ini, eft_error )
        call self%OL_Omega%init_func_from_file  ( Ini, eft_error )
        call self%OL_Gamma1%init_func_from_file ( Ini, eft_error )
        call self%OL_Gamma2%init_func_from_file ( Ini, eft_error )
        call self%OL_Gamma3%init_func_from_file ( Ini, eft_error )

    end subroutine EFTCAMBOLAllocateModelSelection

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that initializes the model parameters based on the values found in an input array.
    subroutine EFTCAMBOLInitModelParameters( self, array )

        implicit none

        class(EFTCAMB_OmegaLambda_gamma)                             :: self   !< the base class
        real(dl), dimension(self%parameter_number), intent(in) :: array  !< input array with the values of the parameters of the model.

        real(dl), allocatable, dimension(:)                    :: temp
        integer :: num_params_function, num_params_temp, i

        num_params_temp     = 1

        ! first elements are w_DE parameters:
        num_params_function = self%OL_Lambda%parameter_number
        allocate( temp(num_params_function) )
        do i = 1, num_params_function
            temp(i)         = array(num_params_temp)
            num_params_temp = num_params_temp +1
        end do
        call self%OL_Lambda%init_parameters(temp)
        deallocate( temp )
        ! then OL_Omega parameters:
        num_params_function = self%OL_Omega%parameter_number
        allocate( temp(num_params_function) )
        do i = 1, num_params_function
            temp(i)         = array(num_params_temp)
            num_params_temp = num_params_temp +1
        end do
        call self%OL_Omega%init_parameters(temp)
        deallocate(temp)
        ! then OL_Gamma1 parameters:
        num_params_function = self%OL_Gamma1%parameter_number
        allocate( temp(num_params_function) )
        do i = 1, num_params_function
            temp(i)         = array(num_params_temp)
            num_params_temp = num_params_temp +1
        end do
        call self%OL_Gamma1%init_parameters(temp)
        deallocate(temp)
        ! then OL_Gamma2 parameters:
        num_params_function = self%OL_Gamma2%parameter_number
        allocate( temp(num_params_function) )
        do i = 1, num_params_function
            temp(i)         = array(num_params_temp)
            num_params_temp = num_params_temp +1
        end do
        call self%OL_Gamma2%init_parameters(temp)
        deallocate(temp)
        ! then OL_Gamma3 parameters:
        num_params_function = self%OL_Gamma3%parameter_number
        allocate( temp(num_params_function) )
        do i = 1, num_params_function
            temp(i)         = array(num_params_temp)
            num_params_temp = num_params_temp +1
        end do
        call self%OL_Gamma3%init_parameters(temp)
        deallocate(temp)

        ! now check the length of the parameters:
        if ( num_params_temp-1 /= self%parameter_number ) then
            write(*,*) 'In EFTCAMBOLInitModelParameters:'
            write(*,*) 'Length of num_params_temp and self%parameter_number do not coincide.'
            write(*,*) 'num_params_temp:', num_params_temp-1
            write(*,*) 'self%parameter_number:', self%parameter_number
            call MpiStop('EFTCAMB error')
        end if

    end subroutine EFTCAMBOLInitModelParameters

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that reads the parameters of the model from file.
    subroutine EFTCAMBOLInitModelParametersFromFile( self, Ini, eft_error )

        implicit none

        class(EFTCAMB_OmegaLambda_gamma)  :: self      !< the base class
        type(TIniFile)                    :: Ini       !< Input ini file
        integer                           :: eft_error !< error code: 0 all fine, 1 initialization failed

        call self%OL_Lambda%init_from_file ( Ini, eft_error )
        call self%OL_Omega%init_from_file  ( Ini, eft_error )
        call self%OL_Gamma1%init_from_file ( Ini, eft_error )
        call self%OL_Gamma2%init_from_file ( Ini, eft_error )
        call self%OL_Gamma3%init_from_file ( Ini, eft_error )

    end subroutine EFTCAMBOLInitModelParametersFromFile

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes the number of parameters of the model.
    subroutine EFTCAMBOLComputeParametersNumber( self )

        implicit none

        class(EFTCAMB_OmegaLambda_gamma)  :: self   !< the base class

        self%parameter_number = 0
        self%parameter_number = self%parameter_number +self%OL_Lambda%parameter_number
        self%parameter_number = self%parameter_number +self%OL_Omega%parameter_number
        self%parameter_number = self%parameter_number +self%OL_Gamma1%parameter_number
        self%parameter_number = self%parameter_number +self%OL_Gamma2%parameter_number
        self%parameter_number = self%parameter_number +self%OL_Gamma3%parameter_number

    end subroutine EFTCAMBOLComputeParametersNumber

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that prints on the screen feedback information about the model.
    subroutine EFTCAMBOLFeedback( self, print_params )

        implicit none

        class(EFTCAMB_OmegaLambda_gamma)  :: self         !< the base class
        logical, optional           :: print_params !< optional flag that decised whether to print numerical values
                                                    !! of the parameters.

        ! print general model informations:
        write(*,*)
        write(*,'(a,a)')    '   Model               =  ', self%name
        write(*,'(a,I3)')   '   Number of params    ='  , self%parameter_number

        ! print model functions informations:
        write(*,*)
        if ( self%OLLambdamodel     /= 0 ) write(*,'(a,I3)') '   OLLambdamodel     =', self%OLLambdamodel
        if ( self%OLOmegamodel      /= 0 ) write(*,'(a,I3)') '   OLOmegamodel      =', self%OLOmegamodel
        if ( self%OLGamma1model     /= 0 ) write(*,'(a,I3)') '   OLGamma1model     =', self%OLGamma1model
        if ( self%OLGamma2model     /= 0 ) write(*,'(a,I3)') '   OLGamma2model     =', self%OLGamma2model
        if ( self%OLGamma3model     /= 0 ) write(*,'(a,I3)') '   OLGamma3model     =', self%OLGamma3model

        write(*,*)
        ! print functions informations:
        call self%OL_Lambda%feedback ( print_params )
        call self%OL_Omega%feedback  ( print_params )
        call self%OL_Gamma1%feedback ( print_params )
        call self%OL_Gamma2%feedback ( print_params )
        call self%OL_Gamma3%feedback ( print_params )

    end subroutine EFTCAMBOLFeedback

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the i-th parameter name of the model
    subroutine EFTCAMBOLParameterNames( self, i, name )

        implicit none

        class(EFTCAMB_OmegaLambda_gamma)  :: self   !< the base class
        integer     , intent(in)    :: i      !< the index of the parameter
        character(*), intent(out)   :: name   !< the output name of the i-th parameter

        integer  :: Nw, NM, NK, NB, NT
        integer  :: j

        ! compute the incremental number of parameters:
        Nw = self%OL_Lambda%parameter_number
        NM = Nw + self%OL_Omega%parameter_number
        NK = NM + self%OL_Gamma1%parameter_number
        NB = NK + self%OL_Gamma2%parameter_number
        NT = NB + self%OL_Gamma3%parameter_number

        ! check validity of input:
        if ( i > self%parameter_number .or. i <= 0 ) then
            write(*,'(a,I3)') 'No parameter corresponding to: ', i
            write(*,'(a,I3)') 'Total number of parameters is: ', self%parameter_number
            call MpiStop('EFTCAMB error')

        ! parameter from Lambda function
        else if ( i <= Nw ) then
            do j = 1, self%OL_Lambda%parameter_number
                if ( i == j ) call self%OL_Lambda%parameter_names( j, name )
            end do
            return

        ! parameter from OL_Omega function
        else if ( i <= NM ) then
            do j = 1, self%OL_Omega%parameter_number
                if ( i-Nw == j ) call self%OL_Omega%parameter_names( j, name )
            end do
            return

        ! parameter from OL_Gamma1 function
        else if ( i <= NK) then
            do j = 1, self%OL_Gamma1%parameter_number
                if ( i-NM == j ) call self%OL_Gamma1%parameter_names( j, name )
            end do
            return

        !parameter from OL_Gamma2 function
        else if ( i <= NB ) then
            do j = 1, self%OL_Gamma2%parameter_number
                if ( i-NK == j ) call self%OL_Gamma2%parameter_names( j, name )
            end do
            return

        !parameter from OL_Gamma3 function
        else if ( i <= NT ) then
            do j = 1, self%OL_Gamma3%parameter_number
                if ( i-NB == j ) call self%OL_Gamma3%parameter_names( j, name )
            end do
            return

        end if

    end subroutine EFTCAMBOLParameterNames

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the i-th parameter name of the model
    subroutine EFTCAMBOLParameterNamesLatex( self, i, latexname )

        implicit none

        class(EFTCAMB_OmegaLambda_gamma)  :: self       !< the base class
        integer     , intent(in)    :: i          !< The index of the parameter
        character(*), intent(out)   :: latexname  !< the output latex name of the i-th parameter

        integer  :: Nw, NM, NK, NB, NT
        integer  :: j

        ! compute the incremental number of parameters:
        Nw = self%OL_Lambda%parameter_number
        NM = Nw + self%OL_Omega%parameter_number
        NK = NM + self%OL_Gamma1%parameter_number
        NB = NK + self%OL_Gamma2%parameter_number
        NT = NB + self%OL_Gamma3%parameter_number

        ! check validity of input:
        if ( i > self%parameter_number .or. i <= 0 ) then
            write(*,'(a,I3)') 'No parameter corresponding to: ', i
            write(*,'(a,I3)') 'Total number of parameters is: ', self%parameter_number
            call MpiStop('EFTCAMB error')

        ! parameter from Lambda function
        else if ( i <= Nw ) then
            do j = 1, self%OL_Lambda%parameter_number
                if ( i == j ) call self%OL_Lambda%parameter_names_latex( j, latexname )
            end do
            return

        ! parameter from OL_Omega function
        else if ( i <= NM ) then
            do j = 1, self%OL_Omega%parameter_number
                if ( i-Nw == j ) call self%OL_Omega%parameter_names_latex( j, latexname )
            end do
            return

        ! parameter from OL_Gamma1 function
        else if ( i <= NK) then
            do j = 1, self%OL_Gamma1%parameter_number
                if ( i-NM == j ) call self%OL_Gamma1%parameter_names_latex( j, latexname )
            end do
            return

        !parameter from OL_Gamma2 function
        else if ( i <= NB ) then
            do j = 1, self%OL_Gamma2%parameter_number
                if ( i-NK == j ) call self%OL_Gamma2%parameter_names_latex( j, latexname )
            end do
            return

        !parameter from OL_Gamma3 function
        else if ( i <= NT ) then

            do j = 1, self%OL_Gamma3%parameter_number
                if ( i-NB == j ) call self%OL_Gamma3%parameter_names_latex( j, latexname )
            end do
            return

        end if

    end subroutine EFTCAMBOLParameterNamesLatex

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the i-th parameter name of the model
    subroutine EFTCAMBOLParameterValues( self, i, value )

        implicit none

        class(EFTCAMB_OmegaLambda_gamma)  :: self   !< the base class
        integer , intent(in)        :: i      !< The index of the parameter
        real(dl), intent(out)       :: value  !< the output value of the i-th parameter

        integer  :: Nw, NM, NK, NB, NT
        integer  :: j

        ! compute the incremental number of parameters:
        Nw = self%OL_Lambda%parameter_number
        NM = Nw + self%OL_Omega%parameter_number
        NK = NM + self%OL_Gamma1%parameter_number
        NB = NK + self%OL_Gamma2%parameter_number
        NT = NB + self%OL_Gamma3%parameter_number

        ! check validity of input:
        if ( i > self%parameter_number .or. i <= 0 ) then
            write(*,'(a,I3)') 'No parameter corresponding to: ', i
            write(*,'(a,I3)') 'Total number of parameters is: ', self%parameter_number
            call MpiStop('EFTCAMB error')

        ! parameter from wDE function
        else if ( i <= Nw ) then
            do j = 1, self%OL_Lambda%parameter_number
                if ( i == j ) call self%OL_Lambda%parameter_value( j, value )
            end do
            return

        ! parameter from OL_Omega function
        else if ( i <= NM ) then
            do j = 1, self%OL_Omega%parameter_number
                if ( i-Nw == j ) call self%OL_Omega%parameter_value( j, value )
            end do
            return

        ! parameter from OL_Gamma1 function
        else if ( i <= NK) then
            do j = 1, self%OL_Gamma1%parameter_number
                if ( i-NM == j ) call self%OL_Gamma1%parameter_value( j, value )
            end do
            return

        !parameter from OL_Gamma2 function
        else if ( i <= NB ) then
            do j = 1, self%OL_Gamma2%parameter_number
                if ( i-NK == j ) call self%OL_Gamma2%parameter_value( j, value )
            end do
            return

        !parameter from OL_Gamma3 function
        else if ( i <= NT ) then

            do j = 1, self%OL_Gamma3%parameter_number
                if ( i-NB == j ) call self%OL_Gamma3%parameter_value( j, value )
            end do
            return

        end if

    end subroutine EFTCAMBOLParameterValues

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes the value of the background EFT functions at a given time.
    subroutine EFTCAMBOLBackgroundEFTFunctions( self, a, eft_par_cache, eft_cache )

        implicit none

        class(EFTCAMB_OmegaLambda_gamma)              :: self          !< the base class
        real(dl), intent(in)                          :: a             !< the input scale factor
        type(TEFTCAMB_parameter_cache), intent(inout) :: eft_par_cache !< the EFTCAMB parameter cache that contains all the physical parameters.
        type(TEFTCAMB_timestep_cache ), intent(inout) :: eft_cache     !< the EFTCAMB timestep cache that contains all the physical values.

        real(dl) :: x, mu, omega_m
        integer  :: ind

        ! compute the EFT functions:
        eft_cache%EFTOmegaV    = self%OL_Omega%value(a)
        eft_cache%EFTOmegaP    = self%OL_Omega%first_derivative(a)
        eft_cache%EFTOmegaPP   = self%OL_Omega%second_derivative(a)
        eft_cache%EFTOmegaPPP  = self%OL_Omega%third_derivative(a)
        eft_cache%EFTOmegaPPPP = self%OL_Omega%fourth_derivative(a)

        ! compute interpolated c and Lambda:
        ! protect against calling zero and convert to log:
        if ( a < exp( self%x_initial ) ) then
            x = log( exp( self%x_initial ) )
        else
            x = log(a)
        end if
        call self%EFTc%precompute(x, ind, mu )

        eft_cache%EFTc         = self%EFTc%value( x, index=ind, coeff=mu )
        eft_cache%EFTcdot      = self%EFTc%first_derivative( x, index=ind, coeff=mu )

        eft_cache%EFTLambda    = self%EFTLambda%value( x, index=ind, coeff=mu )
        eft_cache%EFTLambdadot = self%EFTLambda%first_derivative( x, index=ind, coeff=mu )

    end subroutine EFTCAMBOLBackgroundEFTFunctions

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes the value of the background EFT functions at a given time.
    subroutine EFTCAMBOLSecondOrderEFTFunctions( self, a, eft_par_cache, eft_cache )

        implicit none

        class(EFTCAMB_OmegaLambda_gamma)              :: self          !< the base class
        real(dl), intent(in)                          :: a             !< the input scale factor
        type(TEFTCAMB_parameter_cache), intent(inout) :: eft_par_cache !< the EFTCAMB parameter cache that contains all the physical parameters.
        type(TEFTCAMB_timestep_cache ), intent(inout) :: eft_cache     !< the EFTCAMB timestep cache that contains all the physical values.

        ! compute the EFT functions:
        eft_cache%EFTGamma1V  = self%OL_Gamma1%value(a)
        eft_cache%EFTGamma1P  = self%OL_Gamma1%first_derivative(a)
        eft_cache%EFTGamma2V  = self%OL_Gamma2%value(a)
        eft_cache%EFTGamma2P  = self%OL_Gamma2%first_derivative(a)
        eft_cache%EFTGamma2PP = self%OL_Gamma2%second_derivative(a)
        eft_cache%EFTGamma2PPP = self%OL_Gamma2%third_derivative(a)
        eft_cache%EFTGamma3V  = self%OL_Gamma3%value(a)
        eft_cache%EFTGamma3P  = self%OL_Gamma3%first_derivative(a)
        eft_cache%EFTGamma3PP = self%OL_Gamma3%second_derivative(a)
        eft_cache%EFTGamma3PPP = self%OL_Gamma3%third_derivative(a)
        eft_cache%EFTGamma3PPPP = self%OL_Gamma3%fourth_derivative(a)
        eft_cache%EFTGamma4V  = -eft_cache%EFTGamma3V
        eft_cache%EFTGamma4P  = -eft_cache%EFTGamma3P
        eft_cache%EFTGamma4PP = -self%OL_Gamma3%second_derivative(a)
        eft_cache%EFTGamma5V  = +0.5_dl*eft_cache%EFTGamma3V
        eft_cache%EFTGamma5P  = +0.5_dl*eft_cache%EFTGamma3P
        eft_cache%EFTGamma6V  = 0._dl
        eft_cache%EFTGamma6P  = 0._dl


    end subroutine EFTCAMBOLSecondOrderEFTFunctions

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that initializes the background of the Omega Lambda model.
    subroutine EFTCAMBOLInitBackground( self, params_cache, feedback_level, success, outroot )

        implicit none

        class(EFTCAMB_OmegaLambda_gamma)              :: self           !< the base class
        type(TEFTCAMB_parameter_cache), intent(inout) :: params_cache   !< a EFTCAMB parameter cache containing cosmological parameters
        integer                       , intent(in)    :: feedback_level !< level of feedback from the background code. 0=none; 1=some; 2=chatty.
        logical                       , intent(out)   :: success        !< wether the background initialization succeded or not
        character(LEN=*), optional    , intent(in)    :: outroot        !< the output root for the debug files

        ! some feedback:
        if ( feedback_level>1 ) then
            write(*,'(a)') "***************************************************************"
            write(*,'(a)') ' EFTCAMB Omega Lambda Pure EFT background solver'
            write(*,'(a)')
        end if

        if ( DebugEFTCAMB .or. feedback_level>2 ) then
            call params_cache%print()
        end if

        ! initialize interpolating functions:
        call self%EFTc%initialize ( self%interpolation_num_points, self%x_initial, self%x_final )
        call self%EFTLambda%initialize ( self%interpolation_num_points, self%x_initial, self%x_final )

        ! debug code to print all background quantities:
        if ( DebugEFTCAMB ) then
            print*, 'EFTCAMB DEBUG ( OL ): Printing background results'
            file_1%unit  = 111
            file_2%unit  = 222
            file_3%unit  = 333
            call file_1%CreateFile(  TRIM(outroot)//'background_OL_solution_1.dat'  )
            call file_2%CreateFile(  TRIM(outroot)//'background_OL_solution_2.dat'  )
            call file_3%CreateFile(  TRIM(outroot)//'background_OL_solution_3.dat'  )
            write (file_1%unit,'(a)') '# x a y ydot h2 h2_lcdm delta_h2 w_eff w_DE p_DE rho_DE'
            write (file_2%unit,'(a)') '# x a h2 c cdot Lambda LambdaDot LambdaLCDM DeltaLambda LambdaInput'
            write (file_3%unit,'(a)') '# x a h2 omegaDE omegam omegar omeganu'
        end if

        ! solve the background equations and store the solution:
        call self%solve_background_equations( params_cache, success=success, feedback_level=feedback_level )

        ! close debug files:
        if ( DebugEFTCAMB ) then
            call file_1%close()
            call file_2%close()
            call file_3%close()
        end if

    end subroutine EFTCAMBOLInitBackground

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that solves the Omega Lambda background equations.
    subroutine EFTCAMBOLSolveBackgroundEquations( self, params_cache, success, feedback_level )

        implicit none

        class(EFTCAMB_OmegaLambda_gamma)             :: self             !< the base class.
        type(TEFTCAMB_parameter_cache), intent(in)   :: params_cache     !< a EFTCAMB parameter cache containing cosmological parameters.
        logical , intent(out)                        :: success          !< whether the calculation ended correctly or not
        integer , intent(in)                         :: feedback_level   !< whether the calculation ended correctly or not

        integer, parameter :: num_eq = 1   !<  Number of equations

        real(dl) :: y(num_eq), ydot(num_eq)

        ! odepack quantities:
        integer  :: itol, itask, istate, iopt, LRN, LRS, LRW, LIS, LIN, LIW, JacobianMode, i
        real(dl) :: rtol, atol, t1, t2, B
        real(dl), allocatable :: rwork(:)
        integer,  allocatable :: iwork(:)

        ! background quantities that are shared among the derivs and output routines:
        real(dl) :: a, a2, grhob_t, grhoc_t, grhor_t, grhog_t, grhonu_tot, gpinu_tot
        real(dl) :: grhonu, gpinu, grhormass_t, grho_matter, gpres_matter, grho_m_temp
        real(dl) :: H2, Hdot, omega_m
        real(dl) :: EFTOmega, EFTOmegaP, EFTOmegaPP
        real(dl) :: omega_r_t, omega_m_t, omega_phi_t, omega_nu_t, omega_tot_t
        integer  :: nu_i

        ! other quantities:
        real(dl) :: y1, y2

        ! set omega_m, notice compatibility with CAMB neglecting radiation
        omega_m = params_cache%omegac +params_cache%omegab +params_cache%omegan +params_cache%omegag +params_cache%omegar
        !omega_m = params_cache%omegac +params_cache%omegab +params_cache%omegan

        ! Initialize DLSODA:
        ! set-up the relative and absolute tollerances:
        itol = 1
        rtol = 1.d-14
        atol = 1.d-16

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
        IWORK(6) = 10000  ! maximum number of (internally defined) steps allowed during one call to the solver. The default value is 500.
        IWORK(7) = 0      ! maximum number of messages printed (per problem) warning that T + H = T on a step (H = step size). This must be positive to result in a non-default value.  The default value is 10.
        IWORK(8) = 0      ! the maximum order to be allowed for the nonstiff (Adams) method.  the default value is 12.
        IWORK(9) = 0      ! the maximum order to be allowed for the stiff (BDF) method.  The default value is 5.
        ! additional lsoda stuff:
        CALL XSETF(0) ! suppress odepack printing
        ! Jacobian mode: 1=fullJacobian, 2=not provided
        JacobianMode = 1

        ! solve only backward in time
        if ( .False. ) then
            if ( feedback_level>1 ) write(*,'(a)') '   Backward solver'
            ! Set initial conditions:
            y(1) = params_cache%h0_Mpc**2
            ! store the first value of the EFT functions:
            t1  = self%EFTc%x(self%EFTc%num_points)
            call output( num_eq, self%EFTc%num_points, t1, y )
            ! solve the equations:
            do i=self%EFTc%num_points, 2, -1
                ! set the time step:
                t1 = self%EFTc%x(i)
                t2 = self%EFTc%x(i-1)
                ! solve the system:
                call DLSODA ( derivs, num_eq, y, t1, t2, itol, rtol, atol, itask, istate, iopt, RWORK, LRW, IWORK, LIW, jacobian, JacobianMode)
                ! check istate for LSODA good completion:
                if ( istate < 0 ) then
                    if ( istate == -1 ) then
                        write(*,*) 'DLSODA excessive work'
                        istate = 1
                    else
                        success = .False.
                        write(*,*) 'DLSODA failed with code:', istate
                        deallocate(rwork,iwork)
                        return
                    end if
                end if
                ! compute output EFT functions if needed:
                call output( num_eq, i-1, t2, y )
            end do
        end if
        ! solve forward:
        if ( .True. ) then
            if ( feedback_level>1 ) write(*,'(a)') '   Forward solver'
            ! first guess for initial conditions in the GR limit:
            t1 = self%EFTc%x(1)
            y  = 0._dl
            call output( num_eq, 1, t1, y )
            y  = ( +grho_matter/3._dl +params_cache%h0_Mpc**2*(1._dl-omega_m)*( 1._dl+self%OL_Lambda%value(a) )*a2 )/(1._dl +EFTOmega +a*EFTOmegaP )
            ! try to bracket the root for the initial condition search:
            y1 = a**2*y(1)*0.9_dl/params_cache%h0_Mpc**2
            y2 = a**2*y(1)*1.1_dl/params_cache%h0_Mpc**2
            call zbrac( dummy_solve, y1, y2, success, 1._dl )
            if (.not.success) then
                if ( feedback_level>1 ) write(*,'(a)') '   FAILURE in solution bracketing'
                return
            end if
            ! get the initial condition:
            y1 = zbrent(dummy_solve, y1, y2, 1.d-14, 1._dl, success)
            if (.not.success) then
                if ( feedback_level>1 ) write(*,'(a)') '   FAILURE in initial condition finding'
                return
            end if
            ! now store solution:
            istate = 1
            t1     = self%EFTc%x(1)
            y      = y1
            call output( num_eq, 1, t1, y )
            ! solve the equations:
            do i=2, self%EFTc%num_points
                ! set the time step:
                t1 = self%EFTc%x(i-1)
                t2 = self%EFTc%x(i)
                ! solve the system:
                call DLSODA ( derivs, num_eq, y, t1, t2, itol, rtol, atol, itask, istate, iopt, RWORK, LRW, IWORK, LIW, jacobian, JacobianMode)
                ! check istate for LSODA good completion:
                if ( istate < 0 ) then
                    if ( istate == -1 ) then
                        write(*,*) 'DLSODA excessive work'
                        istate = 1
                    else
                        success = .False.
                        write(*,*) 'DLSODA failed with code:', istate
                        return
                    end if
                end if
                ! compute output EFT functions if needed:
                call output( num_eq, i, t2, y )
            end do
        end if

        ! return:
        success = .True.
        deallocate(rwork,iwork)

        return

    contains

        ! ---------------------------------------------------------------------------------------------
        !> Subroutine is wrapping the solver for the initial condition search:
        function dummy_solve( y_ini )

            implicit none

            real(dl), intent(in) :: y_ini        !< initial condition (in this case the initial condition is a scalar)
            real(dl)             :: dummy_solve  !< the value of y(0)=H0^2

            real(dl) :: y_temp(num_eq)

            ! reset the solver (warning, this is working in shared memory with the main subroutine, the other solver options are there)
            istate = 1
            t1     = self%EFTc%x(1)
            y_temp = y_ini
            ! solve the equations:
            do i=2, self%EFTc%num_points
                ! set the time step:
                t1 = min(self%EFTc%x(i-1), 0._dl)
                t2 = min(self%EFTc%x(i), 0._dl)
                ! solve the system:
                call DLSODA ( derivs, num_eq, y_temp, t1, t2, itol, rtol, atol, itask, istate, iopt, RWORK, LRW, IWORK, LIW, jacobian, JacobianMode)
                ! check istate for LSODA good completion:
                if ( istate < 0 ) then
                    if ( istate == -1 ) then
                        write(*,*) 'DLSODA excessive work'
                        istate = 1
                    else
                        success = .False.
                        write(*,*) 'DLSODA failed with code:', istate
                        return
                    end if
                end if
            end do
            dummy_solve = y_temp(1)

        end function dummy_solve

        ! ---------------------------------------------------------------------------------------------
        !> Subroutine that computes y' given y for the background of 5e
        subroutine derivs( num_eq, x, y, ydot )

            implicit none

            integer , intent(in)                     :: num_eq !< number of equations in the ODE system
            real(dl), intent(in)                     :: x      !< time at which the derivatives of the system are computed
            real(dl), intent(in) , dimension(num_eq) :: y      !< input status of the system
            real(dl), intent(out), dimension(num_eq) :: ydot   !< value of the derivative at time x

            real(dl) :: OL_PM_V, OL_AT_V, OL_PM_P, OL_AT_P, OL_PM_PP, OL_AT_PP
            real(dl) :: coeff_1, coeff_2, coeff_3

            ! 1) convert x in a:
            a  = Exp(x)
            a2 = a*a

            ! 2) compute matter densities:
            grhob_t = params_cache%grhob/a         ! 8\pi G_N \rho_b a^2: baryons background density
            grhoc_t = params_cache%grhoc/a         ! 8\pi G_N \rho_{cdm} a^2: cold dark matter background density
            grhor_t = params_cache%grhornomass/a2  ! 8\pi G_N \rho_{\nu} a^2: massless neutrinos background density
            grhog_t = params_cache%grhog/a2        ! 8\pi G_N \rho_{\gamma} a^2: radiation background density

            grhonu_tot    = 0._dl
            gpinu_tot     = 0._dl
            if ( params_cache%Num_Nu_Massive /= 0 ) then
                do nu_i = 1, params_cache%Nu_mass_eigenstates
                    grhonu      = 0._dl
                    gpinu       = 0._dl
                    grhormass_t = params_cache%grhormass(nu_i)/a2
                    call  ThermalNuBack%rho_P(a*params_cache%nu_masses(nu_i), grhonu, gpinu)
                    grhonu_tot = grhonu_tot + grhormass_t*grhonu ! 8\pi G_N \rho_mnu a^2: massive neutrino background density
                    gpinu_tot  = gpinu_tot  + grhormass_t*gpinu  ! 8\pi G_N P_mnu a^2   : massive neutrino background pressure
                end do
            end if

            grho_matter  = grhob_t +grhoc_t +grhog_t +grhor_t +grhonu_tot
            gpres_matter = gpinu_tot + ( +grhog_t +grhor_t )/3._dl

            ! 3) compute Omega and its time derivatives:
            ! compute the EFT functions:
            EFTOmega   = self%OL_Omega%value(a)
            EFTOmegaP  = self%OL_Omega%first_derivative(a)
            EFTOmegaPP = self%OL_Omega%second_derivative(a)

            ! 4) Compute the equations of motion:
            coeff_1 = 1._dl +EFTOmega +0.5_dl*a*EFTOmegaP
            coeff_2 = 1._dl +EFTOmega +2._dl*a*EFTOmegaP +a2*EFTOmegaPP
            coeff_3 = gpres_matter -3._dl*params_cache%h0_Mpc**2*(1._dl-omega_m)*( 1._dl+self%OL_Lambda%value(a) )*a2

            ! 5) Compute the equations of motion:
            ydot(1) = ( -coeff_2*y(1) -coeff_3*a**2/params_cache%h0_Mpc**2 )/coeff_1 +2._dl*y(1)

        end subroutine derivs

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
            real(dl) :: coeff_1, coeff_2

            ! call derivs to compute everything at the appropriate time step:
            call derivs( num_eq, x, y, ydot )
            ! initial computations:
            a  = exp(x)
            a2 = a*a

            coeff_1 = 1._dl +EFTOmega +0.5_dl*a*EFTOmegaP
            coeff_2 = 1._dl +EFTOmega +2._dl*a*EFTOmegaP +a2*EFTOmegaPP

            pd(1,1) = -coeff_2/coeff_1 +2._dl

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

            logical  :: is_open
            real(dl) :: adotoa, adotoa2, Hdot, w_DE, w_eff

            ! 1) call derivs to make sure everything is initialized at the correct time step:
            call derivs( num_eq, x, y, ydot )

            ! 2) compute c and Lambda:
            a                     = exp(x)
            adotoa2               = y(1)/a**2*params_cache%h0_Mpc**2
            self%EFTc%y(ind)      = +1.5_dl*(1._dl +EFTOmega +a*EFTOmegaP )*adotoa2 -0.5_dl*(grho_matter) -1.5_dl*params_cache%h0_Mpc**2*(1._dl-omega_m)*( 1._dl+self%OL_Lambda%value(a) )*a2
            self%EFTLambda%y(ind) = -3._dl*params_cache%h0_Mpc**2*(1._dl-omega_m)*( 1._dl+self%OL_Lambda%value(a) )*a2

            ! 3) compute cdot:
            if ( adotoa2<0._dl ) then
                self%EFTc%yp(ind) = double_NaN
                return
            end if

            adotoa  = sqrt(adotoa2)
            Hdot    = 0.5_dl*(ydot(1)/a**2*params_cache%h0_Mpc**2 -2._dl*adotoa2)

            self%EFTc%yp(ind)  = +3._dl*adotoa*( &
                & -( 1._dl +EFTOmega -0.5_dl*a2*EFTOmegaPP )*adotoa2 &
                & +( 1._dl +EFTOmega +a*EFTOmegaP )*Hdot &
                & +0.5_dl*( grho_matter + gpres_matter ) &
                & -0.5_dl*a**3*params_cache%h0_Mpc**2*(1._dl-omega_m)*self%OL_Lambda%first_derivative(a) )
            self%EFTLambda%yp(ind) = -3._dl*params_cache%h0_Mpc**2*adotoa*(1._dl-omega_m)*( self%OL_Lambda%first_derivative(a) )*a**3

            ! 4) compute auxiliary quantities:
            w_DE       = ( -2._dl*Hdot -adotoa2 -gpres_matter )/( 3._dl*adotoa2 -grho_matter )
            w_eff      = -1._dl/3._dl -2._dl/3._dl*Hdot/adotoa2

            ! 5) debug code:
            if ( DebugEFTCAMB ) then
                 write (file_1%unit,'(20ES15.4E3)') x, a, y, ydot, adotoa**2, ( grho_matter +params_cache%grhov*a2 )/3._dl, (adotoa**2-( grho_matter +params_cache%grhov*a2 )/3._dl)/adotoa**2, w_eff, w_DE, ( -2._dl*Hdot -adotoa2 -gpres_matter ), ( 3._dl*adotoa2 -grho_matter )
                 write (file_2%unit,'(20ES15.4E3)') x, a, adotoa**2, self%EFTc%y(ind), self%EFTc%yp(ind), self%EFTLambda%y(ind), self%EFTLambda%yp(ind), -params_cache%grhov*a2, (self%EFTLambda%y(ind)+params_cache%grhov*a2)/(params_cache%grhov*a2), self%OL_Lambda%value(a)
                 write (file_3%unit,'(20ES15.4E3)') x, a, adotoa**2, ( 3._dl*adotoa2 -grho_matter )/3._dl/adotoa**2, (grhob_t +grhoc_t)/3._dl/adotoa**2, (grhog_t +grhor_t)/3._dl/adotoa**2, (grhonu_tot)/3._dl/adotoa**2
            end if

        end subroutine output

    end subroutine EFTCAMBOLSolveBackgroundEquations

    ! ---------------------------------------------------------------------------------------------

end module EFTCAMB_Omega_Lambda_gamma
