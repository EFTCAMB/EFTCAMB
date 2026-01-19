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

!> @file 007p1_Pure_EFT_std.f90
!! This file contains the definition of the Pure EFT model in which the EFT is described
!! by six functions of time and w_DE. Please refer to the numerical notes for details.


!----------------------------------------------------------------------------------------
!> This module contains the definition of the Pure EFT model in which the EFT is described
!! by six functions of time and w_DE. Please refer to the numerical notes for details.

!> @author Bin Hu, Marco Raveri, Simone Peirone, Gen Ye

module EFTCAMB_pure_EFT_std

    use precision
    use IniObjects
    use MpiUtils
    use EFTCAMB_cache
    use EFT_def
    use EFTCAMB_abstract_parametrizations_1D
    use EFTCAMB_parametrizations_1D
    use EFTCAMB_abstract_model_designer
    use equispaced_linear_interpolation_1D
    use MassiveNu

    implicit none

    private

    public EFTCAMB_std_pure_EFT

    !----------------------------------------------------------------------------------------
    !> This is the pure EFT model with the six functions of time and w_DE.
    type, extends ( EFTCAMB_designer_model ) :: EFTCAMB_std_pure_EFT

        ! the pure EFT functions model selection flags:
        integer  :: PureEFTmodelOmega   !< Model selection flag for Pure EFT Omega.
        integer  :: EFTwDE              !< Model selection flag for Pure EFT w DE.
        integer  :: PureEFTmodelGamma1  !< Model selection flag for Pure EFT Gamma1.
        integer  :: PureEFTmodelGamma2  !< Model selection flag for Pure EFT Gamma2.
        integer  :: PureEFTmodelGamma3  !< Model selection flag for Pure EFT Gamma3.
        integer  :: PureEFTmodelGamma4  !< Model selection flag for Pure EFT Gamma4.
        integer  :: PureEFTmodelGamma5  !< Model selection flag for Pure EFT Gamma5.
        integer  :: PureEFTmodelGamma6  !< Model selection flag for Pure EFT Gamma6.

        ! the pure EFT functions model selection flags:
        integer  :: PureEFTmodelOmega_ODE   !< Model selection flag for Pure EFT Omega with Omega_DE argument.
        integer  :: PureEFTmodelGamma1_ODE  !< Model selection flag for Pure EFT Gamma1 with Omega_DE argument.
        integer  :: PureEFTmodelGamma2_ODE  !< Model selection flag for Pure EFT Gamma2 with Omega_DE argument.
        integer  :: PureEFTmodelGamma3_ODE  !< Model selection flag for Pure EFT Gamma3 with Omega_DE argument.
        integer  :: PureEFTmodelGamma4_ODE  !< Model selection flag for Pure EFT Gamma4 with Omega_DE argument.
        integer  :: PureEFTmodelGamma5_ODE  !< Model selection flag for Pure EFT Gamma5 with Omega_DE argument.
        integer  :: PureEFTmodelGamma6_ODE  !< Model selection flag for Pure EFT Gamma6 with Omega_DE argument.

        ! Each function X is the sum of a function of scale factor, defined by Xmodel, and a function of Omega_DE, defined by Xmodel_ODE. E.g. Gamma1 = F(a) + G(Omega_DE) 

        ! selection flag for Horndeski:
        logical  :: PureEFTHorndeski    !< Selects wether to use the Horndeski bound on EFT functions.

        ! flag deciding whether computes Omega_DE
        logical  :: PureEFTOmegaDE

        ! DE energy density from wDE integration
        type( equispaced_linear_interpolate_function_1D ) :: rhov_t

        ! the pure EFT functions:
        class( parametrized_function_1D ), allocatable :: PureEFTOmega    !< The pure EFT function Omega.
        class( parametrized_function_1D ), allocatable :: PureEFTwDE      !< The pure EFT function w_DE.
        class( parametrized_function_1D ), allocatable :: PureEFTGamma1   !< The pure EFT function Gamma1.
        class( parametrized_function_1D ), allocatable :: PureEFTGamma2   !< The pure EFT function Gamma2.
        class( parametrized_function_1D ), allocatable :: PureEFTGamma3   !< The pure EFT function Gamma3.
        class( parametrized_function_1D ), allocatable :: PureEFTGamma4   !< The pure EFT function Gamma4.
        class( parametrized_function_1D ), allocatable :: PureEFTGamma5   !< The pure EFT function Gamma5.
        class( parametrized_function_1D ), allocatable :: PureEFTGamma6   !< The pure EFT function Gamma6.

        ! the pure EFT functions with OmegaDE as argument
        class( parametrized_function_1D ), allocatable :: PureEFTOmega_ODE    !< The pure EFT function Omega.
        class( parametrized_function_1D ), allocatable :: PureEFTGamma1_ODE   !< The pure EFT function Gamma1.
        class( parametrized_function_1D ), allocatable :: PureEFTGamma2_ODE   !< The pure EFT function Gamma2.
        class( parametrized_function_1D ), allocatable :: PureEFTGamma3_ODE   !< The pure EFT function Gamma3.
        class( parametrized_function_1D ), allocatable :: PureEFTGamma4_ODE   !< The pure EFT function Gamma4.
        class( parametrized_function_1D ), allocatable :: PureEFTGamma5_ODE   !< The pure EFT function Gamma5.
        class( parametrized_function_1D ), allocatable :: PureEFTGamma6_ODE   !< The pure EFT function Gamma6.

        ! background solver parameters:
        integer  :: designer_num_points = 6000                            !< Number of points sampled by the designer code.
        real(dl) :: x_initial           = log(10._dl**(-8))               !< log(a start)
        real(dl) :: x_final             = 0.0_dl                          !< log(a final)

    contains

        ! initialization of the model:
        procedure :: read_model_selection            => EFTCAMBPureEFTstdReadModelSelectionFromFile  !< subroutine that reads the parameters of the model from file
        procedure :: allocate_model_selection        => EFTCAMBPureEFTstdAllocateModelSelection      !< subroutine that allocates the model selection.
        procedure :: init_model_parameters           => EFTCAMBPureEFTstdInitModelParameters         !< subroutine that initializes the model parameters based on the values found in an input array.
        procedure :: init_model_parameters_from_file => EFTCAMBPureEFTstdInitModelParametersFromFile !< subroutine that reads the parameters of the model from file.

        ! background solver:
        procedure :: initialize_background           => EFTCAMBPureEFTstdInitBackground               !< subroutine that initializes the background.
        procedure :: solve_designer_equations        => EFTCAMBPureEFTstdSolveDesignerEquations       !< subroutine that solves the for Mp.

        ! utility functions:
        procedure :: compute_param_number  => EFTCAMBPureEFTstdComputeParametersNumber    !< subroutine that computes the number of parameters of the model.
        procedure :: feedback              => EFTCAMBPureEFTstdFeedback                   !< subroutine that prints on the screen feedback information about the model.
        procedure :: parameter_names       => EFTCAMBPureEFTstdParameterNames             !< subroutine that returns the i-th parameter name of the model.
        procedure :: parameter_names_latex => EFTCAMBPureEFTstdParameterNamesLatex        !< subroutine that returns the i-th parameter name of the model.
        procedure :: parameter_values      => EFTCAMBPureEFTstdParameterValues            !< subroutine that returns the i-th parameter value.
        ! CAMB related procedures:
        procedure :: compute_background_EFT_functions  => EFTCAMBPureEFTstdBackgroundEFTFunctions   !< subroutine that computes the value of the background EFT functions at a given time.
        procedure :: compute_secondorder_EFT_functions => EFTCAMBPureEFTstdSecondOrderEFTFunctions  !< subroutine that computes the value of the second order EFT functions at a given time.
        procedure :: compute_dtauda                    => EFTCAMBPureEFTstdComputeDtauda            !< function that computes dtauda = 1/sqrt(a^2H^2).
        procedure :: compute_adotoa                    => EFTCAMBPureEFTstdComputeAdotoa            !< subroutine that computes adotoa = H and its two derivatives wrt conformal time.
        procedure :: compute_H_derivs                  => EFTCAMBPureEFTstdComputeHubbleDer         !< subroutine that computes the three derivatives wrt conformal time of H.
        ! stability procedures:
        procedure :: additional_model_stability        => EFTCAMBPureEFTstdAdditionalModelStability !< function that computes model specific stability requirements.

    end type EFTCAMB_std_pure_EFT

    ! ---------------------------------------------------------------------------------------------

contains

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that reads the parameters of the model from file.
    subroutine EFTCAMBPureEFTstdReadModelSelectionFromFile( self, Ini, eft_error )

        implicit none

        class(EFTCAMB_std_pure_EFT) :: self      !< the base class
        type(TIniFile)              :: Ini       !< Input ini file
        integer                     :: eft_error !< error code: 0 all fine, 1 initialization failed

        ! read model selection flags:
        self%PureEFTmodelOmega  = Ini%Read_Int( 'PureEFTmodelOmega'  , 0 )
        self%EFTwDE             = Ini%Read_Int( 'EFTwDE'             , 0 )
        self%PureEFTmodelGamma1 = Ini%Read_Int( 'PureEFTmodelGamma1' , 0 )
        self%PureEFTmodelGamma2 = Ini%Read_Int( 'PureEFTmodelGamma2' , 0 )
        self%PureEFTmodelGamma3 = Ini%Read_Int( 'PureEFTmodelGamma3' , 0 )
        self%PureEFTmodelGamma4 = Ini%Read_Int( 'PureEFTmodelGamma4' , 0 )
        self%PureEFTmodelGamma5 = Ini%Read_Int( 'PureEFTmodelGamma5' , 0 )
        self%PureEFTmodelGamma6 = Ini%Read_Int( 'PureEFTmodelGamma6' , 0 )
        self%PureEFTmodelOmega_ODE  = Ini%Read_Int( 'PureEFTmodelOmega_ODE'  , 0 )
        self%PureEFTmodelGamma1_ODE = Ini%Read_Int( 'PureEFTmodelGamma1_ODE' , 0 )
        self%PureEFTmodelGamma2_ODE = Ini%Read_Int( 'PureEFTmodelGamma2_ODE' , 0 )
        self%PureEFTmodelGamma3_ODE = Ini%Read_Int( 'PureEFTmodelGamma3_ODE' , 0 )
        self%PureEFTmodelGamma4_ODE = Ini%Read_Int( 'PureEFTmodelGamma4_ODE' , 0 )
        self%PureEFTmodelGamma5_ODE = Ini%Read_Int( 'PureEFTmodelGamma5_ODE' , 0 )
        self%PureEFTmodelGamma6_ODE = Ini%Read_Int( 'PureEFTmodelGamma6_ODE' , 0 )
        ! read the Horndeski flag:
        self%PureEFTHorndeski   = Ini%Read_Logical( 'PureEFTHorndeski' , .false. )
        ! read precision parameters
        self%designer_num_points = Ini%Read_Int( 'model_background_num_points', 6000 )
        self%x_initial = Log( Ini%Read_Double( 'model_background_a_ini', 1d-8 ) )
        self%x_final = Log( Ini%Read_Double( 'model_background_a_final', 1._dl ) )

    end subroutine EFTCAMBPureEFTstdReadModelSelectionFromFile

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that allocates the model selection.
    subroutine EFTCAMBPureEFTstdAllocateModelSelection( self, Ini, eft_error )

        implicit none

        class(EFTCAMB_std_pure_EFT) :: self      !< the base class
        type(TIniFile)              :: Ini       !< Input ini file
        integer                     :: eft_error !< error code: 0 all fine, 1 initialization failed

        integer                     :: temp_feedback 

        ! get feedback flag:
        temp_feedback = Ini%Read_Int('feedback_level', 0)

        ! allocate Omega:
        call allocate_parametrized_1D_function( self%PureEFTOmega, self%PureEFTmodelOmega, 'EFTOmega', '\Omega', eft_error, temp_feedback )
        if ( eft_error == 1 ) return
        call allocate_parametrized_1D_function( self%PureEFTOmega_ODE, self%PureEFTmodelOmega_ODE, 'EFTOmega_ODE', '\Omega_{ODE}', eft_error, temp_feedback )
        if ( eft_error == 1 ) return
        ! allocate wDE:
        if ( allocated(self%PureEFTwDE) ) deallocate(self%PureEFTwDE)
        select case ( self%EFTwDE )
            case(0)
                allocate( wDE_LCDM_parametrization_1D::self%PureEFTwDE )
                call self%PureEFTwDE%set_name( 'EFTw', 'w' )
            case(1)
                allocate( constant_parametrization_1D::self%PureEFTwDE )
                call self%PureEFTwDE%set_name( 'EFTw', 'w' )
            case(2)
                allocate( CPL_parametrization_1D::self%PureEFTwDE )
                call self%PureEFTwDE%set_param_names( ['EFTw0', 'EFTwa'], ['w_0', 'w_a'] )
                call self%PureEFTwDE%set_name( 'EFTw', 'w' )
            case(3)
                allocate( JBP_parametrization_1D::self%PureEFTwDE )
                call self%PureEFTwDE%set_param_names( ['EFTw0', 'EFTwa', 'EFTwn'], [ 'w_0', 'w_a', 'n  ' ] )
                call self%PureEFTwDE%set_name( 'EFTw', 'w' )
            case(4)
                allocate( turning_point_parametrization_1D::self%PureEFTwDE )
                call self%PureEFTwDE%set_param_names( ['EFTw0 ', 'EFTwa ', 'EFTwat'], ['w_0', 'w_a', 'a_t'] )
                call self%PureEFTwDE%set_name( 'EFTw', 'w' )
            case(5)
                allocate( taylor_parametrization_1D::self%PureEFTwDE )
                call self%PureEFTwDE%set_param_names( ['EFTw0', 'EFTwa', 'EFTw2', 'EFTw3'], ['w_0', 'w_a', 'w_2', 'w_3'] )
                call self%PureEFTwDE%set_name( 'EFTw', 'w' )
            case default
                call allocate_parametrized_1D_function( self%PureEFTwDE, self%EFTwDE, 'EFTw', 'w', eft_error, temp_feedback )
                if ( eft_error == 1 ) return
        end select
        ! allocate Gamma1:
        call allocate_parametrized_1D_function( self%PureEFTGamma1, self%PureEFTmodelGamma1, 'EFTGamma1', '{\gamma_1}', eft_error, temp_feedback )
        if ( eft_error == 1 ) return
        call allocate_parametrized_1D_function( self%PureEFTGamma1_ODE, self%PureEFTmodelGamma1_ODE, 'EFTGamma1_ODE', '{\gamma_{1,ODE}}', eft_error, temp_feedback )
        if ( eft_error == 1 ) return
        ! allocate Gamma2:
        call allocate_parametrized_1D_function( self%PureEFTGamma2, self%PureEFTmodelGamma2, 'EFTGamma2', '{\gamma_2}', eft_error, temp_feedback )
        if ( eft_error == 1 ) return
        call allocate_parametrized_1D_function( self%PureEFTGamma2_ODE, self%PureEFTmodelGamma2_ODE, 'EFTGamma2_ODE', '{\gamma_{2,ODE}}', eft_error, temp_feedback )
        if ( eft_error == 1 ) return
        ! allocate Gamma3:
        call allocate_parametrized_1D_function( self%PureEFTGamma3, self%PureEFTmodelGamma3, 'EFTGamma3', '{\gamma_3}', eft_error, temp_feedback )
        if ( eft_error == 1 ) return
        call allocate_parametrized_1D_function( self%PureEFTGamma3_ODE, self%PureEFTmodelGamma3_ODE, 'EFTGamma3_ODE', '{\gamma_{3,ODE}}', eft_error, temp_feedback )
        if ( eft_error == 1 ) return

        ! allocate the other functions only if not Horndeski:
        if ( .not. self%PureEFTHorndeski ) then
            ! allocate Gamma4:
            call allocate_parametrized_1D_function( self%PureEFTGamma4, self%PureEFTmodelGamma4, 'EFTGamma4', '{\gamma_4}', eft_error, temp_feedback )
            if ( eft_error == 1 ) return
            call allocate_parametrized_1D_function( self%PureEFTGamma4_ODE, self%PureEFTmodelGamma4_ODE, 'EFTGamma4_ODE', '{\gamma_{4,ODE}}', eft_error, temp_feedback )
            if ( eft_error == 1 ) return
            ! allocate Gamma5:
            call allocate_parametrized_1D_function( self%PureEFTGamma5, self%PureEFTmodelGamma5, 'EFTGamma5', '{\gamma_5}', eft_error, temp_feedback )
            if ( eft_error == 1 ) return
            call allocate_parametrized_1D_function( self%PureEFTGamma5_ODE, self%PureEFTmodelGamma5_ODE, 'EFTGamma5_ODE', '{\gamma_{5,ODE}}', eft_error, temp_feedback )
            if ( eft_error == 1 ) return
            ! allocate Gamma6:
            call allocate_parametrized_1D_function( self%PureEFTGamma6, self%PureEFTmodelGamma6, 'EFTGamma6', '{\gamma_6}', eft_error, temp_feedback )
            if ( eft_error == 1 ) return
            call allocate_parametrized_1D_function( self%PureEFTGamma6_ODE, self%PureEFTmodelGamma6_ODE, 'EFTGamma6_ODE', '{\gamma_{6,ODE}}', eft_error, temp_feedback )
            if ( eft_error == 1 ) return
        end if

        ! additional initialization of the function:
        call self%PureEFTOmega%init_func_from_file ( Ini, eft_error )
        call self%PureEFTwDE%init_func_from_file   ( Ini, eft_error )
        call self%PureEFTGamma1%init_func_from_file( Ini, eft_error )
        call self%PureEFTGamma2%init_func_from_file( Ini, eft_error )
        call self%PureEFTGamma3%init_func_from_file( Ini, eft_error )

        call self%PureEFTOmega_ODE%init_func_from_file ( Ini, eft_error )
        call self%PureEFTGamma1_ODE%init_func_from_file( Ini, eft_error )
        call self%PureEFTGamma2_ODE%init_func_from_file( Ini, eft_error )
        call self%PureEFTGamma3_ODE%init_func_from_file( Ini, eft_error )

        if ( .not. self%PureEFTHorndeski ) then
            call self%PureEFTGamma4%init_func_from_file( Ini, eft_error )
            call self%PureEFTGamma5%init_func_from_file( Ini, eft_error )
            call self%PureEFTGamma6%init_func_from_file( Ini, eft_error )

            call self%PureEFTGamma4_ODE%init_func_from_file( Ini, eft_error )
            call self%PureEFTGamma5_ODE%init_func_from_file( Ini, eft_error )
            call self%PureEFTGamma6_ODE%init_func_from_file( Ini, eft_error )
        end if

        if (self%PureEFTmodelOmega_ODE>1 .or. self%PureEFTmodelGamma1_ODE>1 .or. self%PureEFTmodelGamma2_ODE>1 .or. self%PureEFTmodelGamma3_ODE>1 .or. self%PureEFTmodelGamma4_ODE>1 .or. self%PureEFTmodelGamma5_ODE>1 .or. self%PureEFTmodelGamma6_ODE>1) then
            self%PureEFTOmegaDE = .true.
        else
            self%PureEFTOmegaDE = .false.
        end if

    end subroutine EFTCAMBPureEFTstdAllocateModelSelection

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that initializes the model parameters based on the values found in an input array.
    subroutine EFTCAMBPureEFTstdInitModelParameters( self, array )

        implicit none

        class(EFTCAMB_std_pure_EFT)                            :: self   !< the base class
        real(dl), dimension(self%parameter_number), intent(in) :: array  !< input array with the values of the parameters of the model.

        real(dl), allocatable, dimension(:)                    :: temp
        integer :: num_params_function, num_params_temp, i
        
        num_params_temp     = 1

        ! first elements are Omega parameters:
        num_params_function = self%PureEFTOmega%parameter_number
        allocate( temp(num_params_function) )
        do i = 1, num_params_function
            temp(i)         = array(num_params_temp)
            num_params_temp = num_params_temp +1
        end do
        call self%PureEFTOmega%init_parameters(temp)
        deallocate( temp )
        ! then w_DE parameters:
        num_params_function = self%PureEFTwDE%parameter_number
        allocate( temp(num_params_function) )
        do i = 1, num_params_function
            temp(i)         = array(num_params_temp)
            num_params_temp = num_params_temp +1
        end do
        call self%PureEFTwDE%init_parameters(temp)
        deallocate(temp)
        ! then gamma1:
        num_params_function = self%PureEFTGamma1%parameter_number
        allocate( temp(num_params_function) )
        do i = 1, num_params_function
            temp(i)         = array(num_params_temp)
            num_params_temp = num_params_temp +1
        end do
        call self%PureEFTGamma1%init_parameters(temp)
        deallocate(temp)
        ! then gamma2:
        num_params_function = self%PureEFTGamma2%parameter_number
        allocate( temp(num_params_function) )
        do i = 1, num_params_function
            temp(i)         = array(num_params_temp)
            num_params_temp = num_params_temp +1
        end do
        call self%PureEFTGamma2%init_parameters(temp)
        deallocate(temp)
        ! then gamma3:
        num_params_function = self%PureEFTGamma3%parameter_number
        allocate( temp(num_params_function) )
        do i = 1, num_params_function
            temp(i)         = array(num_params_temp)
            num_params_temp = num_params_temp +1
        end do
        call self%PureEFTGamma3%init_parameters(temp)
        deallocate(temp)

        ! The Omega_DE parameters
        ! first elements are Omega parameters:
        num_params_function = self%PureEFTOmega_ODE%parameter_number
        allocate( temp(num_params_function) )
        do i = 1, num_params_function
            temp(i)         = array(num_params_temp)
            num_params_temp = num_params_temp +1
        end do
        call self%PureEFTOmega_ODE%init_parameters(temp)
        deallocate( temp )
        ! then gamma1:
        num_params_function = self%PureEFTGamma1_ODE%parameter_number
        allocate( temp(num_params_function) )
        do i = 1, num_params_function
            temp(i)         = array(num_params_temp)
            num_params_temp = num_params_temp +1
        end do
        call self%PureEFTGamma1_ODE%init_parameters(temp)
        deallocate(temp)
        ! then gamma2:
        num_params_function = self%PureEFTGamma2_ODE%parameter_number
        allocate( temp(num_params_function) )
        do i = 1, num_params_function
            temp(i)         = array(num_params_temp)
            num_params_temp = num_params_temp +1
        end do
        call self%PureEFTGamma2_ODE%init_parameters(temp)
        deallocate(temp)
        ! then gamma3:
        num_params_function = self%PureEFTGamma3_ODE%parameter_number
        allocate( temp(num_params_function) )
        do i = 1, num_params_function
            temp(i)         = array(num_params_temp)
            num_params_temp = num_params_temp +1
        end do
        call self%PureEFTGamma3_ODE%init_parameters(temp)
        deallocate(temp)

        ! then beyond Horndeski parameters:
        if ( .not. self%PureEFTHorndeski ) then

            ! gamma4:
            num_params_function = self%PureEFTGamma4%parameter_number
            allocate( temp(num_params_function) )
            do i = 1, num_params_function
                temp(i)         = array(num_params_temp)
                num_params_temp = num_params_temp +1
            end do
            call self%PureEFTGamma4%init_parameters(temp)
            deallocate(temp)
            ! gamma5:
            num_params_function = self%PureEFTGamma5%parameter_number
            allocate( temp(num_params_function) )
            do i = 1, num_params_function
                temp(i)         = array(num_params_temp)
                num_params_temp = num_params_temp +1
            end do
            call self%PureEFTGamma5%init_parameters(temp)
            deallocate(temp)
            ! gamma6:
            num_params_function = self%PureEFTGamma6%parameter_number
            allocate( temp(num_params_function) )
            do i = 1, num_params_function
                temp(i)         = array(num_params_temp)
                num_params_temp = num_params_temp +1
            end do
            call self%PureEFTGamma6%init_parameters(temp)
            deallocate(temp)

            ! gamma4:
            num_params_function = self%PureEFTGamma4_ODE%parameter_number
            allocate( temp(num_params_function) )
            do i = 1, num_params_function
                temp(i)         = array(num_params_temp)
                num_params_temp = num_params_temp +1
            end do
            call self%PureEFTGamma4_ODE%init_parameters(temp)
            deallocate(temp)
            ! gamma5:
            num_params_function = self%PureEFTGamma5_ODE%parameter_number
            allocate( temp(num_params_function) )
            do i = 1, num_params_function
                temp(i)         = array(num_params_temp)
                num_params_temp = num_params_temp +1
            end do
            call self%PureEFTGamma5_ODE%init_parameters(temp)
            deallocate(temp)
            ! gamma6:
            num_params_function = self%PureEFTGamma6_ODE%parameter_number
            allocate( temp(num_params_function) )
            do i = 1, num_params_function
                temp(i)         = array(num_params_temp)
                num_params_temp = num_params_temp +1
            end do
            call self%PureEFTGamma6_ODE%init_parameters(temp)
            deallocate(temp)

        end if

        ! now check the length of the parameters:
        if ( num_params_temp-1 /= self%parameter_number ) then
            write(*,*) 'In EFTCAMBPureEFTstdInitModelParameters:'
            write(*,*) 'Length of num_params_temp and self%parameter_number do not coincide.'
            write(*,*) 'num_params_temp:', num_params_temp-1
            write(*,*) 'self%parameter_number:', self%parameter_number
            call MpiStop('EFTCAMB error')
        end if

    end subroutine EFTCAMBPureEFTstdInitModelParameters

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that reads the parameters of the model from file.
    subroutine EFTCAMBPureEFTstdInitModelParametersFromFile( self, Ini, eft_error )

        implicit none

        class(EFTCAMB_std_pure_EFT) :: self      !< the base class
        type(TIniFile)              :: Ini       !< Input ini file
        integer                     :: eft_error !< error code: 0 all fine, 1 initialization failed
        
        call self%PureEFTOmega%init_from_file ( Ini, eft_error )
        call self%PureEFTwDE%init_from_file   ( Ini, eft_error )
        call self%PureEFTGamma1%init_from_file( Ini, eft_error )
        call self%PureEFTGamma2%init_from_file( Ini, eft_error )
        call self%PureEFTGamma3%init_from_file( Ini, eft_error )

        call self%PureEFTOmega_ODE%init_from_file ( Ini, eft_error )
        call self%PureEFTGamma1_ODE%init_from_file( Ini, eft_error )
        call self%PureEFTGamma2_ODE%init_from_file( Ini, eft_error )
        call self%PureEFTGamma3_ODE%init_from_file( Ini, eft_error )

        if ( .not. self%PureEFTHorndeski ) then
            call self%PureEFTGamma4%init_from_file( Ini, eft_error )
            call self%PureEFTGamma5%init_from_file( Ini, eft_error )
            call self%PureEFTGamma6%init_from_file( Ini, eft_error )

            call self%PureEFTGamma4_ODE%init_from_file( Ini, eft_error )
            call self%PureEFTGamma5_ODE%init_from_file( Ini, eft_error )
            call self%PureEFTGamma6_ODE%init_from_file( Ini, eft_error )
        end if

    end subroutine EFTCAMBPureEFTstdInitModelParametersFromFile

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes the number of parameters of the model.
    subroutine EFTCAMBPureEFTstdComputeParametersNumber( self )

        implicit none

        class(EFTCAMB_std_pure_EFT)  :: self   !< the base class

        self%parameter_number = 0
        self%parameter_number = self%parameter_number +self%PureEFTOmega%parameter_number
        self%parameter_number = self%parameter_number +self%PureEFTwDE%parameter_number
        self%parameter_number = self%parameter_number +self%PureEFTGamma1%parameter_number
        self%parameter_number = self%parameter_number +self%PureEFTGamma2%parameter_number
        self%parameter_number = self%parameter_number +self%PureEFTGamma3%parameter_number

        self%parameter_number = self%parameter_number +self%PureEFTOmega_ODE%parameter_number
        self%parameter_number = self%parameter_number +self%PureEFTGamma1_ODE%parameter_number
        self%parameter_number = self%parameter_number +self%PureEFTGamma2_ODE%parameter_number
        self%parameter_number = self%parameter_number +self%PureEFTGamma3_ODE%parameter_number
        if ( .not. self%PureEFTHorndeski ) then
            self%parameter_number = self%parameter_number +self%PureEFTGamma4%parameter_number
            self%parameter_number = self%parameter_number +self%PureEFTGamma5%parameter_number
            self%parameter_number = self%parameter_number +self%PureEFTGamma6%parameter_number

            self%parameter_number = self%parameter_number +self%PureEFTGamma4_ODE%parameter_number
            self%parameter_number = self%parameter_number +self%PureEFTGamma5_ODE%parameter_number
            self%parameter_number = self%parameter_number +self%PureEFTGamma6_ODE%parameter_number
        end if

    end subroutine EFTCAMBPureEFTstdComputeParametersNumber

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that prints on the screen feedback information about the model.
    subroutine EFTCAMBPureEFTstdFeedback( self, print_params )

        implicit none

        class(EFTCAMB_std_pure_EFT)  :: self         !< the base class
        logical, optional            :: print_params !< optional flag that decised whether to print numerical values
                                                     !! of the parameters.

        ! print general model informations:
        write(*,*)
        write(*,'(a,a)')    '   Model               =  ', self%name
        if ( self%PureEFTHorndeski ) then
            write(*,"(a)")  '   Pure EFT Horndeski'
        end if
        write(*,'(a,I3)')   '   Number of params    ='  , self%parameter_number

        ! print model functions informations:
        write(*,*)
        if ( self%PureEFTmodelOmega  /= 0 ) write(*,'(a,I3)') '   PureEFTmodelOmega   =', self%PureEFTmodelOmega
        if ( self%EFTwDE             /= 0 ) write(*,'(a,I3)') '   EFTwDE              =', self%EFTwDE
        if ( self%PureEFTmodelGamma1 /= 0 ) write(*,'(a,I3)') '   PureEFTmodelGamma1  =', self%PureEFTmodelGamma1
        if ( self%PureEFTmodelGamma2 /= 0 ) write(*,'(a,I3)') '   PureEFTmodelGamma2  =', self%PureEFTmodelGamma2
        if ( self%PureEFTmodelGamma3 /= 0 ) write(*,'(a,I3)') '   PureEFTmodelGamma3  =', self%PureEFTmodelGamma3

        if ( self%PureEFTmodelOmega_ODE  /= 0 ) write(*,'(a,I3)') '   PureEFTmodelOmega_ODE   =', self%PureEFTmodelOmega
        if ( self%PureEFTmodelGamma1_ODE /= 0 ) write(*,'(a,I3)') '   PureEFTmodelGamma1_ODE  =', self%PureEFTmodelGamma1_ODE
        if ( self%PureEFTmodelGamma2_ODE /= 0 ) write(*,'(a,I3)') '   PureEFTmodelGamma2_ODE  =', self%PureEFTmodelGamma2_ODE
        if ( self%PureEFTmodelGamma3_ODE /= 0 ) write(*,'(a,I3)') '   PureEFTmodelGamma3_ODE  =', self%PureEFTmodelGamma3_ODE
        if ( .not. self%PureEFTHorndeski ) then
            if ( self%PureEFTmodelGamma4 /= 0 ) write(*,'(a,I3)') '   PureEFTmodelGamma4  =', self%PureEFTmodelGamma4
            if ( self%PureEFTmodelGamma5 /= 0 ) write(*,'(a,I3)') '   PureEFTmodelGamma5  =', self%PureEFTmodelGamma5
            if ( self%PureEFTmodelGamma6 /= 0 ) write(*,'(a,I3)') '   PureEFTmodelGamma6  =', self%PureEFTmodelGamma6

            if ( self%PureEFTmodelGamma4_ODE /= 0 ) write(*,'(a,I3)') '   PureEFTmodelGamma4_ODE  =', self%PureEFTmodelGamma4
            if ( self%PureEFTmodelGamma5_ODE /= 0 ) write(*,'(a,I3)') '   PureEFTmodelGamma5_ODE  =', self%PureEFTmodelGamma5
            if ( self%PureEFTmodelGamma6_ODE /= 0 ) write(*,'(a,I3)') '   PureEFTmodelGamma6_ODE  =', self%PureEFTmodelGamma6
        end if

        write(*,*)
        ! print functions informations:
        call self%PureEFTOmega%feedback  ( print_params )
        call self%PureEFTwDE%feedback    ( print_params )
        call self%PureEFTGamma1%feedback ( print_params )
        call self%PureEFTGamma2%feedback ( print_params )
        call self%PureEFTGamma3%feedback ( print_params )

        call self%PureEFTOmega_ODE%feedback  ( print_params )
        call self%PureEFTGamma1_ODE%feedback ( print_params )
        call self%PureEFTGamma2_ODE%feedback ( print_params )
        call self%PureEFTGamma3_ODE%feedback ( print_params )
        if ( .not. self%PureEFTHorndeski ) then
            call self%PureEFTGamma4%feedback( print_params )
            call self%PureEFTGamma5%feedback( print_params )
            call self%PureEFTGamma6%feedback( print_params )

            call self%PureEFTGamma4_ODE%feedback( print_params )
            call self%PureEFTGamma5_ODE%feedback( print_params )
            call self%PureEFTGamma6_ODE%feedback( print_params )
        end if

    end subroutine EFTCAMBPureEFTstdFeedback

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the i-th parameter name of the model
    subroutine EFTCAMBPureEFTstdParameterNames( self, i, name )

        implicit none

        class(EFTCAMB_std_pure_EFT) :: self   !< the base class
        integer     , intent(in)    :: i      !< the index of the parameter
        character(*), intent(out)   :: name   !< the output name of the i-th parameter

        integer  :: Nc, NOmega, Nw, N1, N2, N3, N4, N5, N6, NOmegao, N1o, N2o, N3o, N4o, N5o, N6o
        integer  :: j

        ! compute the incremental number of parameters:
        Nc     = 0
        NOmega = Nc + self%PureEFTOmega%parameter_number
        Nw     = NOmega + self%PureEFTwDE%parameter_number
        N1     = Nw + self%PureEFTGamma1%parameter_number
        N2     = N1 + self%PureEFTGamma2%parameter_number
        N3     = N2 + self%PureEFTGamma3%parameter_number
        
        NOmegao = N3 + self%PureEFTOmega_ODE%parameter_number
        N1o     = NOmegao + self%PureEFTGamma1_ODE%parameter_number
        N2o     = N1o + self%PureEFTGamma2_ODE%parameter_number
        N3o     = N2o + self%PureEFTGamma3_ODE%parameter_number
        if ( .not. self%PureEFTHorndeski ) then
            N4 = N3o + self%PureEFTGamma4%parameter_number
            N5 = N4 + self%PureEFTGamma5%parameter_number
            N6 = N5 + self%PureEFTGamma6%parameter_number

            N4o = N6 + self%PureEFTGamma4_ODE%parameter_number
            N5o = N4o + self%PureEFTGamma5_ODE%parameter_number
            N6o = N5o + self%PureEFTGamma6_ODE%parameter_number
        end if

        ! check validity of input:
        if ( i > self%parameter_number .or. i <= 0 ) then
            write(*,'(a,I3)') 'No parameter corresponding to: ', i
            write(*,'(a,I3)') 'Total number of parameters is: ', self%parameter_number
            call MpiStop('EFTCAMB error')
        
        ! parameter from Omega function
        else if ( i <= NOmega ) then
            do j = 1, self%PureEFTOmega%parameter_number
                if ( i-Nc == j ) call self%PureEFTOmega%parameter_names( j, name )
            end do
            return

        ! parameter from wDE function
        else if ( i <= Nw ) then
            do j = 1, self%PureEFTwDE%parameter_number
                if ( i-NOmega == j ) call self%PureEFTwDE%parameter_names( j, name )
            end do
            return

        ! parameter from Gamma1 function
        else if ( i <= N1) then
            do j = 1, self%PureEFTGamma1%parameter_number
                if ( i-Nw == j ) call self%PureEFTGamma1%parameter_names( j, name )
            end do
            return

        !parameter from Gamma2 function
        else if ( i <= N2 ) then
            do j = 1, self%PureEFTGamma2%parameter_number
                if ( i-N1 == j ) call self%PureEFTGamma2%parameter_names( j, name )
            end do
            return

        !parameter from Gamma3 function
        else if ( i <= N3 ) then

            do j = 1, self%PureEFTGamma3%parameter_number
                if ( i-N2 == j ) call self%PureEFTGamma3%parameter_names( j, name )
            end do
            return

        ! parameter from Omega_ODE function
        else if ( i <= NOmegao ) then
            do j = 1, self%PureEFTOmega_ODE%parameter_number
                if ( i-N3 == j ) call self%PureEFTOmega_ODE%parameter_names( j, name )
            end do
            return

        ! parameter from Gamma1_ODE function
        else if ( i <= N1o) then
            do j = 1, self%PureEFTGamma1_ODE%parameter_number
                if ( i-NOmegao == j ) call self%PureEFTGamma1_ODE%parameter_names( j, name )
            end do
            return

        !parameter from Gamma2_ODE function
        else if ( i <= N2o ) then
            do j = 1, self%PureEFTGamma2_ODE%parameter_number
                if ( i-N1o == j ) call self%PureEFTGamma2_ODE%parameter_names( j, name )
            end do
            return

        !parameter from Gamma3_ODE function
        else if ( i <= N3o ) then

            do j = 1, self%PureEFTGamma3_ODE%parameter_number
                if ( i-N2o == j ) call self%PureEFTGamma3_ODE%parameter_names( j, name )
            end do
            return

        !parameter from Gamma4 function
        else if ( .not. self%PureEFTHorndeski .and. i <= N4 ) then

            do j = 1, self%PureEFTGamma4%parameter_number
                if ( i-N3o == j ) call self%PureEFTGamma4%parameter_names( j, name )
            end do
            return

        !parameter from Gamma5 function
        else if ( .not. self%PureEFTHorndeski .and. i <= N5 ) then

            do j = 1, self%PureEFTGamma5%parameter_number
                if ( i-N4 == j ) call self%PureEFTGamma5%parameter_names( j, name )
            end do
            return

        !parameter from Gamma6 function
        else if ( .not. self%PureEFTHorndeski .and. i <= N6 ) then

            do j = 1, self%PureEFTGamma6%parameter_number
                if ( i-N5 == j ) call self%PureEFTGamma6%parameter_names( j, name )
            end do
            return
        
        !parameter from Gamma4_ODE function
        else if ( .not. self%PureEFTHorndeski .and. i <= N4o ) then

            do j = 1, self%PureEFTGamma4_ODE%parameter_number
                if ( i-N6 == j ) call self%PureEFTGamma4_ODE%parameter_names( j, name )
            end do
            return

        !parameter from Gamma5_ODE function
        else if ( .not. self%PureEFTHorndeski .and. i <= N5o ) then

            do j = 1, self%PureEFTGamma5_ODE%parameter_number
                if ( i-N4o == j ) call self%PureEFTGamma5_ODE%parameter_names( j, name )
            end do
            return

        !parameter from Gamma6_ODE function
        else if ( .not. self%PureEFTHorndeski .and. i <= N6o ) then

            do j = 1, self%PureEFTGamma6_ODE%parameter_number
                if ( i-N5o == j ) call self%PureEFTGamma6_ODE%parameter_names( j, name )
            end do
            return

        end if

    end subroutine EFTCAMBPureEFTstdParameterNames

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the i-th parameter name of the model
    subroutine EFTCAMBPureEFTstdParameterNamesLatex( self, i, latexname )

        implicit none

        class(EFTCAMB_std_pure_EFT) :: self       !< the base class
        integer     , intent(in)    :: i          !< The index of the parameter
        character(*), intent(out)   :: latexname  !< the output latex name of the i-th parameter

        integer  :: Nc, NOmega, Nw, N1, N2, N3, N4, N5, N6, NOmegao, N1o, N2o, N3o, N4o, N5o, N6o
        integer  :: j

        ! compute the incremental number of parameters:
        Nc     = 0
        NOmega = Nc + self%PureEFTOmega%parameter_number
        Nw     = NOmega + self%PureEFTwDE%parameter_number
        N1     = Nw + self%PureEFTGamma1%parameter_number
        N2     = N1 + self%PureEFTGamma2%parameter_number
        N3     = N2 + self%PureEFTGamma3%parameter_number
        
        NOmegao = N3 + self%PureEFTOmega_ODE%parameter_number
        N1o     = NOmegao + self%PureEFTGamma1_ODE%parameter_number
        N2o     = N1o + self%PureEFTGamma2_ODE%parameter_number
        N3o     = N2o + self%PureEFTGamma3_ODE%parameter_number
        if ( .not. self%PureEFTHorndeski ) then
            N4 = N3o + self%PureEFTGamma4%parameter_number
            N5 = N4 + self%PureEFTGamma5%parameter_number
            N6 = N5 + self%PureEFTGamma6%parameter_number

            N4o = N6 + self%PureEFTGamma4%parameter_number
            N5o = N4o + self%PureEFTGamma5%parameter_number
            N6o = N5o + self%PureEFTGamma6%parameter_number
        end if

        ! check validity of input:
        if ( i > self%parameter_number .or. i <= 0 ) then
            write(*,'(a,I3)') 'No parameter corresponding to: ', i
            write(*,'(a,I3)') 'Total number of parameters is: ', self%parameter_number
            call MpiStop('EFTCAMB error')
        
        ! parameter from Omega function
        else if ( i <= NOmega ) then
            do j = 1, self%PureEFTOmega%parameter_number
                if ( i-Nc == j ) call self%PureEFTOmega%parameter_names_latex( j, latexname )
            end do
            return

        ! parameter from wDE function
        else if ( i <= Nw ) then
            do j = 1, self%PureEFTwDE%parameter_number
                if ( i-NOmega == j ) call self%PureEFTwDE%parameter_names_latex( j, latexname )
            end do
            return

        ! parameter from Gamma1 function
        else if ( i <= N1) then
            do j = 1, self%PureEFTGamma1%parameter_number
                if ( i-Nw == j ) call self%PureEFTGamma1%parameter_names_latex( j, latexname )
            end do
            return

        !parameter from Gamma2 function
        else if ( i <= N2 ) then
            do j = 1, self%PureEFTGamma2%parameter_number
                if ( i-N1 == j ) call self%PureEFTGamma2%parameter_names_latex( j, latexname )
            end do
            return

        !parameter from Gamma3 function
        else if ( i <= N3 ) then

            do j = 1, self%PureEFTGamma3%parameter_number
                if ( i-N2 == j ) call self%PureEFTGamma3%parameter_names_latex( j, latexname )
            end do
            return

        ! parameter from Omega_ODE function
        else if ( i <= NOmegao ) then
            do j = 1, self%PureEFTOmega_ODE%parameter_number
                if ( i-N3 == j ) call self%PureEFTOmega_ODE%parameter_names_latex( j, latexname )
            end do
            return

        ! parameter from Gamma1_ODE function
        else if ( i <= N1o) then
            do j = 1, self%PureEFTGamma1_ODE%parameter_number
                if ( i-NOmegao == j ) call self%PureEFTGamma1_ODE%parameter_names_latex( j, latexname )
            end do
            return

        !parameter from Gamma2_ODE function
        else if ( i <= N2o ) then
            do j = 1, self%PureEFTGamma2_ODE%parameter_number
                if ( i-N1o == j ) call self%PureEFTGamma2_ODE%parameter_names_latex( j, latexname )
            end do
            return

        !parameter from Gamma3_ODE function
        else if ( i <= N3o ) then

            do j = 1, self%PureEFTGamma3_ODE%parameter_number
                if ( i-N2o == j ) call self%PureEFTGamma3_ODE%parameter_names_latex( j, latexname )
            end do
            return

        !parameter from Gamma4 function
        else if ( .not. self%PureEFTHorndeski .and. i <= N4 ) then

            do j = 1, self%PureEFTGamma4%parameter_number
                if ( i-N3o == j ) call self%PureEFTGamma4%parameter_names_latex( j, latexname )
            end do
            return

        !parameter from Gamma5 function
        else if ( .not. self%PureEFTHorndeski .and. i <= N5 ) then

            do j = 1, self%PureEFTGamma5%parameter_number
                if ( i-N4 == j ) call self%PureEFTGamma5%parameter_names_latex( j, latexname )
            end do
            return

        !parameter from Gamma6 function
        else if ( .not. self%PureEFTHorndeski .and. i <= N6 ) then

            do j = 1, self%PureEFTGamma6%parameter_number
                if ( i-N5 == j ) call self%PureEFTGamma6%parameter_names_latex( j, latexname )
            end do
            return
        
        !parameter from Gamma4_ODE function
        else if ( .not. self%PureEFTHorndeski .and. i <= N4o ) then

            do j = 1, self%PureEFTGamma4_ODE%parameter_number
                if ( i-N6 == j ) call self%PureEFTGamma4_ODE%parameter_names_latex( j, latexname )
            end do
            return

        !parameter from Gamma5_ODE function
        else if ( .not. self%PureEFTHorndeski .and. i <= N5o ) then

            do j = 1, self%PureEFTGamma5_ODE%parameter_number
                if ( i-N4o == j ) call self%PureEFTGamma5_ODE%parameter_names_latex( j, latexname )
            end do
            return

        !parameter from Gamma6_ODE function
        else if ( .not. self%PureEFTHorndeski .and. i <= N6o ) then

            do j = 1, self%PureEFTGamma6_ODE%parameter_number
                if ( i-N5o == j ) call self%PureEFTGamma6_ODE%parameter_names_latex( j, latexname )
            end do
            return

        end if

    end subroutine EFTCAMBPureEFTstdParameterNamesLatex

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the i-th parameter name of the model
    subroutine EFTCAMBPureEFTstdParameterValues( self, i, value )

        implicit none

        class(EFTCAMB_std_pure_EFT) :: self   !< the base class
        integer , intent(in)        :: i      !< The index of the parameter
        real(dl), intent(out)       :: value  !< the output value of the i-th parameter

        integer  :: Nc, NOmega, Nw, N1, N2, N3, N4, N5, N6, NOmegao, N1o, N2o, N3o, N4o, N5o, N6o
        integer  :: j

        ! compute the incremental number of parameters:
        Nc     = 0
        NOmega = Nc + self%PureEFTOmega%parameter_number
        Nw     = NOmega + self%PureEFTwDE%parameter_number
        N1     = Nw + self%PureEFTGamma1%parameter_number
        N2     = N1 + self%PureEFTGamma2%parameter_number
        N3     = N2 + self%PureEFTGamma3%parameter_number
        
        NOmegao = N3 + self%PureEFTOmega_ODE%parameter_number
        N1o     = NOmegao + self%PureEFTGamma1_ODE%parameter_number
        N2o     = N1o + self%PureEFTGamma2_ODE%parameter_number
        N3o     = N2o + self%PureEFTGamma3_ODE%parameter_number
        if ( .not. self%PureEFTHorndeski ) then
            N4 = N3o + self%PureEFTGamma4%parameter_number
            N5 = N4 + self%PureEFTGamma5%parameter_number
            N6 = N5 + self%PureEFTGamma6%parameter_number

            N4o = N6 + self%PureEFTGamma4%parameter_number
            N5o = N4o + self%PureEFTGamma5%parameter_number
            N6o = N5o + self%PureEFTGamma6%parameter_number
        end if

        ! check validity of input:
        if ( i > self%parameter_number .or. i <= 0 ) then
            write(*,'(a,I3)') 'No parameter corresponding to: ', i
            write(*,'(a,I3)') 'Total number of parameters is: ', self%parameter_number
            call MpiStop('EFTCAMB error')
        
        ! parameter from Omega function
        else if ( i <= NOmega ) then
            do j = 1, self%PureEFTOmega%parameter_number
                if ( i-Nc == j ) call self%PureEFTOmega%parameter_value( j, value )
            end do
            return

        ! parameter from wDE function
        else if ( i <= Nw ) then
            do j = 1, self%PureEFTwDE%parameter_number
                if ( i-NOmega == j ) call self%PureEFTwDE%parameter_value( j, value )
            end do
            return

        ! parameter from Gamma1 function
        else if ( i <= N1) then
            do j = 1, self%PureEFTGamma1%parameter_number
                if ( i-Nw == j ) call self%PureEFTGamma1%parameter_value( j, value )
            end do
            return

        !parameter from Gamma2 function
        else if ( i <= N2 ) then
            do j = 1, self%PureEFTGamma2%parameter_number
                if ( i-N1 == j ) call self%PureEFTGamma2%parameter_value( j, value )
            end do
            return

        !parameter from Gamma3 function
        else if ( i <= N3 ) then

            do j = 1, self%PureEFTGamma3%parameter_number
                if ( i-N2 == j ) call self%PureEFTGamma3%parameter_value( j, value )
            end do
            return

        ! parameter from Omega_ODE function
        else if ( i <= NOmegao ) then
            do j = 1, self%PureEFTOmega_ODE%parameter_number
                if ( i-N3 == j ) call self%PureEFTOmega_ODE%parameter_value( j, value )
            end do
            return

        ! parameter from Gamma1_ODE function
        else if ( i <= N1o) then
            do j = 1, self%PureEFTGamma1_ODE%parameter_number
                if ( i-NOmegao == j ) call self%PureEFTGamma1_ODE%parameter_value( j, value )
            end do
            return

        !parameter from Gamma2_ODE function
        else if ( i <= N2o ) then
            do j = 1, self%PureEFTGamma2_ODE%parameter_number
                if ( i-N1o == j ) call self%PureEFTGamma2_ODE%parameter_value( j, value )
            end do
            return

        !parameter from Gamma3_ODE function
        else if ( i <= N3o ) then

            do j = 1, self%PureEFTGamma3_ODE%parameter_number
                if ( i-N2o == j ) call self%PureEFTGamma3_ODE%parameter_value( j, value )
            end do
            return

        !parameter from Gamma4 function
        else if ( .not. self%PureEFTHorndeski .and. i <= N4 ) then

            do j = 1, self%PureEFTGamma4%parameter_number
                if ( i-N3o == j ) call self%PureEFTGamma4%parameter_value( j, value )
            end do
            return

        !parameter from Gamma5 function
        else if ( .not. self%PureEFTHorndeski .and. i <= N5 ) then

            do j = 1, self%PureEFTGamma5%parameter_number
                if ( i-N4 == j ) call self%PureEFTGamma5%parameter_value( j, value )
            end do
            return

        !parameter from Gamma6 function
        else if ( .not. self%PureEFTHorndeski .and. i <= N6 ) then

            do j = 1, self%PureEFTGamma6%parameter_number
                if ( i-N5 == j ) call self%PureEFTGamma6%parameter_value( j, value )
            end do
            return
        
        !parameter from Gamma4_ODE function
        else if ( .not. self%PureEFTHorndeski .and. i <= N4o ) then

            do j = 1, self%PureEFTGamma4_ODE%parameter_number
                if ( i-N6 == j ) call self%PureEFTGamma4_ODE%parameter_value( j, value )
            end do
            return

        !parameter from Gamma5_ODE function
        else if ( .not. self%PureEFTHorndeski .and. i <= N5o ) then

            do j = 1, self%PureEFTGamma5_ODE%parameter_number
                if ( i-N4o == j ) call self%PureEFTGamma5_ODE%parameter_value( j, value )
            end do
            return

        !parameter from Gamma6_ODE function
        else if ( .not. self%PureEFTHorndeski .and. i <= N6o ) then

            do j = 1, self%PureEFTGamma6_ODE%parameter_number
                if ( i-N5o == j ) call self%PureEFTGamma6_ODE%parameter_value( j, value )
            end do
            return

        end if

    end subroutine EFTCAMBPureEFTstdParameterValues

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that initializes the background.
    subroutine EFTCAMBPureEFTstdInitBackground( self, params_cache, feedback_level, success, outroot )

        implicit none

        class(EFTCAMB_std_pure_EFT)                    :: self           !< the base class
        type(TEFTCAMB_parameter_cache), intent(inout) :: params_cache   !< a EFTCAMB parameter cache containing cosmological parameters
        integer                       , intent(in)    :: feedback_level !< level of feedback from the background code. 0=none; 1=some; 2=chatty.
        logical                       , intent(out)   :: success        !< wether the background initialization succeded or not
        character(LEN=*), optional    , intent(in)    :: outroot        !< the output root for the debug files

        ! some feedback:
        if ( feedback_level>0 ) then
        write(*,'(a)') "***************************************************************"
        write(*,'(a)') ' EFTCAMB designer pure EFT background solver'
        write(*,'(a)')
        end if

        call self%rhov_t%initialize( self%designer_num_points, self%x_initial, self%x_final )

        success = .True.
        ! solve the background equations and store the solution:
        call self%solve_designer_equations( params_cache, success=success, feedback_level=feedback_level )

    end subroutine EFTCAMBPureEFTstdInitBackground

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that solves for M by integrating apha_M
    subroutine EFTCAMBPureEFTstdSolveDesignerEquations( self, params_cache, success, feedback_level )

        implicit none

        class(EFTCAMB_std_pure_EFT)                   :: self          !< the base class.
        type(TEFTCAMB_parameter_cache), intent(in)   :: params_cache  !< a EFTCAMB parameter cache containing cosmological parameters.
        logical , intent(out)                        :: success       !< whether the calculation ended correctly or not
        integer , intent(in)                         :: feedback_level   !< whether be noisy

        integer, parameter :: num_eq = 1   !<  Number of equations
        real(dl) :: y(num_eq), ydot(num_eq)

        real(dl) :: grhom0, grhor0, grhov0

        integer  :: itol, itask, istate, iopt, LRN, LRS, LRW, LIS, LIN, LIW, JacobianMode, i
        real(dl) :: rtol, atol, t1, t2
        real(dl), allocatable :: rwork(:)
        integer,  allocatable :: iwork(:)

        ! 1) Cosmological parameters:
        grhom0 = params_cache%grhob + params_cache%grhoc
        grhor0 = 3._dl * params_cache%h0_Mpc**2 * (params_cache%omegag + params_cache%omegar)
        grhov0 = params_cache%grhov

        ! 2) Initial values
        y(1) = 1._dl
        success = .True.

        ! 3) Initialize DLSODA:
        ! set-up the relative and absolute tollerances:
        itol = 1
        rtol = 1.d-12
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
        IWORK(6) = 100    ! maximum number of (internally defined) steps allowed during one call to the solver. The default value is 500.
        IWORK(7) = 0      ! maximum number of messages printed (per problem) warning that T + H = T on a step (H = step size). This must be positive to result in a non-default value.  The default value is 10.
        IWORK(8) = 0      ! the maximum order to be allowed for the nonstiff (Adams) method.  the default value is 12.
        IWORK(9) = 0      ! the maximum order to be allowed for the stiff (BDF) method.  The default value is 5.
        ! additional lsoda stuff:
        CALL XSETF(0) ! suppress odepack printing
        ! Jacobian mode: 1=fullJacobian, 2=not provided
        JacobianMode = 2

        ! initial time:
        t1  = self%rhov_t%x( self%designer_num_points )
        ! store initial step:
        call output( num_eq, self%designer_num_points, t1, y)
        ! solve the equations:
        do i=self%designer_num_points, 2, -1
            ! set the time step:
            t1 = self%rhov_t%x(i)
            t2 = self%rhov_t%x(i-1)
            ! solve the system:
            call DLSODA ( derivs, num_eq, y, t1, t2, itol, rtol, atol, itask, istate, iopt, RWORK, LRW, IWORK, LIW, jacobian, JacobianMode)
            ! check istate for LSODA good completion:
            if ( istate < 0 ) then
                if ( istate == -1 ) then
                    if ( feedback_level>1 ) write(*,*) 'DLSODA excessive work'
                    istate = 1
                else
                    success = .False.
                    write(*,*) 'DLSODA failed with code:', istate
                    return
                end if
            end if
            ! compute output functions if needed:
            call output( num_eq, i-1, t2, y )
        end do

        return

    contains
        ! ---------------------------------------------------------------------------------------------
        !> Subroutine that computes y' given y
        subroutine derivs( num_eq, x, y, ydot )

        implicit none

        integer , intent(in)                     :: num_eq !< number of equations in the ODE system
        real(dl), intent(in)                     :: x      !< time at which the derivatives of the system are computed
        real(dl), intent(in) , dimension(num_eq) :: y      !< input status of the system
        real(dl), intent(out), dimension(num_eq) :: ydot   !< value of the derivative at time x

        real(dl) :: a
        real(dl) :: grhor, grhom, grhov, grhon
        real(dl) :: grhonu, gpinu, grhormass_t
        integer  :: nu_i
        
        a = Exp(x)

        grhom = grhom0 / a
        grhor = grhor0 / a**2
        grhov = grhov0 * y(1) * a**2
        grhon = 0._dl
        if ( params_cache%Num_Nu_Massive /= 0 ) then
            do nu_i = 1, params_cache%Nu_mass_eigenstates
                grhonu      = 0._dl
                gpinu       = 0._dl
                grhormass_t = params_cache%grhormass(nu_i)/a**2
                call ThermalNuBack%rho_P(a*params_cache%nu_masses(nu_i), grhonu, gpinu )
                grhon = grhon +grhormass_t*grhonu
            end do
        end if

        ydot(1) = -3._dl*(1._dl + self%PureEFTwDE%value(a))*y(1)

        end subroutine derivs

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
        !> Subroutine that takes the solution of the background DGP equations and stores the values of the EFT functions.
        subroutine output( num_eq, ind, x, y)

        implicit none

        integer , intent(in)                     :: num_eq !< number of equations in the ODE system.
        integer , intent(in)                     :: ind    !< index of the EFT functions interpolation tables to fill.
        real(dl), intent(in)                     :: x      !< time at which the derivatives of the system are computed.
        real(dl), intent(in) , dimension(num_eq) :: y      !< input status of the system.

        real(dl) :: a

        self%rhov_t%y( ind ) = grhov0 * y(1)
        
        end subroutine output

    end subroutine EFTCAMBPureEFTstdSolveDesignerEquations

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes the value of the background EFT functions at a given time.
    subroutine EFTCAMBPureEFTstdBackgroundEFTFunctions( self, a, eft_par_cache, eft_cache )

        implicit none

        class(EFTCAMB_std_pure_EFT)                  :: self          !< the base class
        real(dl), intent(in)                         :: a             !< the input scale factor
        type(TEFTCAMB_parameter_cache), intent(inout) :: eft_par_cache !< the EFTCAMB parameter cache that contains all the physical parameters.
        type(TEFTCAMB_timestep_cache ), intent(inout) :: eft_cache     !< the EFTCAMB timestep cache that contains all the physical values.

        real(dl) :: H, Ht, Ht2, Ht3, Ht4
        real(dl) :: grhor, grhom, grhov, grhon, gpn, gdpn, gddpn, gdddpn
        real(dl) :: grhonu, gpinu, gpinudot, gpinudotdot, gpinudotdotdot, grhormass_t
        integer  :: nu_i
        real(dl) :: wDE, dwDE, d2wDE, d3wDE
        real(dl) :: Ode, dOde, d2Ode, d3Ode, d4Ode

        grhom = eft_cache%grhob_t + eft_cache%grhoc_t
        grhor = eft_cache%grhor_t + eft_cache%grhog_t
        grhov = eft_cache%grhov_t

        wDE   = self%PureEFTwDE%value(a)
        dwDE  = self%PureEFTwDE%first_derivative(a)
        d2wDE = self%PureEFTwDE%second_derivative(a)
        d3wDE = self%PureEFTwDE%third_derivative(a)
        
        H   = eft_cache%adotoa
        Ht  = eft_cache%Hdot
        Ht2 = eft_cache%Hdotdot
        Ht3 = eft_cache%Hdotdotdot
        if ( self%PureEFTOmegaDE ) then
            grhon = 0._dl
            gpn = 0._dl
            gdpn = 0._dl
            gddpn = 0._dl
            gdddpn = 0._dl
            if ( eft_par_cache%Num_Nu_Massive /= 0 ) then
                do nu_i = 1, eft_par_cache%Nu_mass_eigenstates
                    grhonu         = 0._dl
                    gpinu          = 0._dl
                    gpinudot       = 0._dl
                    gpinudotdot    = 0._dl
                    gpinudotdotdot = 0._dl
                    grhormass_t = eft_par_cache%grhormass(nu_i)/a**2
                    call ThermalNuBack%rho_P(a*eft_par_cache%nu_masses(nu_i), grhonu, gpinu )
                    grhon          = grhon +grhormass_t*grhonu
                    gpn            = gpn  +grhormass_t*gpinu
                    gpinudot       = ThermalNuBack%pidot( a*eft_par_cache%nu_masses(nu_i), H, gpinu )
                    gdpn           = gdpn  + grhormass_t*( gpinudot - 4._dl*H*gpinu)
                    gpinudotdot    = ThermalNuBack%pidotdot( a*eft_par_cache%nu_masses(nu_i), H, Ht, gpinu, gpinudot)
                    gddpn          = gddpn - 4._dl*H*grhormass_t*( gpinudot - 4._dl*H*gpinu )+ grhormass_t*( gpinudotdot - 4._dl*Ht*gpinu - 4._dl*H*gpinudot )
                    gpinudotdotdot = ThermalNuBack%pidotdotdot(a*eft_par_cache%nu_masses(nu_i),H,Ht,Ht2,gpinu,gpinudot,gpinudotdot)
                    gdddpn = gdddpn + grhormass_t*(gpinudotdotdot &
                        & - 12._dl*H*gpinudotdot &
                        & + (48._dl*H**2 - 12._dl*Ht)*gpinudot &
                        & + (-64._dl*H**3 + 48._dl*H*Ht - 4._dl*Ht2)*gpinu )
                end do
            end if
            Ht4 = (-3._dl*gdddpn - 27._dl*gdpn*H**2 - 12._dl*gdpn*Ht - 3._dl*gpn*Ht2 + grhom*Ht2 + grhon*Ht2 + 4._dl*grhor*Ht2 + grhov*Ht2 - 3._dl*a*dwDE*grhov*Ht2 + 6._dl*grhov*Ht2*wDE + 9._dl*grhov*Ht2*wDE**2 - 3._dl*H*(5._dl*gddpn + 9._dl*gpn*Ht + grhom*Ht + grhon*Ht + 8._dl*grhor*Ht + grhov*Ht*(3._dl*a**2*d2wDE + (1._dl + 3._dl*wDE)**3 - 3._dl*a*dwDE*(2._dl + 9._dl*wDE))) + H**3*(-15._dl*gpn + grhom + grhon + 16._dl*grhor + grhov*(-3._dl*a**3*d3wDE + (1._dl + 3._dl*wDE)**4 + 3._dl*a**2*(d2wDE + 9._dl*dwDE**2 + 12._dl*d2wDE*wDE) - 9._dl*a*dwDE*(1._dl + 8._dl*wDE + 18._dl*wDE**2))))/6._dl

            Ode    = grhov / 3._dl / H**2
            dOde   = -(grhov*(2._dl*Ht + H**2*(1._dl + 3._dl*wDE)))/(3._dl*a*H**4)
            d2Ode  = (grhov*(8._dl*Ht**2 - 2._dl*H*Ht2 + 6._dl*H**2*(Ht + 2._dl*Ht*wDE) + H**4*(2._dl - 3._dl*a*dwDE + 9._dl*wDE + 9._dl*wDE**2)))/(3._dl*a**2*H**6)
            d3Ode  = -(grhov*(48._dl*Ht**3 - 26._dl*H*Ht*Ht2 - 6._dl*H**3*Ht2*(2._dl + 3._dl*wDE) + 2._dl*H**4*Ht*(11._dl - 9._dl*a*dwDE + 36._dl*wDE + 27._dl*wDE**2) + 3._dl*H**6*(2._dl + a**2*d2wDE - 5._dl*a*dwDE + (11._dl - 9._dl*a*dwDE)*wDE + 18._dl*wDE**2 + 9._dl*wDE**3) + 2._dl*H**2*(Ht3 + 12._dl*Ht**2*(2._dl + 3._dl*wDE))))/(3._dl*a**3*H**8)
            d4Ode  = (grhov*(384._dl*Ht**4 - 326._dl*H*Ht**2*Ht2 - 2._dl*H**5*Ht2*(35._dl - 18._dl*a*dwDE + 90._dl*wDE + 54._dl*wDE**2) + H**2*(26._dl*Ht2**2 + 38._dl*Ht*Ht3 + 96._dl*Ht**3*(5._dl + 6._dl*wDE)) - 2._dl*H**3*(Ht4 + 26*Ht*Ht2*(5 + 6*wDE)) + 4*H**6*Ht*(25 + 6*a**2*d2wDE + 105*wDE + 135*wDE**2 + 54*wDE**3 - 3*a*dwDE*(13._dl + 18._dl*wDE)) + 4._dl*H**4*(Ht3*(5._dl + 6._dl*wDE) + 2._dl*Ht**2*(35._dl - 18._dl*a*dwDE + 90._dl*wDE + 54._dl*wDE**2)) + 3._dl*H**8*(8._dl - a**3*d3wDE + 50._dl*wDE + 105._dl*wDE**2 + 90._dl*wDE**3 + 27._dl*wDE**4 - 2._dl*a*dwDE*(13._dl + 39._dl*wDE + 27._dl*wDE**2) + a**2*(9._dl*dwDE**2 + d2wDE*(7._dl + 12._dl*wDE)))))/(3._dl*a**4*H**10)
        else
            Ode = 0._dl
            dOde = 0._dl
            d2Ode = 0._dl
            d3Ode = 0._dl
            d4Ode = 0._dl
        end if

        eft_cache%EFTOmegaV    = self%PureEFTOmega%value(a) + self%PureEFTOmega_ODE%value(Ode)
        eft_cache%EFTOmegaP    = self%PureEFTOmega%first_derivative(a) + (dOde*self%PureEFTOmega_ODE%first_derivative(Ode))
        eft_cache%EFTOmegaPP   = self%PureEFTOmega%second_derivative(a) + (d2Ode*self%PureEFTOmega_ODE%first_derivative(Ode) + dOde**2*self%PureEFTOmega_ODE%second_derivative(Ode))
        eft_cache%EFTOmegaPPP  = self%PureEFTOmega%third_derivative(a) + (d3Ode*self%PureEFTOmega_ODE%first_derivative(Ode) + 3*d2Ode*dOde*self%PureEFTOmega_ODE%second_derivative(Ode) + dOde**3*self%PureEFTOmega_ODE%third_derivative(Ode))
        eft_cache%EFTOmegaPPPP = self%PureEFTOmega%fourth_derivative(a) + (d4Ode*self%PureEFTOmega_ODE%first_derivative(Ode) + 3._dl*d2Ode**2*self%PureEFTOmega_ODE%second_derivative(Ode) + 4._dl*d3Ode*dOde*self%PureEFTOmega_ODE%second_derivative(Ode) + 6._dl*d2Ode*dOde**2*self%PureEFTOmega_ODE%third_derivative(Ode) + dOde**4*self%PureEFTOmega_ODE%fourth_derivative(Ode))
        eft_cache%EFTc         = ( eft_cache%adotoa**2 - eft_cache%Hdot )*( eft_cache%EFTOmegaV + 0.5_dl*a*eft_cache%EFTOmegaP ) &
            & -0.5_dl*( a*eft_cache%adotoa )**2*eft_cache%EFTOmegaPP&
            & +0.5_dl*eft_cache%grhov_t*( 1._dl+ self%PureEFTwDE%value(a) )
        eft_cache%EFTLambda    = -eft_cache%EFTOmegaV*( 2._dl*eft_cache%Hdot +eft_cache%adotoa**2 ) &
            & -a*eft_cache%EFTOmegaP*( 2._dl*eft_cache%adotoa**2 + eft_cache%Hdot ) &
            & -( a*eft_cache%adotoa )**2*eft_cache%EFTOmegaPP &
            & +self%PureEFTwDE%value(a)*eft_cache%grhov_t
        eft_cache%EFTcdot      = +0.5_dl*eft_cache%adotoa*eft_cache%grhov_t*( -3._dl*(1._dl +self%PureEFTwDE%value(a))**2 + a*self%PureEFTwDE%first_derivative(a) ) &
            & -eft_cache%EFTOmegaV*( eft_cache%Hdotdot -4._dl*eft_cache%adotoa*eft_cache%Hdot +2._dl*eft_cache%adotoa**3 ) &
            & +0.5_dl*a*eft_cache%EFTOmegaP*( -eft_cache%Hdotdot +eft_cache%adotoa*eft_cache%Hdot +eft_cache%adotoa**3) &
            & +0.5_dl*a**2*eft_cache%adotoa*eft_cache%EFTOmegaPP*( eft_cache%adotoa**2 -3._dl*eft_cache%Hdot ) &
            & -0.5_dl*(a*eft_cache%adotoa)**3*eft_cache%EFTOmegaPPP
            eft_cache%EFTcdotdot   = (-(a**4*eft_cache%EFTOmegaPPPP*eft_cache%adotoa**4) - 6._dl*a**3*eft_cache%EFTOmegaPPP*eft_cache%adotoa**2*eft_cache%Hdot + a**2*eft_cache%EFTOmegaPP*(eft_cache%adotoa**4 + 4._dl*eft_cache%adotoa**2*eft_cache%Hdot - 3._dl*eft_cache%Hdot**2 - 4._dl*eft_cache%adotoa*eft_cache%Hdotdot) + a*eft_cache%EFTOmegaP*(-5._dl*eft_cache%adotoa**4 + 10._dl*eft_cache%adotoa**2*eft_cache%Hdot + eft_cache%Hdot**2 - eft_cache%Hdotdotdot) + 2._dl*(4._dl*eft_cache%adotoa**4 - 14._dl*eft_cache%adotoa**2*eft_cache%Hdot + 4._dl*eft_cache%Hdot**2 + 6._dl*eft_cache%adotoa*eft_cache%Hdotdot - eft_cache%Hdotdotdot)*eft_cache%EFTOmegaV + grhov*eft_cache%Hdot*(a*self%PureEFTwDE%first_derivative(a) - 3._dl*(1 + self%PureEFTwDE%value(a))**2) + grhov*eft_cache%adotoa**2*(a**2*self%PureEFTwDE%second_derivative(a) + 9._dl*(1 + self%PureEFTwDE%value(a))**3 - a*self%PureEFTwDE%first_derivative(a)*(8._dl + 9._dl*self%PureEFTwDE%value(a))))/2._dl
        eft_cache%EFTLambdadot = -2._dl*eft_cache%EFTOmegaV*( eft_cache%Hdotdot -eft_cache%adotoa*eft_cache%Hdot -eft_cache%adotoa**3 ) &
            & -a*eft_cache%EFTOmegaP*( +eft_cache%Hdotdot +5._dl*eft_cache%adotoa*eft_cache%Hdot -eft_cache%adotoa**3  ) &
            & -a**2*eft_cache%EFTOmegaPP*eft_cache%adotoa*( +2._dl*eft_cache%adotoa**2 +3._dl*eft_cache%Hdot )&
            & -(a*eft_cache%adotoa)**3*eft_cache%EFTOmegaPPP &
            & +eft_cache%grhov_t*eft_cache%adotoa*( a*self%PureEFTwDE%first_derivative(a) -3._dl*self%PureEFTwDE%value(a)*(1._dl +self%PureEFTwDE%value(a) ))
        eft_cache%EFTLambdadotdot = -(a**4*eft_cache%EFTOmegaPPPP*eft_cache%adotoa**4) + eft_cache%EFTOmegaP*(-(a*eft_cache%Hdotdotdot) - 5._dl*a*eft_cache%Hdot**2 - 6._dl*a*eft_cache%Hdotdot*eft_cache%adotoa + 10._dl*a*eft_cache%Hdot*eft_cache%adotoa**2 + a*eft_cache%adotoa**4) + eft_cache%EFTOmegaPP*(-3._dl*a**2*eft_cache%Hdot**2 - 4._dl*a**2*eft_cache%Hdotdot*eft_cache%adotoa - 11._dl*a**2*eft_cache%Hdot*eft_cache%adotoa**2 + a**2*eft_cache%adotoa**4) + eft_cache%EFTOmegaPPP*(-6._dl*a**3*eft_cache%Hdot*eft_cache%adotoa**2 - 3._dl*a**3*eft_cache%adotoa**4) + (-2._dl*eft_cache%Hdotdotdot + 2._dl*eft_cache%Hdot**2 + 6._dl*eft_cache%Hdotdot*eft_cache%adotoa + 2._dl*eft_cache%Hdot*eft_cache%adotoa**2 - 4._dl*eft_cache%adotoa**4)*eft_cache%EFTOmegaV + eft_cache%Hdot*eft_cache%grhov_t*(a*self%PureEFTwDE%first_derivative(a) - 3._dl*self%PureEFTwDE%value(a) - 3._dl*self%PureEFTwDE%value(a)**2) + eft_cache%adotoa**2*eft_cache%grhov_t*(a**2*self%PureEFTwDE%second_derivative(a) + 9._dl*self%PureEFTwDE%value(a) + 18._dl*self%PureEFTwDE%value(a)**2 + 9._dl*self%PureEFTwDE%value(a)**3 - a*self%PureEFTwDE%first_derivative(a)*(5._dl + 9._dl*self%PureEFTwDE%value(a)))

    end subroutine EFTCAMBPureEFTstdBackgroundEFTFunctions

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes the value of the background EFT functions at a given time.
    subroutine EFTCAMBPureEFTstdSecondOrderEFTFunctions( self, a, eft_par_cache, eft_cache )

        implicit none

        class(EFTCAMB_std_pure_EFT)                  :: self          !< the base class
        real(dl), intent(in)                         :: a             !< the input scale factor
        type(TEFTCAMB_parameter_cache), intent(inout) :: eft_par_cache !< the EFTCAMB parameter cache that contains all the physical parameters.
        type(TEFTCAMB_timestep_cache ), intent(inout) :: eft_cache     !< the EFTCAMB timestep cache that contains all the physical values.

        real(dl) :: H, Ht, Ht2, Ht3, Ht4
        real(dl) :: grhor, grhom, grhov, grhon, gpn, gdpn, gddpn, gdddpn
        real(dl) :: grhonu, gpinu, gpinudot, gpinudotdot, gpinudotdotdot, grhormass_t
        integer  :: nu_i
        real(dl) :: wDE, dwDE, d2wDE, d3wDE
        real(dl) :: Ode, dOde, d2Ode, d3Ode, d4Ode

        grhom = eft_cache%grhob_t + eft_cache%grhoc_t
        grhor = eft_cache%grhor_t + eft_cache%grhog_t
        grhov = eft_cache%grhov_t

        wDE   = self%PureEFTwDE%value(a)
        dwDE  = self%PureEFTwDE%first_derivative(a)
        d2wDE = self%PureEFTwDE%second_derivative(a)
        d3wDE = self%PureEFTwDE%third_derivative(a)
        
        H   = eft_cache%adotoa
        Ht  = eft_cache%Hdot
        Ht2 = eft_cache%Hdotdot
        Ht3 = eft_cache%Hdotdotdot
        if ( self%PureEFTOmegaDE ) then
            grhon = 0._dl
            gpn = 0._dl
            gdpn = 0._dl
            gddpn = 0._dl
            gdddpn = 0._dl
            if ( eft_par_cache%Num_Nu_Massive /= 0 ) then
                do nu_i = 1, eft_par_cache%Nu_mass_eigenstates
                    grhonu         = 0._dl
                    gpinu          = 0._dl
                    gpinudot       = 0._dl
                    gpinudotdot    = 0._dl
                    gpinudotdotdot = 0._dl
                    grhormass_t = eft_par_cache%grhormass(nu_i)/a**2
                    call ThermalNuBack%rho_P(a*eft_par_cache%nu_masses(nu_i), grhonu, gpinu )
                    grhon          = grhon +grhormass_t*grhonu
                    gpn            = gpn  +grhormass_t*gpinu
                    gpinudot       = ThermalNuBack%pidot( a*eft_par_cache%nu_masses(nu_i), H, gpinu )
                    gdpn           = gdpn  + grhormass_t*( gpinudot - 4._dl*H*gpinu)
                    gpinudotdot    = ThermalNuBack%pidotdot( a*eft_par_cache%nu_masses(nu_i), H, Ht, gpinu, gpinudot)
                    gddpn          = gddpn - 4._dl*H*grhormass_t*( gpinudot - 4._dl*H*gpinu )+ grhormass_t*( gpinudotdot - 4._dl*Ht*gpinu - 4._dl*H*gpinudot )
                    gpinudotdotdot = ThermalNuBack%pidotdotdot(a*eft_par_cache%nu_masses(nu_i),H,Ht,Ht2,gpinu,gpinudot,gpinudotdot)
                    gdddpn = gdddpn + grhormass_t*(gpinudotdotdot &
                        & - 12._dl*H*gpinudotdot &
                        & + (48._dl*H**2 - 12._dl*Ht)*gpinudot &
                        & + (-64._dl*H**3 + 48._dl*H*Ht - 4._dl*Ht2)*gpinu )
                end do
            end if
            Ht4 = (-3._dl*gdddpn - 27._dl*gdpn*H**2 - 12._dl*gdpn*Ht - 3._dl*gpn*Ht2 + grhom*Ht2 + grhon*Ht2 + 4._dl*grhor*Ht2 + grhov*Ht2 - 3._dl*a*dwDE*grhov*Ht2 + 6._dl*grhov*Ht2*wDE + 9._dl*grhov*Ht2*wDE**2 - 3._dl*H*(5._dl*gddpn + 9._dl*gpn*Ht + grhom*Ht + grhon*Ht + 8._dl*grhor*Ht + grhov*Ht*(3._dl*a**2*d2wDE + (1._dl + 3._dl*wDE)**3 - 3._dl*a*dwDE*(2._dl + 9._dl*wDE))) + H**3*(-15._dl*gpn + grhom + grhon + 16._dl*grhor + grhov*(-3._dl*a**3*d3wDE + (1._dl + 3._dl*wDE)**4 + 3._dl*a**2*(d2wDE + 9._dl*dwDE**2 + 12._dl*d2wDE*wDE) - 9._dl*a*dwDE*(1._dl + 8._dl*wDE + 18._dl*wDE**2))))/6._dl

            Ode    = grhov / 3._dl / H**2
            dOde   = -(grhov*(2._dl*Ht + H**2*(1._dl + 3._dl*wDE)))/(3._dl*a*H**4)
            d2Ode  = (grhov*(8._dl*Ht**2 - 2._dl*H*Ht2 + 6._dl*H**2*(Ht + 2._dl*Ht*wDE) + H**4*(2._dl - 3._dl*a*dwDE + 9._dl*wDE + 9._dl*wDE**2)))/(3._dl*a**2*H**6)
            d3Ode  = -(grhov*(48._dl*Ht**3 - 26._dl*H*Ht*Ht2 - 6._dl*H**3*Ht2*(2._dl + 3._dl*wDE) + 2._dl*H**4*Ht*(11._dl - 9._dl*a*dwDE + 36._dl*wDE + 27._dl*wDE**2) + 3._dl*H**6*(2._dl + a**2*d2wDE - 5._dl*a*dwDE + (11._dl - 9._dl*a*dwDE)*wDE + 18._dl*wDE**2 + 9._dl*wDE**3) + 2._dl*H**2*(Ht3 + 12._dl*Ht**2*(2._dl + 3._dl*wDE))))/(3._dl*a**3*H**8)
            d4Ode  = (grhov*(384._dl*Ht**4 - 326._dl*H*Ht**2*Ht2 - 2._dl*H**5*Ht2*(35._dl - 18._dl*a*dwDE + 90._dl*wDE + 54._dl*wDE**2) + H**2*(26._dl*Ht2**2 + 38._dl*Ht*Ht3 + 96._dl*Ht**3*(5._dl + 6._dl*wDE)) - 2._dl*H**3*(Ht4 + 26*Ht*Ht2*(5 + 6*wDE)) + 4*H**6*Ht*(25 + 6*a**2*d2wDE + 105*wDE + 135*wDE**2 + 54*wDE**3 - 3*a*dwDE*(13._dl + 18._dl*wDE)) + 4._dl*H**4*(Ht3*(5._dl + 6._dl*wDE) + 2._dl*Ht**2*(35._dl - 18._dl*a*dwDE + 90._dl*wDE + 54._dl*wDE**2)) + 3._dl*H**8*(8._dl - a**3*d3wDE + 50._dl*wDE + 105._dl*wDE**2 + 90._dl*wDE**3 + 27._dl*wDE**4 - 2._dl*a*dwDE*(13._dl + 39._dl*wDE + 27._dl*wDE**2) + a**2*(9._dl*dwDE**2 + d2wDE*(7._dl + 12._dl*wDE)))))/(3._dl*a**4*H**10)
        else
            Ode = 0._dl
            dOde = 0._dl
            d2Ode = 0._dl
            d3Ode = 0._dl
            d4Ode = 0._dl
        end if

        eft_cache%EFTGamma1V  = self%PureEFTGamma1%value(a) + self%PureEFTGamma1_ODE%value(Ode)
        eft_cache%EFTGamma1P  = self%PureEFTGamma1%first_derivative(a) + (dOde*self%PureEFTGamma1_ODE%first_derivative(Ode))
        eft_cache%EFTGamma1PP  = self%PureEFTGamma1%second_derivative(a) + (d2Ode*self%PureEFTGamma1_ODE%first_derivative(Ode) + dOde**2*self%PureEFTGamma1_ODE%second_derivative(Ode))
        eft_cache%EFTGamma2V  = self%PureEFTGamma2%value(a) + self%PureEFTGamma2_ODE%value(Ode)
        eft_cache%EFTGamma2P  = self%PureEFTGamma2%first_derivative(a) + (dOde*self%PureEFTGamma2_ODE%first_derivative(Ode))
        eft_cache%EFTGamma2PP = self%PureEFTGamma2%second_derivative(a) + (d2Ode*self%PureEFTGamma2_ODE%first_derivative(Ode) + dOde**2*self%PureEFTGamma2_ODE%second_derivative(Ode))
        eft_cache%EFTGamma2PPP = self%PureEFTGamma2%third_derivative(a) + (d3Ode*self%PureEFTGamma2_ODE%first_derivative(Ode) + 3._dl*d2Ode*dOde*self%PureEFTGamma2_ODE%second_derivative(Ode) + dOde**3*self%PureEFTGamma2_ODE%third_derivative(Ode))
        eft_cache%EFTGamma3V  = self%PureEFTGamma3%value(a) + self%PureEFTGamma3_ODE%value(Ode)
        eft_cache%EFTGamma3P  = self%PureEFTGamma3%first_derivative(a) + (dOde*self%PureEFTGamma3_ODE%first_derivative(Ode))
        eft_cache%EFTGamma3PP = self%PureEFTGamma3%second_derivative(a) + (d2Ode*self%PureEFTGamma3_ODE%first_derivative(Ode) + dOde**2*self%PureEFTGamma3_ODE%second_derivative(Ode))
        eft_cache%EFTGamma3PPP = self%PureEFTGamma3%third_derivative(a) + (d3Ode*self%PureEFTGamma3_ODE%first_derivative(Ode) + 3._dl*d2Ode*dOde*self%PureEFTGamma3_ODE%second_derivative(Ode) + dOde**3*self%PureEFTGamma3_ODE%third_derivative(Ode))
        eft_cache%EFTGamma3PPPP = self%PureEFTGamma3%fourth_derivative(a) + (d4Ode*self%PureEFTGamma3_ODE%first_derivative(Ode) + 3._dl*d2Ode**2*self%PureEFTGamma3_ODE%second_derivative(Ode) + 4._dl*d3Ode*dOde*self%PureEFTGamma3_ODE%second_derivative(Ode) + 6._dl*d2Ode*dOde**2*self%PureEFTGamma3_ODE%third_derivative(Ode) + dOde**4*self%PureEFTGamma3_ODE%fourth_derivative(Ode))
        if ( self%PureEFTHorndeski ) then
            eft_cache%EFTGamma4V  = -eft_cache%EFTGamma3V
            eft_cache%EFTGamma4P  = -eft_cache%EFTGamma3P
            eft_cache%EFTGamma4PP = -eft_cache%EFTGamma3PP
            eft_cache%EFTGamma5V  = +0.5_dl*eft_cache%EFTGamma3V
            eft_cache%EFTGamma5P  = +0.5_dl*eft_cache%EFTGamma3P
            eft_cache%EFTGamma6V  = 0._dl
            eft_cache%EFTGamma6P  = 0._dl
        else
            eft_cache%EFTGamma4V  = self%PureEFTGamma4%value(a) + self%PureEFTGamma4_ODE%value(Ode)
            eft_cache%EFTGamma4P  = self%PureEFTGamma4%first_derivative(a) + (dOde*self%PureEFTGamma4_ODE%first_derivative(Ode))
            eft_cache%EFTGamma4PP = self%PureEFTGamma4%second_derivative(a) + (d2Ode*self%PureEFTGamma4_ODE%first_derivative(Ode) + dOde**2*self%PureEFTGamma4_ODE%second_derivative(Ode))
            eft_cache%EFTGamma5V  = self%PureEFTGamma5%value(a) + self%PureEFTGamma5_ODE%value(Ode)
            eft_cache%EFTGamma5P  = self%PureEFTGamma5%first_derivative(a) + (dOde*self%PureEFTGamma5_ODE%first_derivative(Ode))
            eft_cache%EFTGamma6V  = self%PureEFTGamma6%value(a) + self%PureEFTGamma6_ODE%value(Ode)
            eft_cache%EFTGamma6P  = self%PureEFTGamma6%first_derivative(a) + (dOde*self%PureEFTGamma6_ODE%first_derivative(Ode))
        end if

    end subroutine EFTCAMBPureEFTstdSecondOrderEFTFunctions

    ! ---------------------------------------------------------------------------------------------
    !> Function that computes dtauda = 1/sqrt(a^2H^2).For pure EFT std this has to be overridden
    !! for performance reasons.
    function EFTCAMBPureEFTstdComputeDtauda( self, a, eft_par_cache, eft_cache )

        implicit none

        class(EFTCAMB_std_pure_EFT)                  :: self          !< the base class
        real(dl), intent(in)                         :: a             !< the input scale factor
        type(TEFTCAMB_parameter_cache), intent(inout) :: eft_par_cache !< the EFTCAMB parameter cache that contains all the physical parameters.
        type(TEFTCAMB_timestep_cache ), intent(inout) :: eft_cache     !< the EFTCAMB timestep cache that contains all the physical values.

        real(dl) :: EFTCAMBPureEFTstdComputeDtauda               !< the output dtauda

        real(dl) :: temp

        real(dl) :: x, mu
        integer  :: ind

        x = log(a)
        if ( x < self%rhov_t%x_initial ) then
            temp = eft_cache%grhoa2 + a**4*self%rhov_t%y(1)
        else
            call self%rhov_t%precompute( x, ind, mu )
            temp = eft_cache%grhoa2 + a**4*self%rhov_t%value( x, index=ind, coeff=mu )
        end if
        EFTCAMBPureEFTstdComputeDtauda = sqrt(3._dl/temp)

    end function EFTCAMBPureEFTstdComputeDtauda

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes adotoa = H and its two derivatives wrt conformal time.
    subroutine EFTCAMBPureEFTstdComputeAdotoa( self, a, eft_par_cache, eft_cache )

        implicit none

        class(EFTCAMB_std_pure_EFT)                  :: self          !< the base class
        real(dl), intent(in)                         :: a             !< the input scale factor
        type(TEFTCAMB_parameter_cache), intent(inout) :: eft_par_cache !< the EFTCAMB parameter cache that contains all the physical parameters.
        type(TEFTCAMB_timestep_cache ), intent(inout) :: eft_cache     !< the EFTCAMB timestep cache that contains all the physical values.

        real(dl) :: x, mu
        integer  :: ind

        x = log(a)
        if ( x < self%rhov_t%x_initial ) then
            eft_cache%grhov_t = a**2*self%rhov_t%y(1)
        else
            call self%rhov_t%precompute( x, ind, mu )
            eft_cache%grhov_t = a**2*self%rhov_t%value( x, index=ind, coeff=mu )
        end if
        eft_cache%adotoa  = sqrt( ( eft_cache%grhom_t +eft_cache%grhov_t )/3._dl )

    end subroutine EFTCAMBPureEFTstdComputeAdotoa

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes the three derivatives wrt conformal time of H.
    subroutine EFTCAMBPureEFTstdComputeHubbleDer( self, a, eft_par_cache, eft_cache )

        implicit none

        class(EFTCAMB_std_pure_EFT)                  :: self          !< the base class
        real(dl), intent(in)                         :: a             !< the input scale factor
        type(TEFTCAMB_parameter_cache), intent(inout) :: eft_par_cache !< the EFTCAMB parameter cache that contains all the physical parameters.
        type(TEFTCAMB_timestep_cache ), intent(inout) :: eft_cache     !< the EFTCAMB timestep cache that contains all the physical values.

        eft_cache%gpiv_t  = self%PureEFTwDE%value(a)*eft_cache%grhov_t
        eft_cache%Hdot    = -0.5_dl*( eft_cache%adotoa**2 +eft_cache%gpresm_t +eft_cache%gpiv_t )
        eft_cache%Hdotdot = eft_cache%adotoa*( ( eft_cache%grhob_t +eft_cache%grhoc_t)/6._dl +2._dl*( eft_cache%grhor_t +eft_cache%grhog_t)/3._dl ) &
                          & +eft_cache%adotoa*eft_cache%grhov_t*( 1._dl/6._dl +self%PureEFTwDE%value(a) +1.5_dl*self%PureEFTwDE%value(a)**2 -0.5_dl*a*self%PureEFTwDE%first_derivative(a) ) &
                          & +eft_cache%adotoa*eft_cache%grhonu_tot/6._dl -0.5_dl*eft_cache%adotoa*eft_cache%gpinu_tot -0.5_dl*eft_cache%gpinudot_tot
        eft_cache%Hdotdotdot = (2._dl*(eft_cache%grhob_t + eft_cache%grhoc_t)*eft_cache%adotoa**2 + eft_cache%Hdot &
                          &* (eft_cache%grhob_t + eft_cache%grhoc_t) - 3._dl*eft_cache%adotoa**2 *(eft_cache%grhob_t + eft_cache%grhoc_t))/(6._dl) &
                          &+ (2._dl/3._dl)*(2._dl*eft_cache%adotoa**2 * (eft_cache%grhor_t + eft_cache%grhog_t) + eft_cache%Hdot &
                          &* (eft_cache%grhor_t + eft_cache%grhog_t) - 4._dl*eft_cache%adotoa**2 * (eft_cache%grhor_t + eft_cache%grhog_t)) &
                          &+ (2._dl*eft_cache%adotoa**2 * eft_cache%grhov_t - 3._dl*eft_cache%adotoa**2 * eft_cache%grhov_t &
                          &* (1._dl + self%PureEFTwDE%value(a)) + eft_cache%Hdot*eft_cache%grhov_t)*(1._dl/6._dl + self%PureEFTwDE%value(a) &
                          &+ (3._dl/2._dl)*self%PureEFTwDE%value(a)**2 - (1._dl/2._dl)*a*self%PureEFTwDE%first_derivative(a)) &
                          &+ a*eft_cache%adotoa**2 * eft_cache%grhov_t * ((1._dl/2._dl)*self%PureEFTwDE%first_derivative(a) &
                          &+ 3._dl*self%PureEFTwDE%value(a)*self%PureEFTwDE%first_derivative(a) - (a/2._dl)*self%PureEFTwDE%second_derivative(a)) &
                          &+ (eft_cache%adotoa**2/3._dl)*eft_cache%grhonu_tot - eft_cache%adotoa**2 * eft_cache%gpinu_tot &
                          &- (3._dl/2._dl)*eft_cache%adotoa*eft_cache%gpinudot_tot + (eft_cache%Hdot/6._dl)*eft_cache%grhonu_tot &
                          &+ (eft_cache%adotoa/6._dl)*eft_cache%grhonudot_tot - (eft_cache%Hdot/2._dl)*eft_cache%gpinu_tot &
                          &- (1._dl/2._dl)*eft_cache%gpinudotdot_tot
    end subroutine EFTCAMBPureEFTstdComputeHubbleDer

    ! ---------------------------------------------------------------------------------------------
    !> Function that computes model specific stability requirements.
    function EFTCAMBPureEFTstdAdditionalModelStability( self, a, eft_par_cache, eft_cache )

        implicit none

        class(EFTCAMB_std_pure_EFT)                  :: self          !< the base class
        real(dl), intent(in)                         :: a             !< the input scale factor.
        type(TEFTCAMB_parameter_cache), intent(inout) :: eft_par_cache !< the EFTCAMB parameter cache that contains all the physical parameters.
        type(TEFTCAMB_timestep_cache ), intent(inout) :: eft_cache     !< the EFTCAMB timestep cache that contains all the physical values.

        logical :: EFTCAMBPureEFTstdAdditionalModelStability          !< the return value of the stability computation. True if the model specific stability criteria are met, false otherwise.

        EFTCAMBPureEFTstdAdditionalModelStability = .True.
        if ( self%PureEFTwDE%value(a) > -1._dl/3._dl ) EFTCAMBPureEFTstdAdditionalModelStability = .False.

    end function EFTCAMBPureEFTstdAdditionalModelStability

    ! ---------------------------------------------------------------------------------------------

end module EFTCAMB_pure_EFT_std

!----------------------------------------------------------------------------------------
