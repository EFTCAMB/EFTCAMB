!----------------------------------------------------------------------------------------
!
! This file is part of EFTCAMB.
!
! Copyright (C) 2013-2016 by the EFTCAMB authors
!
! The EFTCAMB code is free software;
! You can use it, redistribute it, and/or modify it under the terms
! of the GNU General Public License as published by the Free Software Foundation;
! either version 3 of the License, or (at your option) any later version.
! The full text of the license can be found in the file eftcamb/LICENSE at
! the top level of the EFTCAMB distribution.
!
!----------------------------------------------------------------------------------------

!> @file 07p1_Pure_EFT_std.f90
!! This file contains the definition of the Pure EFT model in which the EFT is described
!! by six functions of time and w_DE. Please refer to the numerical notes for details.


!----------------------------------------------------------------------------------------
!> This module contains the definition of the Pure EFT model in which the EFT is described
!! by six functions of time and w_DE. Please refer to the numerical notes for details.

!> @author Bin Hu, Marco Raveri, Simone Peirone

module EFTCAMB_pure_EFT_std

    use precision
    use IniFile
    use AMLutils
    use EFTCAMB_cache
    use EFT_def
    use EFTCAMB_abstract_parametrizations_1D
    use EFTCAMB_neutral_parametrization_1D
    use EFTCAMB_constant_parametrization_1D
    use EFTCAMB_linear_parametrizations_1D
    use EFTCAMB_power_law_parametrizations_1D
    use EFTCAMB_exponential_parametrizations_1D
    use EFTCAMB_CPL_parametrizations_1D
    use EFTCAMB_JBP_parametrizations_1D
    use EFTCAMB_turning_point_parametrizations_1D
    use EFTCAMB_taylor_parametrizations_1D
    use EFTCAMB_abstract_model_designer

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

        ! selection flag for Horndeski:
        logical  :: PureEFTHorndeski    !< Selects wether to use the Horndeski bound on EFT functions.

        ! the pure EFT functions:
        class( parametrized_function_1D ), allocatable :: PureEFTOmega    !< The pure EFT function Omega.
        class( parametrized_function_1D ), allocatable :: PureEFTwDE      !< The pure EFT function w_DE.
        class( parametrized_function_1D ), allocatable :: PureEFTGamma1   !< The pure EFT function Gamma1.
        class( parametrized_function_1D ), allocatable :: PureEFTGamma2   !< The pure EFT function Gamma2.
        class( parametrized_function_1D ), allocatable :: PureEFTGamma3   !< The pure EFT function Gamma3.
        class( parametrized_function_1D ), allocatable :: PureEFTGamma4   !< The pure EFT function Gamma4.
        class( parametrized_function_1D ), allocatable :: PureEFTGamma5   !< The pure EFT function Gamma5.
        class( parametrized_function_1D ), allocatable :: PureEFTGamma6   !< The pure EFT function Gamma6.

    contains

        ! initialization of the model:
        procedure :: read_model_selection            => EFTCAMBPureEFTstdReadModelSelectionFromFile  !< subroutine that reads the parameters of the model from file
        procedure :: allocate_model_selection        => EFTCAMBPureEFTstdAllocateModelSelection      !< subroutine that allocates the model selection.
        procedure :: init_model_parameters           => EFTCAMBPureEFTstdInitModelParameters         !< subroutine that initializes the model parameters based on the values found in an input array.
        procedure :: init_model_parameters_from_file => EFTCAMBPureEFTstdInitModelParametersFromFile !< subroutine that reads the parameters of the model from file.
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
        procedure :: compute_H_derivs                  => EFTCAMBPureEFTstdComputeHubbleDer         !< subroutine that computes the two derivatives wrt conformal time of H.
        ! stability procedures:
        procedure :: additional_model_stability        => EFTCAMBPureEFTstdAdditionalModelStability !< function that computes model specific stability requirements.

    end type EFTCAMB_std_pure_EFT

    ! ---------------------------------------------------------------------------------------------

contains

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that reads the parameters of the model from file.
    subroutine EFTCAMBPureEFTstdReadModelSelectionFromFile( self, Ini )

        implicit none

        class(EFTCAMB_std_pure_EFT) :: self   !< the base class
        type(TIniFile)              :: Ini    !< Input ini file

        ! read model selection flags:
        self%PureEFTmodelOmega  = Ini_Read_Int_File( Ini, 'PureEFTmodelOmega'  , 0 )
        self%EFTwDE             = Ini_Read_Int_File( Ini, 'EFTwDE'             , 0 )
        self%PureEFTmodelGamma1 = Ini_Read_Int_File( Ini, 'PureEFTmodelGamma1' , 0 )
        self%PureEFTmodelGamma2 = Ini_Read_Int_File( Ini, 'PureEFTmodelGamma2' , 0 )
        self%PureEFTmodelGamma3 = Ini_Read_Int_File( Ini, 'PureEFTmodelGamma3' , 0 )
        self%PureEFTmodelGamma4 = Ini_Read_Int_File( Ini, 'PureEFTmodelGamma4' , 0 )
        self%PureEFTmodelGamma5 = Ini_Read_Int_File( Ini, 'PureEFTmodelGamma5' , 0 )
        self%PureEFTmodelGamma6 = Ini_Read_Int_File( Ini, 'PureEFTmodelGamma6' , 0 )
        ! read the Horndeski flag:
        self%PureEFTHorndeski   = Ini_Read_Logical_File( Ini, 'PureEFTHorndeski' , .false. )

    end subroutine EFTCAMBPureEFTstdReadModelSelectionFromFile

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that allocates the model selection.
    subroutine EFTCAMBPureEFTstdAllocateModelSelection( self )

        implicit none

        class(EFTCAMB_std_pure_EFT)                       :: self              !< the base class

        ! allocate Omega:
        if ( allocated(self%PureEFTOmega) ) deallocate(self%PureEFTOmega)
        select case ( self%PureEFTmodelOmega )
            case(0)
                allocate( zero_parametrization_1D::self%PureEFTOmega )
            case(1)
                allocate( constant_parametrization_1D::self%PureEFTOmega )
            case(2)
                allocate( linear_parametrization_1D::self%PureEFTOmega )
            case(3)
                allocate( power_law_parametrization_1D::self%PureEFTOmega )
                call self%PureEFTOmega%set_param_names( ['EFTOmega0  ', 'EFTOmegaExp'], ['\Omega_0^{\rm EFT}', 'n^{\rm EFT}       '] )
            case(4)
                allocate( exponential_parametrization_1D::self%PureEFTOmega )
                call self%PureEFTOmega%set_param_names( ['EFTOmega0  ', 'EFTOmegaExp'], ['\Omega_0^{\rm EFT}', 'n^{\rm EFT}       '] )
            case default
                write(*,'(a,I3)') 'No model corresponding to PureEFTmodelOmega =', self%PureEFTmodelOmega
                write(*,'(a)')    'Please select an appropriate model.'
        end select
        ! allocate wDE:
        if ( allocated(self%PureEFTwDE) ) deallocate(self%PureEFTwDE)
        select case ( self%EFTwDE )
            case(0)
                allocate( wDE_LCDM_parametrization_1D::self%PureEFTwDE )
            case(1)
                allocate( constant_parametrization_1D::self%PureEFTwDE )
            case(2)
                allocate( CPL_parametrization_1D::self%PureEFTwDE )
                call self%PureEFTwDE%set_param_names( ['EFTw0', 'EFTwa'], ['w_0', 'w_a'] )
            case(3)
                allocate( JBP_parametrization_1D::self%PureEFTwDE )
                call self%PureEFTwDE%set_param_names( ['EFTw0', 'EFTwa', 'EFTwn'], [ 'w_0', 'w_a', 'n  ' ] )
            case(4)
                allocate( turning_point_parametrization_1D::self%PureEFTwDE )
                call self%PureEFTwDE%set_param_names( ['EFTw0 ', 'EFTwa ', 'EFTwat'], ['w_0', 'w_a', 'a_t'] )
            case(5)
                allocate( taylor_parametrization_1D::self%PureEFTwDE )
                call self%PureEFTwDE%set_param_names( ['EFTw0', 'EFTwa', 'EFTw2', 'EFTw3'], ['w_0', 'w_a', 'w_2', 'w_3'] )
            case default
                write(*,'(a,I3)') 'No model corresponding to EFTwDE =', self%EFTwDE
                write(*,'(a)')    'Please select an appropriate model.'
        end select
        ! allocate Gamma1:
        if ( allocated(self%PureEFTGamma1) ) deallocate(self%PureEFTGamma1)
        select case ( self%PureEFTmodelGamma1 )
            case(0)
                allocate( zero_parametrization_1D::self%PureEFTGamma1 )
            case(1)
                allocate( constant_parametrization_1D::self%PureEFTGamma1 )
            case(2)
                allocate( linear_parametrization_1D::self%PureEFTGamma1 )
            case(3)
                allocate( power_law_parametrization_1D::self%PureEFTGamma1 )
                call self%PureEFTGamma1%set_param_names( ['EFTGamma10  ', 'EFTGamma1Exp'], ['\gamma_0^{(1) {\rm EFT}}        ', '\gamma_{\rm exp}^{(1) {\rm EFT}}'] )
            case(4)
                allocate( exponential_parametrization_1D::self%PureEFTGamma1 )
                call self%PureEFTGamma1%set_param_names( ['EFTGamma10  ', 'EFTGamma1Exp'], ['\gamma_0^{(1) {\rm EFT}}        ', '\gamma_{\rm exp}^{(1) {\rm EFT}}'] )
            case default
                write(*,'(a,I3)') 'No model corresponding to PureEFTmodelGamma1 =', self%PureEFTmodelGamma1
                write(*,'(a)')    'Please select an appropriate model.'
        end select
        ! allocate Gamma2:
        if ( allocated(self%PureEFTGamma2) ) deallocate(self%PureEFTGamma2)
        select case ( self%PureEFTmodelGamma2 )
            case(0)
                allocate( zero_parametrization_1D::self%PureEFTGamma2 )
            case(1)
                allocate( constant_parametrization_1D::self%PureEFTGamma2 )
            case(2)
                allocate( linear_parametrization_1D::self%PureEFTGamma2 )
            case(3)
                allocate( power_law_parametrization_1D::self%PureEFTGamma2 )
                call self%PureEFTGamma2%set_param_names( ['EFTGamma20  ', 'EFTGamma2Exp'], ['\gamma_0^{(2) {\rm EFT}}        ', '\gamma_{\rm exp}^{(2) {\rm EFT}}'] )
            case(4)
                allocate( exponential_parametrization_1D::self%PureEFTGamma2 )
                call self%PureEFTGamma2%set_param_names( ['EFTGamma20  ', 'EFTGamma2Exp'], ['\gamma_0^{(2) {\rm EFT}}        ', '\gamma_{\rm exp}^{(2) {\rm EFT}}'] )
            case default
                write(*,'(a,I3)') 'No model corresponding to PureEFTmodelGamma2 =', self%PureEFTmodelGamma2
                write(*,'(a)')    'Please select an appropriate model.'
        end select
        ! allocate Gamma3:
        if ( allocated(self%PureEFTGamma3) ) deallocate(self%PureEFTGamma3)
        select case ( self%PureEFTmodelGamma3 )
            case(0)
                allocate( zero_parametrization_1D::self%PureEFTGamma3 )
            case(1)
                allocate( constant_parametrization_1D::self%PureEFTGamma3 )
            case(2)
                allocate( linear_parametrization_1D::self%PureEFTGamma3 )
            case(3)
                allocate( power_law_parametrization_1D::self%PureEFTGamma3 )
                call self%PureEFTGamma3%set_param_names( ['EFTGamma30  ', 'EFTGamma3Exp'], ['\gamma_0^{(3) {\rm EFT}}        ', '\gamma_{\rm exp}^{(3) {\rm EFT}}'] )
            case(4)
                allocate( exponential_parametrization_1D::self%PureEFTGamma3 )
                call self%PureEFTGamma3%set_param_names( ['EFTGamma30  ', 'EFTGamma3Exp'], ['\gamma_0^{(3) {\rm EFT}}        ', '\gamma_{\rm exp}^{(3) {\rm EFT}}'] )
            case default
                write(*,'(a,I3)') 'No model corresponding to PureEFTmodelGamma3 =', self%PureEFTmodelGamma3
                write(*,'(a)')    'Please select an appropriate model.'
        end select
        ! allocate the other functions only if not Horndeski:
        if ( .not. self%PureEFTHorndeski ) then
            ! allocate Gamma4:
            if ( allocated(self%PureEFTGamma4) ) deallocate(self%PureEFTGamma4)
            select case ( self%PureEFTmodelGamma4 )
                case(0)
                    allocate( zero_parametrization_1D::self%PureEFTGamma4 )
                case(1)
                    allocate( constant_parametrization_1D::self%PureEFTGamma4 )
                case(2)
                    allocate( linear_parametrization_1D::self%PureEFTGamma4 )
                case(3)
                    allocate( power_law_parametrization_1D::self%PureEFTGamma4 )
                    call self%PureEFTGamma4%set_param_names( ['EFTGamma40  ', 'EFTGamma4Exp'], ['\gamma_0^{(4) {\rm EFT}}        ', '\gamma_{\rm exp}^{(4) {\rm EFT}}'] )
                case(4)
                    allocate( exponential_parametrization_1D::self%PureEFTGamma4 )
                    call self%PureEFTGamma4%set_param_names( ['EFTGamma40  ', 'EFTGamma4Exp'], ['\gamma_0^{(4) {\rm EFT}}        ', '\gamma_{\rm exp}^{(4) {\rm EFT}}'] )
                case default
                    write(*,'(a,I3)') 'No model corresponding to PureEFTmodelGamma4 =', self%PureEFTmodelGamma4
                    write(*,'(a)')    'Please select an appropriate model.'
            end select
            ! allocate Gamma5:
            if ( allocated(self%PureEFTGamma5) ) deallocate(self%PureEFTGamma5)
            select case ( self%PureEFTmodelGamma5 )
                case(0)
                    allocate( zero_parametrization_1D::self%PureEFTGamma5 )
                case(1)
                    allocate( constant_parametrization_1D::self%PureEFTGamma5 )
                case(2)
                    allocate( linear_parametrization_1D::self%PureEFTGamma5 )
                case(3)
                    allocate( power_law_parametrization_1D::self%PureEFTGamma5 )
                    call self%PureEFTGamma5%set_param_names( ['EFTGamma50  ', 'EFTGamma5Exp'], ['\gamma_0^{(5) {\rm EFT}}        ', '\gamma_{\rm exp}^{(5) {\rm EFT}}'] )
                case(4)
                    allocate( exponential_parametrization_1D::self%PureEFTGamma5 )
                    call self%PureEFTGamma5%set_param_names( ['EFTGamma50  ', 'EFTGamma5Exp'], ['\gamma_0^{(5) {\rm EFT}}        ', '\gamma_{\rm exp}^{(5) {\rm EFT}}'] )
                case default
                    write(*,'(a,I3)') 'No model corresponding to PureEFTmodelGamma5 =', self%PureEFTmodelGamma5
                    write(*,'(a)')    'Please select an appropriate model.'
            end select
            ! allocate Gamma6:
            if ( allocated(self%PureEFTGamma6) ) deallocate(self%PureEFTGamma6)
            select case ( self%PureEFTmodelGamma6 )
                case(0)
                    allocate( zero_parametrization_1D::self%PureEFTGamma6 )
                case(1)
                    allocate( constant_parametrization_1D::self%PureEFTGamma6 )
                case(2)
                    allocate( linear_parametrization_1D::self%PureEFTGamma6 )
                case(3)
                    allocate( power_law_parametrization_1D::self%PureEFTGamma6 )
                    call self%PureEFTGamma6%set_param_names( ['EFTGamma60  ', 'EFTGamma6Exp'], ['\gamma_0^{(6) {\rm EFT}}        ', '\gamma_{\rm exp}^{(6) {\rm EFT}}'] )
                case(4)
                    allocate( exponential_parametrization_1D::self%PureEFTGamma6 )
                    call self%PureEFTGamma6%set_param_names( ['EFTGamma60  ', 'EFTGamma6Exp'], ['\gamma_0^{(6) {\rm EFT}}        ', '\gamma_{\rm exp}^{(6) {\rm EFT}}'] )
                case default
                    write(*,'(a,I3)') 'No model corresponding to PureEFTmodelGamma6 =', self%PureEFTmodelGamma6
                    write(*,'(a)')    'Please select an appropriate model.'
            end select
        end if

        ! initialize the names:
        call self%PureEFTOmega%set_name ( 'EFTOmega' , '\Omega'       )
        call self%PureEFTwDE%set_name   ( 'EFTw'     , 'w'            )
        call self%PureEFTGamma1%set_name( 'EFTGamma1', '\gamma^{(1)}' )
        call self%PureEFTGamma2%set_name( 'EFTGamma2', '\gamma^{(2)}' )
        call self%PureEFTGamma3%set_name( 'EFTGamma3', '\gamma^{(3)}' )

        if ( .not. self%PureEFTHorndeski ) then
            call self%PureEFTGamma4%set_name( 'EFTGamma4', '\gamma^{(4)}' )
            call self%PureEFTGamma5%set_name( 'EFTGamma5', '\gamma^{(5)}' )
            call self%PureEFTGamma6%set_name( 'EFTGamma6', '\gamma^{(6)}' )
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
    subroutine EFTCAMBPureEFTstdInitModelParametersFromFile( self, Ini )

        implicit none

        class(EFTCAMB_std_pure_EFT) :: self   !< the base class
        type(TIniFile)              :: Ini    !< Input ini file

        call self%PureEFTOmega%init_from_file ( Ini )
        call self%PureEFTwDE%init_from_file   ( Ini )
        call self%PureEFTGamma1%init_from_file( Ini )
        call self%PureEFTGamma2%init_from_file( Ini )
        call self%PureEFTGamma3%init_from_file( Ini )

        if ( .not. self%PureEFTHorndeski ) then
            call self%PureEFTGamma4%init_from_file( Ini )
            call self%PureEFTGamma5%init_from_file( Ini )
            call self%PureEFTGamma6%init_from_file( Ini )
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
        if ( .not. self%PureEFTHorndeski ) then
            self%parameter_number = self%parameter_number +self%PureEFTGamma4%parameter_number
            self%parameter_number = self%parameter_number +self%PureEFTGamma5%parameter_number
            self%parameter_number = self%parameter_number +self%PureEFTGamma6%parameter_number
        end if

    end subroutine EFTCAMBPureEFTstdComputeParametersNumber

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that prints on the screen feedback information about the model.
    subroutine EFTCAMBPureEFTstdFeedback( self )

        implicit none

        class(EFTCAMB_std_pure_EFT)  :: self   !< the base class

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
        if ( .not. self%PureEFTHorndeski ) then
            if ( self%PureEFTmodelGamma4 /= 0 ) write(*,'(a,I3)') '   PureEFTmodelGamma4  =', self%PureEFTmodelGamma4
            if ( self%PureEFTmodelGamma5 /= 0 ) write(*,'(a,I3)') '   PureEFTmodelGamma5  =', self%PureEFTmodelGamma5
            if ( self%PureEFTmodelGamma6 /= 0 ) write(*,'(a,I3)') '   PureEFTmodelGamma6  =', self%PureEFTmodelGamma6
        end if

        write(*,*)
        ! print functions informations:
        call self%PureEFTOmega%feedback()
        call self%PureEFTwDE%feedback()
        call self%PureEFTGamma1%feedback()
        call self%PureEFTGamma2%feedback()
        call self%PureEFTGamma3%feedback()
        if ( .not. self%PureEFTHorndeski ) then
            call self%PureEFTGamma4%feedback()
            call self%PureEFTGamma5%feedback()
            call self%PureEFTGamma6%feedback()
        end if

    end subroutine EFTCAMBPureEFTstdFeedback

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the i-th parameter name of the model
    subroutine EFTCAMBPureEFTstdParameterNames( self, i, name )

        implicit none

        class(EFTCAMB_std_pure_EFT) :: self   !< the base class
        integer     , intent(in)    :: i      !< the index of the parameter
        character(*), intent(out)   :: name   !< the output name of the i-th parameter

        integer  :: NOmega, Nw, N1, N2, N3, N4, N5, N6
        integer  :: j

        ! compute the incremental number of parameters:
        NOmega = self%PureEFTOmega%parameter_number
        Nw     = NOmega + self%PureEFTwDE%parameter_number
        N1     = Nw + self%PureEFTGamma1%parameter_number
        N2     = N1 + self%PureEFTGamma2%parameter_number
        N3     = N2 + self%PureEFTGamma3%parameter_number
        if ( .not. self%PureEFTHorndeski ) then
            N4 = N3 + self%PureEFTGamma4%parameter_number
            N4 = N4 + self%PureEFTGamma5%parameter_number
            N6 = N5 + self%PureEFTGamma6%parameter_number
        end if

        ! check validity of input:
        if ( i > self%parameter_number .or. i <= 0 ) then
            write(*,'(a,I3)') 'No parameter corresponding to: ', i
            write(*,'(a,I3)') 'Total number of parameters is: ', self%parameter_number
            call MpiStop('EFTCAMB error')

        ! parameter from Omega function
        else if ( i <= NOmega ) then
            do j = 1, self%PureEFTOmega%parameter_number
                if ( i == j ) call self%PureEFTOmega%parameter_names( j, name )
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

        !parameter from Gamma4 function
        else if ( .not. self%PureEFTHorndeski .and. i <= N4 ) then

            do j = 1, self%PureEFTGamma4%parameter_number
                if ( i-N3 == j ) call self%PureEFTGamma4%parameter_names( j, name )
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

        end if

    end subroutine EFTCAMBPureEFTstdParameterNames

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the i-th parameter name of the model
    subroutine EFTCAMBPureEFTstdParameterNamesLatex( self, i, latexname )

        implicit none

        class(EFTCAMB_std_pure_EFT) :: self       !< the base class
        integer     , intent(in)    :: i          !< The index of the parameter
        character(*), intent(out)   :: latexname  !< the output latex name of the i-th parameter

        integer  :: NOmega, Nw, N1, N2, N3, N4, N5, N6
        integer  :: j

        ! compute the incremental number of parameters:
        NOmega = self%PureEFTOmega%parameter_number
        Nw = NOmega + self%PureEFTwDE%parameter_number
        N1 = Nw + self%PureEFTGamma1%parameter_number
        N2 = N1 + self%PureEFTGamma2%parameter_number
        N3 = N2 + self%PureEFTGamma3%parameter_number
        if ( .not. self%PureEFTHorndeski ) then
            N4 = N3 + self%PureEFTGamma4%parameter_number
            N4 = N4 + self%PureEFTGamma5%parameter_number
            N6 = N5 + self%PureEFTGamma6%parameter_number
        end if

        ! check validity of input:
        if ( i > self%parameter_number .or. i <= 0) then
            write(*,'(a,I3)') 'No parameter corresponding to: ', i
            write(*,'(a,I3)') 'Total number of parameters is: ', self%parameter_number
            call MpiStop('EFTCAMB error')

        ! parameter from Omega function
        else if ( i <= NOmega ) then
            do j = 1, self%PureEFTOmega%parameter_number
                if ( i == j ) call self%PureEFTOmega%parameter_names_latex( j, latexname )
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

            !parameter from Gamma4 function
        else if ( .not. self%PureEFTHorndeski .and. i <= N4 ) then
            do j = 1, self%PureEFTGamma4%parameter_number
                if ( i-N3 == j ) call self%PureEFTGamma4%parameter_names_latex( j, latexname )
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

        end if

    end subroutine EFTCAMBPureEFTstdParameterNamesLatex

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the i-th parameter name of the model
    subroutine EFTCAMBPureEFTstdParameterValues( self, i, value )

        implicit none

        class(EFTCAMB_std_pure_EFT) :: self   !< the base class
        integer , intent(in)        :: i      !< The index of the parameter
        real(dl), intent(out)       :: value  !< the output value of the i-th parameter

        integer  :: NOmega, Nw, N1, N2, N3, N4, N5, N6
        integer  :: j

        ! compute the incremental number of parameters:
        NOmega = self%PureEFTOmega%parameter_number
        Nw = NOmega + self%PureEFTwDE%parameter_number
        N1 = Nw + self%PureEFTGamma1%parameter_number
        N2 = N1 + self%PureEFTGamma2%parameter_number
        N3 = N2 + self%PureEFTGamma3%parameter_number
        if ( .not. self%PureEFTHorndeski ) then
            N4 = N3 + self%PureEFTGamma4%parameter_number
            N4 = N4 + self%PureEFTGamma5%parameter_number
            N6 = N5 + self%PureEFTGamma6%parameter_number
        end if

        ! check validity of input:
        if ( i > self%parameter_number .or. i <= 0) then
            write(*,'(a,I3)') 'EFTCAMB error: no parameter corresponding to number ',i
            write(*,'(a,I3)') 'Total number of parameters is ', self%parameter_number
            call MpiStop('EFTCAMB error')

        ! parameter from Omega function
        else if ( i <= NOmega ) then
            do j = 1, self%PureEFTOmega%parameter_number
                if ( i == j ) call self%PureEFTOmega%parameter_value( j, value )
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

        !parameter from Gamma4 function
        else if ( .not. self%PureEFTHorndeski .and. i <= N4 ) then
            do j = 1, self%PureEFTGamma4%parameter_number
                if ( i-N3 == j ) call self%PureEFTGamma4%parameter_value( j, value )
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

        end if

    end subroutine EFTCAMBPureEFTstdParameterValues

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes the value of the background EFT functions at a given time.
    subroutine EFTCAMBPureEFTstdBackgroundEFTFunctions( self, a, eft_par_cache, eft_cache )

        implicit none

        class(EFTCAMB_std_pure_EFT)                  :: self          !< the base class
        real(dl), intent(in)                         :: a             !< the input scale factor
        type(EFTCAMB_parameter_cache), intent(inout) :: eft_par_cache !< the EFTCAMB parameter cache that contains all the physical parameters.
        type(EFTCAMB_timestep_cache ), intent(inout) :: eft_cache     !< the EFTCAMB timestep cache that contains all the physical values.

        eft_cache%EFTOmegaV    = self%PureEFTOmega%value(a)
        eft_cache%EFTOmegaP    = self%PureEFTOmega%first_derivative(a)
        eft_cache%EFTOmegaPP   = self%PureEFTOmega%second_derivative(a)
        eft_cache%EFTOmegaPPP  = self%PureEFTOmega%third_derivative(a)
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
        eft_cache%EFTLambdadot = -2._dl*eft_cache%EFTOmegaV*( eft_cache%Hdotdot -eft_cache%adotoa*eft_cache%Hdot -eft_cache%adotoa**3 ) &
            & -a*eft_cache%EFTOmegaP*( 4._dl*eft_cache%adotoa*eft_cache%Hdot +eft_cache%Hdotdot ) &
            & -a**2*eft_cache%EFTOmegaPP*eft_cache%adotoa*( +2._dl*eft_cache%adotoa**2 +3._dl*eft_cache%Hdot )&
            & -(a*eft_cache%adotoa)**3*eft_cache%EFTOmegaPPP &
            & +eft_cache%grhov_t*eft_cache%adotoa*( a*self%PureEFTwDE%first_derivative(a) -3._dl*self%PureEFTwDE%value(a)*(1._dl +self%PureEFTwDE%value(a) ))

    end subroutine EFTCAMBPureEFTstdBackgroundEFTFunctions

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes the value of the background EFT functions at a given time.
    subroutine EFTCAMBPureEFTstdSecondOrderEFTFunctions( self, a, eft_par_cache, eft_cache )

        implicit none

        class(EFTCAMB_std_pure_EFT)                  :: self          !< the base class
        real(dl), intent(in)                         :: a             !< the input scale factor
        type(EFTCAMB_parameter_cache), intent(inout) :: eft_par_cache !< the EFTCAMB parameter cache that contains all the physical parameters.
        type(EFTCAMB_timestep_cache ), intent(inout) :: eft_cache     !< the EFTCAMB timestep cache that contains all the physical values.

        eft_cache%EFTGamma1V  = self%PureEFTGamma1%value(a)
        eft_cache%EFTGamma1P  = self%PureEFTGamma1%first_derivative(a)
        eft_cache%EFTGamma2V  = self%PureEFTGamma2%value(a)
        eft_cache%EFTGamma2P  = self%PureEFTGamma2%first_derivative(a)
        eft_cache%EFTGamma3V  = self%PureEFTGamma3%value(a)
        eft_cache%EFTGamma3P  = self%PureEFTGamma3%first_derivative(a)
        if ( self%PureEFTHorndeski ) then
            eft_cache%EFTGamma4V  = -eft_cache%EFTGamma3V
            eft_cache%EFTGamma4P  = -eft_cache%EFTGamma3P
            eft_cache%EFTGamma4PP = -self%PureEFTGamma3%second_derivative(a)
            eft_cache%EFTGamma5V  = +0.5_dl*eft_cache%EFTGamma3V
            eft_cache%EFTGamma5P  = +0.5_dl*eft_cache%EFTGamma3P
            eft_cache%EFTGamma6V  = 0._dl
            eft_cache%EFTGamma6P  = 0._dl
        else
            eft_cache%EFTGamma4V  = self%PureEFTGamma4%value(a)
            eft_cache%EFTGamma4P  = self%PureEFTGamma4%first_derivative(a)
            eft_cache%EFTGamma4PP = self%PureEFTGamma4%second_derivative(a)
            eft_cache%EFTGamma5V  = self%PureEFTGamma5%value(a)
            eft_cache%EFTGamma5P  = self%PureEFTGamma5%first_derivative(a)
            eft_cache%EFTGamma6V  = self%PureEFTGamma6%value(a)
            eft_cache%EFTGamma6P  = self%PureEFTGamma6%first_derivative(a)
        end if

    end subroutine EFTCAMBPureEFTstdSecondOrderEFTFunctions

    ! ---------------------------------------------------------------------------------------------
    !> Function that computes dtauda = 1/sqrt(a^2H^2).For pure EFT std this has to be overridden
    !! for performance reasons.
    function EFTCAMBPureEFTstdComputeDtauda( self, a, eft_par_cache, eft_cache )

        implicit none

        class(EFTCAMB_std_pure_EFT)                  :: self          !< the base class
        real(dl), intent(in)                         :: a             !< the input scale factor
        type(EFTCAMB_parameter_cache), intent(inout) :: eft_par_cache !< the EFTCAMB parameter cache that contains all the physical parameters.
        type(EFTCAMB_timestep_cache ), intent(inout) :: eft_cache     !< the EFTCAMB timestep cache that contains all the physical values.

        real(dl) :: EFTCAMBPureEFTstdComputeDtauda               !< the output dtauda

        real(dl) :: temp

        temp = eft_cache%grhoa2 +eft_par_cache%grhov*a*a*self%PureEFTwDE%integral(a)
        EFTCAMBPureEFTstdComputeDtauda = sqrt(3/temp)

    end function EFTCAMBPureEFTstdComputeDtauda

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes adotoa = H and its two derivatives wrt conformal time.
    subroutine EFTCAMBPureEFTstdComputeAdotoa( self, a, eft_par_cache, eft_cache )

        implicit none

        class(EFTCAMB_std_pure_EFT)                  :: self          !< the base class
        real(dl), intent(in)                         :: a             !< the input scale factor
        type(EFTCAMB_parameter_cache), intent(inout) :: eft_par_cache !< the EFTCAMB parameter cache that contains all the physical parameters.
        type(EFTCAMB_timestep_cache ), intent(inout) :: eft_cache     !< the EFTCAMB timestep cache that contains all the physical values.

        eft_cache%grhov_t = eft_par_cache%grhov*self%PureEFTwDE%integral(a)
        eft_cache%adotoa  = sqrt( ( eft_cache%grhom_t +eft_cache%grhov_t )/3._dl )

    end subroutine EFTCAMBPureEFTstdComputeAdotoa

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes the two derivatives wrt conformal time of H.
    subroutine EFTCAMBPureEFTstdComputeHubbleDer( self, a, eft_par_cache, eft_cache )

        implicit none

        class(EFTCAMB_std_pure_EFT)                  :: self          !< the base class
        real(dl), intent(in)                         :: a             !< the input scale factor
        type(EFTCAMB_parameter_cache), intent(inout) :: eft_par_cache !< the EFTCAMB parameter cache that contains all the physical parameters.
        type(EFTCAMB_timestep_cache ), intent(inout) :: eft_cache     !< the EFTCAMB timestep cache that contains all the physical values.

        eft_cache%gpiv_t  = self%PureEFTwDE%value(a)*eft_cache%grhov_t
        eft_cache%Hdot    = -0.5_dl*( eft_cache%adotoa**2 +eft_cache%gpresm_t +eft_cache%gpiv_t )

        ! IW POSSIBLE BUG
        !eft_cache%Hdotdot = eft_cache%adotoa*( ( eft_cache%grhob_t +eft_cache%grhoc_t)/6._dl +2._dl*( eft_cache%grhor_t +eft_cache%grhog_t)/3._dl ) &
        !    & +eft_cache%adotoa*eft_cache%grhov_t*( 1._dl/6._dl +self%PureEFTwDE%value(a) +1.5_dl*self%PureEFTwDE%value(a)**2 -0.5_dl*a*self%PureEFTwDE%first_derivative(a) ) &
        !    & +eft_cache%adotoa*eft_cache%grhonu_tot/6._dl -0.5_dl*eft_cache%adotoa*eft_cache%gpinu_tot -0.5_dl*eft_cache%gpinudot_tot

        eft_cache%Hdotdot = 2._dl*eft_cache%adotoa*eft_cache%Hdot &
            & + 0.5_dl*eft_cache%adotoa*(eft_cache%grhob_t + eft_cache%grhoc_t + 8._dl*(eft_cache%grhog_t+eft_cache%grhor_t)/3._dl)&
            & + 0.5_dl*eft_cache%adotoa*eft_cache%grhov_t*( (1._dl+self%PureEFTwDE%value(a) )*(1._dl+3._dl*self%PureEFTwDE%value(a)) -a*self%PureEFTwDE%first_derivative(a))&
            & + eft_cache%adotoa/6._dl*eft_cache%grhonu_tot -0.5_dl*eft_cache%adotoa*eft_cache%gpinu_tot -0.5_dl*eft_cache%gpinudot_tot

    end subroutine EFTCAMBPureEFTstdComputeHubbleDer

    ! ---------------------------------------------------------------------------------------------
    !> Function that computes model specific stability requirements.
    function EFTCAMBPureEFTstdAdditionalModelStability( self, a, eft_par_cache, eft_cache )

        implicit none

        class(EFTCAMB_std_pure_EFT)                  :: self          !< the base class
        real(dl), intent(in)                         :: a             !< the input scale factor.
        type(EFTCAMB_parameter_cache), intent(inout) :: eft_par_cache !< the EFTCAMB parameter cache that contains all the physical parameters.
        type(EFTCAMB_timestep_cache ), intent(inout) :: eft_cache     !< the EFTCAMB timestep cache that contains all the physical values.

        logical :: EFTCAMBPureEFTstdAdditionalModelStability          !< the return value of the stability computation. True if the model specific stability criteria are met, false otherwise.

        EFTCAMBPureEFTstdAdditionalModelStability = .True.
        if ( self%PureEFTwDE%value(a) > -1._dl/3._dl ) EFTCAMBPureEFTstdAdditionalModelStability = .False.

    end function EFTCAMBPureEFTstdAdditionalModelStability

    ! ---------------------------------------------------------------------------------------------

end module EFTCAMB_pure_EFT_std

!----------------------------------------------------------------------------------------
