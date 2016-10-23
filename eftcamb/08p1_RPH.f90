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

!> @file 08p1_RPH.f90
!! This file contains the definition of an EFT reparametrization of Horndeski models.
!! This is described by four functions of time and w_DE.
!! Please refer to the numerical notes for details.

!----------------------------------------------------------------------------------------
!> This module contains the definition of an EFT reparametrization of Horndeski models.
!! This is described by four functions of time and w_DE.
!! Please refer to the numerical notes for details.

!> @author Bin Hu, Marco Raveri

module EFTCAMB_Reparametrized_Horndeski

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

    public EFTCAMB_RPH

    !----------------------------------------------------------------------------------------
    !> This is the type that contains the definition of RPH.
    type, extends ( EFTCAMB_designer_model ) :: EFTCAMB_RPH

        ! the RPH functions model selection flags:
        integer  :: RPHwDE              !< Model selection flag for RPH w DE.
        integer  :: RPHmassPmodel       !< Model selection flag for RPH Planck mass model.
        integer  :: RPHkineticitymodel  !< Model selection flag for RPH kineticity model.
        integer  :: RPHbraidingmodel    !< Model selection flag for RPH braiding model.
        integer  :: RPHtensormodel      !< Model selection flag for RPH tensor model.

        ! the RPH functions:
        class( parametrized_function_1D ), allocatable :: RPH_wDE         !< The RPH function w_DE.
        class( parametrized_function_1D ), allocatable :: RPH_PlanckMass  !< The RPH function Planck Mass.
        class( parametrized_function_1D ), allocatable :: RPH_Kineticity  !< The RPH function Kineticity.
        class( parametrized_function_1D ), allocatable :: RPH_Braiding    !< The RPH function Braiding.
        class( parametrized_function_1D ), allocatable :: RPH_Tensor      !< The RPH function Tensor.

    contains

        ! initialization of the model:
        procedure :: read_model_selection            => EFTCAMBRPHReadModelSelectionFromFile  !< subroutine that reads the parameters of the model from file
        procedure :: allocate_model_selection        => EFTCAMBRPHAllocateModelSelection      !< subroutine that allocates the model selection.
        procedure :: init_model_parameters           => EFTCAMBRPHInitModelParameters         !< subroutine that initializes the model parameters based on the values found in an input array.
        procedure :: init_model_parameters_from_file => EFTCAMBRPHInitModelParametersFromFile !< subroutine that reads the parameters of the model from file.
        ! utility functions:
        procedure :: compute_param_number  => EFTCAMBRPHComputeParametersNumber    !< subroutine that computes the number of parameters of the model.
        procedure :: feedback              => EFTCAMBRPHFeedback                   !< subroutine that prints on the screen feedback information about the model.
        procedure :: parameter_names       => EFTCAMBRPHParameterNames             !< subroutine that returns the i-th parameter name of the model.
        procedure :: parameter_names_latex => EFTCAMBRPHParameterNamesLatex        !< subroutine that returns the i-th parameter name of the model.
        procedure :: parameter_values      => EFTCAMBRPHParameterValues            !< subroutine that returns the i-th parameter value.
        ! CAMB related procedures:
        procedure :: compute_background_EFT_functions  => EFTCAMBRPHBackgroundEFTFunctions   !< subroutine that computes the value of the background EFT functions at a given time.
        procedure :: compute_secondorder_EFT_functions => EFTCAMBRPHSecondOrderEFTFunctions  !< subroutine that computes the value of the second order EFT functions at a given time.
        procedure :: compute_dtauda                    => EFTCAMBRPHComputeDtauda            !< function that computes dtauda = 1/sqrt(a^2H^2).
        procedure :: compute_adotoa                    => EFTCAMBRPHComputeAdotoa            !< subroutine that computes adotoa = H and its two derivatives wrt conformal time.
        procedure :: compute_H_derivs                  => EFTCAMBRPHComputeHubbleDer         !< subroutine that computes the two derivatives wrt conformal time of H.
        ! stability procedures:
        procedure :: additional_model_stability        => EFTCAMBRPHAdditionalModelStability !< function that computes model specific stability requirements.

    end type EFTCAMB_RPH

    ! ---------------------------------------------------------------------------------------------

contains

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that reads the parameters of the model from file.
    subroutine EFTCAMBRPHReadModelSelectionFromFile( self, Ini )

        implicit none

        class(EFTCAMB_RPH) :: self   !< the base class
        type(TIniFile)              :: Ini    !< Input ini file

        ! read model selection flags:
        self%RPHwDE             = Ini_Read_Int_File( Ini, 'RPHwDE'            , 0 )
        self%RPHmassPmodel      = Ini_Read_Int_File( Ini, 'RPHmassPmodel'     , 0 )
        self%RPHkineticitymodel = Ini_Read_Int_File( Ini, 'RPHkineticitymodel', 0 )
        self%RPHbraidingmodel   = Ini_Read_Int_File( Ini, 'RPHbraidingmodel'  , 0 )
        self%RPHtensormodel     = Ini_Read_Int_File( Ini, 'RPHtensormodel'    , 0 )

    end subroutine EFTCAMBRPHReadModelSelectionFromFile

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that allocates the model selection.
    subroutine EFTCAMBRPHAllocateModelSelection( self )

        implicit none

        class(EFTCAMB_RPH)                       :: self              !< the base class

        ! allocate wDE:
        if ( allocated(self%RPH_wDE) ) deallocate(self%RPH_wDE)
        select case ( self%RPHwDE )
            case(0)
                allocate( wDE_LCDM_parametrization_1D::self%RPH_wDE )
            case(1)
                allocate( constant_parametrization_1D::self%RPH_wDE )
            case(2)
                allocate( CPL_parametrization_1D::self%RPH_wDE )
                call self%RPH_wDE%set_param_names( ['EFTw0', 'EFTwa'], ['w_0', 'w_a'] )
            case(3)
                allocate( JBP_parametrization_1D::self%RPH_wDE )
                call self%RPH_wDE%set_param_names( ['EFTw0', 'EFTwa', 'EFTwn'], [ 'w_0', 'w_a', 'n  ' ] )
            case(4)
                allocate( turning_point_parametrization_1D::self%RPH_wDE )
                call self%RPH_wDE%set_param_names( ['EFTw0 ', 'EFTwa ', 'EFTwat'], ['w_0', 'w_a', 'a_t'] )
            case(5)
                allocate( taylor_parametrization_1D::self%RPH_wDE )
                call self%RPH_wDE%set_param_names( ['EFTw0', 'EFTwa', 'EFTw2', 'EFTw3'], ['w_0', 'w_a', 'w_2', 'w_3'] )
            case default
                write(*,'(a,I3)') 'No model corresponding to RPH_wDE =', self%RPHwDE
                write(*,'(a)')    'Please select an appropriate model.'
        end select
        ! allocate RPH_PlanckMass:
        if ( allocated(self%RPH_PlanckMass) ) deallocate(self%RPH_PlanckMass)
        select case ( self%RPHmassPmodel )
            case(0)
                allocate( zero_parametrization_1D::self%RPH_PlanckMass )
            case(1)
                allocate( constant_parametrization_1D::self%RPH_PlanckMass )
            case(2)
                allocate( linear_parametrization_1D::self%RPH_PlanckMass )
            case(3)
                allocate( power_law_parametrization_1D::self%RPH_PlanckMass )
                call self%RPH_PlanckMass%set_param_names( ['RPHmassP0  ', 'RPHmassPExp'], ['\tilde{M}_0                  ', '\tilde{M}_{\rm exp}^{\rm RPH}'] )
            case(4)
                allocate( exponential_parametrization_1D::self%RPH_PlanckMass )
                call self%RPH_PlanckMass%set_param_names( ['RPHmassP0  ', 'RPHmassPExp'], ['\tilde{M}_0                  ', '\tilde{M}_{\rm exp}^{\rm RPH}'] )
            case default
                write(*,'(a,I3)') 'No model corresponding to RPH_PlanckMass =', self%RPHmassPmodel
                write(*,'(a)')    'Please select an appropriate model.'
        end select
        ! allocate RPH_Kineticity:
        if ( allocated(self%RPH_Kineticity) ) deallocate(self%RPH_Kineticity)
        select case ( self%RPHkineticitymodel )
            case(0)
                allocate( zero_parametrization_1D::self%RPH_Kineticity )
            case(1)
                allocate( constant_parametrization_1D::self%RPH_Kineticity )
            case(2)
                allocate( linear_parametrization_1D::self%RPH_Kineticity )
            case(3)
                allocate( power_law_parametrization_1D::self%RPH_Kineticity )
                call self%RPH_Kineticity%set_param_names( ['RPHkineticity0  ', 'RPHkineticityExp'], ['\alpha_0^{\rm K}        ', '\alpha_{\rm exp}^{\rm K}'] )
            case(4)
                allocate( exponential_parametrization_1D::self%RPH_Kineticity )
                call self%RPH_Kineticity%set_param_names( ['RPHkineticity0  ', 'RPHkineticityExp'], ['\alpha_0^{\rm K}        ', '\alpha_{\rm exp}^{\rm K}'] )
            case default
                write(*,'(a,I3)') 'No model corresponding to RPH_Kineticity =', self%RPHkineticitymodel
                write(*,'(a)')    'Please select an appropriate model.'
        end select
        ! allocate RPH_Braiding:
        if ( allocated(self%RPH_Braiding) ) deallocate(self%RPH_Braiding)
        select case ( self%RPHbraidingmodel )
            case(0)
                allocate( zero_parametrization_1D::self%RPH_Braiding )
            case(1)
                allocate( constant_parametrization_1D::self%RPH_Braiding )
            case(2)
                allocate( linear_parametrization_1D::self%RPH_Braiding )
            case(3)
                allocate( power_law_parametrization_1D::self%RPH_Braiding )
                call self%RPH_Braiding%set_param_names( ['RPHbraiding0  ', 'RPHbraidingExp'], ['\alpha_0^{\rm B}        ', '\alpha_{\rm exp}^{\rm B}'] )
            case(4)
                allocate( exponential_parametrization_1D::self%RPH_Braiding )
                call self%RPH_Braiding%set_param_names( ['RPHbraiding0  ', 'RPHbraidingExp'], ['\alpha_0^{\rm B}        ', '\alpha_{\rm exp}^{\rm B}'] )
            case default
                write(*,'(a,I3)') 'No model corresponding to RPH_Braiding =', self%RPHbraidingmodel
                write(*,'(a)')    'Please select an appropriate model.'
        end select
        ! allocate RPH_Braiding:
        if ( allocated(self%RPH_Tensor) ) deallocate(self%RPH_Tensor)
        select case ( self%RPHtensormodel )
            case(0)
                allocate( zero_parametrization_1D::self%RPH_Tensor )
            case(1)
                allocate( constant_parametrization_1D::self%RPH_Tensor )
            case(2)
                allocate( linear_parametrization_1D::self%RPH_Tensor )
            case(3)
                allocate( power_law_parametrization_1D::self%RPH_Tensor )
                call self%RPH_Tensor%set_param_names( ['RPHtensor0  ', 'RPHtensorExp'], ['\alpha_0^{\rm T}        ', '\alpha_{\rm exp}^{\rm T}'] )
            case(4)
                allocate( exponential_parametrization_1D::self%RPH_Tensor )
                call self%RPH_Tensor%set_param_names( ['RPHtensor0  ', 'RPHtensorExp'], ['\alpha_0^{\rm T}        ', '\alpha_{\rm exp}^{\rm T}'] )
            case default
                write(*,'(a,I3)') 'No model corresponding to RPH_Tensor =', self%RPHtensormodel
                write(*,'(a)')    'Please select an appropriate model.'
        end select

        ! initialize the names:
        call self%RPH_wDE%set_name        ( 'RPHw'         , 'w'              )
        call self%RPH_PlanckMass%set_name ( 'RPHmassP'     , '\tilde{M}'      )
        call self%RPH_Kineticity%set_name ( 'RPHkineticity', '\alpha^{\rm K}' )
        call self%RPH_Braiding%set_name   ( 'RPHbraiding'  , '\alpha^{\rm B}' )
        call self%RPH_Tensor%set_name     ( 'RPHtensor'    , '\alpha^{\rm T}' )

    end subroutine EFTCAMBRPHAllocateModelSelection

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that initializes the model parameters based on the values found in an input array.
    subroutine EFTCAMBRPHInitModelParameters( self, array )

        implicit none

        class(EFTCAMB_RPH)                            :: self   !< the base class
        real(dl), dimension(self%parameter_number), intent(in) :: array  !< input array with the values of the parameters of the model.

        real(dl), allocatable, dimension(:)                    :: temp
        integer :: num_params_function, num_params_temp, i

        num_params_temp     = 1

        ! first elements are w_DE parameters:
        num_params_function = self%RPH_wDE%parameter_number
        allocate( temp(num_params_function) )
        do i = 1, num_params_function
            temp(i)         = array(num_params_temp)
            num_params_temp = num_params_temp +1
        end do
        call self%RPH_wDE%init_parameters(temp)
        deallocate( temp )
        ! then RPH_PlanckMass parameters:
        num_params_function = self%RPH_PlanckMass%parameter_number
        allocate( temp(num_params_function) )
        do i = 1, num_params_function
            temp(i)         = array(num_params_temp)
            num_params_temp = num_params_temp +1
        end do
        call self%RPH_PlanckMass%init_parameters(temp)
        deallocate(temp)
        ! then RPH_Kineticity parameters:
        num_params_function = self%RPH_Kineticity%parameter_number
        allocate( temp(num_params_function) )
        do i = 1, num_params_function
            temp(i)         = array(num_params_temp)
            num_params_temp = num_params_temp +1
        end do
        call self%RPH_Kineticity%init_parameters(temp)
        deallocate(temp)
        ! then RPH_Braiding parameters:
        num_params_function = self%RPH_Braiding%parameter_number
        allocate( temp(num_params_function) )
        do i = 1, num_params_function
            temp(i)         = array(num_params_temp)
            num_params_temp = num_params_temp +1
        end do
        call self%RPH_Braiding%init_parameters(temp)
        deallocate(temp)
        ! then RPH_Tensor parameters:
        num_params_function = self%RPH_Tensor%parameter_number
        allocate( temp(num_params_function) )
        do i = 1, num_params_function
            temp(i)         = array(num_params_temp)
            num_params_temp = num_params_temp +1
        end do
        call self%RPH_Tensor%init_parameters(temp)
        deallocate(temp)

        ! now check the length of the parameters:
        if ( num_params_temp-1 /= self%parameter_number ) then
            write(*,*) 'In EFTCAMBRPHInitModelParameters:'
            write(*,*) 'Length of num_params_temp and self%parameter_number do not coincide.'
            write(*,*) 'num_params_temp:', num_params_temp-1
            write(*,*) 'self%parameter_number:', self%parameter_number
            call MpiStop('EFTCAMB error')
        end if

    end subroutine EFTCAMBRPHInitModelParameters

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that reads the parameters of the model from file.
    subroutine EFTCAMBRPHInitModelParametersFromFile( self, Ini )

        implicit none

        class(EFTCAMB_RPH)          :: self   !< the base class
        type(TIniFile)              :: Ini    !< Input ini file

        call self%RPH_wDE%init_from_file( Ini )
        call self%RPH_PlanckMass%init_from_file( Ini )
        call self%RPH_Kineticity%init_from_file( Ini )
        call self%RPH_Braiding%init_from_file( Ini )
        call self%RPH_Tensor%init_from_file( Ini )

    end subroutine EFTCAMBRPHInitModelParametersFromFile

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes the number of parameters of the model.
    subroutine EFTCAMBRPHComputeParametersNumber( self )

        implicit none

        class(EFTCAMB_RPH)  :: self   !< the base class

        self%parameter_number = 0
        self%parameter_number = self%parameter_number +self%RPH_wDE%parameter_number
        self%parameter_number = self%parameter_number +self%RPH_PlanckMass%parameter_number
        self%parameter_number = self%parameter_number +self%RPH_Kineticity%parameter_number
        self%parameter_number = self%parameter_number +self%RPH_Braiding%parameter_number
        self%parameter_number = self%parameter_number +self%RPH_Tensor%parameter_number

    end subroutine EFTCAMBRPHComputeParametersNumber

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that prints on the screen feedback information about the model.
    subroutine EFTCAMBRPHFeedback( self )

        implicit none

        class(EFTCAMB_RPH)  :: self   !< the base class

        ! print general model informations:
        write(*,*)
        write(*,'(a,a)')    '   Model               =  ', self%name
        write(*,'(a,I3)')   '   Number of params    ='  , self%parameter_number

        ! print model functions informations:
        write(*,*)
        if ( self%RPHwDE             /= 0 ) write(*,'(a,I3)') '   RPHwDE             =', self%RPHwDE
        if ( self%RPHmassPmodel      /= 0 ) write(*,'(a,I3)') '   RPHmassPmodel      =', self%RPHmassPmodel
        if ( self%RPHkineticitymodel /= 0 ) write(*,'(a,I3)') '   RPHkineticitymodel =', self%RPHkineticitymodel
        if ( self%RPHbraidingmodel   /= 0 ) write(*,'(a,I3)') '   RPHbraidingmodel   =', self%RPHbraidingmodel
        if ( self%RPHtensormodel     /= 0 ) write(*,'(a,I3)') '   RPHtensormodel     =', self%RPHtensormodel

        write(*,*)
        ! print functions informations:
        call self%RPH_wDE%feedback()
        call self%RPH_PlanckMass%feedback()
        call self%RPH_Kineticity%feedback()
        call self%RPH_Braiding%feedback()
        call self%RPH_Tensor%feedback()

    end subroutine EFTCAMBRPHFeedback

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the i-th parameter name of the model
    subroutine EFTCAMBRPHParameterNames( self, i, name )

        implicit none

        class(EFTCAMB_RPH) :: self   !< the base class
        integer     , intent(in)    :: i      !< the index of the parameter
        character(*), intent(out)   :: name   !< the output name of the i-th parameter

        integer  :: Nw, NM, NK, NB, NT
        integer  :: j

        ! compute the incremental number of parameters:
        Nw = self%RPH_wDE%parameter_number
        NM = Nw + self%RPH_PlanckMass%parameter_number
        NK = NM + self%RPH_Kineticity%parameter_number
        NB = NK + self%RPH_Braiding%parameter_number
        NT = NB + self%RPH_Tensor%parameter_number

        ! check validity of input:
        if ( i > self%parameter_number .or. i <= 0 ) then
            write(*,'(a,I3)') 'No parameter corresponding to: ', i
            write(*,'(a,I3)') 'Total number of parameters is: ', self%parameter_number
            call MpiStop('EFTCAMB error')

        ! parameter from wDE function
        else if ( i <= Nw ) then
            do j = 1, self%RPH_wDE%parameter_number
                if ( i == j ) call self%RPH_wDE%parameter_names( j, name )
            end do
            return

        ! parameter from RPH_PlanckMass function
        else if ( i <= NM ) then
            do j = 1, self%RPH_PlanckMass%parameter_number
                if ( i-Nw == j ) call self%RPH_PlanckMass%parameter_names( j, name )
            end do
            return

        ! parameter from RPH_Kineticity function
        else if ( i <= NK) then
            do j = 1, self%RPH_Kineticity%parameter_number
                if ( i-NM == j ) call self%RPH_Kineticity%parameter_names( j, name )
            end do
            return

        !parameter from RPH_Braiding function
        else if ( i <= NB ) then
            do j = 1, self%RPH_Braiding%parameter_number
                if ( i-NK == j ) call self%RPH_Braiding%parameter_names( j, name )
            end do
            return

        !parameter from RPH_Tensor function
        else if ( i <= NT ) then

            do j = 1, self%RPH_Tensor%parameter_number
                if ( i-NB == j ) call self%RPH_Tensor%parameter_names( j, name )
            end do
            return

        end if

    end subroutine EFTCAMBRPHParameterNames

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the i-th parameter name of the model
    subroutine EFTCAMBRPHParameterNamesLatex( self, i, latexname )

        implicit none

        class(EFTCAMB_RPH) :: self       !< the base class
        integer     , intent(in)    :: i          !< The index of the parameter
        character(*), intent(out)   :: latexname  !< the output latex name of the i-th parameter

        integer  :: Nw, NM, NK, NB, NT
        integer  :: j

        ! compute the incremental number of parameters:
        Nw = self%RPH_wDE%parameter_number
        NM = Nw + self%RPH_PlanckMass%parameter_number
        NK = NM + self%RPH_Kineticity%parameter_number
        NB = NK + self%RPH_Braiding%parameter_number
        NT = NB + self%RPH_Tensor%parameter_number

        ! check validity of input:
        if ( i > self%parameter_number .or. i <= 0 ) then
            write(*,'(a,I3)') 'No parameter corresponding to: ', i
            write(*,'(a,I3)') 'Total number of parameters is: ', self%parameter_number
            call MpiStop('EFTCAMB error')

        ! parameter from wDE function
        else if ( i <= Nw ) then
            do j = 1, self%RPH_wDE%parameter_number
                if ( i == j ) call self%RPH_wDE%parameter_names_latex( j, latexname )
            end do
            return

        ! parameter from RPH_PlanckMass function
        else if ( i <= NM ) then
            do j = 1, self%RPH_PlanckMass%parameter_number
                if ( i-Nw == j ) call self%RPH_PlanckMass%parameter_names_latex( j, latexname )
            end do
            return

        ! parameter from RPH_Kineticity function
        else if ( i <= NK) then
            do j = 1, self%RPH_Kineticity%parameter_number
                if ( i-NM == j ) call self%RPH_Kineticity%parameter_names_latex( j, latexname )
            end do
            return

        !parameter from RPH_Braiding function
        else if ( i <= NB ) then
            do j = 1, self%RPH_Braiding%parameter_number
                if ( i-NK == j ) call self%RPH_Braiding%parameter_names_latex( j, latexname )
            end do
            return

        !parameter from RPH_Tensor function
        else if ( i <= NT ) then

            do j = 1, self%RPH_Tensor%parameter_number
                if ( i-NB == j ) call self%RPH_Tensor%parameter_names_latex( j, latexname )
            end do
            return

        end if

    end subroutine EFTCAMBRPHParameterNamesLatex

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the i-th parameter name of the model
    subroutine EFTCAMBRPHParameterValues( self, i, value )

        implicit none

        class(EFTCAMB_RPH) :: self   !< the base class
        integer , intent(in)        :: i      !< The index of the parameter
        real(dl), intent(out)       :: value  !< the output value of the i-th parameter

        integer  :: Nw, NM, NK, NB, NT
        integer  :: j

        ! compute the incremental number of parameters:
        Nw = self%RPH_wDE%parameter_number
        NM = Nw + self%RPH_PlanckMass%parameter_number
        NK = NM + self%RPH_Kineticity%parameter_number
        NB = NK + self%RPH_Braiding%parameter_number
        NT = NB + self%RPH_Tensor%parameter_number

        ! check validity of input:
        if ( i > self%parameter_number .or. i <= 0 ) then
            write(*,'(a,I3)') 'No parameter corresponding to: ', i
            write(*,'(a,I3)') 'Total number of parameters is: ', self%parameter_number
            call MpiStop('EFTCAMB error')

        ! parameter from wDE function
        else if ( i <= Nw ) then
            do j = 1, self%RPH_wDE%parameter_number
                if ( i == j ) call self%RPH_wDE%parameter_value( j, value )
            end do
            return

        ! parameter from RPH_PlanckMass function
        else if ( i <= NM ) then
            do j = 1, self%RPH_PlanckMass%parameter_number
                if ( i-Nw == j ) call self%RPH_PlanckMass%parameter_value( j, value )
            end do
            return

        ! parameter from RPH_Kineticity function
        else if ( i <= NK) then
            do j = 1, self%RPH_Kineticity%parameter_number
                if ( i-NM == j ) call self%RPH_Kineticity%parameter_value( j, value )
            end do
            return

        !parameter from RPH_Braiding function
        else if ( i <= NB ) then
            do j = 1, self%RPH_Braiding%parameter_number
                if ( i-NK == j ) call self%RPH_Braiding%parameter_value( j, value )
            end do
            return

        !parameter from RPH_Tensor function
        else if ( i <= NT ) then

            do j = 1, self%RPH_Tensor%parameter_number
                if ( i-NB == j ) call self%RPH_Tensor%parameter_value( j, value )
            end do
            return

        end if

    end subroutine EFTCAMBRPHParameterValues

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes the value of the background EFT functions at a given time.
    subroutine EFTCAMBRPHBackgroundEFTFunctions( self, a, eft_par_cache, eft_cache )

        implicit none

        class(EFTCAMB_RPH)                           :: self          !< the base class
        real(dl), intent(in)                         :: a             !< the input scale factor
        type(EFTCAMB_parameter_cache), intent(inout) :: eft_par_cache !< the EFTCAMB parameter cache that contains all the physical parameters.
        type(EFTCAMB_timestep_cache ), intent(inout) :: eft_cache     !< the EFTCAMB timestep cache that contains all the physical values.

        real(dl) :: RPH_PM_V, RPH_AT_V, RPH_PM_P, RPH_AT_P, RPH_PM_PP, RPH_AT_PP, RPH_PM_PPP, RPH_AT_PPP

        ! precompute some functions:
        RPH_PM_V               = self%RPH_PlanckMass%value(a)
        RPH_AT_V               = self%RPH_Tensor%value(a)
        RPH_PM_P               = self%RPH_PlanckMass%first_derivative(a)
        RPH_AT_P               = self%RPH_Tensor%first_derivative(a)
        RPH_PM_PP              = self%RPH_PlanckMass%second_derivative(a)
        RPH_AT_PP              = self%RPH_Tensor%second_derivative(a)
        RPH_PM_PPP             = self%RPH_PlanckMass%third_derivative(a)
        RPH_AT_PPP             = self%RPH_Tensor%third_derivative(a)
        ! compute the EFT functions:
        eft_cache%EFTOmegaV    = RPH_PM_V +RPH_AT_V*(1._dl +RPH_PM_V)
        eft_cache%EFTOmegaP    = RPH_PM_V*RPH_AT_P +(1._dl +RPH_AT_V)*RPH_PM_P
        eft_cache%EFTOmegaPP   = RPH_PM_V*RPH_AT_PP +2._dl*RPH_PM_P*RPH_AT_P +(1._dl +RPH_AT_V)*RPH_PM_PP
        eft_cache%EFTOmegaPPP  = RPH_PM_V*RPH_AT_PPP +3._dl*RPH_PM_PP*RPH_AT_P +3._dl*RPH_PM_P*RPH_AT_PP &
            & +(1._dl +RPH_AT_V)*RPH_PM_PPP
        eft_cache%EFTc         = ( eft_cache%adotoa**2 - eft_cache%Hdot )*( eft_cache%EFTOmegaV + 0.5_dl*a*eft_cache%EFTOmegaP ) &
            & -0.5_dl*( a*eft_cache%adotoa )**2*eft_cache%EFTOmegaPP&
            & +0.5_dl*eft_cache%grhov_t*( 1._dl+ self%RPH_wDE%value(a) )
        eft_cache%EFTLambda    = -eft_cache%EFTOmegaV*( 2._dl*eft_cache%Hdot +eft_cache%adotoa**2 ) &
            & -a*eft_cache%EFTOmegaP*( 2._dl*eft_cache%adotoa**2 + eft_cache%Hdot ) &
            & -( a*eft_cache%adotoa )**2*eft_cache%EFTOmegaPP &
            & +self%RPH_wDE%value(a)*eft_cache%grhov_t
        eft_cache%EFTcdot      = +0.5_dl*eft_cache%adotoa*eft_cache%grhov_t*( -3._dl*(1._dl +self%RPH_wDE%value(a))**2 + a*self%RPH_wDE%first_derivative(a) ) &
            & -eft_cache%EFTOmegaV*( eft_cache%Hdotdot -4._dl*eft_cache%adotoa*eft_cache%Hdot +2._dl*eft_cache%adotoa**3 ) &
            & +0.5_dl*a*eft_cache%EFTOmegaP*( -eft_cache%Hdotdot +eft_cache%adotoa*eft_cache%Hdot +eft_cache%adotoa**3) &
            & +0.5_dl*a**2*eft_cache%adotoa*eft_cache%EFTOmegaPP*( eft_cache%adotoa**2 -3._dl*eft_cache%Hdot ) &
            & -0.5_dl*(a*eft_cache%adotoa)**3*eft_cache%EFTOmegaPPP
        eft_cache%EFTLambdadot = -2._dl*eft_cache%EFTOmegaV*( eft_cache%Hdotdot -eft_cache%adotoa*eft_cache%Hdot -eft_cache%adotoa**3 ) &
            & -a*eft_cache%EFTOmegaP*( 4._dl*eft_cache%adotoa*eft_cache%Hdot +eft_cache%Hdotdot ) &
            & -a**2*eft_cache%EFTOmegaPP*eft_cache%adotoa*( +2._dl*eft_cache%adotoa**2 +3._dl*eft_cache%Hdot )&
            & -(a*eft_cache%adotoa)**3*eft_cache%EFTOmegaPPP &
            & +eft_cache%grhov_t*eft_cache%adotoa*( a*self%RPH_wDE%first_derivative(a) -3._dl*self%RPH_wDE%value(a)*(1._dl +self%RPH_wDE%value(a) ))

    end subroutine EFTCAMBRPHBackgroundEFTFunctions

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes the value of the background EFT functions at a given time.
    subroutine EFTCAMBRPHSecondOrderEFTFunctions( self, a, eft_par_cache, eft_cache )

        implicit none

        class(EFTCAMB_RPH)                           :: self          !< the base class
        real(dl), intent(in)                         :: a             !< the input scale factor
        type(EFTCAMB_parameter_cache), intent(inout) :: eft_par_cache !< the EFTCAMB parameter cache that contains all the physical parameters.
        type(EFTCAMB_timestep_cache ), intent(inout) :: eft_cache     !< the EFTCAMB timestep cache that contains all the physical values.

        real(dl) :: RPH_PM_V, RPH_AT_V, RPH_AK_V, RPH_AB_V, RPH_PM_P, RPH_AT_P, RPH_AK_P, RPH_AB_P
        real(dl) :: RPH_PM_PP, RPH_AT_PP

        ! precompute some functions:
        RPH_PM_V               = self%RPH_PlanckMass%value(a)
        RPH_AT_V               = self%RPH_Tensor%value(a)
        RPH_AK_V               = self%RPH_Kineticity%value(a)
        RPH_AB_V               = self%RPH_Braiding%value(a)
        RPH_PM_P               = self%RPH_PlanckMass%first_derivative(a)
        RPH_AT_P               = self%RPH_Tensor%first_derivative(a)
        RPH_AK_P               = self%RPH_Kineticity%first_derivative(a)
        RPH_AB_P               = self%RPH_Braiding%first_derivative(a)
        RPH_PM_PP              = self%RPH_PlanckMass%second_derivative(a)
        RPH_AT_PP              = self%RPH_Tensor%second_derivative(a)
        ! compute the EFT functions:
        eft_cache%EFTGamma1V  = 0.25_dl*( RPH_AK_V*(1._dl +RPH_PM_V)*eft_cache%adotoa**2 &
            & -2._dl*eft_cache%EFTc )/(eft_par_cache%h0_Mpc**2*a**2)
        eft_cache%EFTGamma1P  = - 0.5_dl*( RPH_AK_V*(1._dl +RPH_PM_V)*eft_cache%adotoa**2 &
            & -2._dl*eft_cache%EFTc )/(eft_par_cache%h0_Mpc**2*a**3) &
            & +0.25_dl*( RPH_AK_P*(1._dl +RPH_PM_V)*eft_cache%adotoa**2 &
            & +RPH_AK_V*RPH_PM_P*eft_cache%adotoa**2 &
            & +2._dl*RPH_AK_V*(1._dl +RPH_PM_V)*eft_cache%Hdot/a &
            & -4._dl*eft_cache%EFTc -2._dl*eft_cache%EFTcdot/a/eft_cache%adotoa )/(eft_par_cache%h0_Mpc**2*a**2)
        eft_cache%EFTGamma2V  = ( +2._dl*RPH_AB_V*(1._dl +RPH_PM_V) &
            & -a*eft_cache%EFTOmegaP )*eft_cache%adotoa/(eft_par_cache%h0_Mpc*a)
        eft_cache%EFTGamma2P  = -0.5_dl*( +2._dl*RPH_AB_V*(1._dl +RPH_PM_V) &
            & -a*eft_cache%EFTOmegaP )*eft_cache%adotoa/(eft_par_cache%h0_Mpc*a) &
            & -( -2._dl*(1._dl +RPH_PM_V)*( RPH_AB_P*eft_cache%adotoa**2 &
            & + RPH_AB_V*eft_cache%Hdot/a) &
            & - 2._dl*RPH_AB_V*eft_cache%adotoa**2*RPH_PM_P &
            & + eft_cache%EFTOmegaP*( eft_cache%adotoa**2 +eft_cache%Hdot ) &
            & + a*eft_cache%adotoa**2*eft_cache%EFTOmegaPP )/(eft_par_cache%h0_Mpc*a*eft_cache%adotoa)
        eft_cache%EFTGamma3V  = -RPH_AT_V*(1._dl +RPH_PM_V)
        eft_cache%EFTGamma3P  = -RPH_PM_P*RPH_AT_V -(1._dl +RPH_PM_V)*RPH_AT_P
        eft_cache%EFTGamma4V  = -eft_cache%EFTGamma3V
        eft_cache%EFTGamma4P  = -eft_cache%EFTGamma3P
        eft_cache%EFTGamma4PP = +(1._dl +RPH_PM_V)*RPH_AT_PP +RPH_PM_PP*RPH_AT_V &
            & +2._dl*RPH_PM_P*RPH_AT_P
        eft_cache%EFTGamma5V  = +0.5_dl*eft_cache%EFTGamma3V
        eft_cache%EFTGamma5P  = +0.5_dl*eft_cache%EFTGamma3P
        eft_cache%EFTGamma6V  = 0._dl
        eft_cache%EFTGamma6P  = 0._dl

    end subroutine EFTCAMBRPHSecondOrderEFTFunctions

    ! ---------------------------------------------------------------------------------------------
    !> Function that computes dtauda = 1/sqrt(a^2H^2).For pure EFT std this has to be overridden
    !! for performance reasons.
    function EFTCAMBRPHComputeDtauda( self, a, eft_par_cache, eft_cache )

        implicit none

        class(EFTCAMB_RPH)                           :: self          !< the base class
        real(dl), intent(in)                         :: a             !< the input scale factor
        type(EFTCAMB_parameter_cache), intent(inout) :: eft_par_cache !< the EFTCAMB parameter cache that contains all the physical parameters.
        type(EFTCAMB_timestep_cache ), intent(inout) :: eft_cache     !< the EFTCAMB timestep cache that contains all the physical values.

        real(dl) :: EFTCAMBRPHComputeDtauda                           !< the output dtauda

        real(dl) :: temp

        temp = eft_cache%grhoa2 +eft_par_cache%grhov*a*a*self%RPH_wDE%integral(a)
        EFTCAMBRPHComputeDtauda = sqrt(3/temp)

    end function EFTCAMBRPHComputeDtauda

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes adotoa = H and its two derivatives wrt conformal time.
    subroutine EFTCAMBRPHComputeAdotoa( self, a, eft_par_cache, eft_cache )

        implicit none

        class(EFTCAMB_RPH)                           :: self          !< the base class
        real(dl), intent(in)                         :: a             !< the input scale factor
        type(EFTCAMB_parameter_cache), intent(inout) :: eft_par_cache !< the EFTCAMB parameter cache that contains all the physical parameters.
        type(EFTCAMB_timestep_cache ), intent(inout) :: eft_cache     !< the EFTCAMB timestep cache that contains all the physical values.

        eft_cache%grhov_t = eft_par_cache%grhov*self%RPH_wDE%integral(a)
        eft_cache%adotoa  = sqrt( ( eft_cache%grhom_t +eft_cache%grhov_t )/3._dl )

    end subroutine EFTCAMBRPHComputeAdotoa

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes the two derivatives wrt conformal time of H.
    subroutine EFTCAMBRPHComputeHubbleDer( self, a, eft_par_cache, eft_cache )

        implicit none

        class(EFTCAMB_RPH)                           :: self          !< the base class
        real(dl), intent(in)                         :: a             !< the input scale factor
        type(EFTCAMB_parameter_cache), intent(inout) :: eft_par_cache !< the EFTCAMB parameter cache that contains all the physical parameters.
        type(EFTCAMB_timestep_cache ), intent(inout) :: eft_cache     !< the EFTCAMB timestep cache that contains all the physical values.

        eft_cache%gpiv_t  = self%RPH_wDE%value(a)*eft_cache%grhov_t
        eft_cache%Hdot    = -0.5_dl*( eft_cache%adotoa**2 +eft_cache%gpresm_t +eft_cache%gpiv_t )

        ! IW POSSIBLE BUG
        !eft_cache%Hdotdot = eft_cache%adotoa*( ( eft_cache%grhob_t +eft_cache%grhoc_t)/6._dl +2._dl*( eft_cache%grhor_t +eft_cache%grhog_t)/3._dl ) &
        !    & +eft_cache%adotoa*eft_cache%grhov_t*( 1._dl/6._dl +self%RPH_wDE%value(a) +1.5_dl*self%RPH_wDE%value(a)**2 -0.5_dl*a*self%RPH_wDE%first_derivative(a) ) &
        !    & +eft_cache%adotoa*eft_cache%grhonu_tot/6._dl -0.5_dl*eft_cache%adotoa*eft_cache%gpinu_tot -0.5_dl*eft_cache%gpinudot_tot

        eft_cache%Hdotdot = 2._dl*eft_cache%adotoa*eft_cache%Hdot &
            & + 0.5_dl*eft_cache%adotoa*(eft_cache%grhob_t + eft_cache%grhoc_t + 8._dl*(eft_cache%grhog_t+eft_cache%grhor_t)/3._dl)&
            & + 0.5_dl*eft_cache%adotoa*eft_cache%grhov_t*( (1._dl+self%RPH_wDE%value(a) )*(1._dl+3._dl*self%RPH_wDE%value(a)) -a*self%RPH_wDE%first_derivative(a))&
            & + eft_cache%adotoa/6._dl*eft_cache%grhonu_tot -0.5_dl*eft_cache%adotoa*eft_cache%gpinu_tot -0.5_dl*eft_cache%gpinudot_tot

    end subroutine EFTCAMBRPHComputeHubbleDer

    ! ---------------------------------------------------------------------------------------------
    !> Function that computes model specific stability requirements.
    function EFTCAMBRPHAdditionalModelStability( self, a, eft_par_cache, eft_cache )

        implicit none

        class(EFTCAMB_RPH)                           :: self          !< the base class
        real(dl), intent(in)                         :: a             !< the input scale factor.
        type(EFTCAMB_parameter_cache), intent(inout) :: eft_par_cache !< the EFTCAMB parameter cache that contains all the physical parameters.
        type(EFTCAMB_timestep_cache ), intent(inout) :: eft_cache     !< the EFTCAMB timestep cache that contains all the physical values.

        logical :: EFTCAMBRPHAdditionalModelStability          !< the return value of the stability computation. True if the model specific stability criteria are met, false otherwise.

        EFTCAMBRPHAdditionalModelStability = .True.
        if ( self%RPH_wDE%value(a) > -1._dl/3._dl ) EFTCAMBRPHAdditionalModelStability = .False.

    end function EFTCAMBRPHAdditionalModelStability

    ! ---------------------------------------------------------------------------------------------

end module EFTCAMB_Reparametrized_Horndeski
