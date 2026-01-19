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

!> @file 007p2_RPH.f90
!! This file contains the definition of an EFT reparametrization of Horndeski models.
!! This is described by four functions of time and w_DE.
!! Please refer to the numerical notes for details.

!----------------------------------------------------------------------------------------
!> This module contains the definition of an EFT reparametrization of Horndeski models.
!! This is described by four functions of time and w_DE.
!! Please refer to the numerical notes for details.

!> @author Bin Hu, Marco Raveri, Gen Ye

module EFTCAMB_Reparametrized_Horndeski

    use precision
    use IniObjects
    use EFTCAMB_cache
    use MpiUtils
    use EFT_def
    use EFTCAMB_abstract_parametrizations_1D
    use EFTCAMB_parametrizations_1D
    use equispaced_linear_interpolation_1D
    use EFTCAMB_abstract_model_designer
    use MassiveNu

    implicit none

    private

    public EFTCAMB_RPH

    !----------------------------------------------------------------------------------------
    !> This is the type that contains the definition of RPH.
    type, extends ( EFTCAMB_designer_model ) :: EFTCAMB_RPH

        ! the RPH functions model selection flags:
        ! Each function X is the sum of a function of scale factor, defined by Xmodel, and a function of Omega_DE, defined by Xmodel_ODE. E.g. alpha_M = F(a) + G(Omega_DE)  
        integer  :: RPHwDE              !< Model selection flag for RPH w DE.
        integer  :: RPHalphaMmodel      !< Model selection flag for RPH alphaM.
        integer  :: RPHmassPmodel       !< Model selection flag for RPH Planck mass model.
        integer  :: RPHkineticitymodel  !< Model selection flag for RPH kineticity model.
        integer  :: RPHbraidingmodel    !< Model selection flag for RPH braiding model.
        integer  :: RPHtensormodel      !< Model selection flag for RPH tensor model.

        integer  :: RPHalphaMmodel_ODE      !< Model selection flag for RPH alphaM.
        integer  :: RPHmassPmodel_ODE       !< Model selection flag for RPH Planck mass model.
        integer  :: RPHkineticitymodel_ODE  !< Model selection flag for RPH kineticity model.
        integer  :: RPHbraidingmodel_ODE    !< Model selection flag for RPH braiding model.
        integer  :: RPHtensormodel_ODE      !< Model selection flag for RPH tensor model.

        logical  :: RPHusealphaM            !< Flag whether using alphaM instead of effective Planck Mass
        logical  :: RPHhasOmegaDE           !< Flag whether using the OmegaDE parametrization
        logical  :: RPHintegratefromtoday   !< Flag whether integrating from today or radiation era
        
        real(dl) :: RPH_M0                  !< Modification to Planck Mass at today, 1 + M0 = Meff^2(a=1)/m_0^2

        ! the RPH functions:
        class( parametrized_function_1D ), allocatable :: RPH_wDE         !< The RPH function w_DE.
        class( parametrized_function_1D ), allocatable :: RPH_alphaM      !< The RPH function alpha Mass.
        class( parametrized_function_1D ), allocatable :: RPH_PlanckMass  !< The RPH function Planck Mass.
        class( parametrized_function_1D ), allocatable :: RPH_Kineticity  !< The RPH function Kineticity.
        class( parametrized_function_1D ), allocatable :: RPH_Braiding    !< The RPH function Braiding.
        class( parametrized_function_1D ), allocatable :: RPH_Tensor      !< The RPH function Tensor.

        ! the RPH functions of Omega_DE
        class( parametrized_function_1D ), allocatable :: RPH_PlanckMass_ODE  !< The RPH function Planck Mass.
        class( parametrized_function_1D ), allocatable :: RPH_alphaM_ODE      !< The RPH function alpha Mass.
        class( parametrized_function_1D ), allocatable :: RPH_Kineticity_ODE  !< The RPH function Kineticity.
        class( parametrized_function_1D ), allocatable :: RPH_Braiding_ODE    !< The RPH function Braiding.
        class( parametrized_function_1D ), allocatable :: RPH_Tensor_ODE      !< The RPH function Tensor.

        ! output of background solver:
        type( equispaced_linear_interpolate_function_1D ) :: rhov_t           !< The DE energy density
        type( equispaced_linear_interpolate_function_1D ) :: adotoa           !< The conformal Hubble
        type( equispaced_linear_interpolate_function_1D ) :: Hdot             !< The conformal derivative of conformal Hubble
        type( equispaced_linear_interpolate_function_1D ) :: Hdot2            !< The second conformal derivative of conformal Hubble
        type( equispaced_linear_interpolate_function_1D ) :: Hdot3            !< The third conformal derivative of conformal Hubble
        type( equispaced_linear_interpolate_function_1D ) :: Hdot4            !< The fourth conformal derivative of conformal Hubble
        type( equispaced_linear_interpolate_function_1D ) :: RPH_M            !< The RPH Planck Mass function M = Meff^2-1
        type( equispaced_linear_interpolate_function_1D ) :: Ode              !< The DE energy fraction
        type( equispaced_linear_interpolate_function_1D ) :: dOde             !< The derivative of DE energy fraction wrt scale factor
        type( equispaced_linear_interpolate_function_1D ) :: d2Ode            !< The second derivative of DE energy fraction wrt scale factor
        type( equispaced_linear_interpolate_function_1D ) :: d3Ode            !< The third derivative of DE energy fraction wrt scale factor
        type( equispaced_linear_interpolate_function_1D ) :: d4Ode             !< The fourth derivative of DE energy fraction wrt scale factor

        

        ! background solver parameters:
        integer  :: designer_num_points = 6000                            !< Number of points sampled by the designer code.
        real(dl) :: x_initial           = log(10._dl**(-8))               !< log(a start)
        real(dl) :: x_final             = 0.0_dl                          !< log(a final)

    contains

        ! initialization of the model:
        procedure :: read_model_selection            => EFTCAMBRPHReadModelSelectionFromFile  !< subroutine that reads the parameters of the model from file
        procedure :: allocate_model_selection        => EFTCAMBRPHAllocateModelSelection      !< subroutine that allocates the model selection.
        procedure :: init_model_parameters           => EFTCAMBRPHInitModelParameters         !< subroutine that initializes the model parameters based on the values found in an input array.
        procedure :: init_model_parameters_from_file => EFTCAMBRPHInitModelParametersFromFile !< subroutine that reads the parameters of the model from file.

        ! background solver:
        procedure :: initialize_background           => EFTCAMBRPHInitBackground               !< subroutine that initializes the background.
        procedure :: solve_designer_equations        => EFTCAMBRPHSolveDesignerEquations       !< subroutine that solves for the Mp.

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
    subroutine EFTCAMBRPHReadModelSelectionFromFile( self, Ini, eft_error )

        implicit none

        class(EFTCAMB_RPH) :: self      !< the base class
        type(TIniFile)     :: Ini       !< Input ini file
        integer            :: eft_error !< error code: 0 all fine, 1 initialization failed

        ! read alphaM flag
        self%RPHusealphaM       = Ini%Read_Logical( 'RPHusealphaM', .false. )
        ! read integration flag
        self%RPHintegratefromtoday = Ini%Read_Logical( 'RPHintegratefromtoday', .false. )
        
        ! read model selection flags:
        self%RPHwDE             = Ini%Read_Int( 'RPHwDE'            , 0 )
        self%RPHalphaMmodel     = Ini%Read_Int( 'RPHalphaMmodel'    , 0 )
        self%RPHmassPmodel      = Ini%Read_Int( 'RPHmassPmodel'     , 0 )
        self%RPHkineticitymodel = Ini%Read_Int( 'RPHkineticitymodel', 0 )
        self%RPHbraidingmodel   = Ini%Read_Int( 'RPHbraidingmodel'  , 0 )
        self%RPHtensormodel     = Ini%Read_Int( 'RPHtensormodel'    , 0 )

        self%RPHalphaMmodel_ODE      = Ini%Read_Int( 'RPHalphaMmodel_ODE'   , 0 )
        self%RPHmassPmodel_ODE       = Ini%Read_Int( 'RPHmassPmodel_ODE'    , 0 )
        self%RPHkineticitymodel_ODE = Ini%Read_Int( 'RPHkineticitymodel_ODE', 0 )
        self%RPHbraidingmodel_ODE   = Ini%Read_Int( 'RPHbraidingmodel_ODE'  , 0 )
        self%RPHtensormodel_ODE     = Ini%Read_Int( 'RPHtensormodel_ODE'    , 0 )

        if ( (self%RPHusealphaM .and. self%RPHalphaMmodel_ODE > 0) .or. ((.not. self%RPHusealphaM) .and. self%RPHmassPmodel_ODE > 0) .or. self%RPHkineticitymodel_ODE > 0 .or. self%RPHbraidingmodel_ODE > 0 .or. self%RPHtensormodel_ODE > 0 ) then
            self%RPHhasOmegaDE = .true.
        else
            self%RPHhasOmegaDE = .false.
        end if

        ! read precision parameters
        self%designer_num_points = Ini%Read_Int( 'model_background_num_points', 6000 )
        self%x_initial = Log( Ini%Read_Double( 'model_background_a_ini', 1d-8 ) )
        self%x_final = Log( Ini%Read_Double( 'model_background_a_final', 1._dl ) )

    end subroutine EFTCAMBRPHReadModelSelectionFromFile

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that allocates the model selection.
    subroutine EFTCAMBRPHAllocateModelSelection( self, Ini, eft_error )

        implicit none

        class(EFTCAMB_RPH)  :: self      !< the base class
        type(TIniFile)      :: Ini       !< Input ini file
        integer             :: eft_error !< error code: 0 all fine, 1 initialization failed

        integer             :: temp_feedback 

        ! get feedback flag:
        temp_feedback = Ini%Read_Int('feedback_level', 0)

        ! allocate wDE:
        if ( allocated(self%RPH_wDE) ) deallocate(self%RPH_wDE)
        select case ( self%RPHwDE )
            case(0)
                allocate( wDE_LCDM_parametrization_1D::self%RPH_wDE )
                call self%RPH_wDE%set_name( 'RPHw', 'w' )
            case(1)
                allocate( constant_parametrization_1D::self%RPH_wDE )
                call self%RPH_wDE%set_name( 'RPHw', 'w' )
            case(2)
                allocate( CPL_parametrization_1D::self%RPH_wDE )
                call self%RPH_wDE%set_param_names( ['RPHw0', 'RPHwa'], ['w_0', 'w_a'] )
                call self%RPH_wDE%set_name( 'RPHw', 'w' )
            case(3)
                allocate( JBP_parametrization_1D::self%RPH_wDE )
                call self%RPH_wDE%set_param_names( ['RPHw0', 'RPHwa', 'RPHwn'], [ 'w_0', 'w_a', 'n  ' ] )
                call self%RPH_wDE%set_name( 'RPHw', 'w' )
            case(4)
                allocate( turning_point_parametrization_1D::self%RPH_wDE )
                call self%RPH_wDE%set_param_names( ['RPHw0 ', 'RPHwa ', 'RPHwat'], ['w_0', 'w_a', 'a_t'] )
                call self%RPH_wDE%set_name( 'RPHw', 'w' )
            case(5)
                allocate( taylor_parametrization_1D::self%RPH_wDE )
                call self%RPH_wDE%set_param_names( ['RPHw0', 'RPHwa', 'RPHw2', 'RPHw3'], ['w_0', 'w_a', 'w_2', 'w_3'] )
                call self%RPH_wDE%set_name( 'RPHw', 'w' )
            case default
                call allocate_parametrized_1D_function( self%RPH_wDE, self%RPHwDE, 'RPHw', 'w', eft_error, temp_feedback )
                if ( eft_error == 1 ) return
        end select
        
        if ( self%RPHusealphaM ) then
            ! allocate RPH_alphaM:
            call allocate_parametrized_1D_function( self%RPH_alphaM, self%RPHalphaMmodel, 'RPHalphaM', '{\alpha^{\rm M}}', eft_error, temp_feedback )
            if ( eft_error == 1 ) return
        else
            ! allocate RPH_PlanckMass:
            call allocate_parametrized_1D_function( self%RPH_PlanckMass, self%RPHmassPmodel, 'RPHmassP', '\tilde{M}', eft_error, temp_feedback )
            if ( eft_error == 1 ) return
        end if
        ! allocate RPH_Kineticity:
        call allocate_parametrized_1D_function( self%RPH_Kineticity, self%RPHkineticitymodel, 'RPHkineticity', '{\alpha^{\rm K}}', eft_error, temp_feedback )
        if ( eft_error == 1 ) return
        ! allocate RPH_Braiding:
        call allocate_parametrized_1D_function( self%RPH_Braiding, self%RPHbraidingmodel, 'RPHbraiding', '{\alpha^{\rm B}}', eft_error, temp_feedback )
        if ( eft_error == 1 ) return
        ! allocate RPH_Tensor:
        call allocate_parametrized_1D_function( self%RPH_Tensor, self%RPHtensormodel, 'RPHtensor', '{\alpha^{\rm T}}', eft_error, temp_feedback )
        if ( eft_error == 1 ) return

        if ( self%RPHusealphaM ) then
            ! allocate RPH_alphaM_ODE:
            call allocate_parametrized_1D_function( self%RPH_alphaM_ODE, self%RPHalphaMmodel_ODE, 'RPHalphaM_ODE', '{\alpha^{\rm M,ODE}}', eft_error, temp_feedback )
            if ( eft_error == 1 ) return
        else
            ! allocate RPH_PlanckMass:
            call allocate_parametrized_1D_function( self%RPH_PlanckMass_ODE, self%RPHmassPmodel_ODE, 'RPHmassP_ODE', '\tilde{M}_{ODE}', eft_error, temp_feedback )
            if ( eft_error == 1 ) return
        end if
        ! allocate RPH_Kineticity_ODE:
        call allocate_parametrized_1D_function( self%RPH_Kineticity_ODE, self%RPHkineticitymodel_ODE, 'RPHkineticity_ODE', '{\alpha^{\rm K,ODE}}', eft_error, temp_feedback )
        if ( eft_error == 1 ) return
        ! allocate RPH_Braiding_ODE:
        call allocate_parametrized_1D_function( self%RPH_Braiding_ODE, self%RPHbraidingmodel_ODE, 'RPHbraiding_ODE', '{\alpha^{\rm B,ODE}}', eft_error, temp_feedback )
        if ( eft_error == 1 ) return
        ! allocate RPH_Tensor_ODE:
        call allocate_parametrized_1D_function( self%RPH_Tensor_ODE, self%RPHtensormodel_ODE, 'RPHtensor_ODE', '{\alpha^{\rm T,ODE}}', eft_error, temp_feedback )
        if ( eft_error == 1 ) return

        ! additional initialization of the function:
        call self%RPH_wDE%init_func_from_file       ( Ini, eft_error )
        if ( self%RPHusealphaM ) then
            call self%RPH_alphaM%init_func_from_file( Ini, eft_error )
        else
            call self%RPH_PlanckMass%init_func_from_file( Ini, eft_error )
        end if
        call self%RPH_Kineticity%init_func_from_file( Ini, eft_error )
        call self%RPH_Braiding%init_func_from_file  ( Ini, eft_error )
        call self%RPH_Tensor%init_func_from_file    ( Ini, eft_error )

        if ( self%RPHusealphaM ) then
            call self%RPH_alphaM_ODE%init_func_from_file( Ini, eft_error )
        else
            call self%RPH_PlanckMass_ODE%init_func_from_file( Ini, eft_error )
        end if
        call self%RPH_Kineticity_ODE%init_func_from_file( Ini, eft_error )
        call self%RPH_Braiding_ODE%init_func_from_file  ( Ini, eft_error )
        call self%RPH_Tensor_ODE%init_func_from_file    ( Ini, eft_error )

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
        if ( self%RPHusealphaM ) then
            ! then RPH_alphaM parameters:
            num_params_function = self%RPH_alphaM%parameter_number
            allocate( temp(num_params_function) )
            do i = 1, num_params_function
                temp(i)         = array(num_params_temp)
                num_params_temp = num_params_temp +1
            end do
            call self%RPH_alphaM%init_parameters(temp)
            self%RPH_M0 = array(num_params_temp)
            num_params_temp = num_params_temp + 1
            deallocate(temp)
        else
            ! then RPH_PlanckMass parameters:
            num_params_function = self%RPH_PlanckMass%parameter_number
            allocate( temp(num_params_function) )
            do i = 1, num_params_function
                temp(i)         = array(num_params_temp)
                num_params_temp = num_params_temp +1
            end do
            call self%RPH_PlanckMass%init_parameters(temp)
            deallocate(temp)
        end if
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

        ! The ODE version
        if ( self%RPHusealphaM ) then
            ! then RPH_alphaM parameters:
            num_params_function = self%RPH_alphaM_ODE%parameter_number
            allocate( temp(num_params_function) )
            do i = 1, num_params_function
                temp(i)         = array(num_params_temp)
                num_params_temp = num_params_temp +1
            end do
            call self%RPH_alphaM_ODE%init_parameters(temp)
            deallocate(temp)
        else
            ! then RPH_PlanckMass parameters:
            num_params_function = self%RPH_PlanckMass_ODE%parameter_number
            allocate( temp(num_params_function) )
            do i = 1, num_params_function
                temp(i)         = array(num_params_temp)
                num_params_temp = num_params_temp +1
            end do
            call self%RPH_PlanckMass_ODE%init_parameters(temp)
            deallocate(temp)
        end if
        ! then RPH_Kineticity parameters:
        num_params_function = self%RPH_Kineticity_ODE%parameter_number
        allocate( temp(num_params_function) )
        do i = 1, num_params_function
            temp(i)         = array(num_params_temp)
            num_params_temp = num_params_temp +1
        end do
        call self%RPH_Kineticity_ODE%init_parameters(temp)
        deallocate(temp)
        ! then RPH_Braiding parameters:
        num_params_function = self%RPH_Braiding_ODE%parameter_number
        allocate( temp(num_params_function) )
        do i = 1, num_params_function
            temp(i)         = array(num_params_temp)
            num_params_temp = num_params_temp +1
        end do
        call self%RPH_Braiding_ODE%init_parameters(temp)
        deallocate(temp)
        ! then RPH_Tensor parameters:
        num_params_function = self%RPH_Tensor_ODE%parameter_number
        allocate( temp(num_params_function) )
        do i = 1, num_params_function
            temp(i)         = array(num_params_temp)
            num_params_temp = num_params_temp +1
        end do
        call self%RPH_Tensor_ODE%init_parameters(temp)
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
    subroutine EFTCAMBRPHInitModelParametersFromFile( self, Ini, eft_error )

        implicit none

        class(EFTCAMB_RPH) :: self      !< the base class
        type(TIniFile)     :: Ini       !< Input ini file
        integer            :: eft_error !< error code: 0 all fine, 1 initialization failed

        call self%RPH_wDE%init_from_file( Ini, eft_error )
        if ( self%RPHusealphaM ) then
            call self%RPH_alphaM%init_from_file( Ini, eft_error )
            self%RPH_M0 = Ini%Read_Double( 'RPH_M0', 0._dl )
        else
            call self%RPH_PlanckMass%init_from_file( Ini, eft_error )
        end if
        call self%RPH_Kineticity%init_from_file( Ini, eft_error )
        call self%RPH_Braiding%init_from_file( Ini, eft_error )
        call self%RPH_Tensor%init_from_file( Ini, eft_error )

        if ( self%RPHusealphaM ) then
            call self%RPH_alphaM_ODE%init_from_file( Ini, eft_error )
        else   
            call self%RPH_PlanckMass_ODE%init_from_file( Ini, eft_error )
        end if
        call self%RPH_Kineticity_ODE%init_from_file( Ini, eft_error )
        call self%RPH_Braiding_ODE%init_from_file( Ini, eft_error )
        call self%RPH_Tensor_ODE%init_from_file( Ini, eft_error )

    end subroutine EFTCAMBRPHInitModelParametersFromFile

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes the number of parameters of the model.
    subroutine EFTCAMBRPHComputeParametersNumber( self )

        implicit none

        class(EFTCAMB_RPH)  :: self   !< the base class

        self%parameter_number = 0
        self%parameter_number = self%parameter_number +self%RPH_wDE%parameter_number
        if ( self%RPHusealphaM ) then
            self%parameter_number = self%parameter_number + self%RPH_alphaM%parameter_number + 1
        else
            self%parameter_number = self%parameter_number +self%RPH_PlanckMass%parameter_number
        end if
        self%parameter_number = self%parameter_number +self%RPH_Kineticity%parameter_number
        self%parameter_number = self%parameter_number +self%RPH_Braiding%parameter_number
        self%parameter_number = self%parameter_number +self%RPH_Tensor%parameter_number

        if ( self%RPHusealphaM ) then
            self%parameter_number = self%parameter_number +self%RPH_alphaM_ODE%parameter_number
        else
            self%parameter_number = self%parameter_number +self%RPH_PlanckMass_ODE%parameter_number
        end if
        self%parameter_number = self%parameter_number +self%RPH_Kineticity_ODE%parameter_number
        self%parameter_number = self%parameter_number +self%RPH_Braiding_ODE%parameter_number
        self%parameter_number = self%parameter_number +self%RPH_Tensor_ODE%parameter_number

    end subroutine EFTCAMBRPHComputeParametersNumber

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that prints on the screen feedback information about the model.
    subroutine EFTCAMBRPHFeedback( self, print_params )

        implicit none

        class(EFTCAMB_RPH)  :: self         !< the base class
        logical, optional   :: print_params !< optional flag that decised whether to print numerical values
                                            !! of the parameters.

        ! print general model informations:
        write(*,*)
        write(*,'(a,a)')    '   Model               =  ', self%name
        write(*,'(a,I3)')   '   Number of params    ='  , self%parameter_number

        ! print model functions informations:
        write(*,*)
        if ( self%RPHwDE             /= 0 ) write(*,'(a,I3)') '   RPHwDE             =', self%RPHwDE
        if ( self%RPHusealphaM ) then
            if ( self%RPHalphaMmodel      /= 0 ) write(*,'(a,I3)') '   RPHalphaMmodel      =', self%RPHalphaMmodel
            if ( self%RPHalphaMmodel      /= 0 ) write(*,'(a,I3)') '   RPH_M0              =', self%RPH_M0
        else
            if ( self%RPHmassPmodel      /= 0 ) write(*,'(a,I3)') '   RPHmassPmodel      =', self%RPHmassPmodel
        end if
        if ( self%RPHkineticitymodel /= 0 ) write(*,'(a,I3)') '   RPHkineticitymodel =', self%RPHkineticitymodel
        if ( self%RPHbraidingmodel   /= 0 ) write(*,'(a,I3)') '   RPHbraidingmodel   =', self%RPHbraidingmodel
        if ( self%RPHtensormodel     /= 0 ) write(*,'(a,I3)') '   RPHtensormodel     =', self%RPHtensormodel

        if ( self%RPHusealphaM ) then
            if ( self%RPHalphaMmodel_ODE      /= 0 ) write(*,'(a,I3)') '   RPHalphaMmodel_ODE      =', self%RPHalphaMmodel_ODE
        else
            if ( self%RPHmassPmodel_ODE      /= 0 ) write(*,'(a,I3)') '   RPHmassPmodel_ODE      =', self%RPHmassPmodel_ODE
        end if
        if ( self%RPHkineticitymodel_ODE /= 0 ) write(*,'(a,I3)') '   RPHkineticitymodel_ODE =', self%RPHkineticitymodel_ODE
        if ( self%RPHbraidingmodel_ODE   /= 0 ) write(*,'(a,I3)') '   RPHbraidingmodel_ODE   =', self%RPHbraidingmodel_ODE
        if ( self%RPHtensormodel_ODE     /= 0 ) write(*,'(a,I3)') '   RPHtensormodel_ODE     =', self%RPHtensormodel_ODE

        write(*,*)
        ! print functions informations:
        call self%RPH_wDE%feedback        ( print_params )
        if ( self%RPHusealphaM ) then
            call self%RPH_alphaM%feedback ( print_params )
        else
            call self%RPH_PlanckMass%feedback ( print_params )
        end if
        call self%RPH_Kineticity%feedback ( print_params )
        call self%RPH_Braiding%feedback   ( print_params )
        call self%RPH_Tensor%feedback     ( print_params )
        
        if ( self%RPHusealphaM ) then
            call self%RPH_alphaM_ODE%feedback ( print_params )
        else
            call self%RPH_PlanckMass_ODE%feedback ( print_params )
        end if
        call self%RPH_Kineticity_ODE%feedback ( print_params )
        call self%RPH_Braiding_ODE%feedback   ( print_params )
        call self%RPH_Tensor_ODE%feedback     ( print_params )

    end subroutine EFTCAMBRPHFeedback

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the i-th parameter name of the model
    subroutine EFTCAMBRPHParameterNames( self, i, name )

        implicit none

        class(EFTCAMB_RPH) :: self   !< the base class
        integer     , intent(in)    :: i      !< the index of the parameter
        character(*), intent(out)   :: name   !< the output name of the i-th parameter

        integer  :: Nw, NM, NK, NB, NT, NMo, NKo, NBo, NTo
        integer  :: j

        ! compute the incremental number of parameters:
        Nw = self%RPH_wDE%parameter_number
        if ( self%RPHusealphaM ) then
            NM = Nw + self%RPH_alphaM%parameter_number + 1
        else
            NM = Nw + self%RPH_PlanckMass%parameter_number
        end if
        NK = NM + self%RPH_Kineticity%parameter_number
        NB = NK + self%RPH_Braiding%parameter_number
        NT = NB + self%RPH_Tensor%parameter_number

        if ( self%RPHusealphaM ) then
            NMo = NT + self%RPH_alphaM_ODE%parameter_number
        else
            NMo = NT + self%RPH_Tensor_ODE%parameter_number
        end if
        NKo = NMo + self%RPH_Kineticity_ODE%parameter_number
        NBo = NKo + self%RPH_Braiding_ODE%parameter_number
        NTo = NBo + self%RPH_Tensor_ODE%parameter_number

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
            if ( self%RPHusealphaM ) then
                do j = 1, self%RPH_alphaM%parameter_number
                    if ( i-Nw == j ) call self%RPH_alphaM%parameter_names( j, name )
                end do
                if ( i == NM ) name = 'RPH_M0'
            else
                do j = 1, self%RPH_PlanckMass%parameter_number
                    if ( i-Nw == j ) call self%RPH_PlanckMass%parameter_names( j, name )
                end do
            end if
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

        ! parameter from RPH_PlanckMass function
        else if ( i <= NMo ) then
            if ( self%RPHusealphaM ) then
                do j = 1, self%RPH_alphaM_ODE%parameter_number
                    if ( i-NT == j ) call self%RPH_alphaM_ODE%parameter_names( j, name )
                end do
            else
                do j = 1, self%RPH_PlanckMass_ODE%parameter_number
                    if ( i-NT == j ) call self%RPH_PlanckMass_ODE%parameter_names( j, name )
                end do
            end if
            return

        ! parameter from RPH_Kineticity function
        else if ( i <= NKo) then
            do j = 1, self%RPH_Kineticity_ODE%parameter_number
                if ( i-NMo == j ) call self%RPH_Kineticity_ODE%parameter_names( j, name )
            end do
            return

        !parameter from RPH_Braiding function
        else if ( i <= NBo ) then
            do j = 1, self%RPH_Braiding_ODE%parameter_number
                if ( i-NKo == j ) call self%RPH_Braiding_ODE%parameter_names( j, name )
            end do
            return

        !parameter from RPH_Tensor function
        else if ( i <= NTo ) then

            do j = 1, self%RPH_Tensor_ODE%parameter_number
                if ( i-NBo == j ) call self%RPH_Tensor_ODE%parameter_names( j, name )
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

        integer  :: Nw, NM, NK, NB, NT, NMo, NKo, NBo, NTo
        integer  :: j

        ! compute the incremental number of parameters:
        Nw = self%RPH_wDE%parameter_number
        if ( self%RPHusealphaM ) then
            NM = Nw + self%RPH_alphaM%parameter_number + 1
        else
            NM = Nw + self%RPH_PlanckMass%parameter_number
        end if
        NK = NM + self%RPH_Kineticity%parameter_number
        NB = NK + self%RPH_Braiding%parameter_number
        NT = NB + self%RPH_Tensor%parameter_number

        if ( self%RPHusealphaM ) then
            NMo = NT + self%RPH_alphaM_ODE%parameter_number
        else
            NMo = NT + self%RPH_Tensor_ODE%parameter_number
        end if
        NKo = NMo + self%RPH_Kineticity_ODE%parameter_number
        NBo = NKo + self%RPH_Braiding_ODE%parameter_number
        NTo = NBo + self%RPH_Tensor_ODE%parameter_number

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
            if ( self%RPHusealphaM ) then
                do j = 1, self%RPH_alphaM%parameter_number
                    if ( i-Nw == j ) call self%RPH_alphaM%parameter_names_latex( j, latexname )
                end do
                if ( i == NM ) latexname = 'RPH_M0'
            else
                do j = 1, self%RPH_PlanckMass%parameter_number
                    if ( i-Nw == j ) call self%RPH_PlanckMass%parameter_names_latex( j, latexname )
                end do
            end if
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

        ! parameter from RPH_PlanckMass function
        else if ( i <= NMo ) then
            if ( self%RPHusealphaM ) then
                do j = 1, self%RPH_alphaM_ODE%parameter_number
                    if ( i-NT == j ) call self%RPH_alphaM_ODE%parameter_names_latex( j, latexname )
                end do
            else
                do j = 1, self%RPH_PlanckMass_ODE%parameter_number
                    if ( i-NT == j ) call self%RPH_PlanckMass_ODE%parameter_names_latex( j, latexname )
                end do
            end if
            return

        ! parameter from RPH_Kineticity function
        else if ( i <= NKo) then
            do j = 1, self%RPH_Kineticity_ODE%parameter_number
                if ( i-NMo == j ) call self%RPH_Kineticity_ODE%parameter_names_latex( j, latexname )
            end do
            return

        !parameter from RPH_Braiding function
        else if ( i <= NBo ) then
            do j = 1, self%RPH_Braiding_ODE%parameter_number
                if ( i-NKo == j ) call self%RPH_Braiding_ODE%parameter_names_latex( j, latexname )
            end do
            return

        !parameter from RPH_Tensor function
        else if ( i <= NTo ) then

            do j = 1, self%RPH_Tensor_ODE%parameter_number
                if ( i-NBo == j ) call self%RPH_Tensor_ODE%parameter_names_latex( j, latexname )
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

        integer  :: Nw, NM, NK, NB, NT, NMo, NKo, NBo, NTo
        integer  :: j

        ! compute the incremental number of parameters:
        Nw = self%RPH_wDE%parameter_number
        if ( self%RPHusealphaM ) then
            NM = Nw + self%RPH_alphaM%parameter_number + 1
        else
            NM = Nw + self%RPH_PlanckMass%parameter_number
        end if
        NK = NM + self%RPH_Kineticity%parameter_number
        NB = NK + self%RPH_Braiding%parameter_number
        NT = NB + self%RPH_Tensor%parameter_number

        if ( self%RPHusealphaM ) then
            NMo = NT + self%RPH_alphaM_ODE%parameter_number
        else
            NMo = NT + self%RPH_Tensor_ODE%parameter_number
        end if
        NKo = NMo + self%RPH_Kineticity_ODE%parameter_number
        NBo = NKo + self%RPH_Braiding_ODE%parameter_number
        NTo = NBo + self%RPH_Tensor_ODE%parameter_number

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
            if ( self%RPHusealphaM ) then
                do j = 1, self%RPH_alphaM%parameter_number
                    if ( i-Nw == j ) call self%RPH_alphaM%parameter_value( j, value )
                end do
                if ( i == NM ) value = self%RPH_M0
            else
                do j = 1, self%RPH_PlanckMass%parameter_number
                    if ( i-Nw == j ) call self%RPH_PlanckMass%parameter_value( j, value )
                end do
            end if
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
        
        ! parameter from RPH_PlanckMass function
        else if ( i <= NMo ) then
            if ( self%RPHusealphaM ) then
                do j = 1, self%RPH_alphaM_ODE%parameter_number
                    if ( i-NT == j ) call self%RPH_alphaM_ODE%parameter_value( j, value )
                end do
            else
                do j = 1, self%RPH_PlanckMass_ODE%parameter_number
                    if ( i-NT == j ) call self%RPH_PlanckMass_ODE%parameter_value( j, value )
                end do
            end if
            return

        ! parameter from RPH_Kineticity function
        else if ( i <= NKo) then
            do j = 1, self%RPH_Kineticity_ODE%parameter_number
                if ( i-NMo == j ) call self%RPH_Kineticity_ODE%parameter_value( j, value )
            end do
            return

        !parameter from RPH_Braiding function
        else if ( i <= NBo ) then
            do j = 1, self%RPH_Braiding_ODE%parameter_number
                if ( i-NKo == j ) call self%RPH_Braiding_ODE%parameter_value( j, value )
            end do
            return

        !parameter from RPH_Tensor function
        else if ( i <= NTo ) then

            do j = 1, self%RPH_Tensor_ODE%parameter_number
                if ( i-NBo == j ) call self%RPH_Tensor_ODE%parameter_value( j, value )
            end do
            return

        end if

    end subroutine EFTCAMBRPHParameterValues

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that initializes the background.
    subroutine EFTCAMBRPHInitBackground( self, params_cache, feedback_level, success, outroot )

        implicit none

        class(EFTCAMB_RPH)                    :: self           !< the base class
        type(TEFTCAMB_parameter_cache), intent(inout) :: params_cache   !< a EFTCAMB parameter cache containing cosmological parameters
        integer                       , intent(in)    :: feedback_level !< level of feedback from the background code. 0=none; 1=some; 2=chatty.
        logical                       , intent(out)   :: success        !< wether the background initialization succeded or not
        character(LEN=*), optional    , intent(in)    :: outroot        !< the output root for the debug files

        ! some feedback:
        if ( feedback_level>0 ) then
        write(*,'(a)') "***************************************************************"
        write(*,'(a)') ' EFTCAMB designer RPH background solver'
        write(*,'(a)')
        end if

        call self%rhov_t%initialize( self%designer_num_points, self%x_initial, self%x_final )
        call self%adotoa%initialize( self%designer_num_points, self%x_initial, self%x_final )
        call self%Hdot%initialize( self%designer_num_points, self%x_initial, self%x_final )
        call self%Hdot2%initialize( self%designer_num_points, self%x_initial, self%x_final )
        call self%Hdot3%initialize( self%designer_num_points, self%x_initial, self%x_final )
        call self%Hdot4%initialize( self%designer_num_points, self%x_initial, self%x_final )
        if ( self%RPHusealphaM ) then
            call self%RPH_M%initialize( self%designer_num_points, self%x_initial, self%x_final )
            end if
        if ( self%RPHhasOmegaDE ) then
            call self%Ode%initialize( self%designer_num_points, self%x_initial, self%x_final )
            call self%dOde%initialize( self%designer_num_points, self%x_initial, self%x_final )
            call self%d2Ode%initialize( self%designer_num_points, self%x_initial, self%x_final )
            call self%d3Ode%initialize( self%designer_num_points, self%x_initial, self%x_final )
            call self%d4Ode%initialize( self%designer_num_points, self%x_initial, self%x_final )
        end if

        success = .True.
        ! solve the background equations and store the solution:
        call self%solve_designer_equations( params_cache, success, feedback_level )

    end subroutine EFTCAMBRPHInitBackground

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that solves for M by integrating alpha_M and rho_DE
    subroutine EFTCAMBRPHSolveDesignerEquations( self, params_cache, success, feedback_level )

        implicit none

        class(EFTCAMB_RPH)                   :: self          !< the base class.
        type(TEFTCAMB_parameter_cache), intent(in)   :: params_cache  !< a EFTCAMB parameter cache containing cosmological parameters.
        logical , intent(out)                        :: success       !< whether the calculation ended correctly or not
        integer , intent(in)                         :: feedback_level   !< whether be noisy

        integer :: num_eq       !<  Number of equations
        real(dl), allocatable :: y(:), ydot(:)

        real(dl) :: grhom0, grhor0, grhov0

        integer  :: itol, itask, istate, iopt, LRN, LRS, LRW, LIS, LIN, LIW, JacobianMode, i
        real(dl) :: rtol, atol, t1, t2
        real(dl), allocatable :: rwork(:)
        integer,  allocatable :: iwork(:)

        ! This routine integrates rho_DE and the \log (1+M) = \log Meff2 = \int d\log a \alpha_M if using alphaM
        if ( self%RPHusealphaM ) then
            num_eq = 2
        else
            num_eq = 1
        end if

        allocate(y(num_eq))
        allocate(ydot(num_eq))

        ! 1) Cosmological parameters:
        grhom0 = params_cache%grhob + params_cache%grhoc
        grhor0 = 3._dl * params_cache%h0_Mpc**2 * (params_cache%omegag + params_cache%omegar)
        grhov0 = params_cache%grhov

        ! 2) Initial values
        y(1) = 1._dl
        if ( self%RPHusealphaM ) y(2) = log(1._dl + self%RPH_M0)
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

        if (self%RPHintegratefromtoday) then
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
        else
            ! first integrate DE to find initial condition
            t1 = self%rhov_t%x(1)
            t2 = self%rhov_t%x(self%designer_num_points)
            call DLSODA ( derivs, num_eq, y, t1, t2, itol, rtol, atol, itask, istate, iopt, RWORK, LRW, IWORK, LIW, jacobian, JacobianMode)
            if ( istate < 0 ) then
                if ( istate == -1 ) then
                    if ( feedback_level>1 ) write(*,*) 'w_DE integration: DLSODA excessive work'
                    istate = 1
                else
                    success = .False.
                    write(*,*) 'w_DE integration: DLSODA failed with code:', istate
                    return
                end if
            end if
            ! set the initial condition again
            y(1) = 1._dl/y(1) ! ensuring that rhov_t/rhov_0 = 1 at a=1
            if ( self%RPHusealphaM ) y(2) = log(1._dl + self%RPH_M0)
            
            ! reinitialize DLSODA
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
            deallocate(rwork)
            deallocate(iwork)
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
            ! Jacobian mode: 1=fullJacobian, 2=not provided
            JacobianMode = 2
            
            ! initial time:
            t1  = self%rhov_t%x( 1 )
            ! store initial step:
            call output( num_eq, 1, t1, y)
            ! solve the equations:
            do i=1, self%designer_num_points-1
                ! set the time step:
                t1 = self%rhov_t%x(i)
                t2 = self%rhov_t%x(i+1)
                ! solve the system:
                call DLSODA ( derivs, num_eq, y, t1, t2, itol, rtol, atol, itask, istate, iopt, RWORK, LRW, IWORK, LIW, jacobian, JacobianMode)
                ! check istate for LSODA good completion:
                if ( istate < 0 ) then
                    if ( istate == -1 ) then
                        if ( feedback_level>1 ) write(*,*) 'DLSODA excessive work'
                        istate = 1
                    else
                        success = .False.
                        write(*,*) t1, t2, i
                        write(*,*) 'DLSODA failed with code:', istate
                        return
                    end if
                end if
                ! compute output functions if needed:
                call output( num_eq, i+1, t2, y )
            end do
        end if

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

            real(dl) :: a, Ode
            real(dl) :: grhor, grhom, grhov, grhon
            real(dl) :: grhonu, gpinu, grhormass_t
            integer  :: nu_i
            
            a = Exp(x)

            ydot(1) = -3._dl*(1._dl + self%RPH_wDE%value(a))*y(1)
            
            if ( self%RPHusealphaM ) then
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

                Ode = grhov/(grhom + grhor + grhon + grhov)
                ! if (Abs(Ode) < 1d-9) Ode = 0._dl

                ydot(2) = self%RPH_alphaM%value(a) + self%RPH_alphaM_ODE%value(Ode)
            end if

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
        !> Subroutine that takes the solution of the background equations and stores the values of the EFT functions.
        subroutine output( num_eq, ind, x, y)

            implicit none

            integer , intent(in)                     :: num_eq !< number of equations in the ODE system.
            integer , intent(in)                     :: ind    !< index of the EFT functions interpolation tables to fill.
            real(dl), intent(in)                     :: x      !< time at which the derivatives of the system are computed.
            real(dl), intent(in) , dimension(num_eq) :: y      !< input status of the system.

            real(dl) :: a
            real(dl) :: wDE, dwDE, d2wDE, d3wDE, d4wDE
            real(dl) :: Ht, Ht2, Ht3, Ht4, H, dH, d2H, d3H, d4H
            real(dl) :: grhor, grhom, grhov, grhon, gpn, gdpn, gddpn, gdddpn
            real(dl) :: grhonu, gpinu, gpinudot, gpinudotdot, gpinudotdotdot, grhormass_t
            real(dl) :: Ode, dOde, d2Ode, d3Ode, d4Ode
            integer  :: nu_i
            
            a = Exp(x)

            ! short vars for wDE
            wDE   = self%RPH_wDE%value(a)
            dwDE  = self%RPH_wDE%first_derivative(a)
            d2wDE = self%RPH_wDE%second_derivative(a)
            d3wDE = self%RPH_wDE%third_derivative(a)

            !1) conformal Hubble and derivatives
            grhom = grhom0 / a
            grhor = grhor0 / a**2
            grhov = grhov0 * y(1) * a**2
            grhon = 0._dl
            gpn = 0._dl
            gdpn = 0._dl
            gddpn = 0._dl
            gdddpn = 0._dl
            if ( params_cache%Num_Nu_Massive /= 0 ) then
                do nu_i = 1, params_cache%Nu_mass_eigenstates
                    grhonu      = 0._dl
                    gpinu       = 0._dl
                    grhormass_t = params_cache%grhormass(nu_i)/a**2
                    call ThermalNuBack%rho_P(a*params_cache%nu_masses(nu_i), grhonu, gpinu )
                    grhon = grhon +grhormass_t*grhonu
                    gpn   = gpn  +grhormass_t*gpinu
                end do
            end if
            H = sqrt( ( grhom + grhor + grhov + grhon )/3._dl )
            Ht = -0.5_dl*( H**2 + grhor/3._dl + grhov*wDE + gpn )
            ! Ht2
            grhon = 0._dl
            gpn = 0._dl
            gdpn = 0._dl
            gddpn = 0._dl
            gdddpn = 0._dl
            if ( params_cache%Num_Nu_Massive /= 0 ) then
                do nu_i = 1, params_cache%Nu_mass_eigenstates
                    grhonu      = 0._dl
                    gpinu       = 0._dl
                    gpinudot    = 0._dl
                    grhormass_t = params_cache%grhormass(nu_i)/a**2
                    call ThermalNuBack%rho_P(a*params_cache%nu_masses(nu_i), grhonu, gpinu )
                    grhon    = grhon +grhormass_t*grhonu
                    gpn      = gpn  +grhormass_t*gpinu
                    gpinudot = ThermalNuBack%pidot( a*params_cache%nu_masses(nu_i), H, gpinu )
                    gdpn     = gdpn  + grhormass_t*( gpinudot - 4._dl*H*gpinu)
                end do
            end if
            Ht2 = H*( grhom/6._dl + 2._dl*grhor/3._dl ) &
            & + H*grhov*( 1._dl/6._dl +wDE +1.5_dl*wDE**2 -0.5_dl*dwDE ) &
            & + H*grhon/6._dl -0.5_dl*H*gpn -0.5_dl*gdpn
            ! Ht3
            grhon = 0._dl
            gpn = 0._dl
            gdpn = 0._dl
            gddpn = 0._dl
            gdddpn = 0._dl
            if ( params_cache%Num_Nu_Massive /= 0 ) then
                do nu_i = 1, params_cache%Nu_mass_eigenstates
                    grhonu      = 0._dl
                    gpinu       = 0._dl
                    gpinudot    = 0._dl
                    gpinudotdot = 0._dl
                    grhormass_t = params_cache%grhormass(nu_i)/a**2
                    call ThermalNuBack%rho_P(a*params_cache%nu_masses(nu_i), grhonu, gpinu )
                    grhon       = grhon +grhormass_t*grhonu
                    gpn         = gpn  +grhormass_t*gpinu
                    gpinudot    = ThermalNuBack%pidot( a*params_cache%nu_masses(nu_i), H, gpinu )
                    gdpn        = gdpn  + grhormass_t*( gpinudot - 4._dl*H*gpinu)
                    gpinudotdot = ThermalNuBack%pidotdot( a*params_cache%nu_masses(nu_i), H, Ht, gpinu, gpinudot)
                    gddpn       = gddpn - 4._dl*H*grhormass_t*( gpinudot - 4._dl*H*gpinu )+ grhormass_t*( gpinudotdot - 4._dl*Ht*gpinu - 4._dl*H*gpinudot )
                end do
            end if
            Ht3 = H**2*( -grhom/6._dl - 4._dl*grhor/3._dl &
                & + grhov*(-1._dl/6._dl - 1.5_dl*wDE - 4.5_dl*wDE**2 - 4.5_dl*wDE**3 + dwDE + 4.5_dl*wDE*dwDE - 0.5_dl*d2wDE) &
                & -1._dl/6._dl*grhon - 1.5_dl*gpn ) &
                & + Ht*( grhom/6._dl + 2._dl/3._dl*grhor &
                & + grhov*(1._dl/6._dl + wDE +1.5_dl*wDE**2 - 0.5_dl*dwDE) &
                & + grhon/6._dl - 0.5_dl*gpn ) &
                & -1.5_dl*H*gdpn - 0.5_dl*gddpn
            ! Ht4
            grhon = 0._dl
            gpn = 0._dl
            gdpn = 0._dl
            gddpn = 0._dl
            gdddpn = 0._dl
            if ( params_cache%Num_Nu_Massive /= 0 ) then
                do nu_i = 1, params_cache%Nu_mass_eigenstates
                    grhonu         = 0._dl
                    gpinu          = 0._dl
                    gpinudot       = 0._dl
                    gpinudotdot    = 0._dl
                    gpinudotdotdot = 0._dl
                    grhormass_t = params_cache%grhormass(nu_i)/a**2
                    call ThermalNuBack%rho_P(a*params_cache%nu_masses(nu_i), grhonu, gpinu )
                    grhon          = grhon +grhormass_t*grhonu
                    gpn            = gpn  +grhormass_t*gpinu
                    gpinudot       = ThermalNuBack%pidot( a*params_cache%nu_masses(nu_i), H, gpinu )
                    gdpn           = gdpn  + grhormass_t*( gpinudot - 4._dl*H*gpinu)
                    gpinudotdot    = ThermalNuBack%pidotdot( a*params_cache%nu_masses(nu_i), H, Ht, gpinu, gpinudot)
                    gddpn          = gddpn - 4._dl*H*grhormass_t*( gpinudot - 4._dl*H*gpinu )+ grhormass_t*( gpinudotdot - 4._dl*Ht*gpinu - 4._dl*H*gpinudot )
                    gpinudotdotdot = ThermalNuBack%pidotdotdot(a*params_cache%nu_masses(nu_i),H,Ht,Ht2,gpinu,gpinudot,gpinudotdot)
                    gdddpn = gdddpn + grhormass_t*(gpinudotdotdot &
                        & - 12._dl*H*gpinudotdot &
                        & + (48._dl*H**2 - 12._dl*Ht)*gpinudot &
                        & + (-64._dl*H**3 + 48._dl*H*Ht - 4._dl*Ht2)*gpinu )
                end do
            end if
            Ht4 = H**3*( &
                & grhom/6._dl + 8._dl/3._dl*grhor &
                & + grhov*( 1._dl/6._dl + 2._dl*wDE + 9._dl*wDE**2 + 18._dl*wDE**3 + 13.5_dl*wDE**4 - 1.5_dl*dwDE - 12._dl*wDE*dwDE - 27._dl*wDE**2*dwDE + 4.5_dl*dwDE**2 + 0.5_dl*d2wDE + 6._dl*wDE*d2wDE - 0.5_dl*d3wDE ) &
                & + grhon/6._dl - 2.5_dl*gpn &
                & ) &
                & + H*Ht*( &
                & -0.5_dl*grhom - 4._dl*grhor &
                & + grhov*( -0.5_dl - 4.5_dl*wDE - 13.5_dl*wDE**2 - 13.5_dl*wDE**3 + 3._dl*dwDE +13.5_dl*wDE*dwDE - 1.5_dl*d2wDE ) &
                & -0.5_dl*grhon - 4.5_dl*gpn &
                & ) &
                & + Ht2*( &
                & grhom/6._dl + 2._dl/3._dl*grhor &
                & + grhov*( 1._dl/6._dl + wDE + 1.5_dl*wDE**2 - 0.5_dl*dwDE ) &
                & + grhon/6._dl - 0.5_dl*gpn &
                & ) &
                & + gdpn*( -2._dl*Ht - 4.5_dl*H**2 ) - 2.5_dl*H*gddpn - 0.5_dl*gdddpn

            if ( self%RPHhasOmegaDE ) then
                ! if (Abs(Ode) > 1d-9) then
                    ! short vars for dnH = d^n H / da^n
                    dH  = Ht/(a*H)
                    d2H = -((H**2*Ht + Ht**2 - H*Ht2)/(a**2*H**3))
                    d3H = (2._dl*H**4*Ht + 3._dl*Ht**3 - 3._dl*H**3*Ht2 - 4._dl*H*Ht*Ht2 + H**2*(3._dl*Ht**2 + Ht3))/(a**3*H**5)
                    d4H = (-6._dl*H**6*Ht - 15._dl*Ht**4 + 11._dl*H**5*Ht2 + 25._dl*H*Ht**2*Ht2 - H**4*(11._dl*Ht**2 + 6._dl*Ht3) - H**2*(18._dl*Ht**3 + 4._dl*Ht2**2 + 7._dl*Ht*Ht3) + H**3*(24._dl*Ht*Ht2 + Ht4))/(a**4*H**7)

                    ! Omega_DE and derivatives wrt a
                    Ode    = grhov / 3._dl / H**2
                    dOde   = -((Ode*(2._dl*a*dH + H + 3._dl*H*wDE))/(a*H))
                    d2Ode  = (-2._dl*(dH*dOde*H - dH**2*Ode + d2H*H*Ode))/H**2 - (dOde + 3._dl*dwDE*Ode + 3._dl*dOde*wDE)/a + (Ode + 3._dl*Ode*wDE)/a**2
                    d3Ode  = (-4._dl*dH**3*Ode)/H**3 + (2._dl*dH*(2._dl*dH*dOde + 3._dl*d2H*Ode))/H**2 - (2._dl*(d2Ode*dH + 2._dl*d2H*dOde + d3H*Ode))/H + (-(a**2*(d2Ode + 6._dl*dOde*dwDE + 3._dl*d2wDE*Ode + 3._dl*d2Ode*wDE)) + 2._dl*a*(dOde + 3._dl*dwDE*Ode + 3._dl*dOde*wDE) - 2._dl*(Ode + 3._dl*Ode*wDE))/a**3
                    d4Ode  = (12._dl*dH**4*Ode)/H**4 - (12._dl*dH**2*(dH*dOde + 2._dl*d2H*Ode))/H**3 - (2._dl*(3._dl*d2H*d2Ode + d3Ode*dH + 3._dl*d3H*dOde + d4H*Ode))/H + (6._dl*d2Ode*dH**2 + 18._dl*d2H*dH*dOde + 6._dl*d2H**2*Ode + 8._dl*d3H*dH*Ode)/H**2 + (3._dl*a**2*(d2Ode + 6._dl*dOde*dwDE + 3._dl*d2wDE*Ode + 3._dl*d2Ode*wDE) - a**3*(d3Ode + 9._dl*d2wDE*dOde + 9._dl*d2Ode*dwDE + 3._dl*d3wDE*Ode + 3._dl*d3Ode*wDE) - 6._dl*a*(dOde + 3._dl*dwDE*Ode + 3._dl*dOde*wDE) + 6._dl*(Ode + 3._dl*Ode*wDE))/a**4
                ! else
                !     Ode   = 0._dl
                !     dOde  = 0._dl
                !     d2Ode = 0._dl
                !     d3Ode = 0._dl
                !     d4Ode = 0._dl
                ! end if
            end if

            ! output
            self%adotoa%y( ind ) = H
            self%Hdot%y( ind )   = Ht
            self%Hdot2%y( ind )  = Ht2
            self%Hdot3%y( ind )  = Ht3
            self%Hdot4%y( ind )  = Ht4
            self%rhov_t%y( ind ) = grhov0 * y(1)
            if ( self%RPHusealphaM ) then
                self%RPH_M%y( ind ) = expm1(y(2))
            end if
            if ( self%RPHhasOmegaDE ) then
                self%Ode%y( ind )   = Ode
                self%dOde%y( ind )  = dOde
                self%d2Ode%y( ind ) = d2Ode
                self%d3Ode%y( ind ) = d3Ode
                self%d4Ode%y( ind ) = d4Ode
            end if
        
        end subroutine output

    end subroutine EFTCAMBRPHSolveDesignerEquations

    ! function from Toshio Fukushima 10.13140/RG.2.2.14468.48001
    function expm1( x )

        implicit none
  
        real(dl), intent(in)  :: x
        real(dl) :: s, c_, e
        real(dl) :: expm1
  
        if (abs(x) .le. 0.6931471805599453094172321214581765680755d0) then
           !
           !   degree 13 minimax polynomial approximation
           !
           expm1=x*(0.99999999999999999301d0 &
                 +x*(0.49999999999999957237d0 &
                 +x*(0.16666666666666769423d0 &
                 +x*(0.041666666666691966524d0 &
                 +x*(0.0083333333333092372395d0 &
                 +x*(0.0013888888884629134308d0 &
                 +x*(0.00019841269861419530922d0 &
                 +x*(0.000024801590367084381006d0 &
                 +x*(2.7557312010044843920d-6 &
                 +x*(2.7556249037054640535d-7 &
                 +x*(2.5053112146542602098d-8 &
                 +x*(2.1055737595306817199d-9 &
                 +x*(1.6058956553927382096d-10)))))))))))))
        else
           !
           !   standard approach using "sinh"
           !
           s=sinh(0.5d0*x)
           c_=sqrt(s*s+1.d0)
           e=s+c_
           expm1=2.d0*e*s
        end if
  
    end function expm1

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes the value of the background EFT functions at a given time.
    subroutine EFTCAMBRPHBackgroundEFTFunctions( self, a, eft_par_cache, eft_cache )

        implicit none

        class(EFTCAMB_RPH)                           :: self          !< the base class
        real(dl), intent(in)                         :: a             !< the input scale factor
        type(TEFTCAMB_parameter_cache), intent(inout) :: eft_par_cache !< the EFTCAMB parameter cache that contains all the physical parameters.
        type(TEFTCAMB_timestep_cache ), intent(inout) :: eft_cache     !< the EFTCAMB timestep cache that contains all the physical values.

        real(dl) :: RPH_PM_V, RPH_AT_V, RPH_PM_P, RPH_AT_P, RPH_PM_PP, RPH_AT_PP, RPH_PM_PPP, RPH_AT_PPP, RPH_PM_PPPP, RPH_AT_PPPP
        real(dl) :: aM, daM, d2aM, d3aM, Mp1
        real(dl) :: Ode, dOde, d2Ode, d3Ode, d4Ode

        real(dl) :: x, mu
        integer  :: ind

        x = log(a)
        call self%adotoa%precompute(x, ind, mu)
        
        ! Omega_DE
        if ( self%RPHhasOmegaDE ) then
            Ode   = self%Ode%value( x, ind, mu )
            dOde  = self%dOde%value( x, ind, mu )
            d2Ode = self%d2Ode%value( x, ind, mu )
            d3Ode = self%d3Ode%value( x, ind, mu )
            d4Ode = self%d4Ode%value( x, ind, mu )
        else
            Ode   = 0._dl
            dOde  = 0._dl
            d2Ode = 0._dl
            d3Ode = 0._dl
            d4Ode = 0._dl
        end if

        ! alphaM and derivatives:
        if ( self%RPHusealphaM ) then
            aM    = self%RPH_alphaM%value(a) + self%RPH_alphaM_ODE%value(Ode)
            daM   = self%RPH_alphaM%first_derivative(a) + (dOde*self%RPH_alphaM_ODE%first_derivative(Ode))
            d2aM  = self%RPH_alphaM%second_derivative(a) + (d2Ode*self%RPH_alphaM_ODE%first_derivative(Ode) + dOde**2*self%RPH_alphaM_ODE%second_derivative(Ode))
            d3aM  = self%RPH_alphaM%third_derivative(a) + (d3Ode*self%RPH_alphaM_ODE%first_derivative(Ode) + 3._dl*d2Ode*dOde*self%RPH_alphaM_ODE%second_derivative(Ode) + dOde**3*self%RPH_alphaM_ODE%third_derivative(Ode))
            
            RPH_PM_V    = self%RPH_M%value(x, ind, mu)
            Mp1         = RPH_PM_V + 1._dl
            RPH_PM_P    = (aM*Mp1)/a
            RPH_PM_PP   = ((-aM + aM**2 + a*daM)*Mp1)/a**2
            RPH_PM_PPP  = ((-3._dl*aM**2 + aM**3 + a*(a*d2aM - 2._dl*daM) + aM*(2._dl + 3._dl*a*daM))*Mp1)/a**3
            RPH_PM_PPPP = ((-6._dl*aM**3 + aM**4 + 2._dl*aM*(-3._dl + 2._dl*a**2*d2aM - 7._dl*a*daM) + aM**2*(11._dl + 6._dl*a*daM) + a*(-3._dl*a*d2aM + a**2*d3aM + 6._dl*daM + 3._dl*a*daM**2))*Mp1)/a**4
        else
            RPH_PM_V    = self%RPH_PlanckMass%value(a) + self%RPH_PlanckMass_ODE%value(Ode)
            Mp1         = RPH_PM_V + 1._dl
            RPH_PM_P    = self%RPH_PlanckMass%first_derivative(a) + (dOde*self%RPH_PlanckMass_ODE%first_derivative(Ode))
            RPH_PM_PP   = self%RPH_PlanckMass%second_derivative(a) + (d2Ode*self%RPH_PlanckMass_ODE%first_derivative(Ode) + dOde**2*self%RPH_PlanckMass_ODE%second_derivative(Ode))
            RPH_PM_PPP  = self%RPH_PlanckMass%third_derivative(a) + (d3Ode*self%RPH_PlanckMass_ODE%first_derivative(Ode) + 3._dl*d2Ode*dOde*self%RPH_PlanckMass_ODE%second_derivative(Ode) + dOde**3*self%RPH_PlanckMass_ODE%third_derivative(Ode))
            RPH_PM_PPPP = self%RPH_PlanckMass%fourth_derivative(a) + (d4Ode*self%RPH_PlanckMass_ODE%first_derivative(Ode) + 3._dl*d2Ode**2*self%RPH_PlanckMass_ODE%second_derivative(Ode) + 4._dl*d3Ode*dOde*self%RPH_PlanckMass_ODE%second_derivative(Ode) + 6._dl*d2Ode*dOde**2*self%RPH_PlanckMass_ODE%third_derivative(Ode) + dOde**4*self%RPH_PlanckMass_ODE%fourth_derivative(Ode))
        end if

        ! alphaKBT and derivatives
        RPH_AT_V               = self%RPH_Tensor%value(a) + self%RPH_Tensor_ODE%value(Ode)
        RPH_AT_P               = self%RPH_Tensor%first_derivative(a) + (dOde*self%RPH_Tensor_ODE%first_derivative(Ode))
        RPH_AT_PP              = self%RPH_Tensor%second_derivative(a) + (d2Ode*self%RPH_Tensor_ODE%first_derivative(Ode) + dOde**2*self%RPH_Tensor_ODE%second_derivative(Ode))
        RPH_AT_PPP             = self%RPH_Tensor%third_derivative(a) + (d3Ode*self%RPH_Tensor_ODE%first_derivative(Ode) + 3._dl*d2Ode*dOde*self%RPH_Tensor_ODE%second_derivative(Ode) + dOde**3*self%RPH_Tensor_ODE%third_derivative(Ode))
        RPH_AT_PPPP            = self%RPH_Tensor%fourth_derivative(a) + (d4Ode*self%RPH_Tensor_ODE%first_derivative(Ode) + 3._dl*d2Ode**2*self%RPH_Tensor_ODE%second_derivative(Ode) + 4._dl*d3Ode*dOde*self%RPH_Tensor_ODE%second_derivative(Ode) + 6._dl*d2Ode*dOde**2*self%RPH_Tensor_ODE%third_derivative(Ode) + dOde**4*self%RPH_Tensor_ODE%fourth_derivative(Ode))
        
        ! compute the EFT functions:
        eft_cache%EFTOmegaV    = RPH_PM_V + RPH_AT_V*Mp1
        eft_cache%EFTOmegaP    = RPH_PM_P + RPH_AT_V*RPH_PM_P + RPH_AT_P*Mp1
        eft_cache%EFTOmegaPP   = RPH_PM_PP + RPH_AT_V*RPH_PM_PP + 2._dl*RPH_AT_P*RPH_PM_P + RPH_AT_PP*Mp1
        eft_cache%EFTOmegaPPP  = RPH_PM_PPP + RPH_AT_V*RPH_PM_PPP + 3._dl*RPH_PM_PP*RPH_AT_P + 3._dl*RPH_AT_PP*RPH_PM_P + RPH_AT_PPP*Mp1
        eft_cache%EFTOmegaPPPP = 6._dl*RPH_AT_PP*RPH_PM_PP + RPH_PM_PPPP + RPH_AT_V*RPH_PM_PPPP + 4._dl*RPH_PM_PPP*RPH_AT_P + 4._dl*RPH_AT_PPP*RPH_PM_P + RPH_AT_PPPP*Mp1
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
        eft_cache%EFTcdotdot   = (-(a**4*eft_cache%EFTOmegaPPPP*eft_cache%adotoa**4) - 6._dl*a**3*eft_cache%EFTOmegaPPP*eft_cache%adotoa**2*eft_cache%Hdot + a**2*eft_cache%EFTOmegaPP*(eft_cache%adotoa**4 + 4._dl*eft_cache%adotoa**2*eft_cache%Hdot - 3._dl*eft_cache%Hdot**2 - 4._dl*eft_cache%adotoa*eft_cache%Hdotdot) + a*eft_cache%EFTOmegaP*(-5._dl*eft_cache%adotoa**4 + 10._dl*eft_cache%adotoa**2*eft_cache%Hdot + eft_cache%Hdot**2 - eft_cache%Hdotdotdot) + 2._dl*(4._dl*eft_cache%adotoa**4 - 14._dl*eft_cache%adotoa**2*eft_cache%Hdot + 4._dl*eft_cache%Hdot**2 + 6._dl*eft_cache%adotoa*eft_cache%Hdotdot - eft_cache%Hdotdotdot)*eft_cache%EFTOmegaV + eft_cache%grhov_t*eft_cache%Hdot*(a*self%RPH_wDE%first_derivative(a) - 3._dl*(1 + self%RPH_wDE%value(a))**2) + eft_cache%grhov_t*eft_cache%adotoa**2*(a**2*self%RPH_wDE%second_derivative(a) + 9._dl*(1 + self%RPH_wDE%value(a))**3 - a*self%RPH_wDE%first_derivative(a)*(8._dl + 9._dl*self%RPH_wDE%value(a))))/2._dl
        eft_cache%EFTLambdadot = -2._dl*eft_cache%EFTOmegaV*( eft_cache%Hdotdot -eft_cache%adotoa*eft_cache%Hdot -eft_cache%adotoa**3 ) &
            & -a*eft_cache%EFTOmegaP*( +eft_cache%Hdotdot +5._dl*eft_cache%adotoa*eft_cache%Hdot -eft_cache%adotoa**3  ) &
            & -a**2*eft_cache%EFTOmegaPP*eft_cache%adotoa*( +2._dl*eft_cache%adotoa**2 +3._dl*eft_cache%Hdot )&
            & -(a*eft_cache%adotoa)**3*eft_cache%EFTOmegaPPP &
            & +eft_cache%grhov_t*eft_cache%adotoa*( a*self%RPH_wDE%first_derivative(a) -3._dl*self%RPH_wDE%value(a)*(1._dl +self%RPH_wDE%value(a) ))
        eft_cache%EFTLambdadotdot = -(a**4*eft_cache%EFTOmegaPPPP*eft_cache%adotoa**4) + eft_cache%EFTOmegaP*(-(a*eft_cache%Hdotdotdot) - 5._dl*a*eft_cache%Hdot**2 - 6._dl*a*eft_cache%Hdotdot*eft_cache%adotoa + 10._dl*a*eft_cache%Hdot*eft_cache%adotoa**2 + a*eft_cache%adotoa**4) + eft_cache%EFTOmegaPP*(-3._dl*a**2*eft_cache%Hdot**2 - 4._dl*a**2*eft_cache%Hdotdot*eft_cache%adotoa - 11._dl*a**2*eft_cache%Hdot*eft_cache%adotoa**2 + a**2*eft_cache%adotoa**4) + eft_cache%EFTOmegaPPP*(-6._dl*a**3*eft_cache%Hdot*eft_cache%adotoa**2 - 3._dl*a**3*eft_cache%adotoa**4) + (-2._dl*eft_cache%Hdotdotdot + 2._dl*eft_cache%Hdot**2 + 6._dl*eft_cache%Hdotdot*eft_cache%adotoa + 2._dl*eft_cache%Hdot*eft_cache%adotoa**2 - 4._dl*eft_cache%adotoa**4)*eft_cache%EFTOmegaV + eft_cache%Hdot*eft_cache%grhov_t*(a*self%RPH_wDE%first_derivative(a) - 3._dl*self%RPH_wDE%value(a) - 3._dl*self%RPH_wDE%value(a)**2) + eft_cache%adotoa**2*eft_cache%grhov_t*(a**2*self%RPH_wDE%second_derivative(a) + 9._dl*self%RPH_wDE%value(a) + 18._dl*self%RPH_wDE%value(a)**2 + 9._dl*self%RPH_wDE%value(a)**3 - a*self%RPH_wDE%first_derivative(a)*(5._dl + 9._dl*self%RPH_wDE%value(a)))

    end subroutine EFTCAMBRPHBackgroundEFTFunctions

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes the value of the background EFT functions at a given time.
    subroutine EFTCAMBRPHSecondOrderEFTFunctions( self, a, eft_par_cache, eft_cache )

        implicit none

        class(EFTCAMB_RPH)                           :: self          !< the base class
        real(dl), intent(in)                         :: a             !< the input scale factor
        type(TEFTCAMB_parameter_cache), intent(inout) :: eft_par_cache !< the EFTCAMB parameter cache that contains all the physical parameters.
        type(TEFTCAMB_timestep_cache ), intent(inout) :: eft_cache     !< the EFTCAMB timestep cache that contains all the physical values.

        real(dl) :: RPH_PM_V, RPH_PM_P, RPH_PM_PP, RPH_PM_PPP, RPH_PM_PPPP, RPH_AT_V, RPH_AT_P, RPH_AT_PP, RPH_AT_PPP, RPH_AT_PPPP, RPH_AK_V, RPH_AK_P, RPH_AK_PP, RPH_AB_V, RPH_AB_P, RPH_AB_PP, RPH_AB_PPP
        real(dl) :: aM, daM, d2aM, d3aM, Mp1
        real(dl) :: Ode, dOde, d2Ode, d3Ode, d4Ode

        real(dl) :: H, Ht, Ht2, Ht3, dH, d2H, d3H

        real(dl) :: x, mu
        integer  :: ind

        x = log(a)
        call self%adotoa%precompute(x, ind, mu)

        ! short var for conformal Hubble
        H = eft_cache%adotoa
        Ht = eft_cache%Hdot
        Ht2 = eft_cache%Hdotdot
        Ht3 = eft_cache%Hdotdotdot
        dH = Ht/(a*H)
        d2H = -((H**2*Ht + Ht**2 - H*Ht2)/(a**2*H**3))
        d3H = (2._dl*H**4*Ht + 3._dl*Ht**3 - 3._dl*H**3*Ht2 - 4._dl*H*Ht*Ht2 + H**2*(3._dl*Ht**2 + Ht3))/(a**3*H**5)
        
        ! Omega_DE
        if ( self%RPHhasOmegaDE ) then
            Ode   = self%Ode%value( x, ind, mu )
            dOde  = self%dOde%value( x, ind, mu )
            d2Ode = self%d2Ode%value( x, ind, mu )
            d3Ode = self%d3Ode%value( x, ind, mu )
            d4Ode = self%d4Ode%value( x, ind, mu )
        else
            Ode   = 0._dl
            dOde  = 0._dl
            d2Ode = 0._dl
            d3Ode = 0._dl
            d4Ode = 0._dl
        end if

        ! alphaM and derivatives:
        if ( self%RPHusealphaM ) then
            aM    = self%RPH_alphaM%value(a) + self%RPH_alphaM_ODE%value(Ode)
            daM   = self%RPH_alphaM%first_derivative(a) + (dOde*self%RPH_alphaM_ODE%first_derivative(Ode))
            d2aM  = self%RPH_alphaM%second_derivative(a) + (d2Ode*self%RPH_alphaM_ODE%first_derivative(Ode) + dOde**2*self%RPH_alphaM_ODE%second_derivative(Ode))
            d3aM  = self%RPH_alphaM%third_derivative(a) + (d3Ode*self%RPH_alphaM_ODE%first_derivative(Ode) + 3._dl*d2Ode*dOde*self%RPH_alphaM_ODE%second_derivative(Ode) + dOde**3*self%RPH_alphaM_ODE%third_derivative(Ode))
            
            RPH_PM_V    = self%RPH_M%value(x, ind, mu)
            Mp1         = RPH_PM_V + 1._dl
            RPH_PM_P    = (aM*Mp1)/a
            RPH_PM_PP   = ((-aM + aM**2 + a*daM)*Mp1)/a**2
            RPH_PM_PPP  = ((-3._dl*aM**2 + aM**3 + a*(a*d2aM - 2._dl*daM) + aM*(2._dl + 3._dl*a*daM))*Mp1)/a**3
            RPH_PM_PPPP = ((-6._dl*aM**3 + aM**4 + 2._dl*aM*(-3._dl + 2._dl*a**2*d2aM - 7._dl*a*daM) + aM**2*(11._dl + 6._dl*a*daM) + a*(-3._dl*a*d2aM + a**2*d3aM + 6._dl*daM + 3._dl*a*daM**2))*Mp1)/a**4
        else
            RPH_PM_V    = self%RPH_PlanckMass%value(a) + self%RPH_PlanckMass_ODE%value(Ode)
            Mp1         = RPH_PM_V + 1._dl
            RPH_PM_P    = self%RPH_PlanckMass%first_derivative(a) + (dOde*self%RPH_PlanckMass_ODE%first_derivative(Ode))
            RPH_PM_PP   = self%RPH_PlanckMass%second_derivative(a) + (d2Ode*self%RPH_PlanckMass_ODE%first_derivative(Ode) + dOde**2*self%RPH_PlanckMass_ODE%second_derivative(Ode))
            RPH_PM_PPP  = self%RPH_PlanckMass%third_derivative(a) + (d3Ode*self%RPH_PlanckMass_ODE%first_derivative(Ode) + 3._dl*d2Ode*dOde*self%RPH_PlanckMass_ODE%second_derivative(Ode) + dOde**3*self%RPH_PlanckMass_ODE%third_derivative(Ode))
            RPH_PM_PPPP = self%RPH_PlanckMass%fourth_derivative(a) + (d4Ode*self%RPH_PlanckMass_ODE%first_derivative(Ode) + 3._dl*d2Ode**2*self%RPH_PlanckMass_ODE%second_derivative(Ode) + 4._dl*d3Ode*dOde*self%RPH_PlanckMass_ODE%second_derivative(Ode) + 6._dl*d2Ode*dOde**2*self%RPH_PlanckMass_ODE%third_derivative(Ode) + dOde**4*self%RPH_PlanckMass_ODE%fourth_derivative(Ode))
        end if

        ! alphaKBT and derivatives
        RPH_AT_V    = self%RPH_Tensor%value(a) + self%RPH_Tensor_ODE%value(Ode)
        RPH_AT_P    = self%RPH_Tensor%first_derivative(a) + (dOde*self%RPH_Tensor_ODE%first_derivative(Ode))
        RPH_AT_PP   = self%RPH_Tensor%second_derivative(a) + (d2Ode*self%RPH_Tensor_ODE%first_derivative(Ode) + dOde**2*self%RPH_Tensor_ODE%second_derivative(Ode))
        RPH_AT_PPP  = self%RPH_Tensor%third_derivative(a) + (d3Ode*self%RPH_Tensor_ODE%first_derivative(Ode) + 3._dl*d2Ode*dOde*self%RPH_Tensor_ODE%second_derivative(Ode) + dOde**3*self%RPH_Tensor_ODE%third_derivative(Ode))
        RPH_AT_PPPP = self%RPH_Tensor%fourth_derivative(a) + (d4Ode*self%RPH_Tensor_ODE%first_derivative(Ode) + 3._dl*d2Ode**2*self%RPH_Tensor_ODE%second_derivative(Ode) + 4._dl*d3Ode*dOde*self%RPH_Tensor_ODE%second_derivative(Ode) + 6._dl*d2Ode*dOde**2*self%RPH_Tensor_ODE%third_derivative(Ode) + dOde**4*self%RPH_Tensor_ODE%fourth_derivative(Ode))
        RPH_AK_V    = self%RPH_Kineticity%value(a) + self%RPH_Kineticity_ODE%value(Ode)
        RPH_AK_P    = self%RPH_Kineticity%first_derivative(a)  + (dOde*self%RPH_Kineticity_ODE%first_derivative(Ode))
        RPH_AK_PP   = self%RPH_Kineticity%second_derivative(a) + (d2Ode*self%RPH_Kineticity_ODE%first_derivative(Ode) + dOde**2*self%RPH_Kineticity_ODE%second_derivative(Ode))
        RPH_AB_V    = self%RPH_Braiding%value(a) + self%RPH_Braiding_ODE%value(Ode)
        RPH_AB_P    = self%RPH_Braiding%first_derivative(a) + (dOde*self%RPH_Braiding_ODE%first_derivative(Ode))
        RPH_AB_PP   = self%RPH_Braiding%second_derivative(a) + (d2Ode*self%RPH_Braiding_ODE%first_derivative(Ode) + dOde**2*self%RPH_Braiding_ODE%second_derivative(Ode))
        RPH_AB_PPP  = self%RPH_Braiding%third_derivative(a) + (d3Ode*self%RPH_Braiding_ODE%first_derivative(Ode) + 3._dl*d2Ode*dOde*self%RPH_Braiding_ODE%second_derivative(Ode) + dOde**3*self%RPH_Braiding_ODE%third_derivative(Ode))
        
        ! output the alternative EFT functions:
        eft_cache%Meff2     = Mp1
        eft_cache%alphaM    = a*RPH_PM_P/Mp1
        eft_cache%alphaMdot = H*(a**2*RPH_PM_PP/Mp1 + eft_cache%alphaM - eft_cache%alphaM**2)
        eft_cache%alphaB    = RPH_AB_V
        eft_cache%alphaBdot = a*H*RPH_AB_P
        eft_cache%alphaK    = RPH_AK_V
        eft_cache%alphaKdot = a*H*RPH_AK_P
        eft_cache%alphaT    = RPH_AT_V
        eft_cache%alphaTdot = a*H*RPH_AT_P

        ! compute the EFT functions:
        eft_cache%EFTGamma1V    = 0.25_dl*( RPH_AK_V*(1._dl +RPH_PM_V)*eft_cache%adotoa**2 &
            & -2._dl*eft_cache%EFTc )/(eft_par_cache%h0_Mpc**2*a**2)
        eft_cache%EFTGamma1P    = - 0.5_dl*( RPH_AK_V*(1._dl +RPH_PM_V)*eft_cache%adotoa**2 &
            & -2._dl*eft_cache%EFTc )/(eft_par_cache%h0_Mpc**2*a**3) &
            & +0.25_dl*( RPH_AK_P*(1._dl +RPH_PM_V)*eft_cache%adotoa**2 &
            & +RPH_AK_V*RPH_PM_P*eft_cache%adotoa**2 &
            & +2._dl*RPH_AK_V*(1._dl +RPH_PM_V)*eft_cache%Hdot/a &
            & -4._dl*eft_cache%EFTc/a -2._dl*eft_cache%EFTcdot/a/eft_cache%adotoa )/(eft_par_cache%h0_Mpc**2*a**2)
        eft_cache%EFTGamma1PP   = (-2._dl*eft_cache%EFTcdotdot*H + 2._dl*eft_cache%EFTcdot*(H**2 + Ht) + H**2*(2._dl*RPH_AK_V*(3*H**3 - 5._dl*H*Ht + Ht2)*Mp1 + a**2*H**3*(RPH_AK_V*RPH_PM_PP + 2._dl*RPH_AK_P*RPH_PM_P + RPH_AK_PP*Mp1) - 4._dl*a*H*(H**2 - Ht)*(RPH_AK_V*RPH_PM_P + RPH_AK_P*Mp1)))/(4._dl*a**4*H**3*eft_par_cache%h0_Mpc**2)
        eft_cache%EFTGamma2V    = ( +2._dl*RPH_AB_V*(1._dl +RPH_PM_V) &
            & -a*eft_cache%EFTOmegaP )*eft_cache%adotoa/(eft_par_cache%h0_Mpc*a)
        eft_cache%EFTGamma2P    = -( +2._dl*RPH_AB_V*(1._dl +RPH_PM_V) &
            & -a*eft_cache%EFTOmegaP )*eft_cache%adotoa/(eft_par_cache%h0_Mpc*a**2) &
            & -( -2._dl*(1._dl +RPH_PM_V)*( RPH_AB_P*eft_cache%adotoa**2 &
            & + RPH_AB_V*eft_cache%Hdot/a) &
            & - 2._dl*RPH_AB_V*eft_cache%adotoa**2*RPH_PM_P &
            & + eft_cache%EFTOmegaP*( eft_cache%adotoa**2 +eft_cache%Hdot ) &
            & + a*eft_cache%adotoa**2*eft_cache%EFTOmegaPP )/(eft_par_cache%h0_Mpc*a*eft_cache%adotoa)
        eft_cache%EFTGamma2PP   = (2._dl*RPH_AB_V*(a**2*RPH_PM_PP*H**4 + 2._dl*a*RPH_PM_P*H**2*(-H**2 + Ht) + (2._dl*H**4 - 3._dl*H**2*Ht - Ht**2 + H*Ht2)*Mp1) + a*(eft_cache%EFTOmegaP*Ht**2 - eft_cache%EFTOmegaP*H*Ht2 + H**4*(-(a**2*eft_cache%EFTOmegaPPP) + 4._dl*a*RPH_AB_P*RPH_PM_P + 2._dl*a*RPH_AB_PP*Mp1 - 4._dl*RPH_AB_P*Mp1) + H**2*Ht*(-2._dl*a*eft_cache%EFTOmegaPP + eft_cache%EFTOmegaP + 4._dl*RPH_AB_P*Mp1)))/(a**3*H**3*eft_par_cache%h0_Mpc)
        eft_cache%EFTGamma2PPP  = (-(a**4*eft_cache%EFTOmegaPPPP*H**6) + 2._dl*RPH_AB_V*(-6._dl*H**6 + 11._dl*H**4*Ht + 3._dl*Ht**3 - 6._dl*H**3*Ht2 - 4._dl*H*Ht*Ht2 + H**2*(6._dl*Ht**2 + Ht3))*Mp1 + a**3*(2._dl*RPH_AB_V*RPH_PM_PPP*H**6 + 6._dl*RPH_PM_PP*RPH_AB_P*H**6 + 6._dl*RPH_AB_PP*RPH_PM_P*H**6 - 3._dl*eft_cache%EFTOmegaPPP*H**4*Ht + 2._dl*RPH_AB_PPP*H**6*Mp1) - 3._dl*a**2*H**2*(2._dl*RPH_AB_V*RPH_PM_PP*H**2*(H**2 - Ht) + 4._dl*RPH_AB_P*RPH_PM_P*H**2*(H**2 - Ht) - eft_cache%EFTOmegaPP*H**2*Ht - eft_cache%EFTOmegaPP*Ht**2 + eft_cache%EFTOmegaPP*H*Ht2 + 2._dl*RPH_AB_PP*H**4*Mp1 - 2._dl*RPH_AB_PP*H**2*Ht*Mp1) + a*(6._dl*RPH_AB_V*RPH_PM_P*H**2*(2._dl*H**4 - 3._dl*H**2*Ht - Ht**2 + H*Ht2) - eft_cache%EFTOmegaP*(2._dl*H**4*Ht + 3._dl*Ht**3 - 3._dl*H**3*Ht2 - 4._dl*H*Ht*Ht2 + H**2*(3._dl*Ht**2 + Ht3)) + 6._dl*RPH_AB_P*H**2*(2._dl*H**4 - 3._dl*H**2*Ht - Ht**2 + H*Ht2)*Mp1))/(a**4*H**5*eft_par_cache%h0_Mpc)
        eft_cache%EFTGamma3V    = -RPH_AT_V*(1._dl +RPH_PM_V)
        eft_cache%EFTGamma3P    = -RPH_PM_P*RPH_AT_V -(1._dl +RPH_PM_V)*RPH_AT_P
        eft_cache%EFTGamma3PP   = -(1._dl + RPH_PM_V)*RPH_AT_PP - RPH_PM_PP*RPH_AT_V - 2._dl*RPH_PM_P*RPH_AT_P
        eft_cache%EFTGamma3PPP  = -(RPH_AT_V*RPH_PM_PPP) - 3._dl*RPH_PM_PP*RPH_AT_P - 3._dl*RPH_AT_PP*RPH_PM_P - RPH_AT_PPP*(1._dl + RPH_PM_V)
        eft_cache%EFTGamma3PPPP = -6._dl*RPH_AT_PP*RPH_PM_PP - RPH_AT_V*RPH_PM_PPPP - 4._dl*RPH_PM_PPP*RPH_AT_P - 4._dl*RPH_AT_PPP*RPH_PM_P - RPH_AT_PPPP*(1._dl + RPH_PM_V)
        eft_cache%EFTGamma4V    = -eft_cache%EFTGamma3V
        eft_cache%EFTGamma4P    = -eft_cache%EFTGamma3P
        eft_cache%EFTGamma4PP   = -eft_cache%EFTGamma3PP
        eft_cache%EFTGamma5V    = +0.5_dl*eft_cache%EFTGamma3V
        eft_cache%EFTGamma5P    = +0.5_dl*eft_cache%EFTGamma3P
        eft_cache%EFTGamma6V    = 0._dl
        eft_cache%EFTGamma6P    = 0._dl

    end subroutine EFTCAMBRPHSecondOrderEFTFunctions

    ! ---------------------------------------------------------------------------------------------
    !> Function that computes dtauda = 1/sqrt(a^2H^2)
    function EFTCAMBRPHComputeDtauda( self, a, eft_par_cache, eft_cache )

        implicit none

        class(EFTCAMB_RPH)                           :: self          !< the base class
        real(dl), intent(in)                         :: a             !< the input scale factor
        type(TEFTCAMB_parameter_cache), intent(inout) :: eft_par_cache !< the EFTCAMB parameter cache that contains all the physical parameters.
        type(TEFTCAMB_timestep_cache ), intent(inout) :: eft_cache     !< the EFTCAMB timestep cache that contains all the physical values.

        real(dl) :: EFTCAMBRPHComputeDtauda                           !< the output dtauda

        real(dl) :: adotoa

        real(dl) :: x

        x = log(a)
        if ( x < self%rhov_t%x_initial ) then
            adotoa = self%adotoa%y(1)/a*Exp(self%adotoa%x_initial)
        else
            adotoa = self%adotoa%value( x )
        end if
        EFTCAMBRPHComputeDtauda = 1._dl/(a*adotoa)

    end function EFTCAMBRPHComputeDtauda

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes adotoa = H.
    subroutine EFTCAMBRPHComputeAdotoa( self, a, eft_par_cache, eft_cache )

        implicit none

        class(EFTCAMB_RPH)                           :: self          !< the base class
        real(dl), intent(in)                         :: a             !< the input scale factor
        type(TEFTCAMB_parameter_cache), intent(inout) :: eft_par_cache !< the EFTCAMB parameter cache that contains all the physical parameters.
        type(TEFTCAMB_timestep_cache ), intent(inout) :: eft_cache     !< the EFTCAMB timestep cache that contains all the physical values.

        real(dl) :: x, mu
        integer  :: ind

        x = log(a)
        if ( x < self%rhov_t%x_initial ) then
            eft_cache%grhov_t = a**2*self%rhov_t%y(1)
            eft_cache%adotoa = self%adotoa%y(1)/a*Exp(self%adotoa%x_initial)
        else
            call self%rhov_t%precompute( x, ind, mu )
            eft_cache%grhov_t = a**2*self%rhov_t%value( x, index=ind, coeff=mu )
            eft_cache%adotoa = self%adotoa%value( x, index=ind, coeff=mu )
        end if

    end subroutine EFTCAMBRPHComputeAdotoa

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes the two derivatives wrt conformal time of H.
    subroutine EFTCAMBRPHComputeHubbleDer( self, a, eft_par_cache, eft_cache )

        implicit none

        class(EFTCAMB_RPH)                           :: self          !< the base class
        real(dl), intent(in)                         :: a             !< the input scale factor
        type(TEFTCAMB_parameter_cache), intent(inout) :: eft_par_cache !< the EFTCAMB parameter cache that contains all the physical parameters.
        type(TEFTCAMB_timestep_cache ), intent(inout) :: eft_cache     !< the EFTCAMB timestep cache that contains all the physical values.

        real(dl) :: x, mu
        integer  :: ind

        x = log(a)

        if ( x < self%rhov_t%x_initial ) then
            eft_cache%gpiv_t     = -eft_cache%grhov_t
            eft_cache%Hdot       = self%Hdot%y(1)
            eft_cache%Hdotdot    = 0._dl
            eft_cache%Hdotdotdot = 0._dl
        else
            call self%adotoa%precompute( x, ind, mu )
            eft_cache%gpiv_t     = self%RPH_wDE%value(a)*eft_cache%grhov_t
            eft_cache%Hdot       = self%Hdot%value( x, ind, mu )
            eft_cache%Hdotdot    = self%Hdot2%value( x, ind, mu )
            eft_cache%Hdotdotdot = self%Hdot3%value( x, ind, mu )
        end if

    end subroutine EFTCAMBRPHComputeHubbleDer

    ! ---------------------------------------------------------------------------------------------
    !> Function that computes model specific stability requirements.
    function EFTCAMBRPHAdditionalModelStability( self, a, eft_par_cache, eft_cache )

        implicit none

        class(EFTCAMB_RPH)                           :: self          !< the base class
        real(dl), intent(in)                         :: a             !< the input scale factor.
        type(TEFTCAMB_parameter_cache), intent(inout) :: eft_par_cache !< the EFTCAMB parameter cache that contains all the physical parameters.
        type(TEFTCAMB_timestep_cache ), intent(inout) :: eft_cache     !< the EFTCAMB timestep cache that contains all the physical values.

        logical :: EFTCAMBRPHAdditionalModelStability          !< the return value of the stability computation. True if the model specific stability criteria are met, false otherwise.

        EFTCAMBRPHAdditionalModelStability = .True.
        if ( self%RPH_wDE%value(a) > -1._dl/3._dl ) EFTCAMBRPHAdditionalModelStability = .False.

    end function EFTCAMBRPHAdditionalModelStability

    ! ---------------------------------------------------------------------------------------------

end module EFTCAMB_Reparametrized_Horndeski
