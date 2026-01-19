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

!> @file 007p8_3p1_wDE_designer.f90
!! This file contains the definition of the EFT corresponding to a full 1+3 theory that 
!! adopts arbitrary wDE.


!----------------------------------------------------------------------------------------
!> This file contains the definition of the EFT corresponding to a full 1+3 theory that 
!! adopts arbitrary wDE.

!> @author Gen Ye

module EFTCAMB_designer_full3p1_wDE

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

    public EFTCAMB_full3p1_wDE

    !----------------------------------------------------------------------------------------
    !> This file contains the definition of the EFT corresponding to a full 1+3 theory that adopts arbitrary wDE.
    type, extends ( EFTCAMB_designer_model ) :: EFTCAMB_full3p1_wDE

        ! the pure EFT functions model selection flags:
        integer  :: Designer3p1modelsigma1   !< Model selection flag for sigma1.
        integer  :: EFTwDE                   !< Model selection flag for w DE.
        integer  :: Designer3p1modelrho1     !< Model selection flag for rho1.
        integer  :: Designer3p1modelrho1N    !< Model selection flag for \partial_N rho1.
        integer  :: Designer3p1modelrho2     !< Model selection flag for rho2.

        ! the pure EFT functions model selection flags:
        integer  :: Designer3p1modelsigma1_ODE   !< Model selection flag for sigma1 with Omega_DE argument.
        integer  :: Designer3p1modelrho1_ODE     !< Model selection flag for rho1 with Omega_DE argument.
        integer  :: Designer3p1modelrho1N_ODE    !< Model selection flag for \partial_N rho1 with Omega_DE argument.
        integer  :: Designer3p1modelrho2_ODE     !< Model selection flag for rho2 with Omega_DE argument.

        ! selection flag for Horndeski:
        logical  :: Designer3p1LuminalGW         !< Selects wether to enforce c_T = 1.

        ! flag deciding whether computes Omega_DE
        logical  :: Designer3p1OmegaDE

        ! flag deciding whether impose ghost stability of this model (beyond GLPV)
        logical  :: Designer3p1GhostStability
        ! flag deciding whether impose gradient stability of this model (beyond GLPV)
        logical  :: Designer3p1GradientStability
        ! stability is checked for kmin < k < kmax with logstep dlgk
        real(dl)   :: Designer3p1Stability_kmin
        real(dl)   :: Designer3p1Stability_kmax
        real(dl)   :: Designer3p1Stability_dlgk

        ! DE energy density from wDE integration
        type( equispaced_linear_interpolate_function_1D ) :: rhov_t

        ! the dark energy equation of state
        class( parametrized_function_1D ), allocatable :: Designer3p1wDE         !< The dark energy equation of state w_DE.
        
        ! the Lagrangian functions:
        class( parametrized_function_1D ), allocatable :: Designer3p1sigma1      !< The lagrangian function sigma1 as a function of scale factor.
        class( parametrized_function_1D ), allocatable :: Designer3p1rho1        !< The lagrangian function rho1 as a function of scale factor.
        class( parametrized_function_1D ), allocatable :: Designer3p1rho1N       !< Paritial derivative wrt N of the lagrangian function rho1 as a function of scale factor.
        class( parametrized_function_1D ), allocatable :: Designer3p1rho2        !< The lagrangian function rho2 as a function of scale factor.

        class( parametrized_function_1D ), allocatable :: Designer3p1sigma1_ODE  !< The lagrangian function sigma1 as a function of dark energy fraction.
        class( parametrized_function_1D ), allocatable :: Designer3p1rho1_ODE    !< The lagrangian function rho1 as a function of dark energy fraction.
        class( parametrized_function_1D ), allocatable :: Designer3p1rho1N_ODE   !< Paritial derivative wrt N of the lagrangian function rho1 as a function of dark energy fraction.
        class( parametrized_function_1D ), allocatable :: Designer3p1rho2_ODE    !< The lagrangian function rho2 as a function of dark energy fraction.

        ! the pure EFT functions:
        type( equispaced_linear_interpolate_function_1D ) :: Designer3p1Omega    !< The pure EFT function Omega.
        type( equispaced_linear_interpolate_function_1D ) :: Designer3p1Gamma1   !< The pure EFT function Gamma1.
        type( equispaced_linear_interpolate_function_1D ) :: Designer3p1Gamma2   !< The pure EFT function Gamma2.
        type( equispaced_linear_interpolate_function_1D ) :: Designer3p1Gamma3   !< The pure EFT function Gamma3.
        type( equispaced_linear_interpolate_function_1D ) :: Designer3p1Gamma5   !< The pure EFT function Gamma5.
        type( equispaced_linear_interpolate_function_1D ) :: Designer3p1Gamma6   !< The pure EFT function Gamma6.

        ! background solver parameters:
        integer  :: designer_num_points = 6000                            !< Number of points sampled by the designer code.
        real(dl) :: x_initial           = log(10._dl**(-8))               !< log(a start)
        real(dl) :: x_final             = 0.0_dl                          !< log(a final)

    contains

        ! initialization of the model:
        procedure :: read_model_selection            => EFTCAMBDesigner3p1ReadModelSelectionFromFile  !< subroutine that reads the parameters of the model from file
        procedure :: allocate_model_selection        => EFTCAMBDesigner3p1AllocateModelSelection      !< subroutine that allocates the model selection.
        procedure :: init_model_parameters           => EFTCAMBDesigner3p1InitModelParameters         !< subroutine that initializes the model parameters based on the values found in an input array.
        procedure :: init_model_parameters_from_file => EFTCAMBDesigner3p1InitModelParametersFromFile !< subroutine that reads the parameters of the model from file.

        ! background solver:
        procedure :: initialize_background           => EFTCAMBDesigner3p1InitBackground               !< subroutine that initializes the background.
        procedure :: solve_designer_equations        => EFTCAMBDesigner3p1SolveDesignerEquations       !< subroutine that solves the for Mp.

        ! utility functions:
        procedure :: compute_param_number  => EFTCAMBDesigner3p1ComputeParametersNumber    !< subroutine that computes the number of parameters of the model.
        procedure :: feedback              => EFTCAMBDesigner3p1Feedback                   !< subroutine that prints on the screen feedback information about the model.
        procedure :: parameter_names       => EFTCAMBDesigner3p1ParameterNames             !< subroutine that returns the i-th parameter name of the model.
        procedure :: parameter_names_latex => EFTCAMBDesigner3p1ParameterNamesLatex        !< subroutine that returns the i-th parameter name of the model.
        procedure :: parameter_values      => EFTCAMBDesigner3p1ParameterValues            !< subroutine that returns the i-th parameter value.
        ! CAMB related procedures:
        procedure :: compute_background_EFT_functions  => EFTCAMBDesigner3p1BackgroundEFTFunctions   !< subroutine that computes the value of the background EFT functions at a given time.
        procedure :: compute_secondorder_EFT_functions => EFTCAMBDesigner3p1SecondOrderEFTFunctions  !< subroutine that computes the value of the second order EFT functions at a given time.
        procedure :: compute_dtauda                    => EFTCAMBDesigner3p1ComputeDtauda            !< function that computes dtauda = 1/sqrt(a^2H^2).
        procedure :: compute_adotoa                    => EFTCAMBDesigner3p1ComputeAdotoa            !< subroutine that computes adotoa = H and its two derivatives wrt conformal time.
        procedure :: compute_H_derivs                  => EFTCAMBDesigner3p1ComputeHubbleDer         !< subroutine that computes the three derivatives wrt conformal time of H.
        ! stability procedures:
        procedure :: additional_model_stability        => EFTCAMBDesigner3p1AdditionalModelStability !< function that computes model specific stability requirements.

    end type EFTCAMB_full3p1_wDE

    ! ---------------------------------------------------------------------------------------------

contains

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that reads the parameters of the model from file.
    subroutine EFTCAMBDesigner3p1ReadModelSelectionFromFile( self, Ini, eft_error )

        implicit none

        class(EFTCAMB_full3p1_wDE) :: self      !< the base class
        type(TIniFile)              :: Ini       !< Input ini file
        integer                     :: eft_error !< error code: 0 all fine, 1 initialization failed

        ! read model selection flags:
        self%Designer3p1modelsigma1  = Ini%Read_Int( 'Designer3p1modelsigma1'  , 0 )
        self%EFTwDE             = Ini%Read_Int( 'EFTwDE'             , 0 )
        self%Designer3p1modelrho1N = Ini%Read_Int( 'Designer3p1modelrho1N' , 0 )
        self%Designer3p1modelrho2 = Ini%Read_Int( 'Designer3p1modelrho2' , 0 )
        
        self%Designer3p1modelsigma1_ODE  = Ini%Read_Int( 'Designer3p1modelsigma1_ODE'  , 0 )
        self%Designer3p1modelrho1N_ODE = Ini%Read_Int( 'Designer3p1modelrho1N_ODE' , 0 )
        self%Designer3p1modelrho2_ODE = Ini%Read_Int( 'Designer3p1modelrho2_ODE' , 0 )
        
        ! read the c_T=1 flag:
        self%Designer3p1LuminalGW   = Ini%Read_Logical( 'Designer3p1LuminalGW' , .false. )
        if (.not. self%Designer3p1LuminalGW) then
            self%Designer3p1modelrho1 = Ini%Read_Int( 'Designer3p1modelrho1' , 0 )
            self%Designer3p1modelrho1_ODE = Ini%Read_Int( 'Designer3p1modelrho1_ODE' , 0 )
        end if

        ! read the stability flags:
        self%Designer3p1GhostStability = Ini%Read_Logical( 'Designer3p1GhostStability', .true. )
        self%Designer3p1GradientStability = Ini%Read_Logical( 'Designer3p1GradientStability', .true. )
        self%Designer3p1Stability_kmin = Ini%Read_Double( 'Designer3p1Stability_kmin', 1d-4 )
        self%Designer3p1Stability_kmax = Ini%Read_Double( 'Designer3p1Stability_kmax', 1._dl )
        self%Designer3p1Stability_dlgk = Ini%Read_Double( 'Designer3p1Stability_dlgk', 0.1_dl )
        
        ! read precision parameters
        self%designer_num_points = Ini%Read_Int( 'model_background_num_points', 6000 )
        self%x_initial = Log( Ini%Read_Double( 'model_background_a_ini', 1d-8 ) )
        self%x_final = Log( Ini%Read_Double( 'model_background_a_final', 1._dl ) )


    end subroutine EFTCAMBDesigner3p1ReadModelSelectionFromFile

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that allocates the model selection.
    subroutine EFTCAMBDesigner3p1AllocateModelSelection( self, Ini, eft_error )

        implicit none

        class(EFTCAMB_full3p1_wDE) :: self      !< the base class
        type(TIniFile)              :: Ini       !< Input ini file
        integer                     :: eft_error !< error code: 0 all fine, 1 initialization failed

        integer                     :: temp_feedback 

        ! get feedback flag:
        temp_feedback = Ini%Read_Int('feedback_level', 0)

        ! allocate wDE:
        if ( allocated(self%Designer3p1wDE) ) deallocate(self%Designer3p1wDE)
        select case ( self%EFTwDE )
            case(0)
                allocate( wDE_LCDM_parametrization_1D::self%Designer3p1wDE )
                call self%Designer3p1wDE%set_name( 'EFTw', 'w' )
            case(1)
                allocate( constant_parametrization_1D::self%Designer3p1wDE )
                call self%Designer3p1wDE%set_name( 'EFTw', 'w' )
            case(2)
                allocate( CPL_parametrization_1D::self%Designer3p1wDE )
                call self%Designer3p1wDE%set_param_names( ['EFTw0', 'EFTwa'], ['w_0', 'w_a'] )
                call self%Designer3p1wDE%set_name( 'EFTw', 'w' )
            case(3)
                allocate( JBP_parametrization_1D::self%Designer3p1wDE )
                call self%Designer3p1wDE%set_param_names( ['EFTw0', 'EFTwa', 'EFTwn'], [ 'w_0', 'w_a', 'n  ' ] )
                call self%Designer3p1wDE%set_name( 'EFTw', 'w' )
            case(4)
                allocate( turning_point_parametrization_1D::self%Designer3p1wDE )
                call self%Designer3p1wDE%set_param_names( ['EFTw0 ', 'EFTwa ', 'EFTwat'], ['w_0', 'w_a', 'a_t'] )
                call self%Designer3p1wDE%set_name( 'EFTw', 'w' )
            case(5)
                allocate( taylor_parametrization_1D::self%Designer3p1wDE )
                call self%Designer3p1wDE%set_param_names( ['EFTw0', 'EFTwa', 'EFTw2', 'EFTw3'], ['w_0', 'w_a', 'w_2', 'w_3'] )
                call self%Designer3p1wDE%set_name( 'EFTw', 'w' )
            case default
                call allocate_parametrized_1D_function( self%Designer3p1wDE, self%EFTwDE, 'EFTw', 'w', eft_error, temp_feedback )
                if ( eft_error == 1 ) return
        end select

        ! allocate the lagrangian functions:
        ! sigma1
        call allocate_parametrized_1D_function( self%Designer3p1sigma1, self%Designer3p1modelsigma1, 'Designer3p1sigma1', '\sigma_1', eft_error, temp_feedback )
        if ( eft_error == 1 ) return
        call allocate_parametrized_1D_function( self%Designer3p1sigma1_ODE, self%Designer3p1modelsigma1_ODE, 'Designer3p1sigma1_ODE', '\sigma_{1,ODE}', eft_error, temp_feedback )
        if ( eft_error == 1 ) return
        ! rho1N
        call allocate_parametrized_1D_function( self%Designer3p1rho1N, self%Designer3p1modelrho1N, 'Designer3p1rho1N', '\partial_N\rho_1', eft_error, temp_feedback )
        if ( eft_error == 1 ) return
        call allocate_parametrized_1D_function( self%Designer3p1rho1N_ODE, self%Designer3p1modelrho1N_ODE, 'Designer3p1rho1N_ODE', '\partial_N\rho_{1,ODE}', eft_error, temp_feedback )
        if ( eft_error == 1 ) return
        ! rho2
        call allocate_parametrized_1D_function( self%Designer3p1rho2, self%Designer3p1modelrho2, 'Designer3p1rho2', '\rho_2', eft_error, temp_feedback )
        if ( eft_error == 1 ) return
        call allocate_parametrized_1D_function( self%Designer3p1rho2_ODE, self%Designer3p1modelrho2_ODE, 'Designer3p1rho2_ODE', '\rho_{2,ODE}', eft_error, temp_feedback )
        if ( eft_error == 1 ) return
        ! rho1 if not enforcing c_T=1, otherwise rho1 = sigma1
        if ( .not. self%Designer3p1LuminalGW ) then
            ! rho1
            call allocate_parametrized_1D_function( self%Designer3p1rho1, self%Designer3p1modelrho1, 'Designer3p1rho1', '\rho_1', eft_error, temp_feedback )
            if ( eft_error == 1 ) return
            call allocate_parametrized_1D_function( self%Designer3p1rho1_ODE, self%Designer3p1modelrho1_ODE, 'Designer3p1rho1_ODE', '\rho_{1,ODE}', eft_error, temp_feedback )
            if ( eft_error == 1 ) return
        end if

        ! additional initialization of the function:
        call self%Designer3p1wDE%init_func_from_file   ( Ini, eft_error )
        
        call self%Designer3p1sigma1%init_func_from_file( Ini, eft_error )
        call self%Designer3p1rho1N%init_func_from_file  ( Ini, eft_error )
        call self%Designer3p1rho2%init_func_from_file  ( Ini, eft_error )

        call self%Designer3p1sigma1_ODE%init_func_from_file( Ini, eft_error )
        call self%Designer3p1rho1N_ODE%init_func_from_file  ( Ini, eft_error )
        call self%Designer3p1rho2_ODE%init_func_from_file  ( Ini, eft_error )

        if ( .not. self%Designer3p1LuminalGW ) then
            call self%Designer3p1rho1%init_func_from_file  ( Ini, eft_error )
            call self%Designer3p1rho1_ODE%init_func_from_file  ( Ini, eft_error )
        end if

        if (self%Designer3p1modelsigma1_ODE>1 .or. self%Designer3p1modelrho1N_ODE>1 .or. self%Designer3p1modelrho2_ODE>1 .or. (self%Designer3p1modelrho1_ODE>1 .and. .not. self%Designer3p1LuminalGW)) then
            self%Designer3p1OmegaDE = .true.
        else
            self%Designer3p1OmegaDE = .false.
        end if

    end subroutine EFTCAMBDesigner3p1AllocateModelSelection

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that initializes the model parameters based on the values found in an input array.
    subroutine EFTCAMBDesigner3p1InitModelParameters( self, array )

        implicit none

        class(EFTCAMB_full3p1_wDE)                            :: self   !< the base class
        real(dl), dimension(self%parameter_number), intent(in) :: array  !< input array with the values of the parameters of the model.

        real(dl), allocatable, dimension(:)                    :: temp
        integer :: num_params_function, num_params_temp, i
        
        num_params_temp     = 1
        
        ! first elements are w_DE parameters:
        num_params_function = self%Designer3p1wDE%parameter_number
        allocate( temp(num_params_function) )
        do i = 1, num_params_function
            temp(i)         = array(num_params_temp)
            num_params_temp = num_params_temp +1
        end do
        call self%Designer3p1wDE%init_parameters(temp)
        deallocate(temp)
        ! then sigma1:
        num_params_function = self%Designer3p1sigma1%parameter_number
        allocate( temp(num_params_function) )
        do i = 1, num_params_function
            temp(i)         = array(num_params_temp)
            num_params_temp = num_params_temp +1
        end do
        call self%Designer3p1sigma1%init_parameters(temp)
        deallocate( temp )
        ! then rho1N:
        num_params_function = self%Designer3p1rho1N%parameter_number
        allocate( temp(num_params_function) )
        do i = 1, num_params_function
            temp(i)         = array(num_params_temp)
            num_params_temp = num_params_temp +1
        end do
        call self%Designer3p1rho1N%init_parameters(temp)
        deallocate( temp )
        ! then rho2:
        num_params_function = self%Designer3p1rho2%parameter_number
        allocate( temp(num_params_function) )
        do i = 1, num_params_function
            temp(i)         = array(num_params_temp)
            num_params_temp = num_params_temp +1
        end do
        call self%Designer3p1rho2%init_parameters(temp)
        deallocate( temp )

        ! The Omega_DE parameters
        ! sigma1:
        num_params_function = self%Designer3p1sigma1_ODE%parameter_number
        allocate( temp(num_params_function) )
        do i = 1, num_params_function
            temp(i)         = array(num_params_temp)
            num_params_temp = num_params_temp +1
        end do
        call self%Designer3p1sigma1_ODE%init_parameters(temp)
        deallocate( temp )
        ! then rho1N:
        num_params_function = self%Designer3p1rho1N_ODE%parameter_number
        allocate( temp(num_params_function) )
        do i = 1, num_params_function
            temp(i)         = array(num_params_temp)
            num_params_temp = num_params_temp +1
        end do
        call self%Designer3p1rho1N_ODE%init_parameters(temp)
        deallocate( temp )
        ! then rho2:
        num_params_function = self%Designer3p1rho2_ODE%parameter_number
        allocate( temp(num_params_function) )
        do i = 1, num_params_function
            temp(i)         = array(num_params_temp)
            num_params_temp = num_params_temp +1
        end do
        call self%Designer3p1rho2_ODE%init_parameters(temp)
        deallocate( temp )

        ! rho1 if not requiring c_T=1:
        if ( .not. self%Designer3p1LuminalGW ) then

            ! rho1:
            num_params_function = self%Designer3p1rho1%parameter_number
            allocate( temp(num_params_function) )
            do i = 1, num_params_function
                temp(i)         = array(num_params_temp)
                num_params_temp = num_params_temp +1
            end do
            call self%Designer3p1rho1%init_parameters(temp)
            deallocate( temp )
            ! rho1:
            num_params_function = self%Designer3p1rho1_ODE%parameter_number
            allocate( temp(num_params_function) )
            do i = 1, num_params_function
                temp(i)         = array(num_params_temp)
                num_params_temp = num_params_temp +1
            end do
            call self%Designer3p1rho1_ODE%init_parameters(temp)
            deallocate( temp )

        end if

        ! now check the length of the parameters:
        if ( num_params_temp-1 /= self%parameter_number ) then
            write(*,*) 'In EFTCAMBDesigner3p1InitModelParameters:'
            write(*,*) 'Length of num_params_temp and self%parameter_number do not coincide.'
            write(*,*) 'num_params_temp:', num_params_temp-1
            write(*,*) 'self%parameter_number:', self%parameter_number
            call MpiStop('EFTCAMB error')
        end if

    end subroutine EFTCAMBDesigner3p1InitModelParameters

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that reads the parameters of the model from file.
    subroutine EFTCAMBDesigner3p1InitModelParametersFromFile( self, Ini, eft_error )

        implicit none

        class(EFTCAMB_full3p1_wDE) :: self      !< the base class
        type(TIniFile)              :: Ini       !< Input ini file
        integer                     :: eft_error !< error code: 0 all fine, 1 initialization failed
        
        call self%Designer3p1wDE%init_from_file   ( Ini, eft_error )

        call self%Designer3p1sigma1%init_from_file ( Ini, eft_error )
        call self%Designer3p1rho1N%init_from_file ( Ini, eft_error )
        call self%Designer3p1rho2%init_from_file ( Ini, eft_error )

        call self%Designer3p1sigma1_ODE%init_from_file ( Ini, eft_error )
        call self%Designer3p1rho1N_ODE%init_from_file ( Ini, eft_error )
        call self%Designer3p1rho2_ODE%init_from_file ( Ini, eft_error )

        if ( .not. self%Designer3p1LuminalGW ) then
            call self%Designer3p1rho1%init_from_file ( Ini, eft_error )

            call self%Designer3p1rho1_ODE%init_from_file ( Ini, eft_error )
        end if

    end subroutine EFTCAMBDesigner3p1InitModelParametersFromFile

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes the number of parameters of the model.
    subroutine EFTCAMBDesigner3p1ComputeParametersNumber( self )

        implicit none

        class(EFTCAMB_full3p1_wDE)  :: self   !< the base class

        self%parameter_number = 0
        self%parameter_number = self%parameter_number +self%Designer3p1wDE%parameter_number
        self%parameter_number = self%parameter_number +self%Designer3p1sigma1%parameter_number
        self%parameter_number = self%parameter_number +self%Designer3p1rho1N%parameter_number
        self%parameter_number = self%parameter_number +self%Designer3p1rho2%parameter_number

        self%parameter_number = self%parameter_number +self%Designer3p1sigma1_ODE%parameter_number
        self%parameter_number = self%parameter_number +self%Designer3p1rho1N_ODE%parameter_number
        self%parameter_number = self%parameter_number +self%Designer3p1rho2_ODE%parameter_number

        if ( .not. self%Designer3p1LuminalGW ) then
            self%parameter_number = self%parameter_number +self%Designer3p1rho1%parameter_number

            self%parameter_number = self%parameter_number +self%Designer3p1rho1_ODE%parameter_number
        end if

    end subroutine EFTCAMBDesigner3p1ComputeParametersNumber

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that prints on the screen feedback information about the model.
    subroutine EFTCAMBDesigner3p1Feedback( self, print_params )

        implicit none

        class(EFTCAMB_full3p1_wDE)  :: self         !< the base class
        logical, optional            :: print_params !< optional flag that decised whether to print numerical values
                                                     !! of the parameters.

        ! print general model informations:
        write(*,*)
        write(*,'(a,a)')    '   Model               =  ', self%name
        if ( self%Designer3p1LuminalGW ) then
            write(*,"(a)")  '   with luminal GW'
        end if
        write(*,'(a,I3)')   '   Number of params    ='  , self%parameter_number

        ! print model functions informations:
        write(*,*)
        if ( self%EFTwDE                  /= 0 ) write(*,'(a,I3)') '   EFTwDE              =', self%EFTwDE
        if ( self%Designer3p1modelsigma1  /= 0 ) write(*,'(a,I3)') '   Designer3p1modelsigma1  =', self%Designer3p1modelsigma1
        if ( self%Designer3p1modelrho1N  /= 0 ) write(*,'(a,I3)') '   Designer3p1modelrho1N  =', self%Designer3p1modelrho1N
        if ( self%Designer3p1modelrho2    /= 0 ) write(*,'(a,I3)') '   Designer3p1modelrho2  =', self%Designer3p1modelrho2

        if ( self%Designer3p1modelsigma1_ODE  /= 0 ) write(*,'(a,I3)') '   Designer3p1modelsigma1_ODE  =', self%Designer3p1modelsigma1_ODE
        if ( self%Designer3p1modelrho1N_ODE    /= 0 ) write(*,'(a,I3)') '   Designer3p1modelrho1N_ODE  =', self%Designer3p1modelrho1N_ODE
        if ( self%Designer3p1modelrho2_ODE    /= 0 ) write(*,'(a,I3)') '   Designer3p1modelrho2_ODE  =', self%Designer3p1modelrho2_ODE
        if ( .not. self%Designer3p1LuminalGW ) then
            if ( self%Designer3p1modelrho1    /= 0 ) write(*,'(a,I3)') '   Designer3p1modelrho1  =', self%Designer3p1modelrho1

            if ( self%Designer3p1modelrho1_ODE    /= 0 ) write(*,'(a,I3)') '   Designer3p1modelrho1_ODE  =', self%Designer3p1modelrho1_ODE
        end if

        write(*,*)
        ! print functions informations:
        call self%Designer3p1wDE%feedback    ( print_params )
        call self%Designer3p1sigma1%feedback ( print_params )
        call self%Designer3p1rho1N%feedback ( print_params )
        call self%Designer3p1rho2%feedback ( print_params )

        call self%Designer3p1sigma1_ODE%feedback ( print_params )
        call self%Designer3p1rho1N_ODE%feedback ( print_params )
        call self%Designer3p1rho2_ODE%feedback ( print_params )
        if ( .not. self%Designer3p1LuminalGW ) then
            call self%Designer3p1rho1%feedback ( print_params )

            call self%Designer3p1rho1_ODE%feedback ( print_params )
        end if

    end subroutine EFTCAMBDesigner3p1Feedback

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the i-th parameter name of the model
    subroutine EFTCAMBDesigner3p1ParameterNames( self, i, name )

        implicit none

        class(EFTCAMB_full3p1_wDE) :: self   !< the base class
        integer     , intent(in)    :: i      !< the index of the parameter
        character(*), intent(out)   :: name   !< the output name of the i-th parameter

        integer  :: Nc, Nw, N1, N1n, N2, N3, N1o, N1no, N2o, N3o
        integer  :: j

        ! compute the incremental number of parameters:
        Nc     = 0
        Nw     = Nc + self%Designer3p1wDE%parameter_number
        
        N1     = Nw + self%Designer3p1sigma1%parameter_number
        N1n    = N1 + self%Designer3p1rho1N%parameter_number
        N2     = N1n + self%Designer3p1rho2%parameter_number
        
        N1o     = N2 + self%Designer3p1sigma1_ODE%parameter_number
        N1no    = N1o + self%Designer3p1rho1N_ODE%parameter_number
        N2o     = N1no + self%Designer3p1rho2_ODE%parameter_number

        if ( .not. self%Designer3p1LuminalGW ) then
            N3 = N2o + self%Designer3p1rho1%parameter_number

            N3o = N3 + self%Designer3p1rho1_ODE%parameter_number
        end if

        ! check validity of input:
        if ( i > self%parameter_number .or. i <= 0 ) then
            write(*,'(a,I3)') 'No parameter corresponding to: ', i
            write(*,'(a,I3)') 'Total number of parameters is: ', self%parameter_number
            call MpiStop('EFTCAMB error')

        ! parameter from wDE function
        else if ( i <= Nw ) then
            do j = 1, self%Designer3p1wDE%parameter_number
                if ( i-Nc == j ) call self%Designer3p1wDE%parameter_names( j, name )
            end do
            return

        ! parameter from sigma function
        else if ( i <= N1) then
            do j = 1, self%Designer3p1sigma1%parameter_number
                if ( i-Nw == j ) call self%Designer3p1sigma1%parameter_names( j, name )
            end do
            return

        ! parameter from rho1N function
        else if ( i <= N1n ) then
            do j = 1, self%Designer3p1rho1N%parameter_number
                if ( i-N1 == j ) call self%Designer3p1rho1N%parameter_names( j, name )
            end do
            return
        
        ! parameter from rho2 function
        else if ( i <= N2 ) then
            do j = 1, self%Designer3p1rho2%parameter_number
                if ( i-N1n == j ) call self%Designer3p1rho2%parameter_names( j, name )
            end do
            return

        ! parameter from sigma1_ODE function
        else if ( i <= N1o) then
            do j = 1, self%Designer3p1sigma1_ODE%parameter_number
                if ( i-N2 == j ) call self%Designer3p1sigma1_ODE%parameter_names( j, name )
            end do
            return

        ! parameter from rho2_ODE function
        else if ( i <= N1no ) then
            do j = 1, self%Designer3p1rho1N_ODE%parameter_number
                if ( i-N1o == j ) call self%Designer3p1rho1N_ODE%parameter_names( j, name )
            end do
            return

        ! parameter from rho2_ODE function
        else if ( i <= N2o ) then
            do j = 1, self%Designer3p1rho2_ODE%parameter_number
                if ( i-N1no == j ) call self%Designer3p1rho2_ODE%parameter_names( j, name )
            end do
            return

        ! parameter from rho1 function
        else if ( .not. self%Designer3p1LuminalGW .and. i <= N3 ) then

            do j = 1, self%Designer3p1rho1%parameter_number
                if ( i-N2o == j ) call self%Designer3p1rho1%parameter_names( j, name )
            end do
            return

        ! parameter from rho1_ODE function
        else if ( .not. self%Designer3p1LuminalGW .and. i <= N3o ) then

            do j = 1, self%Designer3p1rho1_ODE%parameter_number
                if ( i-N3 == j ) call self%Designer3p1rho1_ODE%parameter_names( j, name )
            end do
            return

        end if

    end subroutine EFTCAMBDesigner3p1ParameterNames

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the i-th parameter name of the model
    subroutine EFTCAMBDesigner3p1ParameterNamesLatex( self, i, latexname )

        implicit none

        class(EFTCAMB_full3p1_wDE) :: self       !< the base class
        integer     , intent(in)    :: i          !< The index of the parameter
        character(*), intent(out)   :: latexname  !< the output latex name of the i-th parameter

        integer  :: Nc, Nw, N1, N1n, N2, N3, N1o, N1no, N2o, N3o
        integer  :: j

        ! compute the incremental number of parameters:
        Nc     = 0
        Nw     = Nc + self%Designer3p1wDE%parameter_number
        
        N1     = Nw + self%Designer3p1sigma1%parameter_number
        N1n    = N1 + self%Designer3p1rho1N%parameter_number
        N2     = N1n + self%Designer3p1rho2%parameter_number
        
        N1o     = N2 + self%Designer3p1sigma1_ODE%parameter_number
        N1no    = N1o + self%Designer3p1rho1N_ODE%parameter_number
        N2o     = N1no + self%Designer3p1rho2_ODE%parameter_number

        ! check validity of input:
        if ( i > self%parameter_number .or. i <= 0 ) then
            write(*,'(a,I3)') 'No parameter corresponding to: ', i
            write(*,'(a,I3)') 'Total number of parameters is: ', self%parameter_number
            call MpiStop('EFTCAMB error')

        ! parameter from wDE function
        else if ( i <= Nw ) then
            do j = 1, self%Designer3p1wDE%parameter_number
                if ( i-Nc == j ) call self%Designer3p1wDE%parameter_names_latex( j, latexname )
            end do
            return

        ! parameter from sigma function
        else if ( i <= N1) then
            do j = 1, self%Designer3p1sigma1%parameter_number
                if ( i-Nw == j ) call self%Designer3p1sigma1%parameter_names_latex( j, latexname )
            end do
            return

        ! parameter from rho1N function
        else if ( i <= N1n ) then
            do j = 1, self%Designer3p1rho1N%parameter_number
                if ( i-N1 == j ) call self%Designer3p1rho1N%parameter_names_latex( j, latexname )
            end do
            return

        ! parameter from rho2 function
        else if ( i <= N2 ) then
            do j = 1, self%Designer3p1rho2%parameter_number
                if ( i-N1n == j ) call self%Designer3p1rho2%parameter_names_latex( j, latexname )
            end do
            return

        ! parameter from sigma1_ODE function
        else if ( i <= N1o) then
            do j = 1, self%Designer3p1sigma1_ODE%parameter_number
                if ( i-N2 == j ) call self%Designer3p1sigma1_ODE%parameter_names_latex( j, latexname )
            end do
            return

        ! parameter from rho2_ODE function
        else if ( i <= N1no ) then
            do j = 1, self%Designer3p1rho1N_ODE%parameter_number
                if ( i-N1o == j ) call self%Designer3p1rho1N_ODE%parameter_names_latex( j, latexname )
            end do
            return

        ! parameter from rho2_ODE function
        else if ( i <= N2o ) then
            do j = 1, self%Designer3p1rho2_ODE%parameter_number
                if ( i-N1no == j ) call self%Designer3p1rho2_ODE%parameter_names_latex( j, latexname )
            end do
            return

        ! parameter from rho1 function
        else if ( .not. self%Designer3p1LuminalGW .and. i <= N3 ) then

            do j = 1, self%Designer3p1rho1%parameter_number
                if ( i-N2o == j ) call self%Designer3p1rho1%parameter_names_latex( j, latexname )
            end do
            return

        ! parameter from rho1_ODE function
        else if ( .not. self%Designer3p1LuminalGW .and. i <= N3o ) then

            do j = 1, self%Designer3p1rho1_ODE%parameter_number
                if ( i-N3 == j ) call self%Designer3p1rho1_ODE%parameter_names_latex( j, latexname )
            end do
            return

        end if

    end subroutine EFTCAMBDesigner3p1ParameterNamesLatex

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the i-th parameter name of the model
    subroutine EFTCAMBDesigner3p1ParameterValues( self, i, value )

        implicit none

        class(EFTCAMB_full3p1_wDE) :: self   !< the base class
        integer , intent(in)        :: i      !< The index of the parameter
        real(dl), intent(out)       :: value  !< the output value of the i-th parameter

        integer  :: Nc, Nw, N1, N1n, N2, N3, N1o, N1no, N2o, N3o
        integer  :: j

        ! compute the incremental number of parameters:
        Nc     = 0
        Nw     = Nc + self%Designer3p1wDE%parameter_number
        
        N1     = Nw + self%Designer3p1sigma1%parameter_number
        N1n    = N1 + self%Designer3p1rho1N%parameter_number
        N2     = N1n + self%Designer3p1rho2%parameter_number
        
        N1o     = N2 + self%Designer3p1sigma1_ODE%parameter_number
        N1no    = N1o + self%Designer3p1rho1N_ODE%parameter_number
        N2o     = N1no + self%Designer3p1rho2_ODE%parameter_number

        if ( .not. self%Designer3p1LuminalGW ) then
            N3 = N2o + self%Designer3p1rho1%parameter_number

            N3o = N3 + self%Designer3p1rho1_ODE%parameter_number
        end if

        ! check validity of input:
        if ( i > self%parameter_number .or. i <= 0 ) then
            write(*,'(a,I3)') 'No parameter corresponding to: ', i
            write(*,'(a,I3)') 'Total number of parameters is: ', self%parameter_number
            call MpiStop('EFTCAMB error')
        
        ! parameter from wDE function
        else if ( i <= Nw ) then
            do j = 1, self%Designer3p1wDE%parameter_number
                if ( i-Nc == j ) call self%Designer3p1wDE%parameter_value( j, value )
            end do
            return

        ! parameter from sigma function
        else if ( i <= N1) then
            do j = 1, self%Designer3p1sigma1%parameter_number
                if ( i-Nw == j ) call self%Designer3p1sigma1%parameter_value( j, value )
            end do
            return

        ! parameter from rho1N function
        else if ( i <= N1n ) then
            do j = 1, self%Designer3p1rho1N%parameter_number
                if ( i-N1 == j ) call self%Designer3p1rho1N%parameter_value( j, value )
            end do
            return

        ! parameter from rho2 function
        else if ( i <= N2 ) then
            do j = 1, self%Designer3p1rho2%parameter_number
                if ( i-N1n == j ) call self%Designer3p1rho2%parameter_value( j, value )
            end do
            return

        ! parameter from sigma1_ODE function
        else if ( i <= N1o) then
            do j = 1, self%Designer3p1sigma1_ODE%parameter_number
                if ( i-N2 == j ) call self%Designer3p1sigma1_ODE%parameter_value( j, value )
            end do
            return

        ! parameter from rho1N_ODE function
        else if ( i <= N1no ) then
            do j = 1, self%Designer3p1rho1N_ODE%parameter_number
                if ( i-N1o == j ) call self%Designer3p1rho1N_ODE%parameter_value( j, value )
            end do
            return

        ! parameter from rho2_ODE function
        else if ( i <= N2o ) then
            do j = 1, self%Designer3p1rho2_ODE%parameter_number
                if ( i-N1no == j ) call self%Designer3p1rho2_ODE%parameter_value( j, value )
            end do
            return

        ! parameter from rho1 function
        else if ( .not. self%Designer3p1LuminalGW .and. i <= N3 ) then

            do j = 1, self%Designer3p1rho1%parameter_number
                if ( i-N2o == j ) call self%Designer3p1rho1%parameter_value( j, value )
            end do
            return

        ! parameter from rho1_ODE function
        else if ( .not. self%Designer3p1LuminalGW .and. i <= N3o ) then

            do j = 1, self%Designer3p1rho1_ODE%parameter_number
                if ( i-N3 == j ) call self%Designer3p1rho1_ODE%parameter_value( j, value )
            end do
            return

        end if

    end subroutine EFTCAMBDesigner3p1ParameterValues

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that initializes the background.
    subroutine EFTCAMBDesigner3p1InitBackground( self, params_cache, feedback_level, success, outroot )

        implicit none

        class(EFTCAMB_full3p1_wDE)                    :: self           !< the base class
        type(TEFTCAMB_parameter_cache), intent(inout) :: params_cache   !< a EFTCAMB parameter cache containing cosmological parameters
        integer                       , intent(in)    :: feedback_level !< level of feedback from the background code. 0=none; 1=some; 2=chatty.
        logical                       , intent(out)   :: success        !< wether the background initialization succeded or not
        character(LEN=*), optional    , intent(in)    :: outroot        !< the output root for the debug files

        ! some feedback:
        if ( feedback_level>1 ) then
        write(*,'(a)') "***************************************************************"
        write(*,'(a)') ' EFTCAMB 3+1 designer background solver'
        write(*,'(a)')
        end if

        call self%rhov_t%initialize( self%designer_num_points, self%x_initial, self%x_final )

        call self%Designer3p1Omega%initialize( self%designer_num_points, self%x_initial, self%x_final )
        call self%Designer3p1Gamma1%initialize( self%designer_num_points, self%x_initial, self%x_final )
        call self%Designer3p1Gamma2%initialize( self%designer_num_points, self%x_initial, self%x_final )
        call self%Designer3p1Gamma3%initialize( self%designer_num_points, self%x_initial, self%x_final )
        call self%Designer3p1Gamma5%initialize( self%designer_num_points, self%x_initial, self%x_final )
        call self%Designer3p1Gamma6%initialize( self%designer_num_points, self%x_initial, self%x_final )

        success = .True.
        ! solve the background equations and store the solution:
        call self%solve_designer_equations( params_cache, success=success, feedback_level=feedback_level )

    end subroutine EFTCAMBDesigner3p1InitBackground

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that solves for M by integrating apha_M
    subroutine EFTCAMBDesigner3p1SolveDesignerEquations( self, params_cache, success, feedback_level )

        implicit none

        class(EFTCAMB_full3p1_wDE)                   :: self          !< the base class.
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

        ydot(1) = -3._dl*(1._dl + self%Designer3p1wDE%value(a))*y(1)

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

    end subroutine EFTCAMBDesigner3p1SolveDesignerEquations

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes the value of the background EFT functions at a given time.
    subroutine EFTCAMBDesigner3p1BackgroundEFTFunctions( self, a, eft_par_cache, eft_cache )

        implicit none

        class(EFTCAMB_full3p1_wDE)                  :: self          !< the base class
        real(dl), intent(in)                         :: a             !< the input scale factor
        type(TEFTCAMB_parameter_cache), intent(inout) :: eft_par_cache !< the EFTCAMB parameter cache that contains all the physical parameters.
        type(TEFTCAMB_timestep_cache ), intent(inout) :: eft_cache     !< the EFTCAMB timestep cache that contains all the physical values.

        real(dl) :: H, Ht, Ht2, Ht3, Ht4
        real(dl) :: grhor, grhom, grhov, grhon, gpn, gdpn, gddpn, gdddpn
        real(dl) :: grhonu, gpinu, gpinudot, gpinudotdot, gpinudotdotdot, grhormass_t
        integer  :: nu_i
        real(dl) :: wDE, dwDE, d2wDE, d3wDE
        real(dl) :: Ode, dOde, d2Ode, d3Ode, d4Ode
        real(dl) :: rho1, drho1, d2rho1, d3rho1, d4rho1

        grhom = eft_cache%grhob_t + eft_cache%grhoc_t
        grhor = eft_cache%grhor_t + eft_cache%grhog_t
        grhov = eft_cache%grhov_t

        wDE   = self%Designer3p1wDE%value(a)
        dwDE  = self%Designer3p1wDE%first_derivative(a)
        d2wDE = self%Designer3p1wDE%second_derivative(a)
        d3wDE = self%Designer3p1wDE%third_derivative(a)
        
        H   = eft_cache%adotoa
        Ht  = eft_cache%Hdot
        Ht2 = eft_cache%Hdotdot
        Ht3 = eft_cache%Hdotdotdot
        if ( self%Designer3p1OmegaDE ) then
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

        if ( self%Designer3p1LuminalGW ) then
            rho1   = self%Designer3p1sigma1%value(a) + self%Designer3p1sigma1_ODE%value(Ode)
            drho1  = self%Designer3p1sigma1%first_derivative(a) + (dOde*self%Designer3p1sigma1_ODE%first_derivative(Ode))
            d2rho1  = self%Designer3p1sigma1%second_derivative(a) + (d2Ode*self%Designer3p1sigma1_ODE%first_derivative(Ode) + dOde**2*self%Designer3p1sigma1_ODE%second_derivative(Ode))
            d3rho1  = self%Designer3p1sigma1%third_derivative(a) + (d3Ode*self%Designer3p1sigma1_ODE%first_derivative(Ode) + 3._dl*d2Ode*dOde*self%Designer3p1sigma1_ODE%second_derivative(Ode) + dOde**3*self%Designer3p1sigma1_ODE%third_derivative(Ode))
            d4rho1  = self%Designer3p1sigma1%fourth_derivative(a) + (d4Ode*self%Designer3p1sigma1_ODE%first_derivative(Ode) + 3._dl*d2Ode**2*self%Designer3p1sigma1_ODE%second_derivative(Ode) + 4._dl*d3Ode*dOde*self%Designer3p1sigma1_ODE%second_derivative(Ode) + 6._dl*d2Ode*dOde**2*self%Designer3p1sigma1_ODE%third_derivative(Ode) + dOde**4*self%Designer3p1sigma1_ODE%fourth_derivative(Ode))
        else
            rho1   = self%Designer3p1rho1%value(a) + self%Designer3p1rho1_ODE%value(Ode)
            drho1  = self%Designer3p1rho1%first_derivative(a) + (dOde*self%Designer3p1rho1_ODE%first_derivative(Ode))
            d2rho1 = self%Designer3p1rho1%second_derivative(a) + (d2Ode*self%Designer3p1rho1_ODE%first_derivative(Ode) + dOde**2*self%Designer3p1rho1_ODE%second_derivative(Ode))
            d3rho1 = self%Designer3p1rho1%third_derivative(a) + (d3Ode*self%Designer3p1rho1_ODE%first_derivative(Ode) + 3._dl*d2Ode*dOde*self%Designer3p1rho1_ODE%second_derivative(Ode) + dOde**3*self%Designer3p1rho1_ODE%third_derivative(Ode))
            d4rho1 = self%Designer3p1rho1%fourth_derivative(a) + (d4Ode*self%Designer3p1rho1_ODE%first_derivative(Ode) + 3._dl*d2Ode**2*self%Designer3p1rho1_ODE%second_derivative(Ode) + 4._dl*d3Ode*dOde*self%Designer3p1rho1_ODE%second_derivative(Ode) + 6._dl*d2Ode*dOde**2*self%Designer3p1rho1_ODE%third_derivative(Ode) + dOde**4*self%Designer3p1rho1_ODE%fourth_derivative(Ode))
        end if

        eft_cache%EFTOmegaV    = rho1 - 1._dl
        eft_cache%EFTOmegaP    = drho1
        eft_cache%EFTOmegaPP   = d2rho1
        eft_cache%EFTOmegaPPP  = d3rho1
        eft_cache%EFTOmegaPPPP = d4rho1

        eft_cache%EFTc         = ( eft_cache%adotoa**2 - eft_cache%Hdot )*( eft_cache%EFTOmegaV + 0.5_dl*a*eft_cache%EFTOmegaP ) &
            & -0.5_dl*( a*eft_cache%adotoa )**2*eft_cache%EFTOmegaPP&
            & +0.5_dl*eft_cache%grhov_t*( 1._dl+ self%Designer3p1wDE%value(a) )
        eft_cache%EFTLambda    = -eft_cache%EFTOmegaV*( 2._dl*eft_cache%Hdot +eft_cache%adotoa**2 ) &
            & -a*eft_cache%EFTOmegaP*( 2._dl*eft_cache%adotoa**2 + eft_cache%Hdot ) &
            & -( a*eft_cache%adotoa )**2*eft_cache%EFTOmegaPP &
            & +self%Designer3p1wDE%value(a)*eft_cache%grhov_t
        eft_cache%EFTcdot      = +0.5_dl*eft_cache%adotoa*eft_cache%grhov_t*( -3._dl*(1._dl +self%Designer3p1wDE%value(a))**2 + a*self%Designer3p1wDE%first_derivative(a) ) &
            & -eft_cache%EFTOmegaV*( eft_cache%Hdotdot -4._dl*eft_cache%adotoa*eft_cache%Hdot +2._dl*eft_cache%adotoa**3 ) &
            & +0.5_dl*a*eft_cache%EFTOmegaP*( -eft_cache%Hdotdot +eft_cache%adotoa*eft_cache%Hdot +eft_cache%adotoa**3) &
            & +0.5_dl*a**2*eft_cache%adotoa*eft_cache%EFTOmegaPP*( eft_cache%adotoa**2 -3._dl*eft_cache%Hdot ) &
            & -0.5_dl*(a*eft_cache%adotoa)**3*eft_cache%EFTOmegaPPP
            eft_cache%EFTcdotdot   = (-(a**4*eft_cache%EFTOmegaPPPP*eft_cache%adotoa**4) - 6._dl*a**3*eft_cache%EFTOmegaPPP*eft_cache%adotoa**2*eft_cache%Hdot + a**2*eft_cache%EFTOmegaPP*(eft_cache%adotoa**4 + 4._dl*eft_cache%adotoa**2*eft_cache%Hdot - 3._dl*eft_cache%Hdot**2 - 4._dl*eft_cache%adotoa*eft_cache%Hdotdot) + a*eft_cache%EFTOmegaP*(-5._dl*eft_cache%adotoa**4 + 10._dl*eft_cache%adotoa**2*eft_cache%Hdot + eft_cache%Hdot**2 - eft_cache%Hdotdotdot) + 2._dl*(4._dl*eft_cache%adotoa**4 - 14._dl*eft_cache%adotoa**2*eft_cache%Hdot + 4._dl*eft_cache%Hdot**2 + 6._dl*eft_cache%adotoa*eft_cache%Hdotdot - eft_cache%Hdotdotdot)*eft_cache%EFTOmegaV + grhov*eft_cache%Hdot*(a*self%Designer3p1wDE%first_derivative(a) - 3._dl*(1 + self%Designer3p1wDE%value(a))**2) + grhov*eft_cache%adotoa**2*(a**2*self%Designer3p1wDE%second_derivative(a) + 9._dl*(1 + self%Designer3p1wDE%value(a))**3 - a*self%Designer3p1wDE%first_derivative(a)*(8._dl + 9._dl*self%Designer3p1wDE%value(a))))/2._dl
        eft_cache%EFTLambdadot = -2._dl*eft_cache%EFTOmegaV*( eft_cache%Hdotdot -eft_cache%adotoa*eft_cache%Hdot -eft_cache%adotoa**3 ) &
            & -a*eft_cache%EFTOmegaP*( +eft_cache%Hdotdot +5._dl*eft_cache%adotoa*eft_cache%Hdot -eft_cache%adotoa**3  ) &
            & -a**2*eft_cache%EFTOmegaPP*eft_cache%adotoa*( +2._dl*eft_cache%adotoa**2 +3._dl*eft_cache%Hdot )&
            & -(a*eft_cache%adotoa)**3*eft_cache%EFTOmegaPPP &
            & +eft_cache%grhov_t*eft_cache%adotoa*( a*self%Designer3p1wDE%first_derivative(a) -3._dl*self%Designer3p1wDE%value(a)*(1._dl +self%Designer3p1wDE%value(a) ))
        eft_cache%EFTLambdadotdot = -(a**4*eft_cache%EFTOmegaPPPP*eft_cache%adotoa**4) + eft_cache%EFTOmegaP*(-(a*eft_cache%Hdotdotdot) - 5._dl*a*eft_cache%Hdot**2 - 6._dl*a*eft_cache%Hdotdot*eft_cache%adotoa + 10._dl*a*eft_cache%Hdot*eft_cache%adotoa**2 + a*eft_cache%adotoa**4) + eft_cache%EFTOmegaPP*(-3._dl*a**2*eft_cache%Hdot**2 - 4._dl*a**2*eft_cache%Hdotdot*eft_cache%adotoa - 11._dl*a**2*eft_cache%Hdot*eft_cache%adotoa**2 + a**2*eft_cache%adotoa**4) + eft_cache%EFTOmegaPPP*(-6._dl*a**3*eft_cache%Hdot*eft_cache%adotoa**2 - 3._dl*a**3*eft_cache%adotoa**4) + (-2._dl*eft_cache%Hdotdotdot + 2._dl*eft_cache%Hdot**2 + 6._dl*eft_cache%Hdotdot*eft_cache%adotoa + 2._dl*eft_cache%Hdot*eft_cache%adotoa**2 - 4._dl*eft_cache%adotoa**4)*eft_cache%EFTOmegaV + eft_cache%Hdot*eft_cache%grhov_t*(a*self%Designer3p1wDE%first_derivative(a) - 3._dl*self%Designer3p1wDE%value(a) - 3._dl*self%Designer3p1wDE%value(a)**2) + eft_cache%adotoa**2*eft_cache%grhov_t*(a**2*self%Designer3p1wDE%second_derivative(a) + 9._dl*self%Designer3p1wDE%value(a) + 18._dl*self%Designer3p1wDE%value(a)**2 + 9._dl*self%Designer3p1wDE%value(a)**3 - a*self%Designer3p1wDE%first_derivative(a)*(5._dl + 9._dl*self%Designer3p1wDE%value(a)))

    end subroutine EFTCAMBDesigner3p1BackgroundEFTFunctions

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes the value of the background EFT functions at a given time.
    subroutine EFTCAMBDesigner3p1SecondOrderEFTFunctions( self, a, eft_par_cache, eft_cache )

        implicit none

        class(EFTCAMB_full3p1_wDE)                  :: self          !< the base class
        real(dl), intent(in)                         :: a             !< the input scale factor
        type(TEFTCAMB_parameter_cache), intent(inout) :: eft_par_cache !< the EFTCAMB parameter cache that contains all the physical parameters.
        type(TEFTCAMB_timestep_cache ), intent(inout) :: eft_cache     !< the EFTCAMB timestep cache that contains all the physical values.

        real(dl) :: H, Ht, Ht2, Ht3, Ht4
        real(dl) :: grhor, grhom, grhov, grhon, gpn, gdpn, gddpn, gdddpn
        real(dl) :: grhonu, gpinu, gpinudot, gpinudotdot, gpinudotdotdot, grhormass_t
        integer  :: nu_i
        real(dl) :: wDE, dwDE, d2wDE, d3wDE
        real(dl) :: Ode, dOde, d2Ode, d3Ode, d4Ode
        real(dl) :: sigma1, dsigma1, d2sigma1, d3sigma1, d4sigma1, rho1, drho1, d2rho1, d3rho1, d4rho1, rho1N, drho1N, rho2, drho2
        real(dl) :: at, dat, d2at, d3at

        grhom = eft_cache%grhob_t + eft_cache%grhoc_t
        grhor = eft_cache%grhor_t + eft_cache%grhog_t
        grhov = eft_cache%grhov_t

        wDE   = self%Designer3p1wDE%value(a)
        dwDE  = self%Designer3p1wDE%first_derivative(a)
        d2wDE = self%Designer3p1wDE%second_derivative(a)
        d3wDE = self%Designer3p1wDE%third_derivative(a)
        
        H   = eft_cache%adotoa
        Ht  = eft_cache%Hdot
        Ht2 = eft_cache%Hdotdot
        Ht3 = eft_cache%Hdotdotdot
        if ( self%Designer3p1OmegaDE ) then
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

        at   = a*eft_cache%adotoa
        dat  = eft_cache%Hdot/eft_cache%adotoa + eft_cache%adotoa
        d2at = (-eft_cache%Hdot**2 + eft_cache%Hdotdot*eft_cache%adotoa + eft_cache%Hdot*eft_cache%adotoa**2)/(a*eft_cache%adotoa**3)
        d3at = (3._dl*eft_cache%Hdot**3 - 4._dl*eft_cache%Hdotdot*eft_cache%Hdot*eft_cache%adotoa + eft_cache%Hdotdotdot*eft_cache%adotoa**2 - eft_cache%Hdot*eft_cache%adotoa**4)/(a**2*eft_cache%adotoa**5)

        sigma1   = self%Designer3p1sigma1%value(a) + self%Designer3p1sigma1_ODE%value(Ode)
        dsigma1  = self%Designer3p1sigma1%first_derivative(a) + (dOde*self%Designer3p1sigma1_ODE%first_derivative(Ode))
        d2sigma1 = self%Designer3p1sigma1%second_derivative(a) + (d2Ode*self%Designer3p1sigma1_ODE%first_derivative(Ode) + dOde**2*self%Designer3p1sigma1_ODE%second_derivative(Ode))
        d3sigma1 = self%Designer3p1sigma1%third_derivative(a) + (d3Ode*self%Designer3p1sigma1_ODE%first_derivative(Ode) + 3._dl*d2Ode*dOde*self%Designer3p1sigma1_ODE%second_derivative(Ode) + dOde**3*self%Designer3p1sigma1_ODE%third_derivative(Ode))
        d4sigma1 = self%Designer3p1sigma1%fourth_derivative(a) + (d4Ode*self%Designer3p1sigma1_ODE%first_derivative(Ode) + 3._dl*d2Ode**2*self%Designer3p1sigma1_ODE%second_derivative(Ode) + 4._dl*d3Ode*dOde*self%Designer3p1sigma1_ODE%second_derivative(Ode) + 6._dl*d2Ode*dOde**2*self%Designer3p1sigma1_ODE%third_derivative(Ode) + dOde**4*self%Designer3p1sigma1_ODE%fourth_derivative(Ode))

        rho1N  = self%Designer3p1rho1N%value(a) + self%Designer3p1rho1N_ODE%value(Ode)
        drho1N = self%Designer3p1rho1N%first_derivative(a) + (dOde*self%Designer3p1rho1N_ODE%first_derivative(Ode))

        rho2   = self%Designer3p1rho2%value(a) + self%Designer3p1rho2_ODE%value(Ode)
        drho2  = self%Designer3p1rho2%first_derivative(a) + (dOde*self%Designer3p1rho2_ODE%first_derivative(Ode))

        if ( self%Designer3p1LuminalGW ) then
            rho1   = sigma1
            drho1  = dsigma1
            d2rho1 = d2sigma1
            d3rho1 = d3sigma1
            d4rho1 = d4sigma1
        else
            rho1   = self%Designer3p1rho1%value(a) + self%Designer3p1rho1_ODE%value(Ode)
            drho1  = self%Designer3p1rho1%first_derivative(a) + (dOde*self%Designer3p1rho1_ODE%first_derivative(Ode))
            d2rho1 = self%Designer3p1rho1%second_derivative(a) + (d2Ode*self%Designer3p1rho1_ODE%first_derivative(Ode) + dOde**2*self%Designer3p1rho1_ODE%second_derivative(Ode))
            d3rho1 = self%Designer3p1rho1%third_derivative(a) + (d3Ode*self%Designer3p1rho1_ODE%first_derivative(Ode) + 3._dl*d2Ode*dOde*self%Designer3p1rho1_ODE%second_derivative(Ode) + dOde**3*self%Designer3p1rho1_ODE%third_derivative(Ode))
            d4rho1 = self%Designer3p1rho1%fourth_derivative(a) + (d4Ode*self%Designer3p1rho1_ODE%first_derivative(Ode) + 3._dl*d2Ode**2*self%Designer3p1rho1_ODE%second_derivative(Ode) + 4._dl*d3Ode*dOde*self%Designer3p1rho1_ODE%second_derivative(Ode) + 6._dl*d2Ode*dOde**2*self%Designer3p1rho1_ODE%third_derivative(Ode) + dOde**4*self%Designer3p1rho1_ODE%fourth_derivative(Ode))
        end if

        eft_cache%EFTGamma1V    = - 0.5_dl*( eft_cache%EFTc )/(a*eft_par_cache%h0_Mpc)**2
        eft_cache%EFTGamma1P    = - 0.5_dl*( eft_cache%EFTcdot/at )/(a*eft_par_cache%h0_Mpc)**2
        eft_cache%EFTGamma1PP   = - 0.5_dl*( eft_cache%EFTcdotdot/at**2 - eft_cache%EFTcdot*dat/at**2 )/(a*eft_par_cache%h0_Mpc)**2
        eft_cache%EFTGamma2V    = - ( drho1*at )/eft_par_cache%h0_Mpc
        eft_cache%EFTGamma2P    = - ( at*d2rho1 + dat*drho1 )/eft_par_cache%h0_Mpc
        eft_cache%EFTGamma2PP   = - ( at*d3rho1 + 2._dl*d2rho1*dat + d2at*drho1 )/eft_par_cache%h0_Mpc
        eft_cache%EFTGamma2PPP  = - ( 3._dl*d2at*d2rho1 + at*d4rho1 + 3._dl*d3rho1*dat + d3at*drho1 )/eft_par_cache%h0_Mpc
        eft_cache%EFTGamma3V    = 1._dl/3._dl*sigma1 - rho1 + 2._dl/3._dl 
        eft_cache%EFTGamma3P    = 1._dl/3._dl*dsigma1 - drho1 
        eft_cache%EFTGamma3PP   = 1._dl/3._dl*d2sigma1 - d2rho1 
        eft_cache%EFTGamma3PPP  = 1._dl/3._dl*d3sigma1 - d3rho1
        eft_cache%EFTGamma3PPPP = 1._dl/3._dl*d4sigma1 - d4rho1
        eft_cache%EFTGamma4V    = rho1 - sigma1
        eft_cache%EFTGamma4P    = drho1 - dsigma1
        eft_cache%EFTGamma4PP   = d2rho1 - d2sigma1
        eft_cache%EFTGamma5V    = 0.5_dl*rho1N
        eft_cache%EFTGamma5P    = 0.5_dl*drho1N
        eft_cache%EFTGamma6V    = rho2/8._dl
        eft_cache%EFTGamma6P    = drho2/8._dl

    end subroutine EFTCAMBDesigner3p1SecondOrderEFTFunctions

    ! ---------------------------------------------------------------------------------------------
    !> Function that computes dtauda = 1/sqrt(a^2H^2).For pure EFT std this has to be overridden
    !! for performance reasons.
    function EFTCAMBDesigner3p1ComputeDtauda( self, a, eft_par_cache, eft_cache )

        implicit none

        class(EFTCAMB_full3p1_wDE)                  :: self          !< the base class
        real(dl), intent(in)                         :: a             !< the input scale factor
        type(TEFTCAMB_parameter_cache), intent(inout) :: eft_par_cache !< the EFTCAMB parameter cache that contains all the physical parameters.
        type(TEFTCAMB_timestep_cache ), intent(inout) :: eft_cache     !< the EFTCAMB timestep cache that contains all the physical values.

        real(dl) :: EFTCAMBDesigner3p1ComputeDtauda               !< the output dtauda

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
        EFTCAMBDesigner3p1ComputeDtauda = sqrt(3._dl/temp)

    end function EFTCAMBDesigner3p1ComputeDtauda

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes adotoa = H and its two derivatives wrt conformal time.
    subroutine EFTCAMBDesigner3p1ComputeAdotoa( self, a, eft_par_cache, eft_cache )

        implicit none

        class(EFTCAMB_full3p1_wDE)                  :: self          !< the base class
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

    end subroutine EFTCAMBDesigner3p1ComputeAdotoa

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes the three derivatives wrt conformal time of H.
    subroutine EFTCAMBDesigner3p1ComputeHubbleDer( self, a, eft_par_cache, eft_cache )

        implicit none

        class(EFTCAMB_full3p1_wDE)                  :: self          !< the base class
        real(dl), intent(in)                         :: a             !< the input scale factor
        type(TEFTCAMB_parameter_cache), intent(inout) :: eft_par_cache !< the EFTCAMB parameter cache that contains all the physical parameters.
        type(TEFTCAMB_timestep_cache ), intent(inout) :: eft_cache     !< the EFTCAMB timestep cache that contains all the physical values.

        eft_cache%gpiv_t  = self%Designer3p1wDE%value(a)*eft_cache%grhov_t
        eft_cache%Hdot    = -0.5_dl*( eft_cache%adotoa**2 +eft_cache%gpresm_t +eft_cache%gpiv_t )
        eft_cache%Hdotdot = eft_cache%adotoa*( ( eft_cache%grhob_t +eft_cache%grhoc_t)/6._dl +2._dl*( eft_cache%grhor_t +eft_cache%grhog_t)/3._dl ) &
                          & +eft_cache%adotoa*eft_cache%grhov_t*( 1._dl/6._dl +self%Designer3p1wDE%value(a) +1.5_dl*self%Designer3p1wDE%value(a)**2 -0.5_dl*a*self%Designer3p1wDE%first_derivative(a) ) &
                          & +eft_cache%adotoa*eft_cache%grhonu_tot/6._dl -0.5_dl*eft_cache%adotoa*eft_cache%gpinu_tot -0.5_dl*eft_cache%gpinudot_tot
        eft_cache%Hdotdotdot = (2._dl*(eft_cache%grhob_t + eft_cache%grhoc_t)*eft_cache%adotoa**2 + eft_cache%Hdot &
                          &* (eft_cache%grhob_t + eft_cache%grhoc_t) - 3._dl*eft_cache%adotoa**2 *(eft_cache%grhob_t + eft_cache%grhoc_t))/(6._dl) &
                          &+ (2._dl/3._dl)*(2._dl*eft_cache%adotoa**2 * (eft_cache%grhor_t + eft_cache%grhog_t) + eft_cache%Hdot &
                          &* (eft_cache%grhor_t + eft_cache%grhog_t) - 4._dl*eft_cache%adotoa**2 * (eft_cache%grhor_t + eft_cache%grhog_t)) &
                          &+ (2._dl*eft_cache%adotoa**2 * eft_cache%grhov_t - 3._dl*eft_cache%adotoa**2 * eft_cache%grhov_t &
                          &* (1._dl + self%Designer3p1wDE%value(a)) + eft_cache%Hdot*eft_cache%grhov_t)*(1._dl/6._dl + self%Designer3p1wDE%value(a) &
                          &+ (3._dl/2._dl)*self%Designer3p1wDE%value(a)**2 - (1._dl/2._dl)*a*self%Designer3p1wDE%first_derivative(a)) &
                          &+ a*eft_cache%adotoa**2 * eft_cache%grhov_t * ((1._dl/2._dl)*self%Designer3p1wDE%first_derivative(a) &
                          &+ 3._dl*self%Designer3p1wDE%value(a)*self%Designer3p1wDE%first_derivative(a) - (a/2._dl)*self%Designer3p1wDE%second_derivative(a)) &
                          &+ (eft_cache%adotoa**2/3._dl)*eft_cache%grhonu_tot - eft_cache%adotoa**2 * eft_cache%gpinu_tot &
                          &- (3._dl/2._dl)*eft_cache%adotoa*eft_cache%gpinudot_tot + (eft_cache%Hdot/6._dl)*eft_cache%grhonu_tot &
                          &+ (eft_cache%adotoa/6._dl)*eft_cache%grhonudot_tot - (eft_cache%Hdot/2._dl)*eft_cache%gpinu_tot &
                          &- (1._dl/2._dl)*eft_cache%gpinudotdot_tot
    end subroutine EFTCAMBDesigner3p1ComputeHubbleDer

    ! ---------------------------------------------------------------------------------------------
    !> Function that computes model specific stability requirements.
    function EFTCAMBDesigner3p1AdditionalModelStability( self, a, eft_par_cache, eft_cache )

        implicit none

        class(EFTCAMB_full3p1_wDE)                  :: self          !< the base class
        real(dl), intent(in)                         :: a             !< the input scale factor.
        type(TEFTCAMB_parameter_cache), intent(inout) :: eft_par_cache !< the EFTCAMB parameter cache that contains all the physical parameters.
        type(TEFTCAMB_timestep_cache ), intent(inout) :: eft_cache     !< the EFTCAMB timestep cache that contains all the physical values.

        logical :: EFTCAMBDesigner3p1AdditionalModelStability          !< the return value of the stability computation. True if the model specific stability criteria are met, false otherwise.
                    
        real(dl) :: lgkmin, lgkmax, lgk, k, k2, k4
        real(dl) :: H, Ht
        real(dl) :: grhov
        real(dl) :: wDE
        real(dl) :: Ode, dOde
        real(dl) :: sigma1, dsigma1, rho1, drho1, rho1N, drho1N, rho2, drho2
        real(dl) :: dadt
        real(dl) :: aT, aL, aLt, aH, aHt, b3, b3t, B0, B01, B02, B03, B2, B21, B22, B23, B4, Azz, Bzz

        EFTCAMBDesigner3p1AdditionalModelStability = .True.
        
        if ( self%Designer3p1wDE%value(a) > -1._dl/3._dl ) then
            EFTCAMBDesigner3p1AdditionalModelStability = .False.
            return
        end if
        
        if (self%Designer3p1GhostStability .or. self%Designer3p1GradientStability) then
            
            lgkmin = log10(self%Designer3p1Stability_kmin/a)
            lgkmax = log10(self%Designer3p1Stability_kmax/a)
            lgk = lgkmin
            
            do while ( lgk < lgkmax )
                k = 10._dl**lgk
                k2 = k**2
                k4 = k2*k2

                grhov = eft_cache%grhov_t
                wDE   = self%Designer3p1wDE%value(a)
                H   = eft_cache%adotoa
                Ht  = eft_cache%Hdot
                if ( self%Designer3p1OmegaDE ) then
                    Ode    = grhov / 3._dl / H**2
                    dOde   = -(grhov*(2._dl*Ht + H**2*(1._dl + 3._dl*wDE)))/(3._dl*a*H**4)
                else
                    Ode = 0._dl
                    dOde = 0._dl
                end if

                dadt   = eft_cache%adotoa

                sigma1   = self%Designer3p1sigma1%value(a) + self%Designer3p1sigma1_ODE%value(Ode)
                dsigma1  = dadt*(self%Designer3p1sigma1%first_derivative(a) + (dOde*self%Designer3p1sigma1_ODE%first_derivative(Ode)))

                rho1N  = self%Designer3p1rho1N%value(a) + self%Designer3p1rho1N_ODE%value(Ode)
                drho1N = dadt*(self%Designer3p1rho1N%first_derivative(a) + (dOde*self%Designer3p1rho1N_ODE%first_derivative(Ode)))

                rho2   = self%Designer3p1rho2%value(a) + self%Designer3p1rho2_ODE%value(Ode)
                drho2  = dadt*(self%Designer3p1rho2%first_derivative(a) + (dOde*self%Designer3p1rho2_ODE%first_derivative(Ode)))

                ! c_T = 1
                if ( self%Designer3p1LuminalGW ) then
                    rho1   = sigma1
                    drho1  = dsigma1
                else
                    rho1   = self%Designer3p1rho1%value(a) + self%Designer3p1rho1_ODE%value(Ode)
                    drho1  = dadt*(self%Designer3p1rho1%first_derivative(a) + (dOde*self%Designer3p1rho1_ODE%first_derivative(Ode)))
                end if
                
                aT = rho1/sigma1 - 1._dl
                aL = 1._dl/sigma1 - 1._dl
                aLt = - dsigma1/sigma1**2
                aH = (rho1 + rho1N)/sigma1 - 1._dl
                aHt = (drho1 + drho1N)/sigma1 - dsigma1*(rho1 + rho1N)/sigma1**2
                b3 = rho2/sigma1
                b3t = drho2/sigma1 - dsigma1*rho2/sigma1**2
                
                ! The following are defined using physical H and proper time
                H = H/a
                Ht = Ht/a/a - H**2

                ! custom ghost stability for this beyond GLPV model
                if (self%Designer3p1GhostStability) then
                    Azz = 6._dl*(1._dl + aL)*b3*k2/(6._dl*H**2*(1._dl + aL) + aL*b3*k2)
                    if (Azz <= 0._dl ) then
                        EFTCAMBDesigner3p1AdditionalModelStability = .false.
                        return
                    end if
                end if
                ! custom gradient stability for this beyond GLPV model
                if (self%Designer3p1GradientStability) then
                    B01 = 72._dl*(1._dl + aL)*(aHt*(1._dl + aL) - aLt*(1._dl + aH))
                    B02 = 72._dl*(1._dl + aL)**2*(aH - aT)
                    B03 = -72._dl*(1._dl + aL)**2*(1._dl + aH)
                    B0 = (B03*Ht + B02*H**2 + B01*H)*H**2
                    B21 = 12._dl*(1._dl + aL)*(aHt*aL*b3 - (1._dl + aH)*(b3*aLt + aL*b3t))
                    B22 = 24._dl*aL*(1._dl + aL)*((1._dl + aH)**2 + b3*(1.5_dl*aH - aT))
                    B23 = 12._dl*b3*aL*(1._dl + aL)*(1._dl + aH)
                    B2 = B23*Ht + B22*H**2 + B21*H
                    B4 = 2._dl*aL**2*b3*(2._dl*(1._dl + aH)**2 - b3*(1._dl + aT))
                    Bzz = (B0 + B2*k2 + B4*k4)/(6._dl*H**2*(1._dl + aL) + aL*b3*k2)**2
                    if (Bzz < 0._dl ) then
                        EFTCAMBDesigner3p1AdditionalModelStability = .false.
                        return
                    end if
                end if
                lgk = lgk + self%Designer3p1Stability_dlgk
            end do

        end if

    end function EFTCAMBDesigner3p1AdditionalModelStability

    ! ---------------------------------------------------------------------------------------------

end module EFTCAMB_designer_full3p1_wDE

!----------------------------------------------------------------------------------------
