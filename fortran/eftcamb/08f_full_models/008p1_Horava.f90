! ----------------------------------------------------------------------------------------
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
! ----------------------------------------------------------------------------------------

!> @file 008p1_Horava.f90
!! This file contains the definition of low energy (LE) Horava gravity.
!! Please refer to the numerical notes for details.


! ----------------------------------------------------------------------------------------
!> This module contains the definition of low energy (LE) Horava gravity.
!! Please refer to the numerical notes for details.

!> @author Bin Hu, Marco Raveri

    module EFTCAMB_LE_Horava

    use precision
    use config, only: GlobalError, error_unsupported_params
    use IniObjects
    use MpiUtils
    use EFTCAMB_cache
    use EFT_def
    use EFTCAMB_abstract_model_designer

    implicit none

    private

    public EFTCAMB_Horava

    ! ----------------------------------------------------------------------------------------
    !> This is the type that contains the definition of low energy Horava gravity.
    type, extends (EFTCAMB_designer_model) :: EFTCAMB_Horava

        ! the model selection flag:
        !< Selects whether to consider the version of Horava gravity that evades solar system constraints.
        logical  :: HoravaSolarSystem
        !! Notice that when this is used the parameter Horava_xi is ignored.

        ! the model parameters:
        real(dl) :: Horava_eta      !< first Horava model parameter.
        real(dl) :: Horava_xi       !< second Horava model parameter.
        real(dl) :: Horava_lambda   !< third Horava model parameter.

    contains

    ! initialization of the model:
    !< subroutine that reads the parameters of the model from file
    procedure :: read_model_selection => EFTCAMBHoravaReadModelSelectionFromFile
    !< subroutine that allocates the model selection. For Horava this is a dummy procedure.
    procedure :: allocate_model_selection => EFTCAMBHoravaAllocateModelSelection
    !< subroutine that initializes the model parameters based on the values found in an input array.
    procedure :: init_model_parameters => EFTCAMBHoravaInitModelParameters
    !< subroutine that reads the parameters of the model from file.
    procedure :: init_model_parameters_from_file => EFTCAMBHoravaInitModelParametersFromFile
    ! utility functions:
    !< subroutine that computes the number of parameters of the model.
    procedure :: compute_param_number => EFTCAMBHoravaComputeParametersNumber
    !< subroutine that prints on the screen feedback information about the model.
    procedure :: feedback => EFTCAMBHoravaFeedback
    !< subroutine that returns the i-th parameter name of the model.
    procedure :: parameter_names => EFTCAMBHoravaParameterNames
    !< subroutine that returns the i-th parameter name of the model.
    procedure :: parameter_names_latex => EFTCAMBHoravaParameterNamesLatex
    !< subroutine that returns the i-th parameter value.
    procedure :: parameter_values => EFTCAMBHoravaParameterValues
    ! CAMB related procedures:
    !< subroutine that computes the value of the background EFT functions at a given time.
    procedure :: compute_background_EFT_functions => EFTCAMBHoravaBackgroundEFTFunctions
    !< subroutine that computes the value of the second order EFT functions at a given time.
    procedure :: compute_secondorder_EFT_functions => EFTCAMBHoravaSecondOrderEFTFunctions
    !< function that computes dtauda = 1/sqrt(a^2H^2).
    procedure :: compute_dtauda => EFTCAMBHoravaComputeDtauda
    !< subroutine that computes adotoa = H and its two derivatives wrt conformal time.
    procedure :: compute_adotoa => EFTCAMBHoravaComputeAdotoa
    !< subroutine that computes the two derivatives wrt conformal time of H.
    procedure :: compute_H_derivs => EFTCAMBHoravaComputeHubbleDer
    !< subroutine that computes the pi field equations factors. For details refer to the numerical notes. Notice that in Horava we override it for reasons of numerical stability.
    procedure :: compute_pi_factors => EFTCAMBHoravaComputePiFactors
    ! stability procedures:
    !< function that computes model specific stability requirements.
    procedure :: additional_model_stability => EFTCAMBHoravaAdditionalModelStability

    end type EFTCAMB_Horava

    ! ---------------------------------------------------------------------------------------------

    contains

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that reads the parameters of the model from file.
    subroutine EFTCAMBHoravaReadModelSelectionFromFile(self, Ini, eft_error)

    implicit none

    class(EFTCAMB_Horava) :: self      !< the base class
    class(TIniFile)       :: Ini       !< Input ini file
    integer               :: eft_error !< error code: 0 all fine, 1 initialization failed

    ! read model selection flags:
    self%HoravaSolarSystem = Ini%Read_Logical('HoravaSolarSystem', .false.)

    end subroutine EFTCAMBHoravaReadModelSelectionFromFile

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that allocates the model selection. For Horava this is a dummy procedure.
    subroutine EFTCAMBHoravaAllocateModelSelection(self, Ini, eft_error)

    implicit none

    class(EFTCAMB_Horava) :: self      !< the base class
    class(TIniFile)       :: Ini       !< Input ini file
    integer               :: eft_error !< error code: 0 all fine, 1 initialization failed

    end subroutine EFTCAMBHoravaAllocateModelSelection

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that initializes the model parameters based on the values found in an input array.
    subroutine EFTCAMBHoravaInitModelParameters(self, array)

    implicit none

    class(EFTCAMB_Horava)                                  :: self   !< the base class
    !< input array with the values of the parameters of the model.
    real(dl), dimension(self%parameter_number), intent(in) :: array

    real(dl), allocatable, dimension(:) :: temp
    integer :: num_params_function, num_params_temp, i

    if (self%HoravaSolarSystem) then
        self%Horava_eta = array(1)
        self%Horava_lambda = array(2)
        self%Horava_xi = 0.5_dl*self%Horava_eta
    else
        self%Horava_eta = array(1)
        self%Horava_lambda = array(2)
        self%Horava_xi = array(3)
    end if

    end subroutine EFTCAMBHoravaInitModelParameters

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that reads the parameters of the model from file.
    subroutine EFTCAMBHoravaInitModelParametersFromFile(self, Ini, eft_error)

    implicit none

    class(EFTCAMB_Horava) :: self      !< the base class
    class(TIniFile)       :: Ini       !< Input ini file
    integer               :: eft_error !< error code: 0 all fine, 1 initialization failed

    self%Horava_eta = Ini%Read_Double('Horava_eta', 0._dl)
    self%Horava_lambda = Ini%Read_Double('Horava_lambda', 0._dl)
    if (self%HoravaSolarSystem) then
        self%Horava_xi = 0.5_dl*self%Horava_eta
    else
        self%Horava_xi = Ini%Read_Double('Horava_xi', 0._dl)
    end if

    end subroutine EFTCAMBHoravaInitModelParametersFromFile

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes the number of parameters of the model.
    subroutine EFTCAMBHoravaComputeParametersNumber(self)

    implicit none

    class(EFTCAMB_Horava) :: self   !< the base class

    if (self%HoravaSolarSystem) then
        self%parameter_number = 2
    else
        self%parameter_number = 3
    end if

    end subroutine EFTCAMBHoravaComputeParametersNumber

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that prints on the screen feedback information about the model.
    subroutine EFTCAMBHoravaFeedback(self, print_params)

    implicit none

    class(EFTCAMB_Horava) :: self         !< the base class
    logical, optional     :: print_params !< optional flag that decised whether to print numerical values
    !! of the parameters.

    logical               :: print_params_temp

    ! print general model informations:
    write(*, *)
    write(*, '(a,a)')    '   Model               =  ', self%name
    write(*, '(a,I3)')   '   Number of params    =', self%parameter_number

    if (self%HoravaSolarSystem) then
        write(*, '(a)')  '   Horava with solar system constraints'
    end if

    ! print the values of the parameters:
    if (present(print_params)) then
        print_params_temp = print_params
    else
        print_params_temp = .true.
    end if

    if (print_params_temp) then
        write(*, *)
        write(*, '(a23,a,F12.6)') '   Horava_eta          ', '=', self%Horava_eta
        write(*, '(a23,a,F12.6)') '   Horava_lambda       ', '=', self%Horava_lambda
        if (.not. self%HoravaSolarSystem) then
            write(*, '(a23,a,F12.6)') '   Horava_xi           ', '=', self%Horava_xi
        end if
    end if

    end subroutine EFTCAMBHoravaFeedback

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the i-th parameter name of the model
    subroutine EFTCAMBHoravaParameterNames(self, i, name)

    implicit none

    class(EFTCAMB_Horava) :: self   !< the base class
    integer, intent(in)   :: i      !< the index of the parameter
    character(*), intent(out) :: name   !< the output name of the i-th parameter

    ! check the input index:
    if (i > self%parameter_number) then
        write(*, *) 'Illegal index for parameter_names.'
        write(*, *) 'Maximum value is:', self%parameter_number
        call GlobalError('EFTCAMB error: illegal parameter index in EFTCAMBHoravaParameterNames', &
            error_unsupported_params)
        name = ''
        return
    end if
    ! return the appropriate name:
    if (i == 1) then
        name = 'Horava_eta'
        return
    end if
    if (i == 2) then
        name = 'Horava_lambda'
        return
    end if
    if (i == 3 .and. .not. self%HoravaSolarSystem) then
        name = 'Horava_xi'
        return
    end if

    end subroutine EFTCAMBHoravaParameterNames

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the i-th parameter name of the model
    subroutine EFTCAMBHoravaParameterNamesLatex(self, i, latexname)

    implicit none

    class(EFTCAMB_Horava) :: self       !< the base class
    integer, intent(in)   :: i          !< The index of the parameter
    character(*), intent(out) :: latexname  !< the output latex name of the i-th parameter

    ! check the input index:
    if (i > self%parameter_number) then
        write(*, *) 'Illegal index for parameter_names_latex.'
        write(*, *) 'Maximum value is:', self%parameter_number
        call GlobalError('EFTCAMB error: illegal parameter index in EFTCAMBHoravaParameterNamesLatex', &
            error_unsupported_params)
        latexname = ''
        return
    end if
    ! return the appropriate name:
    if (i == 1) then
        latexname = '\eta_{\rm Ho\v rava}'
        return
    end if
    if (i == 2) then
        latexname = '\lambda_{\rm Ho\v rava}'
        return
    end if
    if (i == 3 .and. .not. self%HoravaSolarSystem) then
        latexname = '\xi_{\rm Ho\v rava}'
        return
    end if

    end subroutine EFTCAMBHoravaParameterNamesLatex

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the i-th parameter name of the model
    subroutine EFTCAMBHoravaParameterValues(self, i, value)

    implicit none

    class(EFTCAMB_Horava) :: self   !< the base class
    integer, intent(in)   :: i      !< The index of the parameter
    real(dl), intent(out) :: value  !< the output value of the i-th parameter

    ! check the input index:
    if (i > self%parameter_number) then
        write(*, *) 'Illegal index for parameter_value.'
        write(*, *) 'Maximum value is:', self%parameter_number
        call GlobalError('EFTCAMB error: illegal parameter index in EFTCAMBHoravaParameterValues', &
            error_unsupported_params)
        value = 0._dl
        return
    end if
    ! return the appropriate name:
    if (i == 1) then
        value = self%Horava_eta
        return
    end if
    if (i == 2) then
        value = self%Horava_lambda
        return
    end if
    if (i == 3 .and. .not. self%HoravaSolarSystem) then
        value = self%Horava_xi
        return
    end if

    end subroutine EFTCAMBHoravaParameterValues

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes the value of the background EFT functions at a given time.
    subroutine EFTCAMBHoravaBackgroundEFTFunctions(self, a, eft_par_cache, eft_cache)

    implicit none

    class(EFTCAMB_Horava) :: self          !< the base class
    real(dl), intent(in)  :: a             !< the input scale factor
    !< the EFTCAMB parameter cache that contains all the physical parameters.
    type(TEFTCAMB_parameter_cache), intent(inout) :: eft_par_cache
    !< the EFTCAMB timestep cache that contains all the physical values.
    type(TEFTCAMB_timestep_cache), intent(inout) :: eft_cache

    ! compute the EFT functions:
    eft_cache%EFTOmegaV = self%Horava_eta/(2._dl + 2._dl*self%Horava_xi - self%Horava_eta)
    eft_cache%EFTOmegaP = 0._dl
    eft_cache%EFTOmegaPP = 0._dl
    eft_cache%EFTOmegaPPP = 0._dl
    eft_cache%EFTc = -(2._dl*self%Horava_xi - 3._dl*self%Horava_lambda)*(eft_cache%Hdot - eft_cache%adotoa**2)/ &
        (2._dl + 2._dl*self%Horava_xi - self%Horava_eta)
    eft_cache%EFTLambda = -eft_cache%grhov_t + 2._dl*(3._dl*self%Horava_lambda - 2._dl*self%Horava_xi)* &
        (0.5_dl*eft_cache%adotoa**2 + eft_cache%Hdot)/(2._dl + 2._dl*self%Horava_xi - self%Horava_eta)
    eft_cache%EFTcdot = (3._dl*self%Horava_lambda - 2._dl*self%Horava_xi)* &
        (eft_cache%Hdotdot - 4._dl*eft_cache%adotoa*eft_cache%Hdot + 2._dl*eft_cache%adotoa**3)/ &
        (2._dl + 2._dl*self%Horava_xi - self%Horava_eta)
    eft_cache%EFTLambdadot = +2._dl*(3._dl*self%Horava_lambda - 2._dl*self%Horava_xi)* &
        (eft_cache%Hdotdot - eft_cache%adotoa*eft_cache%Hdot - eft_cache%adotoa**3)/ &
        (2._dl + 2._dl*self%Horava_xi - self%Horava_eta)

    end subroutine EFTCAMBHoravaBackgroundEFTFunctions

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes the value of the background EFT functions at a given time.
    subroutine EFTCAMBHoravaSecondOrderEFTFunctions(self, a, eft_par_cache, eft_cache)

    implicit none

    class(EFTCAMB_Horava) :: self          !< the base class
    real(dl), intent(in)  :: a             !< the input scale factor
    !< the EFTCAMB parameter cache that contains all the physical parameters.
    type(TEFTCAMB_parameter_cache), intent(inout) :: eft_par_cache
    !< the EFTCAMB timestep cache that contains all the physical values.
    type(TEFTCAMB_timestep_cache), intent(inout) :: eft_cache

    ! compute the EFT functions:
    eft_cache%EFTGamma1V = (2._dl*self%Horava_xi - 3._dl*self%Horava_lambda)/ &
        (2._dl + 2._dl*self%Horava_xi - self%Horava_eta)/(2._dl*a**2*eft_par_cache%h0_Mpc**2)* &
        (eft_cache%Hdot - eft_cache%adotoa**2)
    eft_cache%EFTGamma1P = (2._dl*self%Horava_xi - 3._dl*self%Horava_lambda)/ &
        (2._dl + 2._dl*self%Horava_xi - self%Horava_eta)/(2._dl*a**3*eft_par_cache%h0_Mpc**2)* &
        (eft_cache%Hdotdot/eft_cache%adotoa - 4._dl*eft_cache%Hdot + 2._dl*eft_cache%adotoa**2)
    eft_cache%EFTGamma2V = 0._dl
    eft_cache%EFTGamma2P = 0._dl
    eft_cache%EFTGamma3V = 2._dl*(self%Horava_lambda - self%Horava_xi)/(2._dl + 2._dl*self%Horava_xi - self%Horava_eta)
    eft_cache%EFTGamma3P = 0._dl
    eft_cache%EFTGamma4V = 2._dl*self%Horava_xi/(2._dl + 2._dl*self%Horava_xi - self%Horava_eta)
    eft_cache%EFTGamma4P = 0._dl
    eft_cache%EFTGamma4PP = 0._dl
    eft_cache%EFTGamma5V = 0._dl
    eft_cache%EFTGamma5P = 0._dl
    eft_cache%EFTGamma6V = self%Horava_eta/4._dl/(2._dl + 2._dl*self%Horava_xi - self%Horava_eta)
    eft_cache%EFTGamma6P = 0._dl

    end subroutine EFTCAMBHoravaSecondOrderEFTFunctions

    ! ---------------------------------------------------------------------------------------------
    !> Function that computes dtauda = 1/sqrt(a^2H^2).For pure EFT std this has to be overridden
    !! for performance reasons.
    function EFTCAMBHoravaComputeDtauda(self, a, eft_par_cache, eft_cache)

    implicit none

    class(EFTCAMB_Horava) :: self          !< the base class
    real(dl), intent(in)  :: a             !< the input scale factor
    !< the EFTCAMB parameter cache that contains all the physical parameters.
    type(TEFTCAMB_parameter_cache), intent(inout) :: eft_par_cache
    !< the EFTCAMB timestep cache that contains all the physical values.
    type(TEFTCAMB_timestep_cache), intent(inout) :: eft_cache

    real(dl) :: EFTCAMBHoravaComputeDtauda                        !< the output dtauda

    real(dl) :: temp

    temp = eft_cache%grhoa2 + a**4* &
        (eft_par_cache%grhov + 3._dl*eft_par_cache%h0_Mpc**2* &
        ((self%Horava_eta + 3._dl*self%Horava_lambda - 2._dl*self%Horava_xi)/ &
        (2._dl*self%Horava_xi + 2._dl - self%Horava_eta)))
    EFTCAMBHoravaComputeDtauda = &
        sqrt(3/temp)/sqrt((2._dl + 2._dl*self%Horava_xi - self%Horava_eta)/(3._dl*self%Horava_lambda + 2._dl))

    end function EFTCAMBHoravaComputeDtauda

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes adotoa = H and its two derivatives wrt conformal time.
    subroutine EFTCAMBHoravaComputeAdotoa(self, a, eft_par_cache, eft_cache)

    implicit none

    class(EFTCAMB_Horava) :: self          !< the base class
    real(dl), intent(in)  :: a             !< the input scale factor
    !< the EFTCAMB parameter cache that contains all the physical parameters.
    type(TEFTCAMB_parameter_cache), intent(inout) :: eft_par_cache
    !< the EFTCAMB timestep cache that contains all the physical values.
    type(TEFTCAMB_timestep_cache), intent(inout) :: eft_cache

    eft_cache%grhov_t = (eft_par_cache%grhov + &
        3._dl*eft_par_cache%h0_Mpc**2* &
        ((self%Horava_eta + 3._dl*self%Horava_lambda - 2._dl*self%Horava_xi)/ &
        (2._dl*self%Horava_xi + 2._dl - self%Horava_eta)))*a**2
    eft_cache%adotoa = sqrt((eft_cache%grhom_t + eft_cache%grhov_t)/3._dl)
    eft_cache%adotoa = &
        eft_cache%adotoa*sqrt((2._dl + 2._dl*self%Horava_xi - self%Horava_eta)/(3._dl*self%Horava_lambda + 2._dl))

    end subroutine EFTCAMBHoravaComputeAdotoa

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes the two derivatives wrt conformal time of H.
    subroutine EFTCAMBHoravaComputeHubbleDer(self, a, eft_par_cache, eft_cache)

    implicit none

    class(EFTCAMB_Horava) :: self          !< the base class
    real(dl), intent(in)  :: a             !< the input scale factor
    !< the EFTCAMB parameter cache that contains all the physical parameters.
    type(TEFTCAMB_parameter_cache), intent(inout) :: eft_par_cache
    !< the EFTCAMB timestep cache that contains all the physical values.
    type(TEFTCAMB_timestep_cache), intent(inout) :: eft_cache

    real(dl) :: F, P, Pdot, Pdotdot

    ! H^2 = F (grhom_t + grhov_t)/3 with F = G_cosmo/G_N (arXiv:1508.01787, eq. Hubble), and the matter fluid obeying
    ! the standard continuity equation. P is the total pressure (a^2-weighted), Pdot and Pdotdot its conformal time
    ! derivatives (gpresdotm_t and gpresdotdotm_t are a^2 times the time derivatives of the matter pressure).
    F = (2._dl + 2._dl*self%Horava_xi - self%Horava_eta)/(3._dl*self%Horava_lambda + 2._dl)
    eft_cache%gpiv_t = -eft_cache%grhov_t
    P = eft_cache%gpresm_t + eft_cache%gpiv_t
    Pdot = eft_cache%gpresdotm_t + 2._dl*eft_cache%adotoa*P
    eft_cache%Hdot = -0.5_dl*(eft_cache%adotoa**2 + F*P)
    eft_cache%Hdotdot = -eft_cache%adotoa*eft_cache%Hdot - 0.5_dl*F*Pdot
    Pdotdot = eft_cache%gpresdotdotm_t + 4._dl*eft_cache%adotoa*eft_cache%gpresdotm_t + &
        (2._dl*eft_cache%Hdot + 4._dl*eft_cache%adotoa**2)*P
    eft_cache%Hdotdotdot = -eft_cache%Hdot**2 - eft_cache%adotoa*eft_cache%Hdotdot - 0.5_dl*F*Pdotdot

    end subroutine EFTCAMBHoravaComputeHubbleDer

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes the pi field equations factors. For details refer to the numerical notes. Notice that in Horava we override it for reasons of numerical stability.
    subroutine EFTCAMBHoravaComputePiFactors(self, a, eft_par_cache, eft_cache)

    implicit none

    class(EFTCAMB_Horava) :: self          !< the base class
    real(dl), intent(in)  :: a             !< the input scale factor.
    !< the EFTCAMB parameter cache that contains all the physical parameters.
    type(TEFTCAMB_parameter_cache), intent(inout) :: eft_par_cache
    !< the EFTCAMB timestep cache that contains all the physical values.
    type(TEFTCAMB_timestep_cache), intent(inout) :: eft_cache

    eft_cache%EFTpiA1 = 0._dl
    eft_cache%EFTpiA2 = 0.5_dl*self%Horava_eta
    eft_cache%EFTpiB1 = 0._dl
    eft_cache%EFTpiB2 = eft_cache%adotoa*self%Horava_eta
    eft_cache%EFTpiC = 0._dl
    eft_cache%EFTpiD1 = 0.5_dl*(3._dl*self%Horava_lambda - 2._dl*self%Horava_xi)* &
        (eft_cache%adotoa**2 - eft_cache%Hdot) + 0.5_dl*self%Horava_eta*(eft_cache%adotoa**2 + eft_cache%Hdot)
    eft_cache%EFTpiD2 = 0.5_dl*self%Horava_lambda*(1._dl + self%Horava_xi)
    eft_cache%EFTpiE = 0.5_dl*eft_cache%k**2*self%Horava_lambda*(1._dl + self%Horava_xi)*eft_cache%k*eft_cache%z + &
        0.25_dl*eft_cache%k*self%Horava_xi*(2._dl*self%Horava_xi + 2._dl - self%Horava_eta)*eft_cache%dgq

    end subroutine EFTCAMBHoravaComputePiFactors

    ! ---------------------------------------------------------------------------------------------
    !> Function that computes model specific stability requirements.
    function EFTCAMBHoravaAdditionalModelStability(self, a, eft_par_cache, eft_cache)

    implicit none

    class(EFTCAMB_Horava) :: self          !< the base class
    real(dl), intent(in)  :: a             !< the input scale factor.
    !< the EFTCAMB parameter cache that contains all the physical parameters.
    type(TEFTCAMB_parameter_cache), intent(inout) :: eft_par_cache
    !< the EFTCAMB timestep cache that contains all the physical values.
    type(TEFTCAMB_timestep_cache), intent(inout) :: eft_cache

    !< the return value of the stability computation. True if the model specific stability criteria are met, false otherwise.
    logical :: EFTCAMBHoravaAdditionalModelStability

    EFTCAMBHoravaAdditionalModelStability = .true.

    if (self%Horava_lambda > -2._dl/3._dl .and. self%Horava_lambda < 0._dl) then
        EFTCAMBHoravaAdditionalModelStability = .false.
    end if
    if (self%Horava_eta < 0._dl .or. self%Horava_eta > 2._dl*self%Horava_xi + 2._dl) then
        EFTCAMBHoravaAdditionalModelStability = .false.
    end if

    end function EFTCAMBHoravaAdditionalModelStability

    ! ---------------------------------------------------------------------------------------------

    end module EFTCAMB_LE_Horava
