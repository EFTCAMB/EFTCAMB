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

!> @file 04p016_spline_parametrizations_1D.f90
!! This file contains the definition of a splined 1D function.
!! Inherits from parametrized_function_1D.

! ----------------------------------------------------------------------------------------
!> This module contains the definition of a splined 1D function.
!! Inherits from parametrized_function_1D.

!> @author Marco Raveri

    module EFTCAMB_spline_parametrizations_1D

    use precision
    use IniObjects
    use MpiUtils
    use EFT_splines
    use EFT_def
    use EFTCAMB_mixed_algorithms
    use linear_interpolation_1D
    use EFTCAMB_cache
    use EFTCAMB_abstract_parametrizations_1D

    implicit none

    private

    public spline_parametrization_1D

    ! ---------------------------------------------------------------------------------------------
    !> Type containing the Pade expansion parametrization. Inherits from parametrized_function_1D.
    type, extends (parametrized_function_1D) :: spline_parametrization_1D

        integer                              :: num_pixels     !< number of reconstruction points
        !< the value to report when function called outside the interpolation range
        real(dl)                             :: null_value
        !< array with the value of the coordinates where the reconstruction is defined
        real(dl), dimension(:), allocatable  :: x_coord
        real(dl), dimension(:), allocatable  :: pixels         !< array with the function value in each pixel

        type(linear_interpolate_function_1D) :: interpolation  !< auxiliary interpolation function

        ! internal variables:
        !< internal coefficients for the spline interpolation. Derivatives are computed from the same coefficients.
        real(dl), dimension(:, :), allocatable :: cy

    contains

    ! utility functions:
    !< subroutine that sets the number of parameters of the reconstruction function.
    procedure :: set_param_number => Spline1DSetParamNumber
    !< subroutine that initializes the function parameters based on the values found in an input array.
    procedure :: init_parameters => Spline1DInitParams
    !< subroutine that reads a Ini file looking for the parameters of the function.
    procedure :: init_from_file => Spline1DInitFromFile
    !< subroutine that reads a Ini file looking for initialization parameters for the function..
    procedure :: init_func_from_file => Spline1DInitFunctionFromFile
    !< subroutine that returns the value of the function i-th parameter.
    procedure :: parameter_value => Spline1DParameterValues
    !< subroutine that prints to screen the informations about the function.
    procedure :: feedback => Spline1DFeedback
    !< subroutine that returns the i-th parameter name of the function.
    procedure :: parameter_names => Spline1DParameterNames
    !< subroutine that returns the i-th parameter name of the function in latex format.
    procedure :: parameter_names_latex => Spline1DParameterNamesLatex

    ! evaluation procedures:
    procedure :: value => Spline1DValue               !< function that returns the value of the reconstruction.
    !< function that returns the first derivative of the reconstruction.
    procedure :: first_derivative => Spline1DFirstDerivative
    !< function that returns the second derivative of the reconstruction.
    procedure :: second_derivative => Spline1DSecondDerivative
    !< function that returns the third derivative of the reconstruction.
    procedure :: third_derivative => Spline1DThirdDerivative
    !< function that returns the fourth derivative of the reconstruction.
    procedure :: fourth_derivative => Spline1DFourthDerivative
    !< function that returns the strange integral that we need for w_DE.
    procedure :: integral => Spline1DIntegral

    end type spline_parametrization_1D

    contains

    ! ---------------------------------------------------------------------------------------------
    ! Implementation of the reconstructioin parametrization.
    ! ---------------------------------------------------------------------------------------------

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that sets the number of parameters of the Pade expansion parametrized function.
    subroutine Spline1DSetParamNumber(self)

    implicit none

    class(spline_parametrization_1D) :: self       !< the base class

    ! initialize the number of parameters:
    self%parameter_number = 2*self%num_pixels

    end subroutine Spline1DSetParamNumber

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that initializes the function parameters based on the values found in an input array.
    subroutine Spline1DInitParams(self, array)

    implicit none

    class(spline_parametrization_1D)                       :: self   !< the base class.
    real(dl), dimension(self%parameter_number), intent(in) :: array  !< input array with the values of the parameters.

    integer :: ind

    self%x_coord = array(1:self%num_pixels)
    self%pixels = array(self%num_pixels + 1:)
    ! check the coordinates:
    do ind = 1, self%num_pixels - 1
        if (self%x_coord(ind) >= self%x_coord(ind + 1)) call MpiStop(trim(self%name) // &
            ' error: spline x coordinate must be increasing. Pass correct vaalues for ' // trim(self%name) // 'xN')
    end do
    ! initialize the interpolation:
    call self%interpolation%initialize(self%x_coord, self%null_value)
    self%interpolation%y = self%pixels
    ! create the spline coefficients (natural cubic spline). The derivatives are the derivatives of this spline:
    if (allocated(self%cy)) deallocate(self%cy)
    allocate(self%cy   (5, self%num_pixels - 1))
    call spline3pars(self%x_coord, self%interpolation%y, [2, 2], [0._dl, 0._dl], self%cy)

    end subroutine Spline1DInitParams

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that reads a Ini file looking for the parameters of the function.
    subroutine Spline1DInitFromFile(self, Ini, eft_error)

    implicit none

    class(spline_parametrization_1D) :: self      !< the base class
    type(TIniFile)                   :: Ini       !< Input ini file
    integer                          :: eft_error !< error code: 0 all fine, 1 initialization failed

    character(len=EFT_names_max_length)        :: param_name
    real(dl), dimension(self%parameter_number) :: parameters

    integer :: i

    ! ensure that the number of parameters is properly associated:
    call self%set_param_number()

    ! read the parameters and store them in a vector:
    do i = 1, self%parameter_number
        call self%parameter_names(i, param_name)
        if (i <= self%num_pixels) then
            parameters(i) = Ini%Read_Double(trim(param_name), 1._dl*real(i - 1)/real(self%num_pixels - 1))
        else
            parameters(i) = Ini%Read_Double(trim(param_name), 0._dl)
        end if
    end do

    ! initialize the function parameters from the vector:
    call self%init_parameters(parameters)

    end subroutine Spline1DInitFromFile

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that reads a Ini file looking for initialization parameters for the function..
    subroutine Spline1DInitFunctionFromFile(self, Ini, eft_error)

    implicit none

    class(spline_parametrization_1D) :: self      !< the base class
    type(TIniFile)                   :: Ini       !< Input ini file
    integer                          :: eft_error !< error code: 0 all fine, 1 initialization failed

    ! read number of points:
    self%num_pixels = Ini%Read_Int(trim(self%name) // '_Spline_Pixels', 6)
    if (self%num_pixels < 3) then
        write(*, *) 'Spline requires ', trim(self%name) // '_Spline_Pixels', ' to be at least 6'
        eft_error = 1
        return
    end if
    self%null_value = Ini%Read_Double(trim(self%name) // '_null_value', 0._dl)
    ! allocate:
    if (allocated(self%x_coord)) deallocate(self%x_coord)
    allocate(self%x_coord(self%num_pixels))
    if (allocated(self%pixels)) deallocate(self%pixels)
    allocate(self%pixels(self%num_pixels))
    ! ensure parameter number remains consistent:
    call self%set_param_number()

    end subroutine Spline1DInitFunctionFromFile

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the value of the function i-th parameter.
    subroutine Spline1DParameterValues(self, i, value)

    implicit none

    class(spline_parametrization_1D) :: self        !< the base class
    integer, intent(in)   :: i           !< The index of the parameter
    real(dl), intent(out) :: value       !< the output value of the i-th parameter

    if (i <= self%num_pixels) then
        value = self%x_coord(i)
    else if (i > self%num_pixels .and. i <= self%parameter_number) then
        value = self%pixels(i - self%num_pixels)
    else
        write(*, *) 'Illegal index for parameter_names.'
        write(*, *) 'Maximum value is:', self%parameter_number
        call MpiStop('EFTCAMB error')
    end if

    end subroutine Spline1DParameterValues

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that prints to screen the informations about the function.
    subroutine Spline1DFeedback(self, print_params)

    implicit none

    class(spline_parametrization_1D) :: self         !< the base class
    logical, optional                :: print_params !< optional flag that decised whether to print numerical values
    !! of the parameters.

    integer                             :: i
    real(dl)                            :: param_value
    character(len=EFT_names_max_length) :: param_name
    logical                             :: print_params_temp

    if (present(print_params)) then
        print_params_temp = print_params
    else
        print_params_temp = .true.
    end if

    write(*, '(a,a)')      'spline parametrization: ', self%name
    write(*, '(a,a,a,I3)') ' number of pixels (', trim(self%name) // '_Spline_Pixels', '): ', self%num_pixels
    write(*, '(a,F12.6)')  ' null value: ', self%null_value
    if (print_params_temp) then
        do i = 1, self%parameter_number
            call self%parameter_names(i, param_name)
            call self%parameter_value(i, param_value)
            write(*, '(a23,a,F12.6)') param_name, '=', param_value
        end do
    end if

    end subroutine Spline1DFeedback

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the i-th parameter name.
    subroutine Spline1DParameterNames(self, i, name)

    implicit none

    class(spline_parametrization_1D) :: self   !< the base class
    integer, intent(in)       :: i      !< the index of the parameter
    character(*), intent(out) :: name   !< the output name of the i-th parameter

    ! check the input index:
    if (i > self%parameter_number) then
        write(*, *) 'In parametrized_function_1D:', self%name
        write(*, *) 'Illegal index for parameter_names.'
        write(*, *) 'Maximum value is:', self%parameter_number
        call MpiStop('EFTCAMB error')
    end if
    ! return the parameter name:
    if (allocated(self%param_names)) then
        name = self%param_names(i)%string
    else
        if (i <= self%num_pixels) then
            name = trim(self%name) // 'x' // integer_to_string(i)
        else if (i > self%num_pixels .and. i <= self%parameter_number) then
            name = trim(self%name) // 'v' // integer_to_string(i - self%num_pixels)
        end if
    end if

    end subroutine Spline1DParameterNames

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the latex version of the i-th parameter name.
    subroutine Spline1DParameterNamesLatex(self, i, latexname)

    implicit none

    class(spline_parametrization_1D) :: self       !< the base class
    integer, intent(in)              :: i          !< the index of the parameter
    character(*), intent(out)        :: latexname  !< the output latex name of the i-th parameter

    ! check the input index:
    if (i > self%parameter_number) then
        write(*, *) 'In parametrized_function_1D:', self%name
        write(*, *) 'Illegal index for parameter_names.'
        write(*, *) 'Maximum value is:', self%parameter_number
        call MpiStop('EFTCAMB error')
    end if
    ! return the parameter name:
    if (allocated(self%param_names_latex)) then
        latexname = self%param_names_latex(i)%string
    else
        if (i <= self%num_pixels) then
            latexname = trim(self%name) // 'x_' // integer_to_string(i)
        else if (i > self%num_pixels .and. i <= self%parameter_number) then
            latexname = trim(self%name) // 'v_' // integer_to_string(i - self%num_pixels)
        end if
    end if

    end subroutine Spline1DParameterNamesLatex

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the value of the function in the scale factor.
    function Spline1DValue(self, x, eft_cache)

    implicit none

    class(spline_parametrization_1D)                    :: self      !< the base class
    real(dl), intent(in)                                :: x         !< the input scale factor
    type(TEFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
    real(dl) :: Spline1DValue                                        !< the output value

    real(dl) :: x_temp
    integer  :: ind, ip

    ! Uncomment this to have linear interpolation:
    ! Spline1DValue = self%interpolation%value(x)
    ! Default is cubic spline:
    Spline1DValue = self%null_value
    if (x < self%interpolation%x_initial .or. x > self%interpolation%x_final) return
    ! if inside return the cubic spline value:
    ip = iixmin(x, self%x_coord, 0)
    Spline1DValue = poly3 (x, self%cy(:, ip))

    end function Spline1DValue

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the value of the first derivative, wrt scale factor, of the function.
    function Spline1DFirstDerivative(self, x, eft_cache)

    implicit none

    class(spline_parametrization_1D)                    :: self      !< the base class
    real(dl), intent(in)                                :: x         !< the input scale factor
    type(TEFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
    real(dl) :: Spline1DFirstDerivative                              !< the output value

    real(dl) :: x_temp
    integer  :: ind, ip

    ! Spline1DFirstDerivative = self%interpolation%first_derivative(x)
    ! Default is cubic spline:
    ! outside the nodes the function is the constant null_value:
    Spline1DFirstDerivative = 0._dl
    if (x < self%interpolation%x_initial .or. x > self%interpolation%x_final) return
    ! if inside return the derivative of the cubic spline:
    ip = iixmin(x, self%x_coord, 0)
    Spline1DFirstDerivative = dpoly3 (x, self%cy(:, ip))

    end function Spline1DFirstDerivative

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the second derivative of the function.
    function Spline1DSecondDerivative(self, x, eft_cache)

    implicit none

    class(spline_parametrization_1D)                    :: self      !< the base class
    real(dl), intent(in)                                :: x         !< the input scale factor
    type(TEFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
    real(dl) :: Spline1DSecondDerivative                             !< the output value

    real(dl) :: x_temp
    integer  :: ind, ip

    ! Spline1DSecondDerivative = self%interpolation%second_derivative(x)
    ! Default is cubic spline:
    ! outside the nodes the function is the constant null_value:
    Spline1DSecondDerivative = 0._dl
    if (x < self%interpolation%x_initial .or. x > self%interpolation%x_final) return
    ! if inside return the derivative of the cubic spline:
    ip = iixmin(x, self%x_coord, 0)
    Spline1DSecondDerivative = d2poly3 (x, self%cy(:, ip))

    end function Spline1DSecondDerivative

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the third derivative of the function.
    function Spline1DThirdDerivative(self, x, eft_cache)

    implicit none

    class(spline_parametrization_1D)                    :: self      !< the base class
    real(dl), intent(in)                                :: x         !< the input scale factor
    type(TEFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
    real(dl) :: Spline1DThirdDerivative                              !< the output value

    real(dl) :: x_temp
    integer  :: ind, ip

    ! outside the nodes the function is the constant null_value:
    Spline1DThirdDerivative = 0._dl
    if (x < self%interpolation%x_initial .or. x > self%interpolation%x_final) return
    ! if inside return the derivative of the cubic spline:
    ip = iixmin(x, self%x_coord, 0)
    Spline1DThirdDerivative = 6._dl*self%cy(5, ip)

    end function Spline1DThirdDerivative

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the fourth derivative of the function.
    function Spline1DFourthDerivative(self, x, eft_cache)

    implicit none

    class(spline_parametrization_1D)                    :: self      !< the base class
    real(dl), intent(in)                                :: x         !< the input scale factor
    type(TEFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
    real(dl) :: Spline1DFourthDerivative                             !< the output value

    ! the fourth derivative of a cubic spline vanishes inside each interval, and outside the nodes:
    Spline1DFourthDerivative = 0._dl

    end function Spline1DFourthDerivative

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the integral of the function, as defined in the notes.
    function Spline1DIntegral(self, x, eft_cache)

    implicit none

    class(spline_parametrization_1D)                    :: self      !< the base class
    real(dl), intent(in)                                :: x         !< the input scale factor
    type(TEFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
    real(dl) :: Spline1DIntegral

    ! no closed form: numerical integral of the spline, split at the nodes.
    Spline1DIntegral = self%numerical_integral(x, self%x_coord, eft_cache)

    end function Spline1DIntegral

    ! ---------------------------------------------------------------------------------------------

    end module EFTCAMB_spline_parametrizations_1D

! ----------------------------------------------------------------------------------------
