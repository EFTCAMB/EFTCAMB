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

!> @file 04p009_taylor_series_parametrizations_1D.f90
!! This file contains the definition of the Taylor series parametrization, in the scale
!! factor, around an arbitrary expansion point. For best performances Horner’s algorithm is used.
!! Inherits from parametrized_function_1D.

! ----------------------------------------------------------------------------------------
!> This module contains the definition of the Taylor series parametrization, in the scale
!! factor, around an arbitrary expansion point. For best performances Horner’s algorithm is used.
!! Inherits from parametrized_function_1D.

!> @author Marco Raveri, Simone Peirone

    module EFTCAMB_taylor_series_1D

    use precision
    use IniObjects
    use EFT_def
    use MpiUtils
    use EFTCAMB_mixed_algorithms
    use EFTCAMB_cache
    use EFTCAMB_abstract_parametrizations_1D

    implicit none

    private

    public taylorseries_parametrization_1D

    ! ---------------------------------------------------------------------------------------------
    !> Type containing the Taylor expansion parametrization. Inherits from parametrized_function_1D.
    type, extends (parametrized_function_1D) :: taylorseries_parametrization_1D

        integer                              :: order   !< order of the Taylor expansion
        real(dl)                             :: a0      !< expansion point
        real(dl), dimension(:), allocatable  :: coeff   !< array with the coefficients of the Taylor expansion

    contains

    ! utility functions:
    !< subroutine that sets the number of parameters of the Taylor expansion parametrized function.
    procedure :: set_param_number => TaylorSeriesParametrized1DSetParamNumber
    !< subroutine that initializes the function parameters based on the values found in an input array.
    procedure :: init_parameters => TaylorSeriesParametrized1DInitParams
    !< subroutine that reads a Ini file looking for initialization parameters for the function..
    procedure :: init_func_from_file => TaylorSeriesParametrized1DInitFunctionFromFile
    !< subroutine that returns the value of the function i-th parameter.
    procedure :: parameter_value => TaylorSeriesParametrized1DParameterValues
    !< subroutine that prints to screen the informations about the function.
    procedure :: feedback => TaylorSeriesParametrized1DFeedback
    !< subroutine that returns the i-th parameter name of the function.
    procedure :: parameter_names => TaylorSeriesParametrized1DParameterNames
    !< subroutine that returns the i-th parameter name of the function in latex format.
    procedure :: parameter_names_latex => TaylorSeriesParametrized1DParameterNamesLatex

    ! evaluation procedures:
    !< function that returns the value of the Taylor expansion.
    procedure :: value => TaylorSeriesParametrized1DValue
    !< function that returns the first derivative of the Taylor expansion.
    procedure :: first_derivative => TaylorSeriesParametrized1DFirstDerivative
    !< function that returns the second derivative of the Taylor expansion.
    procedure :: second_derivative => TaylorSeriesParametrized1DSecondDerivative
    !< function that returns the third derivative of the Taylor expansion.
    procedure :: third_derivative => TaylorSeriesParametrized1DThirdDerivative
    !< function that returns the fourth derivative of the Taylor expansion.
    procedure :: fourth_derivative => TaylorSeriesParametrized1DFourthDerivative
    !< function that returns the strange integral that we need for w_DE.
    procedure :: integral => TaylorSeriesParametrized1DIntegral

    end type taylorseries_parametrization_1D

    contains

    ! ---------------------------------------------------------------------------------------------
    ! Implementation of the Taylor series parametrization.
    ! ---------------------------------------------------------------------------------------------

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that sets the number of parameters of the Taylor expansion parametrized function.
    subroutine TaylorSeriesParametrized1DSetParamNumber(self)

    implicit none

    class(taylorseries_parametrization_1D) :: self       !< the base class

    ! initialize the number of parameters:
    self%parameter_number = self%order + 2

    end subroutine TaylorSeriesParametrized1DSetParamNumber

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that initializes the function parameters based on the values found in an input array.
    subroutine TaylorSeriesParametrized1DInitParams(self, array)

    implicit none

    class(taylorseries_parametrization_1D)                 :: self   !< the base class.
    real(dl), dimension(self%parameter_number), intent(in) :: array  !< input array with the values of the parameters.

    integer :: ind

    self%a0 = array(1)
    self%coeff(:) = array(2:)

    end subroutine TaylorSeriesParametrized1DInitParams

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that reads a Ini file looking for initialization parameters for the function..
    subroutine TaylorSeriesParametrized1DInitFunctionFromFile(self, Ini, eft_error)

    implicit none

    class(taylorseries_parametrization_1D) :: self   !< the base class
    class(TIniFile)                        :: Ini    !< Input ini file
    integer                                :: eft_error !< error code: 0 all fine, 1 initialization failed

    ! read order of the Taylor expansion:
    self%order = Ini%Read_Int(trim(self%name) // '_Taylor_order', 0)
    ! allocate:
    if (allocated(self%coeff)) deallocate(self%coeff)
    allocate(self%coeff(self%order + 1))
    self%coeff = 0._dl
    ! ensure parameter number remains consistent:
    call self%set_param_number()

    end subroutine TaylorSeriesParametrized1DInitFunctionFromFile

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the value of the function i-th parameter.
    subroutine TaylorSeriesParametrized1DParameterValues(self, i, value)

    implicit none

    class(taylorseries_parametrization_1D) :: self        !< the base class
    integer, intent(in)   :: i           !< The index of the parameter
    real(dl), intent(out) :: value       !< the output value of the i-th parameter

    if (i == 1) then
        value = self%a0
    else if (i > 1 .and. i <= self%parameter_number) then
        value = self%coeff(i - 1)
    else
        write(*, *) 'Illegal index for parameter_names.'
        write(*, *) 'Maximum value is:', self%parameter_number
        call MpiStop('EFTCAMB error')
    end if

    end subroutine TaylorSeriesParametrized1DParameterValues

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that prints to screen the informations about the function.
    subroutine TaylorSeriesParametrized1DFeedback(self, print_params)

    implicit none

    class(taylorseries_parametrization_1D) :: self         !< the base class
    !< optional flag that decised whether to print numerical values
    logical, optional                      :: print_params
    !! of the parameters.

    integer                                :: i
    real(dl)                               :: param_value
    character(len=EFT_names_max_length)    :: param_name
    logical                                :: print_params_temp

    if (present(print_params)) then
        print_params_temp = print_params
    else
        print_params_temp = .true.
    end if

    write(*, '(a,a)')      'Taylor series parametrization: ', self%name
    write(*, '(a,a,a,I3)') ' of order (', trim(self%name) // '_Taylor_order', '): ', self%order
    if (print_params_temp) then
        do i = 1, self%parameter_number
            call self%parameter_names(i, param_name)
            call self%parameter_value(i, param_value)
            write(*, '(a23,a,F12.6)') param_name, '=', param_value
        end do
    end if

    end subroutine TaylorSeriesParametrized1DFeedback

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the i-th parameter name.
    subroutine TaylorSeriesParametrized1DParameterNames(self, i, name)

    implicit none

    class(taylorseries_parametrization_1D) :: self   !< the base class
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
        if (i == 1) then
            name = trim(self%name) // 'a0'
        else
            name = trim(self%name) // integer_to_string(i - 2)
        end if
    end if

    end subroutine TaylorSeriesParametrized1DParameterNames

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the latex version of the i-th parameter name.
    subroutine TaylorSeriesParametrized1DParameterNamesLatex(self, i, latexname)

    implicit none

    class(taylorseries_parametrization_1D) :: self       !< the base class
    integer, intent(in)       :: i          !< the index of the parameter
    character(*), intent(out) :: latexname  !< the output latex name of the i-th parameter

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
        if (i == 1) then
            latexname = trim(self%name) // '_' // 'a_0'
        else
            latexname = trim(self%name) // '_' // integer_to_string(i - 2)
        end if
    end if

    end subroutine TaylorSeriesParametrized1DParameterNamesLatex

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the value of the function in the scale factor.
    function TaylorSeriesParametrized1DValue(self, x, eft_cache)

    implicit none

    class(taylorseries_parametrization_1D) :: self      !< the base class
    real(dl), intent(in)                   :: x         !< the input scale factor
    type(TEFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
    real(dl) :: TaylorSeriesParametrized1DValue                     !< the output value

    real(dl) :: x_temp
    integer  :: ind

    TaylorSeriesParametrized1DValue = 0._dl

    if (self%order < 0) then
        return
    end if

    x_temp = x - self%a0
    do ind = self%order + 1, 2, -1
        TaylorSeriesParametrized1DValue = (TaylorSeriesParametrized1DValue + self%coeff(ind))*x_temp
    end do
    TaylorSeriesParametrized1DValue = TaylorSeriesParametrized1DValue + self%coeff(1)

    end function TaylorSeriesParametrized1DValue

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the value of the first derivative, wrt scale factor, of the function.
    function TaylorSeriesParametrized1DFirstDerivative(self, x, eft_cache)

    implicit none

    class(taylorseries_parametrization_1D) :: self      !< the base class
    real(dl), intent(in)                   :: x         !< the input scale factor
    type(TEFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
    real(dl) :: TaylorSeriesParametrized1DFirstDerivative           !< the output value

    real(dl) :: x_temp
    integer  :: ind

    TaylorSeriesParametrized1DFirstDerivative = 0._dl

    if (self%order < 1) then
        return
    end if

    x_temp = x - self%a0
    do ind = self%order + 1, 3, -1
        TaylorSeriesParametrized1DFirstDerivative = (TaylorSeriesParametrized1DFirstDerivative &
            + real(ind - 1)*self%coeff(ind))*x_temp
    end do
    TaylorSeriesParametrized1DFirstDerivative = TaylorSeriesParametrized1DFirstDerivative + self%coeff(2)

    end function TaylorSeriesParametrized1DFirstDerivative

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the second derivative of the function.
    function TaylorSeriesParametrized1DSecondDerivative(self, x, eft_cache)

    implicit none

    class(taylorseries_parametrization_1D) :: self      !< the base class
    real(dl), intent(in)                   :: x         !< the input scale factor
    type(TEFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
    real(dl) :: TaylorSeriesParametrized1DSecondDerivative                !< the output value

    real(dl) :: x_temp
    integer  :: ind

    TaylorSeriesParametrized1DSecondDerivative = 0._dl

    if (self%order < 2) then
        return
    end if

    x_temp = x - self%a0
    do ind = self%order + 1, 4, -1
        TaylorSeriesParametrized1DSecondDerivative = (TaylorSeriesParametrized1DSecondDerivative &
            + real(ind - 2)*real(ind - 1)*self%coeff(ind))*x_temp
    end do
    TaylorSeriesParametrized1DSecondDerivative = TaylorSeriesParametrized1DSecondDerivative + 2._dl*self%coeff(3)

    end function TaylorSeriesParametrized1DSecondDerivative

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the third derivative of the function.
    function TaylorSeriesParametrized1DThirdDerivative(self, x, eft_cache)

    implicit none

    class(taylorseries_parametrization_1D) :: self      !< the base class
    real(dl), intent(in)                   :: x         !< the input scale factor
    type(TEFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
    real(dl) :: TaylorSeriesParametrized1DThirdDerivative           !< the output value

    real(dl) :: x_temp
    integer  :: ind

    TaylorSeriesParametrized1DThirdDerivative = 0._dl

    if (self%order < 3) then
        return
    end if

    x_temp = x - self%a0
    do ind = self%order + 1, 5, -1
        TaylorSeriesParametrized1DThirdDerivative = (TaylorSeriesParametrized1DThirdDerivative &
            + real(ind - 3)*real(ind - 2)*real(ind - 1)*self%coeff(ind))*x_temp
    end do
    TaylorSeriesParametrized1DThirdDerivative = TaylorSeriesParametrized1DThirdDerivative + 6._dl*self%coeff(4)

    end function TaylorSeriesParametrized1DThirdDerivative

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the fourth derivative of the function.
    function TaylorSeriesParametrized1DFourthDerivative(self, x, eft_cache)

    implicit none

    class(taylorseries_parametrization_1D) :: self      !< the base class
    real(dl), intent(in)                   :: x         !< the input scale factor
    type(TEFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
    real(dl) :: TaylorSeriesParametrized1DFourthDerivative           !< the output value

    real(dl) :: x_temp
    integer  :: ind

    TaylorSeriesParametrized1DFourthDerivative = 0._dl

    if (self%order < 4) then
        return
    end if

    x_temp = x - self%a0
    do ind = self%order + 1, 6, -1
        TaylorSeriesParametrized1DFourthDerivative = (TaylorSeriesParametrized1DFourthDerivative &
            + real(ind - 4)*real(ind - 3)*real(ind - 2)*real(ind - 1)*self%coeff(ind))*x_temp
    end do
    TaylorSeriesParametrized1DFourthDerivative = TaylorSeriesParametrized1DFourthDerivative + 24._dl*self%coeff(5)

    end function TaylorSeriesParametrized1DFourthDerivative

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the integral of the function, as defined in the notes.
    function TaylorSeriesParametrized1DIntegral(self, x, eft_cache)

    implicit none

    class(taylorseries_parametrization_1D) :: self      !< the base class
    real(dl), intent(in)                   :: x         !< the input scale factor
    type(TEFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
    real(dl) :: TaylorSeriesParametrized1DIntegral                  !< the output value

    ! no closed form: numerical integral (see parametrized_function_1D%numerical_integral).
    TaylorSeriesParametrized1DIntegral = self%numerical_integral(x, eft_cache=eft_cache)

    end function TaylorSeriesParametrized1DIntegral

    ! ---------------------------------------------------------------------------------------------

    end module EFTCAMB_taylor_series_1D

! ----------------------------------------------------------------------------------------
