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

!> @file 04p013_double_exponential_parametrizations_1D.f90
!! This file contains the definition of the second type of exponential parametrization,
!! inheriting from parametrized_function_1D.


!----------------------------------------------------------------------------------------
!> This module contains the definition of the second type of exponential parametrization,
!! inheriting from parametrized_function_1D.

!> @author Marco Raveri

module EFTCAMB_double_exponential_parametrizations_1D

    use precision
    use EFT_def
    use MpiUtils
    use EFTCAMB_cache
    use EFTCAMB_abstract_parametrizations_1D

    implicit none

    private

    public double_exponential_parametrization_1D

    ! ---------------------------------------------------------------------------------------------
    !> Type containing the exponential function parametrization. Inherits from parametrized_function_1D.
    type, extends ( parametrized_function_1D ) :: double_exponential_parametrization_1D

        real(dl) :: coefficient1, coefficient2
        real(dl) :: exponent1, exponent2

    contains

        ! utility functions:
        procedure :: set_param_number      => DoubleExponentialParametrized1DSetParamNumber      !< subroutine that sets the number of parameters of the exponential parametrized function.
        procedure :: init_parameters       => DoubleExponentialParametrized1DInitParams          !< subroutine that initializes the function parameters based on the values found in an input array.
        procedure :: parameter_value       => DoubleExponentialParametrized1DParameterValues     !< subroutine that returns the value of the function i-th parameter.
        procedure :: feedback              => DoubleExponentialParametrized1DFeedback            !< subroutine that prints to screen the informations about the function.

        ! evaluation procedures:
        procedure :: value                 => DoubleExponentialParametrized1DValue               !< function that returns the value of the exponential function.
        procedure :: first_derivative      => DoubleExponentialParametrized1DFirstDerivative     !< function that returns the first derivative of the exponential function.
        procedure :: second_derivative     => DoubleExponentialParametrized1DSecondDerivative    !< function that returns the second derivative of the exponential function.
        procedure :: third_derivative      => DoubleExponentialParametrized1DThirdDerivative     !< function that returns the third derivative of the exponential function.
        procedure :: fourth_derivative     => DoubleExponentialParametrized1DFourthDerivative    !< function that returns the fourth derivative of the exponential function.
        procedure :: integral              => DoubleExponentialParametrized1DIntegral            !< function that returns the strange integral that we need for w_DE.

    end type double_exponential_parametrization_1D

contains

    ! ---------------------------------------------------------------------------------------------
    ! Implementation of the exponential function.
    ! ---------------------------------------------------------------------------------------------

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that sets the number of parameters of the exponential parametrized function.
    subroutine DoubleExponentialParametrized1DSetParamNumber( self )

        implicit none

        class(double_exponential_parametrization_1D) :: self       !< the base class

        ! initialize the number of parameters:
        self%parameter_number = 4

    end subroutine DoubleExponentialParametrized1DSetParamNumber

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that initializes the function parameters based on the values found in an input array.
    subroutine DoubleExponentialParametrized1DInitParams( self, array )

        implicit none

        class(double_exponential_parametrization_1D)            :: self   !< the base class.
        real(dl), dimension(self%parameter_number), intent(in)  :: array  !< input array with the values of the parameters.

        self%coefficient1 = array(1)
        self%exponent1    = array(2)
        self%coefficient2 = array(3)
        self%exponent2    = array(4)

    end subroutine DoubleExponentialParametrized1DInitParams

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the value of the function i-th parameter.
    subroutine DoubleExponentialParametrized1DParameterValues( self, i, value )

        implicit none

        class(double_exponential_parametrization_1D):: self        !< the base class
        integer     , intent(in)             :: i           !< The index of the parameter
        real(dl)    , intent(out)            :: value       !< the output value of the i-th parameter

        select case (i)
            case(1)
                value = self%coefficient1
            case(2)
                value = self%exponent1
            case(3)
                value = self%coefficient2
            case(4)
                value = self%exponent2
            case default
                write(*,*) 'Illegal index for parameter_names.'
                write(*,*) 'Maximum value is:', self%parameter_number
                call MpiStop('EFTCAMB error')
        end select

    end subroutine DoubleExponentialParametrized1DParameterValues

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that prints to screen the informations about the function.
    subroutine DoubleExponentialParametrized1DFeedback( self, print_params )

        implicit none

        class(double_exponential_parametrization_1D) :: self         !< the base class
        logical, optional                     :: print_params !< optional flag that decised whether to print numerical values
                                                              !! of the parameters.

        integer                               :: i
        real(dl)                              :: param_value
        character(len=EFT_names_max_length)   :: param_name
        logical                               :: print_params_temp

        if ( present(print_params) ) then
            print_params_temp = print_params
        else
            print_params_temp = .True.
        end if

        write(*,*)     'Exponential function: ', self%name
        if ( print_params_temp ) then
            do i=1, self%parameter_number
                call self%parameter_names( i, param_name  )
                call self%parameter_value( i, param_value )
                write(*,'(a23,a,F12.6)') param_name, '=', param_value
            end do
        end if

    end subroutine DoubleExponentialParametrized1DFeedback

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the value of the function in the scale factor.
    function DoubleExponentialParametrized1DValue( self, x, eft_cache )

        implicit none

        class(double_exponential_parametrization_1D)              :: self      !< the base class
        real(dl), intent(in)                               :: x         !< the input scale factor
        type(TEFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: DoubleExponentialParametrized1DValue                      !< the output value

        !DoubleExponentialParametrized1DValue = Exp(self%coefficient*x**self%exponent) -1._dl
        DoubleExponentialParametrized1DValue = self%coefficient1*Exp(-x*self%exponent1) +self%coefficient2*Exp(-x*self%exponent2)

    end function DoubleExponentialParametrized1DValue

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the value of the first derivative, wrt scale factor, of the function.
    function DoubleExponentialParametrized1DFirstDerivative( self, x, eft_cache )

        implicit none

        class(double_exponential_parametrization_1D)              :: self      !< the base class
        real(dl), intent(in)                               :: x         !< the input scale factor
        type(TEFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: DoubleExponentialParametrized1DFirstDerivative            !< the output value

        DoubleExponentialParametrized1DFirstDerivative = -self%exponent1*self%coefficient1*Exp(-x*self%exponent1) &
            & -self%exponent2*self%coefficient2*Exp(-x*self%exponent2)

    end function DoubleExponentialParametrized1DFirstDerivative

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the second derivative of the function.
    function DoubleExponentialParametrized1DSecondDerivative( self, x, eft_cache )

        implicit none

        class(double_exponential_parametrization_1D)              :: self      !< the base class
        real(dl), intent(in)                               :: x         !< the input scale factor
        type(TEFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: DoubleExponentialParametrized1DSecondDerivative           !< the output value

        DoubleExponentialParametrized1DSecondDerivative = self%exponent1**2*self%coefficient1*Exp(-x*self%exponent1) &
            & +self%exponent2**2*self%coefficient2*Exp(-x*self%exponent2)

    end function DoubleExponentialParametrized1DSecondDerivative

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the third derivative of the function.
    function DoubleExponentialParametrized1DThirdDerivative( self, x, eft_cache )

        implicit none

        class(double_exponential_parametrization_1D)              :: self      !< the base class
        real(dl), intent(in)                               :: x         !< the input scale factor
        type(TEFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: DoubleExponentialParametrized1DThirdDerivative            !< the output value

        DoubleExponentialParametrized1DThirdDerivative = -self%exponent1**3*self%coefficient1*Exp(-x*self%exponent1) &
            & -self%exponent2**3*self%coefficient2*Exp(-x*self%exponent2)

    end function DoubleExponentialParametrized1DThirdDerivative

    ! --------------------------------------------------------------------------------------------
    !> Function that returns the fourth derivative of the function.
    function DoubleExponentialParametrized1DFourthDerivative( self, x, eft_cache )

        implicit none

        class(double_exponential_parametrization_1D)              :: self      !< the base class
        real(dl), intent(in)                               :: x         !< the input scale factor
        type(TEFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: DoubleExponentialParametrized1DFourthDerivative           !< the output value

        DoubleExponentialParametrized1DFourthDerivative = self%exponent1**4*self%coefficient1*Exp(-x*self%exponent1) &
            & +self%exponent2**4*self%coefficient2*Exp(-x*self%exponent2)

    end function DoubleExponentialParametrized1DFourthDerivative

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the integral of the function, as defined in the notes.
    function DoubleExponentialParametrized1DIntegral( self, x, eft_cache )

        implicit none

        class(double_exponential_parametrization_1D)              :: self      !< the base class
        real(dl), intent(in)                               :: x         !< the input scale factor
        type(TEFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: DoubleExponentialParametrized1DIntegral                   !< the output value

        DoubleExponentialParametrized1DIntegral = 0._dl

        !< No analytic solution >!
        write(*,*) 'DoubleExponentialParametrized1DIntegral is not implemented.'
        write(*,*) 'Calculations cannot proceed.'
        call MpiStop('EFTCAMB error')

    end function DoubleExponentialParametrized1DIntegral

    ! ---------------------------------------------------------------------------------------------

end module EFTCAMB_double_exponential_parametrizations_1D

!----------------------------------------------------------------------------------------
