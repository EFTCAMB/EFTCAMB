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

!> @file 04p4_exponential_parametrizations_1D.f90
!! This file contains the definition of the exponential parametrization,
!! inheriting from parametrized_function_1D.


!----------------------------------------------------------------------------------------
!> This module contains the definition of the exponential parametrization,
!! inheriting from parametrized_function_1D.

!> @author Bin Hu, Marco Raveri, Simone Peirone

module EFTCAMB_exponential_parametrizations_1D

    use precision
    use AMLutils
    use EFT_def
    use EFTCAMB_cache
    use EFTCAMB_abstract_parametrizations_1D

    implicit none

    private

    public exponential_parametrization_1D

    ! ---------------------------------------------------------------------------------------------
    !> Type containing the exponential function parametrization. Inherits from parametrized_function_1D.
    type, extends ( parametrized_function_1D ) :: exponential_parametrization_1D

        real(dl) :: coefficient
        real(dl) :: exponent

    contains

        ! utility functions:
        procedure :: set_param_number      => ExponentialParametrized1DSetParamNumber      !< subroutine that sets the number of parameters of the exponential parametrized function.
        procedure :: init_parameters       => ExponentialParametrized1DInitParams          !< subroutine that initializes the function parameters based on the values found in an input array.
        procedure :: parameter_value       => ExponentialParametrized1DParameterValues     !< subroutine that returns the value of the function i-th parameter.
        procedure :: feedback              => ExponentialParametrized1DFeedback            !< subroutine that prints to screen the informations about the function.

        ! evaluation procedures:
        procedure :: value                 => ExponentialParametrized1DValue               !< function that returns the value of the exponential function.
        procedure :: first_derivative      => ExponentialParametrized1DFirstDerivative     !< function that returns the first derivative of the exponential function.
        procedure :: second_derivative     => ExponentialParametrized1DSecondDerivative    !< function that returns the second derivative of the exponential function.
        procedure :: third_derivative      => ExponentialParametrized1DThirdDerivative     !< function that returns the third derivative of the exponential function.
        procedure :: integral              => ExponentialParametrized1DIntegral            !< function that returns the strange integral that we need for w_DE.

    end type exponential_parametrization_1D

contains

    ! ---------------------------------------------------------------------------------------------
    ! Implementation of the exponential function.
    ! ---------------------------------------------------------------------------------------------

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that sets the number of parameters of the exponential parametrized function.
    subroutine ExponentialParametrized1DSetParamNumber( self )

        implicit none

        class(exponential_parametrization_1D) :: self       !< the base class

        ! initialize the number of parameters:
        self%parameter_number = 2

    end subroutine ExponentialParametrized1DSetParamNumber

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that initializes the function parameters based on the values found in an input array.
    subroutine ExponentialParametrized1DInitParams( self, array )

        implicit none

        class(exponential_parametrization_1D)                   :: self   !< the base class.
        real(dl), dimension(self%parameter_number), intent(in)  :: array  !< input array with the values of the parameters.

        self%coefficient = array(1)
        self%exponent = array(2)

    end subroutine ExponentialParametrized1DInitParams

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the value of the function i-th parameter.
    subroutine ExponentialParametrized1DParameterValues( self, i, value )

        implicit none

        class(exponential_parametrization_1D):: self        !< the base class
        integer     , intent(in)             :: i           !< The index of the parameter
        real(dl)    , intent(out)            :: value       !< the output value of the i-th parameter

        select case (i)
            case(1)
                value = self%coefficient
            case(2)
                value = self%exponent
            case default
                write(*,*) 'Illegal index for parameter_names.'
                write(*,*) 'Maximum value is:', self%parameter_number
                call MpiStop('EFTCAMB error')
        end select

    end subroutine ExponentialParametrized1DParameterValues

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that prints to screen the informations about the function.
    subroutine ExponentialParametrized1DFeedback( self )

        implicit none

        class(exponential_parametrization_1D) :: self         !< the base class

        integer                               :: i
        real(dl)                              :: param_value
        character(len=EFT_names_max_length)   :: param_name

        write(*,*)     'Exponential function: ', self%name
        do i=1, self%parameter_number
            call self%parameter_names( i, param_name  )
            call self%parameter_value( i, param_value )
            write(*,'(a23,a,F12.6)') param_name, '=', param_value
        end do

    end subroutine ExponentialParametrized1DFeedback

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the value of the function in the scale factor.
    function ExponentialParametrized1DValue( self, x, eft_cache )

        implicit none

        class(exponential_parametrization_1D)              :: self      !< the base class
        real(dl), intent(in)                               :: x         !< the input scale factor
        type(EFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: ExponentialParametrized1DValue                      !< the output value

        ExponentialParametrized1DValue = Exp(self%coefficient*x**self%exponent) -1._dl
    end function ExponentialParametrized1DValue

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the value of the first derivative, wrt scale factor, of the function.
    function ExponentialParametrized1DFirstDerivative( self, x, eft_cache )

        implicit none

        class(exponential_parametrization_1D)              :: self      !< the base class
        real(dl), intent(in)                               :: x         !< the input scale factor
        type(EFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: ExponentialParametrized1DFirstDerivative            !< the output value

        ExponentialParametrized1DFirstDerivative = self%coefficient*self%exponent*x**(self%exponent-1._dl)*Exp(self%coefficient*x**self%exponent)

    end function ExponentialParametrized1DFirstDerivative

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the second derivative of the function.
    function ExponentialParametrized1DSecondDerivative( self, x, eft_cache )

        implicit none

        class(exponential_parametrization_1D)              :: self      !< the base class
        real(dl), intent(in)                               :: x         !< the input scale factor
        type(EFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: ExponentialParametrized1DSecondDerivative           !< the output value

        ExponentialParametrized1DSecondDerivative = self%coefficient*self%exponent*(self%exponent -1._dl +x**self%exponent*self%coefficient*self%exponent)*x**(self%exponent-2._dl)*Exp(self%coefficient*x**self%exponent)

    end function ExponentialParametrized1DSecondDerivative

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the third derivative of the function.
    function ExponentialParametrized1DThirdDerivative( self, x, eft_cache )

        implicit none

        class(exponential_parametrization_1D)              :: self      !< the base class
        real(dl), intent(in)                               :: x         !< the input scale factor
        type(EFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: ExponentialParametrized1DThirdDerivative            !< the output value

        ExponentialParametrized1DThirdDerivative = self%coefficient*self%exponent*Exp(self%coefficient*x**self%exponent)*x**(self%exponent-3._dl)&
            & *(2._dl +self%exponent*(-3._dl +self%exponent +x**self%exponent*self%coefficient*(-3._dl +self%exponent*(3._dl +x**self%exponent*self%coefficient))))

    end function ExponentialParametrized1DThirdDerivative

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the integral of the function, as defined in the notes.
    function ExponentialParametrized1DIntegral( self, x, eft_cache )

        implicit none

        class(exponential_parametrization_1D)              :: self      !< the base class
        real(dl), intent(in)                               :: x         !< the input scale factor
        type(EFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: ExponentialParametrized1DIntegral                   !< the output value

        !< No analytic solution >!
        write(*,*) 'ExponentialParametrized1DIntegral is not implemented.'
        write(*,*) 'Calculations cannot proceed.'
        call MpiStop('EFTCAMB error')

    end function ExponentialParametrized1DIntegral

    ! ---------------------------------------------------------------------------------------------

end module EFTCAMB_exponential_parametrizations_1D

!----------------------------------------------------------------------------------------
