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

!> @file 04p012_exponential_2_parametrizations_1D.f90
!! This file contains the definition of the second type of exponential parametrization,
!! inheriting from parametrized_function_1D.


!----------------------------------------------------------------------------------------
!> This module contains the definition of the second type of exponential parametrization,
!! inheriting from parametrized_function_1D.

!> @author Marco Raveri

module EFTCAMB_exponential_parametrizations_2_1D

    use precision
    use EFT_def
    use MpiUtils
    use EFTCAMB_cache
    use EFTCAMB_abstract_parametrizations_1D

    implicit none

    private

    public exponential_parametrization_2_1D

    ! ---------------------------------------------------------------------------------------------
    !> Type containing the exponential function parametrization. Inherits from parametrized_function_1D.
    type, extends ( parametrized_function_1D ) :: exponential_parametrization_2_1D

        real(dl) :: coefficient
        real(dl) :: exponent

    contains

        ! utility functions:
        procedure :: set_param_number      => Exponential2Parametrized1DSetParamNumber      !< subroutine that sets the number of parameters of the exponential parametrized function.
        procedure :: init_parameters       => Exponential2Parametrized1DInitParams          !< subroutine that initializes the function parameters based on the values found in an input array.
        procedure :: parameter_value       => Exponential2Parametrized1DParameterValues     !< subroutine that returns the value of the function i-th parameter.
        procedure :: feedback              => Exponential2Parametrized1DFeedback            !< subroutine that prints to screen the informations about the function.

        ! evaluation procedures:
        procedure :: value                 => Exponential2Parametrized1DValue               !< function that returns the value of the exponential function.
        procedure :: first_derivative      => Exponential2Parametrized1DFirstDerivative     !< function that returns the first derivative of the exponential function.
        procedure :: second_derivative     => Exponential2Parametrized1DSecondDerivative    !< function that returns the second derivative of the exponential function.
        procedure :: third_derivative      => Exponential2Parametrized1DThirdDerivative     !< function that returns the third derivative of the exponential function.
        procedure :: fourth_derivative     => Exponential2Parametrized1DFourthDerivative    !< function that returns the fourth derivative of the exponential function.
        procedure :: integral              => Exponential2Parametrized1DIntegral            !< function that returns the strange integral that we need for w_DE.

    end type exponential_parametrization_2_1D

contains

    ! ---------------------------------------------------------------------------------------------
    ! Implementation of the exponential function.
    ! ---------------------------------------------------------------------------------------------

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that sets the number of parameters of the exponential parametrized function.
    subroutine Exponential2Parametrized1DSetParamNumber( self )

        implicit none

        class(exponential_parametrization_2_1D) :: self       !< the base class

        ! initialize the number of parameters:
        self%parameter_number = 2

    end subroutine Exponential2Parametrized1DSetParamNumber

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that initializes the function parameters based on the values found in an input array.
    subroutine Exponential2Parametrized1DInitParams( self, array )

        implicit none

        class(exponential_parametrization_2_1D)                   :: self   !< the base class.
        real(dl), dimension(self%parameter_number), intent(in)  :: array  !< input array with the values of the parameters.

        self%coefficient = array(1)
        self%exponent    = array(2)

    end subroutine Exponential2Parametrized1DInitParams

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the value of the function i-th parameter.
    subroutine Exponential2Parametrized1DParameterValues( self, i, value )

        implicit none

        class(exponential_parametrization_2_1D):: self        !< the base class
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

    end subroutine Exponential2Parametrized1DParameterValues

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that prints to screen the informations about the function.
    subroutine Exponential2Parametrized1DFeedback( self, print_params )

        implicit none

        class(exponential_parametrization_2_1D) :: self         !< the base class
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

    end subroutine Exponential2Parametrized1DFeedback

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the value of the function in the scale factor.
    function Exponential2Parametrized1DValue( self, x, eft_cache )

        implicit none

        class(exponential_parametrization_2_1D)              :: self      !< the base class
        real(dl), intent(in)                               :: x         !< the input scale factor
        type(TEFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: Exponential2Parametrized1DValue                      !< the output value

        !Exponential2Parametrized1DValue = Exp(self%coefficient*x**self%exponent) -1._dl
        Exponential2Parametrized1DValue = self%coefficient*Exp(-x*self%exponent)

    end function Exponential2Parametrized1DValue

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the value of the first derivative, wrt scale factor, of the function.
    function Exponential2Parametrized1DFirstDerivative( self, x, eft_cache )

        implicit none

        class(exponential_parametrization_2_1D)              :: self      !< the base class
        real(dl), intent(in)                               :: x         !< the input scale factor
        type(TEFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: Exponential2Parametrized1DFirstDerivative            !< the output value

        Exponential2Parametrized1DFirstDerivative = -self%exponent*self%coefficient*Exp(-x*self%exponent)

    end function Exponential2Parametrized1DFirstDerivative

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the second derivative of the function.
    function Exponential2Parametrized1DSecondDerivative( self, x, eft_cache )

        implicit none

        class(exponential_parametrization_2_1D)              :: self      !< the base class
        real(dl), intent(in)                               :: x         !< the input scale factor
        type(TEFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: Exponential2Parametrized1DSecondDerivative           !< the output value

        Exponential2Parametrized1DSecondDerivative = self%exponent**2*self%coefficient*Exp(-x*self%exponent)

    end function Exponential2Parametrized1DSecondDerivative

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the third derivative of the function.
    function Exponential2Parametrized1DThirdDerivative( self, x, eft_cache )

        implicit none

        class(exponential_parametrization_2_1D)              :: self      !< the base class
        real(dl), intent(in)                               :: x         !< the input scale factor
        type(TEFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: Exponential2Parametrized1DThirdDerivative            !< the output value

        Exponential2Parametrized1DThirdDerivative = -self%exponent**3*self%coefficient*Exp(-x*self%exponent)

    end function Exponential2Parametrized1DThirdDerivative

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the fourth derivative of the function.
    function Exponential2Parametrized1DFourthDerivative( self, x, eft_cache )

        implicit none

        class(exponential_parametrization_2_1D)              :: self      !< the base class
        real(dl), intent(in)                               :: x         !< the input scale factor
        type(TEFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: Exponential2Parametrized1DFourthDerivative            !< the output value

        Exponential2Parametrized1DFourthDerivative = self%exponent**4*self%coefficient*Exp(-x*self%exponent)

    end function Exponential2Parametrized1DFourthDerivative

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the integral of the function, as defined in the notes.
    function Exponential2Parametrized1DIntegral( self, x, eft_cache )

        implicit none

        class(exponential_parametrization_2_1D)              :: self      !< the base class
        real(dl), intent(in)                               :: x         !< the input scale factor
        type(TEFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: Exponential2Parametrized1DIntegral                   !< the output value

        Exponential2Parametrized1DIntegral = 0._dl

        !< No analytic solution >!
        write(*,*) 'Exponential2Parametrized1DIntegral is not implemented.'
        write(*,*) 'Calculations cannot proceed.'
        call MpiStop('EFTCAMB error')

    end function Exponential2Parametrized1DIntegral

    ! ---------------------------------------------------------------------------------------------

end module EFTCAMB_exponential_parametrizations_2_1D

!----------------------------------------------------------------------------------------
