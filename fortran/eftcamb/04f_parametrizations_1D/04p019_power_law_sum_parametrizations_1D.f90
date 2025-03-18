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

!> @file 04p019_power_law_sum_parametrizations_1D.f90
!! This file contains the definition of the power law parametrization,
!! inheriting from parametrized_function_1D.


!----------------------------------------------------------------------------------------
!> This module contains the definition of the power law parametrization,
!! inheriting from parametrized_function_1D.
!! This consists in the sum of three power laws, each with a different exponent and a 
!! constant offset.

!> @author Marco Raveri

module EFTCAMB_power_law_sum_parametrizations_1D

    use precision
    use EFT_def
    use MpiUtils
    use EFTCAMB_cache
    use EFTCAMB_abstract_parametrizations_1D

    implicit none

    private

    public power_law_sum_parametrization_1D

    ! ---------------------------------------------------------------------------------------------
    !> Type containing the power law function parametrization. Inherits from parametrized_function_1D.
    type, extends ( parametrized_function_1D ) :: power_law_sum_parametrization_1D

        real(dl) :: offset
        real(dl) :: coefficient_1
        real(dl) :: coefficient_2
        real(dl) :: coefficient_3
        real(dl) :: exponent_1
        real(dl) :: exponent_2
        real(dl) :: exponent_3

    contains

        ! utility functions:
        procedure :: set_param_number      => PowerLawSumParametrized1DSetParamNumber      !< subroutine that sets the number of parameters of the power law parametrized function.
        procedure :: init_parameters       => PowerLawSumParametrized1DInitParams          !< subroutine that initializes the function parameters based on the values found in an input array.
        procedure :: parameter_value       => PowerLawSumParametrized1DParameterValues     !< subroutine that returns the value of the function i-th parameter.
        procedure :: feedback              => PowerLawSumParametrized1DFeedback            !< subroutine that prints to screen the informations about the function.

        ! evaluation procedures:
        procedure :: value                 => PowerLawSumParametrized1DValue               !< function that returns the value of the power law function.
        procedure :: first_derivative      => PowerLawSumParametrized1DFirstDerivative     !< function that returns the first derivative of the power law function.
        procedure :: second_derivative     => PowerLawSumParametrized1DSecondDerivative    !< function that returns the second derivative of the power law function.
        procedure :: third_derivative      => PowerLawSumParametrized1DThirdDerivative     !< function that returns the third derivative of the power law function.
        procedure :: fourth_derivative     => PowerLawSumParametrized1DFourthDerivative    !< function that returns the fourth derivative of the power law function.
        procedure :: integral              => PowerLawSumParametrized1DIntegral            !< function that returns the strange integral that we need for w_DE.

    end type power_law_sum_parametrization_1D

contains

    ! ---------------------------------------------------------------------------------------------
    ! Implementation of the power law function.
    ! ---------------------------------------------------------------------------------------------

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that sets the number of parameters of the power law parametrized function.
    subroutine PowerLawSumParametrized1DSetParamNumber( self )

        implicit none

        class(power_law_sum_parametrization_1D) :: self       !< the base class

        ! initialize the number of parameters:
        self%parameter_number = 7

    end subroutine PowerLawSumParametrized1DSetParamNumber

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that initializes the function parameters based on the values found in an input array.
    subroutine PowerLawSumParametrized1DInitParams( self, array )

        implicit none

        class(power_law_sum_parametrization_1D)                     :: self   !< the base class.
        real(dl), dimension(self%parameter_number), intent(in)  :: array  !< input array with the values of the parameters.

        self%offset = array(1)
        self%coefficient_1 = array(2)
        self%exponent_1 = array(3)
        self%coefficient_2 = array(4)
        self%exponent_2 = array(5)
        self%coefficient_3 = array(6)
        self%exponent_3 = array(7)

    end subroutine PowerLawSumParametrized1DInitParams

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the value of the function i-th parameter.
    subroutine PowerLawSumParametrized1DParameterValues( self, i, value )

        implicit none

        class(power_law_sum_parametrization_1D) :: self        !< the base class
        integer     , intent(in)            :: i           !< The index of the parameter
        real(dl)    , intent(out)           :: value       !< the output value of the i-th parameter

        select case (i)
            case(1)
                value = self%offset
            case(2)
                value = self%coefficient_1
            case(3)
                value = self%exponent_1
            case(4)
                value = self%coefficient_2
            case(5)
                value = self%exponent_2
            case(6)
                value = self%coefficient_3
            case(7)
                value = self%exponent_3
            case default
                write(*,*) 'Illegal index for parameter_names.'
                write(*,*) 'Maximum value is:', self%parameter_number
                call MpiStop('EFTCAMB error')
        end select

    end subroutine PowerLawSumParametrized1DParameterValues

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that prints to screen the informations about the function.
    subroutine PowerLawSumParametrized1DFeedback( self, print_params )

        implicit none

        class(power_law_sum_parametrization_1D) :: self         !< the base class
        logical, optional                   :: print_params !< optional flag that decised whether to print numerical values
                                                            !! of the parameters.

        integer                             :: i
        real(dl)                            :: param_value
        character(len=EFT_names_max_length) :: param_name
        logical                             :: print_params_temp

        if ( present(print_params) ) then
            print_params_temp = print_params
        else
            print_params_temp = .True.
        end if

        write(*,*)     'Power Law function: ', self%name
        if ( print_params_temp ) then
            do i=1, self%parameter_number
                call self%parameter_names( i, param_name  )
                call self%parameter_value( i, param_value )
                write(*,'(a23,a,F12.6)') param_name, '=', param_value
            end do
        end if

    end subroutine PowerLawSumParametrized1DFeedback

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the value of the function in the scale factor.
    function PowerLawSumParametrized1DValue( self, x, eft_cache )

        implicit none

        class(power_law_sum_parametrization_1D)                :: self      !< the base class
        real(dl), intent(in)                               :: x         !< the input scale factor
        type(TEFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: PowerLawSumParametrized1DValue                         !< the output value

        PowerLawSumParametrized1DValue = self%offset  &
            & + self%coefficient_1 * x**self%exponent_1 &
            & + self%coefficient_2 * x**self%exponent_2 & 
            & + self%coefficient_3 * x**self%exponent_3

    end function PowerLawSumParametrized1DValue

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the value of the first derivative, wrt scale factor, of the function.
    function PowerLawSumParametrized1DFirstDerivative( self, x, eft_cache )

        implicit none

        class(power_law_sum_parametrization_1D)                :: self      !< the base class
        real(dl), intent(in)                               :: x         !< the input scale factor
        type(TEFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: PowerLawSumParametrized1DFirstDerivative               !< the output value

        PowerLawSumParametrized1DFirstDerivative = &
            & + self%coefficient_1*self%exponent_1*x**(self%exponent_1-1._dl) &
            & + self%coefficient_2*self%exponent_2*x**(self%exponent_2-1._dl) &
            & + self%coefficient_3*self%exponent_3*x**(self%exponent_3-1._dl)

    end function PowerLawSumParametrized1DFirstDerivative

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the second derivative of the function.
    function PowerLawSumParametrized1DSecondDerivative( self, x, eft_cache )

        implicit none

        class(power_law_sum_parametrization_1D)                :: self      !< the base class
        real(dl), intent(in)                               :: x         !< the input scale factor
        type(TEFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: PowerLawSumParametrized1DSecondDerivative              !< the output value

        PowerLawSumParametrized1DSecondDerivative = &
            & + self%coefficient_1*self%exponent_1*(self%exponent_1-1._dl)*x**(self%exponent_1-2._dl) &
            & + self%coefficient_2*self%exponent_2*(self%exponent_2-1._dl)*x**(self%exponent_2-2._dl) &
            & + self%coefficient_3*self%exponent_3*(self%exponent_3-1._dl)*x**(self%exponent_3-2._dl)

    end function PowerLawSumParametrized1DSecondDerivative

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the third derivative of the function.
    function PowerLawSumParametrized1DThirdDerivative( self, x, eft_cache )

        implicit none

        class(power_law_sum_parametrization_1D)             :: self      !< the base class
        real(dl), intent(in)                                :: x         !< the input scale factor
        type(TEFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: PowerLawSumParametrized1DThirdDerivative               !< the output value

        PowerLawSumParametrized1DThirdDerivative = &
            & + self%coefficient_1*self%exponent_1*(self%exponent_1-1._dl)*(self%exponent_1-2._dl)*x**(self%exponent_1-3._dl) &
            & + self%coefficient_2*self%exponent_2*(self%exponent_2-1._dl)*(self%exponent_2-2._dl)*x**(self%exponent_2-3._dl) &
            & + self%coefficient_3*self%exponent_3*(self%exponent_3-1._dl)*(self%exponent_3-2._dl)*x**(self%exponent_3-3._dl)
        
    end function PowerLawSumParametrized1DThirdDerivative

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the fourth derivative of the function.
    function PowerLawSumParametrized1DFourthDerivative( self, x, eft_cache )

        implicit none

        class(power_law_sum_parametrization_1D)                :: self      !< the base class
        real(dl), intent(in)                               :: x         !< the input scale factor
        type(TEFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: PowerLawSumParametrized1DFourthDerivative               !< the output value

        PowerLawSumParametrized1DFourthDerivative = &
            & + self%coefficient_1*self%exponent_1*(self%exponent_1-1._dl)*(self%exponent_1-2._dl)*(self%exponent_1-3._dl) * x**(self%exponent_1-4._dl) &
            & + self%coefficient_2*self%exponent_2*(self%exponent_2-1._dl)*(self%exponent_2-2._dl)*(self%exponent_2-3._dl)*x**(self%exponent_2-4._dl) &
            & + self%coefficient_3*self%exponent_3*(self%exponent_3-1._dl)*(self%exponent_3-2._dl)*(self%exponent_3-3._dl)*x**(self%exponent_3-4._dl) 
        
    end function PowerLawSumParametrized1DFourthDerivative

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the integral of the function, as defined in the notes.
    function PowerLawSumParametrized1DIntegral( self, x, eft_cache )

        implicit none

        class(power_law_sum_parametrization_1D)                :: self      !< the base class
        real(dl), intent(in)                               :: x         !< the input scale factor
        type(TEFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: PowerLawSumParametrized1DIntegral                      !< the output value

        PowerLawSumParametrized1DIntegral = 0._dl
        !< No analytic solution >!
        write(*,*) 'PowerLawSumParametrized1DIntegral is not implemented.'
        write(*,*) 'Calculations cannot proceed.'
        call MpiStop('EFTCAMB error')

    end function PowerLawSumParametrized1DIntegral

    ! ---------------------------------------------------------------------------------------------

end module EFTCAMB_power_law_sum_parametrizations_1D

!----------------------------------------------------------------------------------------
