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

!> @file 04p3_power_law_parametrizations_1D.f90
!! This file contains the definition of the power law parametrization,
!! inheriting from parametrized_function_1D.


!----------------------------------------------------------------------------------------
!> This module contains the definition of the power law parametrization,
!! inheriting from parametrized_function_1D.

!> @author Bin Hu, Marco Raveri, Simone Peirone

module EFTCAMB_power_law_parametrizations_1D

    use precision
    use AMLutils
    use EFT_def
    use EFTCAMB_cache
    use EFTCAMB_abstract_parametrizations_1D

    implicit none

    private

    public power_law_parametrization_1D

    ! ---------------------------------------------------------------------------------------------
    !> Type containing the power law function parametrization. Inherits from parametrized_function_1D.
    type, extends ( parametrized_function_1D ) :: power_law_parametrization_1D

        real(dl) :: coefficient
        real(dl) :: exponent

    contains

        ! utility functions:
        procedure :: set_param_number      => PowerLawParametrized1DSetParamNumber      !< subroutine that sets the number of parameters of the power law parametrized function.
        procedure :: init_parameters       => PowerLawParametrized1DInitParams          !< subroutine that initializes the function parameters based on the values found in an input array.
        procedure :: parameter_value       => PowerLawParametrized1DParameterValues     !< subroutine that returns the value of the function i-th parameter.
        procedure :: feedback              => PowerLawParametrized1DFeedback            !< subroutine that prints to screen the informations about the function.

        ! evaluation procedures:
        procedure :: value                 => PowerLawParametrized1DValue               !< function that returns the value of the power law function.
        procedure :: first_derivative      => PowerLawParametrized1DFirstDerivative     !< function that returns the first derivative of the power law function.
        procedure :: second_derivative     => PowerLawParametrized1DSecondDerivative    !< function that returns the second derivative of the power law function.
        procedure :: third_derivative      => PowerLawParametrized1DThirdDerivative     !< function that returns the third derivative of the power law function.
        procedure :: integral              => PowerLawParametrized1DIntegral            !< function that returns the strange integral that we need for w_DE.

    end type power_law_parametrization_1D

contains

    ! ---------------------------------------------------------------------------------------------
    ! Implementation of the power law function.
    ! ---------------------------------------------------------------------------------------------

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that sets the number of parameters of the power law parametrized function.
    subroutine PowerLawParametrized1DSetParamNumber( self )

        implicit none

        class(power_law_parametrization_1D) :: self       !< the base class

        ! initialize the number of parameters:
        self%parameter_number = 2

    end subroutine PowerLawParametrized1DSetParamNumber

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that initializes the function parameters based on the values found in an input array.
    subroutine PowerLawParametrized1DInitParams( self, array )

        implicit none

        class(power_law_parametrization_1D)                     :: self   !< the base class.
        real(dl), dimension(self%parameter_number), intent(in)  :: array  !< input array with the values of the parameters.

        self%coefficient = array(1)
        self%exponent = array(2)

    end subroutine PowerLawParametrized1DInitParams

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the value of the function i-th parameter.
    subroutine PowerLawParametrized1DParameterValues( self, i, value )

        implicit none

        class(power_law_parametrization_1D) :: self        !< the base class
        integer     , intent(in)            :: i           !< The index of the parameter
        real(dl)    , intent(out)           :: value       !< the output value of the i-th parameter

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

    end subroutine PowerLawParametrized1DParameterValues

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that prints to screen the informations about the function.
    subroutine PowerLawParametrized1DFeedback( self )

        implicit none

        class(power_law_parametrization_1D) :: self         !< the base class

        integer                             :: i
        real(dl)                            :: param_value
        character(len=EFT_names_max_length) :: param_name

        write(*,*)     'Power Law function: ', self%name
        do i=1, self%parameter_number
            call self%parameter_names( i, param_name  )
            call self%parameter_value( i, param_value )
            write(*,'(a23,a,F12.6)') param_name, '=', param_value
        end do

    end subroutine PowerLawParametrized1DFeedback

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the value of the function in the scale factor.
    function PowerLawParametrized1DValue( self, x, eft_cache )

        implicit none

        class(power_law_parametrization_1D)                :: self      !< the base class
        real(dl), intent(in)                               :: x         !< the input scale factor
        type(EFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: PowerLawParametrized1DValue                         !< the output value

        PowerLawParametrized1DValue = self%coefficient*x**self%exponent

    end function PowerLawParametrized1DValue

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the value of the first derivative, wrt scale factor, of the function.
    function PowerLawParametrized1DFirstDerivative( self, x, eft_cache )

        implicit none

        class(power_law_parametrization_1D)                :: self      !< the base class
        real(dl), intent(in)                               :: x         !< the input scale factor
        type(EFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: PowerLawParametrized1DFirstDerivative               !< the output value

        PowerLawParametrized1DFirstDerivative = self%coefficient*self%exponent*x**(self%exponent-1._dl)

    end function PowerLawParametrized1DFirstDerivative

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the second derivative of the function.
    function PowerLawParametrized1DSecondDerivative( self, x, eft_cache )

        implicit none

        class(power_law_parametrization_1D)                :: self      !< the base class
        real(dl), intent(in)                               :: x         !< the input scale factor
        type(EFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: PowerLawParametrized1DSecondDerivative              !< the output value

        PowerLawParametrized1DSecondDerivative = self%coefficient*self%exponent*(self%exponent-1._dl)*x**(self%exponent-2._dl)

    end function PowerLawParametrized1DSecondDerivative

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the third derivative of the function.
    function PowerLawParametrized1DThirdDerivative( self, x, eft_cache )

        implicit none

        class(power_law_parametrization_1D)                :: self      !< the base class
        real(dl), intent(in)                               :: x         !< the input scale factor
        type(EFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: PowerLawParametrized1DThirdDerivative               !< the output value

        PowerLawParametrized1DThirdDerivative = self%coefficient*self%exponent*(self%exponent-1._dl)*(self%exponent-2._dl)*x**(self%exponent-3._dl)

    end function PowerLawParametrized1DThirdDerivative

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the integral of the function, as defined in the notes.
    function PowerLawParametrized1DIntegral( self, x, eft_cache )

        implicit none

        class(power_law_parametrization_1D)                :: self      !< the base class
        real(dl), intent(in)                               :: x         !< the input scale factor
        type(EFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: PowerLawParametrized1DIntegral                      !< the output value

        if ( self%exponent == 0. ) then
          PowerLawParametrized1DIntegral = x**(-1._dl -3._dl*self%coefficient)
        else
          PowerLawParametrized1DIntegral = 1._dl/(x*exp((3._dl*(-1._dl + x**self%exponent)*self%coefficient)/self%exponent))
        end if

    end function PowerLawParametrized1DIntegral

    ! ---------------------------------------------------------------------------------------------

end module EFTCAMB_power_law_parametrizations_1D

!----------------------------------------------------------------------------------------
