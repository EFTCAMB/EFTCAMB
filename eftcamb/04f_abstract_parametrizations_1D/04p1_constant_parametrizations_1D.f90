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

!> @file 04p1_constant_parametrizations_1D.f90
!! This file contains the definition of the constant parametrization, inheriting from
!! parametrized_function_1D.


!----------------------------------------------------------------------------------------
!> This module contains the definition of the constant parametrization, inheriting from
!! parametrized_function_1D.

!> @author Bin Hu, Marco Raveri, Simone Peirone

module EFTCAMB_constant_parametrization_1D

    use precision
    use AMLutils
    use EFT_def
    use EFTCAMB_cache
    use EFTCAMB_abstract_parametrizations_1D

    implicit none

    private

    public constant_parametrization_1D

    ! ---------------------------------------------------------------------------------------------
    !> Type containing the constant function parametrization. Inherits from parametrized_function_1D.
    type, extends ( parametrized_function_1D ) :: constant_parametrization_1D

        real(dl) :: constant_value

    contains

        ! utility functions:
        procedure :: set_param_number      => ConstantParametrized1DSetParamNumber      !< subroutine that sets the number of parameters of the constant parametrized function.
        procedure :: init_parameters       => ConstantParametrized1DInitParams          !< subroutine that initializes the function parameters based on the values found in an input array.
        procedure :: parameter_value       => ConstantParametrized1DParameterValues     !< subroutine that returns the value of the function i-th parameter.
        procedure :: feedback              => ConstantParametrized1DFeedback            !< subroutine that prints to screen the informations about the function.

        ! evaluation procedures:
        procedure :: value                 => ConstantParametrized1DValue               !< function that returns the value of the constant function.
        procedure :: first_derivative      => ConstantParametrized1DFirstDerivative     !< function that returns the first derivative of the constant function.
        procedure :: second_derivative     => ConstantParametrized1DSecondDerivative    !< function that returns the second derivative of the constant function.
        procedure :: third_derivative      => ConstantParametrized1DThirdDerivative     !< function that returns the third derivative of the constant function.
        procedure :: integral              => ConstantParametrized1DIntegral            !< function that returns the strange integral that we need for w_DE.

    end type constant_parametrization_1D

contains

    ! ---------------------------------------------------------------------------------------------
    ! Implementation of the constant function.
    ! ---------------------------------------------------------------------------------------------

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that sets the number of parameters of the constant parametrized function.
    subroutine ConstantParametrized1DSetParamNumber( self )

        implicit none

        class(constant_parametrization_1D) :: self       !< the base class

        ! initialize the number of parameters:
        self%parameter_number = 1

    end subroutine ConstantParametrized1DSetParamNumber

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that initializes the function parameters based on the values found in an input array.
    subroutine ConstantParametrized1DInitParams( self, array )

        implicit none

        class(constant_parametrization_1D)                      :: self   !< the base class.
        real(dl), dimension(self%parameter_number), intent(in)  :: array  !< input array with the values of the parameters.

        self%constant_value = array(1)

    end subroutine ConstantParametrized1DInitParams

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that prints to screen the informations about the function.
    subroutine ConstantParametrized1DFeedback( self )

        implicit none

        class(constant_parametrization_1D)  :: self       !< the base class

        integer                             :: i
        real(dl)                            :: param_value
        character(len=EFT_names_max_length) :: param_name

        write(*,*)     'Constant function: ', self%name
        do i=1, self%parameter_number
            call self%parameter_names( i, param_name  )
            call self%parameter_value( i, param_value )
            write(*,'(a23,a,F12.6)') param_name, '=', param_value
        end do

    end subroutine ConstantParametrized1DFeedback

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the value of the function i-th parameter.
    subroutine ConstantParametrized1DParameterValues( self, i, value )

        implicit none

        class(constant_parametrization_1D) :: self        !< the base class
        integer     , intent(in)           :: i           !< The index of the parameter
        real(dl)    , intent(out)          :: value       !< the output value of the i-th parameter

        select case (i)
            case(1)
                value = self%constant_value
            case default
                write(*,*) 'In constant_parametrization_1D:', self%name
                write(*,*) 'Illegal index for parameter_names.'
                write(*,*) 'Maximum value is:', self%parameter_number
                call MpiStop('EFTCAMB error')
        end select

    end subroutine ConstantParametrized1DParameterValues

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the value of the constant function.
    function ConstantParametrized1DValue( self, x, eft_cache )

        implicit none

        class(constant_parametrization_1D)                 :: self      !< the base class
        real(dl), intent(in)                               :: x         !< the input scale factor
        type(EFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: ConstantParametrized1DValue                         !< the output value

        ConstantParametrized1DValue = self%constant_value

    end function ConstantParametrized1DValue

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the first derivative of the constant function.
    function ConstantParametrized1DFirstDerivative( self, x, eft_cache )

        implicit none

        class(constant_parametrization_1D)                 :: self      !< the base class
        real(dl), intent(in)                               :: x         !< the input scale factor
        type(EFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: ConstantParametrized1DFirstDerivative               !< the output value

        ConstantParametrized1DFirstDerivative = 0._dl

    end function ConstantParametrized1DFirstDerivative

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the second derivative of the constant function.
    function ConstantParametrized1DSecondDerivative( self, x, eft_cache )

        implicit none

        class(constant_parametrization_1D)                 :: self      !< the base class
        real(dl), intent(in)                               :: x         !< the input scale factor
        type(EFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: ConstantParametrized1DSecondDerivative              !< the output value

        ConstantParametrized1DSecondDerivative = 0._dl

    end function ConstantParametrized1DSecondDerivative

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the third derivative of the constant function.
    function ConstantParametrized1DThirdDerivative( self, x, eft_cache )

        implicit none

        class(constant_parametrization_1D)                 :: self      !< the base class
        real(dl), intent(in)                               :: x         !< the input scale factor
        type(EFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: ConstantParametrized1DThirdDerivative               !< the output value

        ConstantParametrized1DThirdDerivative = 0._dl

    end function ConstantParametrized1DThirdDerivative

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the strange integral that we need for w_DE.
    function ConstantParametrized1DIntegral( self, x, eft_cache )

        implicit none

        class(constant_parametrization_1D)                 :: self      !< the base class
        real(dl), intent(in)                               :: x         !< the input scale factor
        type(EFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: ConstantParametrized1DIntegral                      !< the output value

        ConstantParametrized1DIntegral = x**(-1._dl-3._dl*self%constant_value)

    end function ConstantParametrized1DIntegral

    ! ---------------------------------------------------------------------------------------------

end module EFTCAMB_constant_parametrization_1D

!----------------------------------------------------------------------------------------
