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

!> @file 05p1_constant_parametrizations_2D.f90
!! This file contains the definition of the constant parametrization, inheriting from
!! parametrized_function_2D.


!----------------------------------------------------------------------------------------
!> This module contains the definition of the constant parametrization, inheriting from
!! parametrized_function_2D.

!> @author Simone Peirone, Bin Hu, Marco Raveri

module EFTCAMB_constant_parametrization_2D

    use precision
    use AMLutils
    use EFT_def
    use EFTCAMB_cache
    use EFTCAMB_abstract_parametrizations_2D

    implicit none

    private

    public constant_parametrization_2D

    ! ---------------------------------------------------------------------------------------------
    !> Type containing the constant function parametrization. Inherits from parametrized_function_2D.
    type, extends ( parametrized_function_2D ) :: constant_parametrization_2D

        real(dl) :: constant_value

    contains

        ! utility functions:
        procedure :: set_param_number      => ConstantParametrized2DSetParamNumber      !< subroutine that sets the number of parameters of the constant parametrized function.
        procedure :: init_parameters       => ConstantParametrized2DInitParams          !< subroutine that initializes the function parameters based on the values found in an input array.
        procedure :: parameter_value       => ConstantParametrized2DParameterValues     !< subroutine that returns the value of the function i-th parameter.
        procedure :: feedback              => ConstantParametrized2DFeedback            !< subroutine that prints to screen the informations about the function.

        ! evaluation procedures:
        procedure :: value                 => ConstantParametrized2DValue               !< function that returns the value of the function
        procedure :: first_derivative_x    => ConstantParametrized2DFirstDerivativeX    !< function that returns the first partial derivative of the function with respect to x.
        procedure :: first_derivative_y    => ConstantParametrized2DFirstDerivativeY    !< function that returns the first partial derivative of the function with respect to y.
        procedure :: second_derivative_x   => ConstantParametrized2DSecondDerivativeX   !< function that returns the second derivative of the function.
        procedure :: second_derivative_y   => ConstantParametrized2DSecondDerivativeY   !< function that returns the second derivative of the function.
        procedure :: second_derivative_xy  => ConstantParametrized2DSecondDerivativeXY  !< function that returns the second derivative of the function.

    end type constant_parametrization_2D

contains

    ! ---------------------------------------------------------------------------------------------
    ! Implementation of the constant function.
    ! ---------------------------------------------------------------------------------------------

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that initializes the constant parametrization
    subroutine ConstantParametrized2DSetParamNumber( self )

        implicit none

        class(constant_parametrization_2D) :: self       !< the base class

        ! initialize the number of parameters:
        self%parameter_number = 1

    end subroutine ConstantParametrized2DSetParamNumber

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that initializes the function parameters based on the values found in an input array.
    subroutine ConstantParametrized2DInitParams( self, array )

        implicit none

        class(constant_parametrization_2D)                      :: self   !< the base class.
        real(dl), dimension(self%parameter_number), intent(in)  :: array  !< input array with the values of the parameters.

        self%constant_value = array(1)

    end subroutine ConstantParametrized2DInitParams

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that prints to screen the informations about the function.
    subroutine ConstantParametrized2DFeedback( self )

        implicit none

        class(constant_parametrization_2D)  :: self       !< the base class

        integer                             :: i
        real(dl)                            :: param_value
        character(len=EFT_names_max_length) :: param_name

        write(*,*)     'Constant function: ', self%name
        do i=1, self%parameter_number
            call self%parameter_names( i, param_name  )
            call self%parameter_value( i, param_value )
            write(*,'(a23,a,F12.6)') param_name, '=', param_value
        end do

    end subroutine ConstantParametrized2DFeedback

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the value of the function i-th parameter.
    subroutine ConstantParametrized2DParameterValues( self, i, value )

        implicit none

        class(constant_parametrization_2D) :: self        !< the base class
        integer     , intent(in)           :: i           !< The index of the parameter
        real(dl)    , intent(out)          :: value       !< the output value of the i-th parameter

        select case (i)
            case(1)
                value = self%constant_value
            case default
                write(*,*) 'In constant_parametrization_2D:', self%name
                write(*,*) 'Illegal index for parameter_names.'
                write(*,*) 'Maximum value is:', self%parameter_number
                call MpiStop('EFTCAMB error')
        end select

    end subroutine ConstantParametrized2DParameterValues

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the value of the constant function.
    function ConstantParametrized2DValue( self, x, y, eft_cache )

        implicit none

        class(constant_parametrization_2D)                 :: self      !< the base class
        real(dl), intent(in)                               :: x         !< the first input variable
        real(dl), intent(in)                               :: y         !< the second input variable
        type(EFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: ConstantParametrized2DValue                         !< the output value

        ConstantParametrized2DValue = self%constant_value

    end function ConstantParametrized2DValue

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the first partial derivative of the constant function with respect to x.
    function ConstantParametrized2DFirstDerivativeX( self, x, y, eft_cache )

        implicit none

        class(constant_parametrization_2D)                 :: self      !< the base class
        real(dl), intent(in)                               :: x         !< the input first variable
        real(dl), intent(in)                               :: y         !< the input second variable
        type(EFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: ConstantParametrized2DFirstDerivativeX              !< the output value

        ConstantParametrized2DFirstDerivativeX = 0._dl

    end function ConstantParametrized2DFirstDerivativeX

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the first partial derivative of the constant function with respect to y.
    function ConstantParametrized2DFirstDerivativeY( self, x, y, eft_cache )

        implicit none

        class(constant_parametrization_2D)                 :: self      !< the base class
        real(dl), intent(in)                               :: x         !< the input first variable
        real(dl), intent(in)                               :: y         !< the input second variable
        type(EFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: ConstantParametrized2DFirstDerivativeY              !< the output value

        ConstantParametrized2DFirstDerivativeY = 0._dl

    end function ConstantParametrized2DFirstDerivativeY

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the second partial derivative of the constant function with respect to x.
    function ConstantParametrized2DSecondDerivativeX( self, x, y, eft_cache )

        implicit none

        class(constant_parametrization_2D)                 :: self      !< the base class
        real(dl), intent(in)                               :: x         !< the input first variable
        real(dl), intent(in)                               :: y         !< the input second variable
        type(EFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: ConstantParametrized2DSecondDerivativeX             !< the output value

        ConstantParametrized2DSecondDerivativeX = 0._dl

    end function ConstantParametrized2DSecondDerivativeX

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the second partial derivative of the constant function with respect to y.
    function ConstantParametrized2DSecondDerivativeY( self, x, y, eft_cache )

        implicit none

        class(constant_parametrization_2D)                 :: self      !< the base class
        real(dl), intent(in)                               :: x         !< the input first variable
        real(dl), intent(in)                               :: y         !< the input second variable
        type(EFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: ConstantParametrized2DSecondDerivativeY             !< the output value

        ConstantParametrized2DSecondDerivativeY = 0._dl

    end function ConstantParametrized2DSecondDerivativeY

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the mixed partial derivative of the constant function with respect to x and y.
    function ConstantParametrized2DSecondDerivativeXY( self, x, y, eft_cache )

        implicit none

        class(constant_parametrization_2D)                 :: self      !< the base class
        real(dl), intent(in)                               :: x         !< the input first variable
        real(dl), intent(in)                               :: y         !< the input second variable
        type(EFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: ConstantParametrized2DSecondDerivativeXY            !< the output value

        ConstantParametrized2DSecondDerivativeXY = 0._dl

    end function ConstantParametrized2DSecondDerivativeXY

    !----------------------------------------------------------------------------------------

end module EFTCAMB_constant_parametrization_2D

!----------------------------------------------------------------------------------------
