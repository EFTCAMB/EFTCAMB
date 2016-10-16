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

!> @file 05p0_neutral_parametrizations_2D.f90
!! This file contains the definition of neutral parametrizations that can be used when
!! parametrized functions should take their GR value.


!----------------------------------------------------------------------------------------
!> This module contains the definition of neutral parametrizations that can be used when
!! parametrized functions should take their GR value.

!> @author Simone Peirone, Bin Hu, Marco Raveri

module EFTCAMB_neutral_parametrization_2D

    use precision
    use EFTCAMB_cache
    use EFTCAMB_abstract_parametrizations_2D

    implicit none

    private

    public zero_parametrization_2D

    ! ---------------------------------------------------------------------------------------------
    !> Type containing the constant function parametrization that returns all zero values.
    !! Inherits from parametrized_function_2D.
    type, extends ( parametrized_function_2D ) :: zero_parametrization_2D

    contains

        ! initialization:
        procedure :: set_param_number      => ZeroParametrized2DSetParamNumber     !< subroutine that sets the number of parameters of the zero function = 0.
        procedure :: init_parameters       => ZeroParametrized2DInitParams         !< subroutine that initializes the parameters based on the values found in an input array.
        procedure :: parameter_value       => ZeroParametrized2DParameterValues    !< subroutine that returns the value of the zero function i-th parameter.

        ! evaluation procedures:
        procedure :: value                 => ZeroParametrized2DValue               !< function that returns the value of the zero function.
        procedure :: first_derivative_x    => ZeroParametrized2DFirstDerivativeX    !< function that returns the first partial derivative of the zero function with respect to x.
        procedure :: first_derivative_y    => ZeroParametrized2DFirstDerivativeY    !< function that returns the first partial derivative of the zero function with respect to y.
        procedure :: second_derivative_x   => ZeroParametrized2DSecondDerivativeX   !< function that returns the second derivative of the zero function.
        procedure :: second_derivative_y   => ZeroParametrized2DSecondDerivativeY   !< function that returns the second derivative of the zero function.
        procedure :: second_derivative_xy  => ZeroParametrized2DSecondDerivativeXY  !< function that returns the second derivative of the zero function.

    end type zero_parametrization_2D

contains

    ! ---------------------------------------------------------------------------------------------
    ! Implementation of the constant zero function.
    ! ---------------------------------------------------------------------------------------------

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that sets the number of parameters of the zero function = 0.
    subroutine ZeroParametrized2DSetParamNumber( self )

        implicit none

        class(zero_parametrization_2D)  :: self       !< the base class

        self%parameter_number = 0

    end subroutine ZeroParametrized2DSetParamNumber


    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that initializes the parameters based on the values found in an input array.
    subroutine ZeroParametrized2DInitParams( self, array )

        implicit none

        class(zero_parametrization_2D)                         :: self   !< the base class.
        real(dl), dimension(self%parameter_number), intent(in) :: array  !< input array with the values of the parameters.

    end subroutine ZeroParametrized2DInitParams

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the value of the zero function i-th parameter.
    subroutine ZeroParametrized2DParameterValues( self, i, value )

        implicit none

        class(zero_parametrization_2D)  :: self   !< the base class.
        integer , intent(in)            :: i      !< input number of the parameter
        real(dl), intent(out)           :: value  !< output value of the parameter

    end subroutine ZeroParametrized2DParameterValues

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the value of the zero function.
    function ZeroParametrized2DValue( self, x, y, eft_cache )

        implicit none

        class(zero_parametrization_2D)                     :: self      !< the base class
        real(dl), intent(in)                               :: x         !< the input first variable
        real(dl), intent(in)                               :: y         !< the input second variable
        type(EFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: ZeroParametrized2DValue                             !< the output value

        ZeroParametrized2DValue = 0._dl

    end function ZeroParametrized2DValue

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the first partial derivative of the zero function with respect to x.
    function ZeroParametrized2DFirstDerivativeX( self, x, y, eft_cache )

        implicit none

        class(zero_parametrization_2D)                     :: self      !< the base class
        real(dl), intent(in)                               :: x         !< the input first variable
        real(dl), intent(in)                               :: y         !< the input second variable
        type(EFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: ZeroParametrized2DFirstDerivativeX                  !< the output value

        ZeroParametrized2DFirstDerivativeX = 0._dl

    end function ZeroParametrized2DFirstDerivativeX

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the first partial derivative of the zero function with respect to y.
    function ZeroParametrized2DFirstDerivativeY( self, x, y, eft_cache )

        implicit none

        class(zero_parametrization_2D)                     :: self      !< the base class
        real(dl), intent(in)                               :: x         !< the input first variable
        real(dl), intent(in)                               :: y         !< the input second variable
        type(EFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: ZeroParametrized2DFirstDerivativeY                  !< the output value

        ZeroParametrized2DFirstDerivativeY = 0._dl

    end function ZeroParametrized2DFirstDerivativeY

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the second partial derivative of the zero function with respect to x.
    function ZeroParametrized2DSecondDerivativeX( self, x, y, eft_cache )

        implicit none

        class(zero_parametrization_2D)                     :: self      !< the base class
        real(dl), intent(in)                               :: x         !< the input first variable
        real(dl), intent(in)                               :: y         !< the input second variable
        type(EFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: ZeroParametrized2DSecondDerivativeX                 !< the output value

        ZeroParametrized2DSecondDerivativeX = 0._dl

    end function ZeroParametrized2DSecondDerivativeX

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the second partial derivative of the zero function with respect to y.
    function ZeroParametrized2DSecondDerivativeY( self, x, y, eft_cache )

        implicit none

        class(zero_parametrization_2D)                     :: self      !< the base class
        real(dl), intent(in)                               :: x         !< the input first variable
        real(dl), intent(in)                               :: y         !< the input second variable
        type(EFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: ZeroParametrized2DSecondDerivativeY                 !< the output value

        ZeroParametrized2DSecondDerivativeY = 0._dl

    end function ZeroParametrized2DSecondDerivativeY

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the mixed partial derivative of the zero function with respect to x and y.
    function ZeroParametrized2DSecondDerivativeXY( self, x, y, eft_cache )

        implicit none

        class(zero_parametrization_2D)                     :: self      !< the base class
        real(dl), intent(in)                               :: x         !< the input first variable
        real(dl), intent(in)                               :: y         !< the input second variable
        type(EFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: ZeroParametrized2DSecondDerivativeXY                !< the output value

        ZeroParametrized2DSecondDerivativeXY = 0._dl

    end function ZeroParametrized2DSecondDerivativeXY

    !----------------------------------------------------------------------------------------

end module EFTCAMB_neutral_parametrization_2D

!----------------------------------------------------------------------------------------
