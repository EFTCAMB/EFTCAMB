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

!> @file 04_abstract_parametrizations_2D.f90
!! This file contains the abstract class for generic parametrizations for 2D functions
!! that are used by several models in EFTCAMB. As a rule, when there is a free function
!! in EFT it should be declared as a class inheriting from parametrized_function.
!! This guarantees maximum performances as well as maximum flexibility.


!----------------------------------------------------------------------------------------
!> This module contains the abstract class for generic parametrizations for 2D functions
!! that are used by several models in EFTCAMB. As a rule, when there is a free function
!! in EFT it should be declared as a class inheriting from parametrized_function

!> @author Simone Peirone, Bin Hu, Marco Raveri

module EFTCAMB_abstract_parametrizations_2D

    use precision
    use IniFile

    implicit none

    public parametrized_function_2D

    !----------------------------------------------------------------------------------------
    !> This is the abstract type for parametrized functions. As a rule, when there is a
    !! free function in EFT it should be declared as a class inheriting from parametrized_function.
    !! This guarantees maximum performances as well as maximum flexibility.
    type parametrized_function_2D

        integer                       :: parameter_number !< number of parameters defining the parametrized function
        character(len=:), allocatable :: name             !< name of the function
        character(len=:), allocatable :: name_latex       !< latex name of the function

    contains

        ! initialization procedures:
        procedure :: init                  => ParametrizedFunction2DInitialize          !< subroutine that initializes the name and latex name of the function.
        procedure :: init_from_file        => ParametrizedFunction2DInitFromFile        !< subroutine that reads a Ini file looking for the parameters of the function.
        procedure :: init_parameters       => ParametrizedFunction2DInit                !< subroutine that initializes the function parameters based on the values found in an input array.
        ! utility functions:
        procedure :: feedback              => ParametrizedFunction2DFeedback            !< subroutine that prints to screen the informations about the function.
        procedure :: parameter_names       => ParametrizedFunction2DParameterNames      !< subroutine that returns the i-th parameter name of the function.
        procedure :: parameter_names_latex => ParametrizedFunction2DParameterNamesLatex !< subroutine that returns the i-th parameter name of the function in latex format.
        procedure :: parameter_value       => ParametrizedFunction2DParameterValues     !< subroutine that returns the value of the function i-th parameter.
        ! evaluation procedures:
        procedure :: value                 => ParametrizedFunction2DValue               !< function that returns the value of the function.
        procedure :: first_derivative_x    => ParametrizedFunction2DFirstDerivativeX    !< function that returns the first partial derivative of the function with respect to x.
        procedure :: first_derivative_y    => ParametrizedFunction2DFirstDerivativeY    !< function that returns the first partial derivative of the function with respect to y.
        procedure :: second_derivative_x   => ParametrizedFunction2DSecondDerivativeX   !< function that returns the second partial derivative of the function with respect to x.
        procedure :: second_derivative_y   => ParametrizedFunction2DSecondDerivativeY   !< function that returns the second partial derivative of the function with respect to y.
        procedure :: second_derivative_xy  => ParametrizedFunction2DSecondDerivativeXY  !< function that returns the mixed partial derivative of the function with respect to x and y.

    end type parametrized_function_2D

    ! ---------------------------------------------------------------------------------------------

contains

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that initializes the name and latex name of the parametrization.
    subroutine ParametrizedFunction2DInitialize( self, name, latexname )

        implicit none

        class(parametrized_function_2D)  :: self       !< the base class
        character(*), intent(in)         :: name       !< the name of the function
        character(*), intent(in)         :: latexname  !< the latex name of the function

    end subroutine ParametrizedFunction2DInitialize

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that reads a Ini file looking for the parameters of the function.
    subroutine ParametrizedFunction2DInitFromFile( self, Ini )

        implicit none

        class(parametrized_function_2D) :: self   !< the base class
        type(TIniFile)                  :: Ini    !< Input ini file

    end subroutine ParametrizedFunction2DInitFromFile

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that reads a Ini file looking for the parameters of the function.
    subroutine ParametrizedFunction2DInit( self, array )

        implicit none

        class(parametrized_function_2D)                        :: self   !< the base class.
        real(dl), dimension(self%parameter_number), intent(in) :: array  !< input array with the values of the parameters.

    end subroutine ParametrizedFunction2DInit

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that prints to screen the informations about the function
    subroutine ParametrizedFunction2DFeedback( self )

        implicit none

        class(parametrized_function_2D)  :: self   !< the base class

    end subroutine ParametrizedFunction2DFeedback

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the i-th parameter name
    subroutine ParametrizedFunction2DParameterNames( self, i, name )

        implicit none

        class(parametrized_function_2D) :: self   !< the base class
        integer     , intent(in)        :: i      !< the index of the parameter
        character(*), intent(out)       :: name   !< the output name of the i-th parameter

    end subroutine ParametrizedFunction2DParameterNames

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the latex version of the i-th parameter name
    subroutine ParametrizedFunction2DParameterNamesLatex( self, i, latexname )

        implicit none

        class(parametrized_function_2D) :: self       !< the base class
        integer     , intent(in)        :: i          !< the index of the parameter
        character(*), intent(out)       :: latexname  !< the output latex name of the i-th parameter

    end subroutine ParametrizedFunction2DParameterNamesLatex

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the value of the function i-th parameter.
    subroutine ParametrizedFunction2DParameterValues( self, i, value )

        implicit none

        class(parametrized_function_2D) :: self       !< the base class
        integer  , intent(in)           :: i          !< the index of the parameter
        real(dl) , intent(out)          :: value      !< the output value of the i-th parameter

    end subroutine ParametrizedFunction2DParameterValues

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the value of the function.
    function ParametrizedFunction2DValue( self, x, y )

        implicit none

        class(parametrized_function_2D) :: self  !< the base class
        real(dl), intent(in)            :: x     !< the first input scale factor
        real(dl), intent(in)            :: y     !< the second input scale factor
        real(dl) :: ParametrizedFunction2DValue  !< the output value

        ParametrizedFunction2DValue = 0._dl

    end function ParametrizedFunction2DValue

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the value of the first partial derivative of the function
    !! with respect to the first scale factor.
    function ParametrizedFunction2DFirstDerivativeX( self, x, y )

        implicit none

        class(parametrized_function_2D) :: self             !< the base class
        real(dl), intent(in)            :: x                !< the first input scale factor
        real(dl), intent(in)            :: y                !< the second input scale factor
        real(dl) :: ParametrizedFunction2DFirstDerivativeX  !< the output value

        ParametrizedFunction2DFirstDerivativeX = 0._dl

    end function ParametrizedFunction2DFirstDerivativeX

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the value of the first partial derivative of the function
    !! with respect to the second scale factor.
    function ParametrizedFunction2DFirstDerivativeY( self, x, y )

        implicit none

        class(parametrized_function_2D) :: self             !< the base class
        real(dl), intent(in)            :: x                !< the first input scale factor
        real(dl), intent(in)            :: y                !< the second input scale factor
        real(dl) :: ParametrizedFunction2DFirstDerivativeY  !< the output value

        ParametrizedFunction2DFirstDerivativeY = 0._dl

    end function ParametrizedFunction2DFirstDerivativeY


    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the value of the second partial derivative of the function
    !! with respect to the first scale factor.
    function ParametrizedFunction2DSecondDerivativeX( self, x, y )

        implicit none

        class(parametrized_function_2D) :: self             !< the base class
        real(dl), intent(in)            :: x                !< the first input scale factor
        real(dl), intent(in)            :: y                !< the second input scale factor
        real(dl) :: ParametrizedFunction2DSecondDerivativeX !< the output value

        ParametrizedFunction2DSecondDerivativeX = 0._dl

    end function ParametrizedFunction2DSecondDerivativeX

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the value of the second partial derivative of the function
    !! with respect to the second scale factor.
    function ParametrizedFunction2DSecondDerivativeY( self, x, y )

        implicit none

        class(parametrized_function_2D) :: self             !< the base class
        real(dl), intent(in)            :: x                !< the first input scale factor
        real(dl), intent(in)            :: y                !< the second input scale factor
        real(dl) :: ParametrizedFunction2DSecondDerivativeY !< the output value

        ParametrizedFunction2DSecondDerivativeY = 0._dl

    end function ParametrizedFunction2DSecondDerivativeY

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the value of the mixed partial derivative of the function
    !! with respect to the first and second scale factors.
    function ParametrizedFunction2DSecondDerivativeXY( self, x, y )

        implicit none

        class(parametrized_function_2D) :: self             !< the base class
        real(dl), intent(in)            :: x                !< the first input scale factor
        real(dl), intent(in)            :: y                !< the second input scale factor
        real(dl) :: ParametrizedFunction2DSecondDerivativeXY!< the output value

        ParametrizedFunction2DSecondDerivativeXY = 0._dl

    end function ParametrizedFunction2DSecondDerivativeXY


   ! ---------------------------------------------------------------------------------------------

end module EFTCAMB_abstract_parametrizations_2D

!----------------------------------------------------------------------------------------
