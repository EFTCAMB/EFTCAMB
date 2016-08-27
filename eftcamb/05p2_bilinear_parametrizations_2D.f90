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

!> @file 05p2_bilinear_parametrizations_2D.f90
!! This file contains the definition of the bilinear parametrization, inheriting from
!! parametrized_function_2D.


!----------------------------------------------------------------------------------------
!> This module contains the definition of the bilinear parametrization, inheriting from
!! parametrized_function_2D.

!> @author Simone Peirone, Bin Hu, Marco Raveri

module EFTCAMB_bilinear_parametrizations_2D

    use precision
    use IniFile
    use AMLutils
    use EFTDef
    use EFTCAMB_abstract_parametrizations_2D

    implicit none

    private

    public bilinear_parametrization_2D

    ! ---------------------------------------------------------------------------------------------
    !> Type containing the bilinear function parametrization. Inherits from parametrized_function_2D.
    type, extends ( parametrized_function_2D ) :: bilinear_parametrization_2D

        real(dl) :: linear_value_1
        real(dl) :: linear_value_2

    contains

        ! initialization:
        procedure :: init                  => BilinearParametrized2DInitialize          !< subroutine that initializes the bilinear parametrization.
        procedure :: init_from_file        => BilinearParametrized2DInitFromFile        !< subroutine that reads a Ini file looking for the parameters of the function.
        procedure :: init_parameters       => BilinearParametrized2DInit                !< subroutine that initializes the function parameters based on the values found in an input array.
        ! utility functions:
        procedure :: feedback              => BilinearParametrized2DFeedback            !< subroutine that prints to screen the informations about the function.
        procedure :: parameter_names       => BilinearParametrized2DParameterNames      !< subroutine that returns the i-th parameter name of the function.
        procedure :: parameter_names_latex => BilinearParametrized2DParameterNamesLatex !< subroutine that returns the i-th parameter name of the function in latex format.
        ! evaluation procedures:
        procedure :: value                 => BilinearParametrized2DValue               !< function that returns the value of the function.
        procedure :: first_derivative_x    => BilinearParametrized2DFirstDerivativeX    !< function that returns the first partial derivative of the function with respect to x.
        procedure :: first_derivative_y    => BilinearParametrized2DFirstDerivativeY    !< function that returns the first partial derivative of the function with respect to y.
        procedure :: second_derivative_x   => BilinearParametrized2DSecondDerivativeX       !< function that returns the second derivative of the function with respect to x.
        procedure :: second_derivative_y   => BilinearParametrized2DSecondDerivativeY       !< function that returns the second derivative of the function with respect to y.
        procedure :: second_derivative_xy  => BilinearParametrized2DSecondDerivativeXY      !< function that returns the second derivative of the function with respect to x and y.

    end type bilinear_parametrization_2D

contains

    ! ---------------------------------------------------------------------------------------------
    ! Implementation of the bilinear parametrization for 2D functions.
    ! ---------------------------------------------------------------------------------------------

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that initializes the bilinear parametrization
    subroutine BilinearParametrized2DInitialize( self, name, latexname )

        implicit none

        class(bilinear_parametrization_2D) :: self         !< the base class
        character(*), intent(in)           :: name         !< the name of the function
        character(*), intent(in)           :: latexname    !< the latex name of the function

        ! store the name of the function:
        self%name             = TRIM( name )
        ! store the latex name of the function:
        self%name_latex       = TRIM( latexname )
        ! initialize the number of parameters:
        self%parameter_number = 2

    end subroutine BilinearParametrized2DInitialize

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that reads a Ini file looking for the parameters of the function.
    subroutine BilinearParametrized2DInitFromFile( self, Ini )

        implicit none

        class(bilinear_parametrization_2D)  :: self   !< the base class
        type(TIniFile)                      :: Ini    !< Input ini file

        character(len=EFT_names_max_length) :: param_name_1
        character(len=EFT_names_max_length) :: param_name_2

        call self%parameter_names( 1, param_name_1 )
        call self%parameter_names( 2, param_name_2 )

        self%linear_value_1 = Ini_Read_Double_File( Ini, TRIM(param_name_1), 0._dl )
        self%linear_value_2 = Ini_Read_Double_File( Ini, TRIM(param_name_2), 0._dl )

    end subroutine BilinearParametrized2DInitFromFile

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that initializes the function parameters based on the values found in an input array.
    subroutine BilinearParametrized2DInit( self, array )

        implicit none

        class(bilinear_parametrization_2D)                        :: self   !< the base class.
        real(dl), dimension(self%parameter_number), intent(in)    :: array  !< input array with the values of the parameters.

        self%linear_value_1 = array(1)
        self%linear_value_2 = array(2)

    end subroutine BilinearParametrized2DInit

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that prints to screen the informations about the function.
    subroutine BilinearParametrized2DFeedback( self )

        implicit none

        class(bilinear_parametrization_2D)  :: self         !< the base class

        integer                             :: i
        real(dl)                            :: param_value
        character(len=EFT_names_max_length) :: param_name

        write(*,*)     'Bilinear function: ', self%name
        do i=1, self%parameter_number
            call self%parameter_names( i, param_name  )
            call self%parameter_value( i, param_value )
            write(*,*) param_name, '=', param_value
        end do

    end subroutine BilinearParametrized2DFeedback

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the i-th parameter name
    subroutine BilinearParametrized2DParameterNames( self, i, name )

        implicit none

        class(bilinear_parametrization_2D):: self   !< the base class
        integer     , intent(in)          :: i      !< the index of the parameter
        character(*), intent(out)         :: name   !< the output name of the i-th parameter

        select case (i)
            case(1)
                name = TRIM(self%name)//'_1'
            case(2)
                name = TRIM(self%name)//'_2'
            case default
                write(*,*) 'Illegal index for parameter_names.'
                write(*,*) 'Maximum value is:', self%parameter_number
                stop
        end select

    end subroutine BilinearParametrized2DParameterNames

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the latex version of the i-th parameter name
    subroutine BilinearParametrized2DParameterNamesLatex( self, i, latexname )

        implicit none

        class(bilinear_parametrization_2D) :: self        !< the base class
        integer     , intent(in)           :: i           !< The index of the parameter
        character(*), intent(out)          :: latexname   !< the output latex name of the i-th parameter

        select case (i)
        case(1)
            latexname = TRIM(self%name_latex)//'_1'
        case(2)
            latexname = TRIM(self%name_latex)//'_2'
            case default
                write(*,*) 'Illegal index for parameter_names.'
                write(*,*) 'Maximum value is:', self%parameter_number
                stop
        end select

    end subroutine BilinearParametrized2DParameterNamesLatex

    ! ---------------------------------------------------------------------------------------------

    !> Function that returns the value of the bilinear function in the two variables.

    function BilinearParametrized2DValue( self, x, y )

        implicit none

        class(bilinear_parametrization_2D) :: self  !< the base class

        real(dl), intent(in)               :: x     !< the first input variable
        real(dl), intent(in)               :: y     !< the second input variable

        real(dl) :: BilinearParametrized2DValue     !< the output value

        BilinearParametrized2DValue = self%linear_value_1*x + self%linear_value_2*y

    end function BilinearParametrized2DValue

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the value of the first partial derivative, with respect to the

    ! first variable, of the bilinear function.

    function BilinearParametrized2DFirstDerivativeX( self, x, y )

        implicit none

        class(bilinear_parametrization_2D) :: self           !< the base class

        real(dl), intent(in)             :: x                !< the first input variable
        real(dl), intent(in)             :: y                !< the second input variable

        real(dl) :: BilinearParametrized2DFirstDerivativeX   !< the output value

        BilinearParametrized2DFirstDerivativeX = self%linear_value_1

    end function BilinearParametrized2DFirstDerivativeX

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the value of the first partial derivative, with respect to the
    ! second variable, of the bilinear function.

    function BilinearParametrized2DFirstDerivativeY( self, x, y )

        implicit none

        class(bilinear_parametrization_2D) :: self           !< the base class

        real(dl), intent(in)             :: x                !< the first input variable
        real(dl), intent(in)             :: y                !< the second input variable

        real(dl) :: BilinearParametrized2DFirstDerivativeY   !< the output value

        BilinearParametrized2DFirstDerivativeY = self%linear_value_2

    end function BilinearParametrized2DFirstDerivativeY

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the second partial derivative of the function with respect to x.
    function BilinearParametrized2DSecondDerivativeX( self, x, y )

        implicit none

        class(bilinear_parametrization_2D) :: self          !< the base class
        real(dl), intent(in)           :: x                 !< the input first variable
        real(dl), intent(in)           :: y                 !< the input second variable
        real(dl) :: BilinearParametrized2DSecondDerivativeX !< the output value

        BilinearParametrized2DSecondDerivativeX = 0._dl

    end function BilinearParametrized2DSecondDerivativeX

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the second partial derivative of the function with respect to y.
    function BilinearParametrized2DSecondDerivativeY( self, x, y )

        implicit none

        class(bilinear_parametrization_2D) :: self          !< the base class
        real(dl), intent(in)           :: x                 !< the input first variable
        real(dl), intent(in)           :: y                 !< the input second variable
        real(dl) :: BilinearParametrized2DSecondDerivativeY !< the output value

        BilinearParametrized2DSecondDerivativeY = 0._dl

    end function BilinearParametrized2DSecondDerivativeY

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the mixed partial derivative of the function with respect to x and y.
    function BilinearParametrized2DSecondDerivativeXY( self, x, y )

        implicit none

        class(bilinear_parametrization_2D) :: self          !< the base class
        real(dl), intent(in)           :: x                 !< the input first variable
        real(dl), intent(in)           :: y                 !< the input second variable
        real(dl) :: BilinearParametrized2DSecondDerivativeXY!< the output value

        BilinearParametrized2DSecondDerivativeXY = 0._dl

    end function BilinearParametrized2DSecondDerivativeXY

    !----------------------------------------------------------------------------------------

end module EFTCAMB_bilinear_parametrizations_2D

!----------------------------------------------------------------------------------------
