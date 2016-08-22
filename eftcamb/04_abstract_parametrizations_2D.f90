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
!! that are used by several models in EFTCAMB. When there is a free function
!! in EFT it should be declared as a class inheriting from parametrized_function_2D.
!! This guarantees maximum performances as well as maximum flexibility.


!----------------------------------------------------------------------------------------
!> This module contains the abstract class for generic parametrizations for 2D functions
!! that are used by several models in EFTCAMB. When there is a free function
!! in EFT it should be declared as a class inheriting from parametrized_function_2D.

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
    type, abstract :: parametrized_function_2D

        integer                       :: parameter_number !< number of parameters defining the parametrized function
        character(len=:), allocatable :: name             !< name of the function
        character(len=:), allocatable :: name_latex       !< latex name of the function

    contains

        ! initialization procedures:
        procedure( ParametrizedFunction2DInitialize       ), deferred :: init           !< subroutine that initializes the name and latex name of the function.
        procedure :: init_from_file        => ParametrizedFunction2DInitFromFile        !< subroutine that reads a Ini file looking for the parameters of the function.
        procedure :: init_parameters       => ParametrizedFunction2DInitParams          !< subroutine that initializes the function parameters based on the values found in an input array.
        ! utility functions:
        procedure :: feedback              => ParametrizedFunction2DFeedback            !< subroutine that prints to screen the informations about the function.
        procedure :: parameter_names       => ParametrizedFunction2DParameterNames      !< subroutine that returns the i-th parameter name of the function.
        procedure :: parameter_names_latex => ParametrizedFunction2DParameterNamesLatex !< subroutine that returns the i-th parameter name of the function in latex format.
        procedure :: parameter_value       => ParametrizedFunction2DParameterValues     !< subroutine that returns the value of the function i-th parameter.
        ! evaluation procedures:
        procedure( ParametrizedFunction2DValue              ), deferred :: value               !< function that returns the value of the function.
        procedure( ParametrizedFunction2DFirstDerivativeX   ), deferred :: first_derivative_x  !< function that returns the first partial derivative of the function with respect to x.
        procedure( ParametrizedFunction2DFirstDerivativeY   ), deferred :: first_derivative_y  !< function that returns the first partial derivative of the function with respect to y.
        procedure( ParametrizedFunction2DSecondDerivativeX  ), deferred :: second_derivative_x !< function that returns the second partial derivative of the function with respect to x.
        procedure( ParametrizedFunction2DSecondDerivativeY  ), deferred :: second_derivative_y !< function that returns the second partial derivative of the function with respect to y.
        procedure( ParametrizedFunction2DSecondDerivativeXY ), deferred :: second_derivative_xy!< function that returns the mixed partial derivative of the function with respect to x and y.

    end type parametrized_function_2D

    ! ---------------------------------------------------------------------------------------------
    ! parametrized_function_2D abstract interfaces: these are all the function procedures
    ! that the user HAS to override when writing its own parametrized 2D function.
    ! ---------------------------------------------------------------------------------------------

    abstract interface

        ! ---------------------------------------------------------------------------------------------
        !> Subroutine that initializes the name and latex name of the parametrization.
        subroutine ParametrizedFunction2DInitialize( self, name, latexname )
            use precision
            import parametrized_function_2D
            implicit none
            class(parametrized_function_2D)  :: self       !< the base class
            character(*), intent(in)         :: name       !< the name of the function
            character(*), intent(in)         :: latexname  !< the latex name of the function
        end subroutine ParametrizedFunction2DInitialize

        ! ---------------------------------------------------------------------------------------------
        !> Function that returns the value of the function.
        function ParametrizedFunction2DValue( self, x, y )
            use precision
            import parametrized_function_2D
            implicit none
            class(parametrized_function_2D) :: self  !< the base class
            real(dl), intent(in)            :: x     !< the input first variable
            real(dl), intent(in)            :: y     !< the input second variable
            real(dl) :: ParametrizedFunction2DValue  !< the output value
        end function ParametrizedFunction2DValue

        ! ---------------------------------------------------------------------------------------------
        !> Function that returns the value of the first partial derivative of the function
        !! with respect to the first variable
        function ParametrizedFunction2DFirstDerivativeX( self, x, y )
            use precision
            import parametrized_function_2D
            implicit none
            class(parametrized_function_2D) :: self             !< the base class
            real(dl), intent(in)            :: x                !< the input first variable
            real(dl), intent(in)            :: y                !< the input second variable
            real(dl) :: ParametrizedFunction2DFirstDerivativeX  !< the output value
        end function ParametrizedFunction2DFirstDerivativeX

        ! ---------------------------------------------------------------------------------------------
        !> Function that returns the value of the first partial derivative of the function
        !! with respect to the second variable
        function ParametrizedFunction2DFirstDerivativeY( self, x, y )
            use precision
            import parametrized_function_2D
            implicit none
            class(parametrized_function_2D) :: self             !< the base class
            real(dl), intent(in)            :: x                !< the input first variable
            real(dl), intent(in)            :: y                !< the input second variable
            real(dl) :: ParametrizedFunction2DFirstDerivativeY  !< the output value
        end function ParametrizedFunction2DFirstDerivativeY

        ! ---------------------------------------------------------------------------------------------
        !> Function that returns the value of the second partial derivative of the function
        !! with respect to the first variable.
        function ParametrizedFunction2DSecondDerivativeX( self, x, y )
            use precision
            import parametrized_function_2D
            implicit none
            class(parametrized_function_2D) :: self             !< the base class
            real(dl), intent(in)            :: x                !< the input first variable
            real(dl), intent(in)            :: y                !< the input second variable
            real(dl) :: ParametrizedFunction2DSecondDerivativeX !< the output value
        end function ParametrizedFunction2DSecondDerivativeX

        ! ---------------------------------------------------------------------------------------------
        !> Function that returns the value of the second partial derivative of the function
        !! with respect to the second variable.
        function ParametrizedFunction2DSecondDerivativeY( self, x, y )
            use precision
            import parametrized_function_2D
            implicit none
            class(parametrized_function_2D) :: self             !< the base class
            real(dl), intent(in)            :: x                !< the input first variable
            real(dl), intent(in)            :: y                !< the input second variable
            real(dl) :: ParametrizedFunction2DSecondDerivativeY !< the output value
        end function ParametrizedFunction2DSecondDerivativeY

        ! ---------------------------------------------------------------------------------------------
        !> Function that returns the value of the mixed partial derivative of the function
        !! with respect to the first and second variable.
        function ParametrizedFunction2DSecondDerivativeXY( self, x, y )
            use precision
            import parametrized_function_2D
            implicit none
            class(parametrized_function_2D) :: self             !< the base class
            real(dl), intent(in)            :: x                !< the input first variable
            real(dl), intent(in)            :: y                !< the input second variable
            real(dl) :: ParametrizedFunction2DSecondDerivativeXY!< the output value
        end function ParametrizedFunction2DSecondDerivativeXY

       ! ---------------------------------------------------------------------------------------------

  end interface

    ! ---------------------------------------------------------------------------------------------

contains

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that reads a Ini file looking for the parameters of the function.
    subroutine ParametrizedFunction2DInitFromFile( self, Ini )

        implicit none

        class(parametrized_function_2D) :: self   !< the base class
        type(TIniFile)                  :: Ini    !< Input ini file

    end subroutine ParametrizedFunction2DInitFromFile

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that reads a Ini file looking for the parameters of the function.
    subroutine ParametrizedFunction2DInitParams( self, array )

        implicit none

        class(parametrized_function_2D)                        :: self   !< the base class.
        real(dl), dimension(self%parameter_number), intent(in) :: array  !< input array with the values of the parameters.

    end subroutine ParametrizedFunction2DInitParams

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

end module EFTCAMB_abstract_parametrizations_2D

!----------------------------------------------------------------------------------------
