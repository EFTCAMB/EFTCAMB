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

!> @file 03_abstract_parametrizations_1D.f90
!! This file contains the abstract class for generic parametrizations for 1D functions
!! that are used by several models in EFTCAMB. When there is a free function
!! in EFT it should be declared as a class inheriting from parametrized_function_1D.
!! This guarantees maximum performances as well as maximum flexibility.


!----------------------------------------------------------------------------------------
!> This module contains the abstract class for generic parametrizations for 1D functions
!! that are used by several models in EFTCAMB. When there is a free function
!! in EFT it should be declared as a class inheriting from parametrized_function_1D.

!> @author Bin Hu, Marco Raveri

module EFTCAMB_abstract_parametrizations_1D

    use precision
    use IniFile

    implicit none

    public parametrized_function_1D

    !----------------------------------------------------------------------------------------
    !> This is the abstract type for parametrized functions. As a rule, when there is a
    !! free function in EFT it should be declared as a class inheriting from parametrized_function_1D.
    !! This guarantees maximum performances as well as maximum flexibility.
    type, abstract :: parametrized_function_1D

        integer                       :: parameter_number !< number of parameters defining the parametrized function.
        character(len=:), allocatable :: name             !< name of the function.
        character(len=:), allocatable :: name_latex       !< latex name of the function.

    contains

        ! initialization procedures:
        procedure( ParametrizedFunction1DInitialize       ), deferred :: init              !< subroutine that initializes the name and latex name of the function.
        procedure :: init_from_file        => ParametrizedFunction1DInitFromFile           !< subroutine that reads a Ini file looking for the parameters of the function.
        procedure :: init_parameters       => ParametrizedFunction1DInitParams             !< subroutine that initializes the function parameters based on the values found in an input array.
        ! utility functions:
        procedure :: feedback              => ParametrizedFunction1DFeedback               !< subroutine that prints to screen the informations about the function.
        procedure :: parameter_names       => ParametrizedFunction1DParameterNames         !< subroutine that returns the i-th parameter name of the function.
        procedure :: parameter_names_latex => ParametrizedFunction1DParameterNamesLatex    !< subroutine that returns the i-th parameter name of the function in latex format.
        procedure :: parameter_value       => ParametrizedFunction1DParameterValues        !< subroutine that returns the value of the function i-th parameter.
        ! evaluation procedures:
        procedure( ParametrizedFunction1DValue            ), deferred :: value             !< function that returns the value of the function.
        procedure( ParametrizedFunction1DFirstDerivative  ), deferred :: first_derivative  !< function that returns the first derivative of the function.
        procedure( ParametrizedFunction1DSecondDerivative ), deferred :: second_derivative !< function that returns the second derivative of the function.
        procedure( ParametrizedFunction1DThirdDerivative  ), deferred :: third_derivative  !< function that returns the third derivative of the function.
        procedure( ParametrizedFunction1DIntegral         ), deferred :: integral          !< function that returns the strange integral that we need for w_DE.

    end type parametrized_function_1D

    ! ---------------------------------------------------------------------------------------------
    ! parametrized_function_1D abstract interfaces: these are all the function procedures
    ! that the user HAS to override when writing its own parametrized 1D function.
    ! ---------------------------------------------------------------------------------------------

    abstract interface

        ! ---------------------------------------------------------------------------------------------
        !> Subroutine that initializes the name and latex name of the parametrization.
        subroutine ParametrizedFunction1DInitialize( self, name, latexname )
            use precision
            import parametrized_function_1D
            implicit none
            class(parametrized_function_1D)  :: self       !< the base class
            character(*), intent(in)         :: name       !< the name of the function
            character(*), intent(in)         :: latexname  !< the latex name of the function
        end subroutine ParametrizedFunction1DInitialize

        ! ---------------------------------------------------------------------------------------------
        !> Function that returns the value of the function.
        function ParametrizedFunction1DValue( self, x )
            use precision
            import parametrized_function_1D
            implicit none
            class(parametrized_function_1D) :: self  !< the base class
            real(dl), intent(in)            :: x     !< the input scale factor
            real(dl) :: ParametrizedFunction1DValue  !< the output value
        end function ParametrizedFunction1DValue

        ! ---------------------------------------------------------------------------------------------
        !> Function that returns the value of the first derivative of the function
        !! with respect to the scale factor.
        function ParametrizedFunction1DFirstDerivative( self, x )
            use precision
            import parametrized_function_1D
            implicit none
            class(parametrized_function_1D) :: self             !< the base class
            real(dl), intent(in)            :: x                !< the input scale factor
            real(dl) :: ParametrizedFunction1DFirstDerivative   !< the output value
        end function ParametrizedFunction1DFirstDerivative

        ! ---------------------------------------------------------------------------------------------
        !> Function that returns the value of the second derivative of the function
        !! with respect to the scale factor.
        function ParametrizedFunction1DSecondDerivative( self, x )
            use precision
            import parametrized_function_1D
            implicit none
            class(parametrized_function_1D) :: self             !< the base class
            real(dl), intent(in)            :: x                !< the input scale factor
            real(dl) :: ParametrizedFunction1DSecondDerivative  !< the output value
        end function ParametrizedFunction1DSecondDerivative

        ! ---------------------------------------------------------------------------------------------
        !> Function that returns the value of the third derivative of the function
        !! with respect to the scale factor.
        function ParametrizedFunction1DThirdDerivative( self, x )
            use precision
            import parametrized_function_1D
            implicit none
            class(parametrized_function_1D) :: self            !< the base class
            real(dl), intent(in)            :: x               !< the input scale factor
            real(dl) :: ParametrizedFunction1DThirdDerivative  !< the output value
        end function ParametrizedFunction1DThirdDerivative

        ! ---------------------------------------------------------------------------------------------
        !> Function that returns the integral of the function, as defined in the notes.
        function ParametrizedFunction1DIntegral( self, x )
            use precision
            import parametrized_function_1D
            implicit none
            class(parametrized_function_1D) :: self     !< the base class
            real(dl), intent(in)            :: x        !< the scale factor at which the integral is wanted
            real(dl) :: ParametrizedFunction1DIntegral  !< the output value
        end function ParametrizedFunction1DIntegral

        ! ---------------------------------------------------------------------------------------------

    end interface

    ! ---------------------------------------------------------------------------------------------

contains

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that reads a Ini file looking for the parameters of the function.
    subroutine ParametrizedFunction1DInitFromFile( self, Ini )

        implicit none

        class(parametrized_function_1D) :: self   !< the base class
        type(TIniFile)                  :: Ini    !< Input ini file

    end subroutine ParametrizedFunction1DInitFromFile

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that reads a Ini file looking for the parameters of the function.
    subroutine ParametrizedFunction1DInitParams( self, array )

        implicit none

        class(parametrized_function_1D)                        :: self   !< the base class.
        real(dl), dimension(self%parameter_number), intent(in) :: array  !< input array with the values of the parameters.

    end subroutine ParametrizedFunction1DInitParams

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that prints to screen the informations about the function.
    subroutine ParametrizedFunction1DFeedback( self )

        implicit none

        class(parametrized_function_1D)  :: self   !< the base class

    end subroutine ParametrizedFunction1DFeedback

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the i-th parameter name.
    subroutine ParametrizedFunction1DParameterNames( self, i, name )

        implicit none

        class(parametrized_function_1D) :: self   !< the base class
        integer     , intent(in)        :: i      !< the index of the parameter
        character(*), intent(out)       :: name   !< the output name of the i-th parameter

    end subroutine ParametrizedFunction1DParameterNames

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the latex version of the i-th parameter name.
    subroutine ParametrizedFunction1DParameterNamesLatex( self, i, latexname )

        implicit none

        class(parametrized_function_1D) :: self       !< the base class
        integer     , intent(in)        :: i          !< the index of the parameter
        character(*), intent(out)       :: latexname  !< the output latex name of the i-th parameter

    end subroutine ParametrizedFunction1DParameterNamesLatex

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the value of the function i-th parameter.
    subroutine ParametrizedFunction1DParameterValues( self, i, value )

        implicit none

        class(parametrized_function_1D) :: self       !< the base class
        integer  , intent(in)           :: i          !< the index of the parameter
        real(dl) , intent(out)          :: value      !< the output value of the i-th parameter

    end subroutine ParametrizedFunction1DParameterValues

    ! ---------------------------------------------------------------------------------------------

end module EFTCAMB_abstract_parametrizations_1D

!----------------------------------------------------------------------------------------
