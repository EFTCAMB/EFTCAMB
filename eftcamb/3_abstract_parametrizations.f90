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

!> @file 3_abstract_parametrizations.f90
!! This file contains the abstract class for generic parametrizations for functions
!! that are used by several models in EFTCAMB. As a rule, when there is a free function
!! in EFT it should be declared as a class inheriting from parametrized_function.
!! This guarantees maximum performances as well as maximum flexibility.


!----------------------------------------------------------------------------------------
!> This module contains the abstract class for generic parametrizations for functions
!! that are used by several models in EFTCAMB. As a rule, when there is a free function
!! in EFT it should be declared as a class inheriting from parametrized_function

!> @author Bin Hu, Marco Raveri

module EFTCAMB_abstract_parametrizations

    use precision
    use IniFile

    implicit none

    public parametrized_function

    !----------------------------------------------------------------------------------------
    !> This is the abstract type for parametrized functions. As a rule, when there is a
    !! free function in EFT it should be declared as a class inheriting from parametrized_function.
    !! This guarantees maximum performances as well as maximum flexibility.
    type parametrized_function

        integer                       :: parameter_number !< number of parameters defining the parametrized function
        character(len=:), allocatable :: name             !< name of the function
        character(len=:), allocatable :: name_latex       !< latex name of the function

    contains

        ! initialization procedures:
        procedure :: init                  => ParametrizedFunctionInitialize          !< subroutine that initializes the name and latex name of the function.
        procedure :: init_from_file        => ParametrizedFunctionInitFromFile        !< subroutine that reads a Ini file looking for the parameters of the function.
        procedure :: init_parameters       => ParametrizedFunctionInit                !< subroutine that initializes the function parameters based on the values found in an input array.
        ! utility functions:
        procedure :: feedback              => ParametrizedFunctionFeedback            !< subroutine that prints to screen the informations about the function.
        procedure :: parameter_names       => ParametrizedFunctionParameterNames      !< subroutine that returns the i-th parameter name of the function.
        procedure :: parameter_names_latex => ParametrizedFunctionParameterNamesLatex !< subroutine that returns the i-th parameter name of the function in latex format.
        procedure :: parameter_value       => ParametrizedFunctionParameterValues     !< subroutine that returns the value of the function i-th parameter.
        ! evaluation procedures:
        procedure :: value                 => ParametrizedFunctionValue               !< function that returns the value of the function.
        procedure :: first_derivative      => ParametrizedFunctionFirstDerivative     !< function that returns the first derivative of the function.
        procedure :: second_derivative     => ParametrizedFunctionSecondDerivative    !< function that returns the second derivative of the function.
        procedure :: third_derivative      => ParametrizedFunctionThirdDerivative     !< function that returns the third derivative of the function.
        procedure :: integral              => ParametrizedFunctionIntegral            !< function that returns the strange integral that we need for w_DE.

    end type parametrized_function

    ! ---------------------------------------------------------------------------------------------

contains

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that initializes the name and latex name of the parametrization.
    subroutine ParametrizedFunctionInitialize( self, name, latexname )

        implicit none

        class(parametrized_function)    :: self       !< the base class
        character(*), intent(in)        :: name       !< the name of the function
        character(*), intent(in)        :: latexname  !< the latex name of the function

    end subroutine ParametrizedFunctionInitialize

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that reads a Ini file looking for the parameters of the function.
    subroutine ParametrizedFunctionInitFromFile( self, Ini )

        implicit none

        class(parametrized_function) :: self   !< the base class
        type(TIniFile)               :: Ini    !< Input ini file

    end subroutine ParametrizedFunctionInitFromFile

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that reads a Ini file looking for the parameters of the function.
    subroutine ParametrizedFunctionInit( self, array )

        implicit none

        class(parametrized_function)                           :: self   !< the base class.
        real(dl), dimension(self%parameter_number), intent(in) :: array  !< input array with the values of the parameters.

    end subroutine ParametrizedFunctionInit

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that prints to screen the informations about the function
    subroutine ParametrizedFunctionFeedback( self )

        implicit none

        class(parametrized_function)  :: self   !< the base class

    end subroutine ParametrizedFunctionFeedback

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the i-th parameter name
    subroutine ParametrizedFunctionParameterNames( self, i, name )

        implicit none

        class(parametrized_function)  :: self   !< the base class
        integer     , intent(in)      :: i      !< the index of the parameter
        character(*), intent(out)     :: name   !< the output name of the i-th parameter

    end subroutine ParametrizedFunctionParameterNames

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the latex version of the i-th parameter name
    subroutine ParametrizedFunctionParameterNamesLatex( self, i, latexname )

        implicit none

        class(parametrized_function) :: self       !< the base class
        integer     , intent(in)     :: i          !< the index of the parameter
        character(*), intent(out)    :: latexname  !< the output latex name of the i-th parameter

    end subroutine ParametrizedFunctionParameterNamesLatex

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the value of the function i-th parameter.
    subroutine ParametrizedFunctionParameterValues( self, i, value )

        implicit none

        class(parametrized_function) :: self       !< the base class
        integer  , intent(in)        :: i          !< the index of the parameter
        real(dl) , intent(out)       :: value      !< the output value of the i-th parameter

    end subroutine ParametrizedFunctionParameterValues

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the value of the function.
    function ParametrizedFunctionValue( self, a )

        implicit none

        class(parametrized_function) :: self  !< the base class
        real(dl), intent(in)         :: a     !< the input scale factor
        real(dl) :: ParametrizedFunctionValue !< the output value

        ParametrizedFunctionValue = 0._dl

    end function ParametrizedFunctionValue

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the value of the first derivative of the function
    !! with respect to the scale factor.
    function ParametrizedFunctionFirstDerivative( self, a )

        implicit none

        class(parametrized_function) :: self             !< the base class
        real(dl), intent(in)         :: a                !< the input scale factor
        real(dl) :: ParametrizedFunctionFirstDerivative  !< the output value

        ParametrizedFunctionFirstDerivative = 0._dl

    end function ParametrizedFunctionFirstDerivative

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the value of the second derivative of the function
    !! with respect to the scale factor.
    function ParametrizedFunctionSecondDerivative( self, a )

        implicit none

        class(parametrized_function) :: self             !< the base class
        real(dl), intent(in)         :: a                !< the input scale factor
        real(dl) :: ParametrizedFunctionSecondDerivative !< the output value

        ParametrizedFunctionSecondDerivative = 0._dl

    end function ParametrizedFunctionSecondDerivative

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the value of the third derivative of the function
    !! with respect to the scale factor.
    function ParametrizedFunctionThirdDerivative( self, a )

        implicit none

        class(parametrized_function) :: self            !< the base class
        real(dl), intent(in)         :: a               !< the input scale factor
        real(dl) :: ParametrizedFunctionThirdDerivative !< the output value

        ParametrizedFunctionThirdDerivative = 0._dl

    end function ParametrizedFunctionThirdDerivative

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the integral of the function, as defined in the notes
    function ParametrizedFunctionIntegral( self, a )

        implicit none

        class(parametrized_function) :: self     !< the base class
        real(dl), intent(in)         :: a        !< the scale factor at which the integral is wanted
        real(dl) :: ParametrizedFunctionIntegral !< the output value

        ParametrizedFunctionIntegral = 0._dl

    end function ParametrizedFunctionIntegral

    ! ---------------------------------------------------------------------------------------------

end module EFTCAMB_abstract_parametrizations

!----------------------------------------------------------------------------------------
