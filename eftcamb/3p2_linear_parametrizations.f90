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

!> @file 3p2_linear_parametrizations.f90
!! This file contains the definition of the linear parametrization, inheriting from
!! parametrized_function.


!----------------------------------------------------------------------------------------
!> This module contains the definition of the linear parametrization, inheriting from
!! parametrized_function.

!> @author Bin Hu, Marco Raveri

module EFTCAMB_linear_parametrizations

    use precision
    use EFTDef
    use EFTCAMB_abstract_parametrizations

    implicit none

    ! ---------------------------------------------------------------------------------------------
    !> Type containing the linear function parametrization. Inherits from parametrized_function.
    !! Notice that the derivatives above the first are not overridden since they are zero identically.
    type, extends ( parametrized_function ) :: linear_parametrization

        real(dl) :: linear_value

    contains

        ! initialization:
        procedure :: init                  => LinearParametrizedInitialize          !< subroutine that initializes the constant parametrization
        procedure :: init_from_file        => LinearParametrizedInitFromFile        !< subroutine that reads a Ini file looking for the parameters of the function.
        procedure :: init_parameters       => LinearParametrizedInit                !< subroutine that initializes the function parameters based on the values found in an input array.
        ! utility functions:
        procedure :: feedback              => LinearParametrizedFeedback            !< subroutine that prints to screen the informations about the function.
        procedure :: parameter_names       => LinearParametrizedParameterNames      !< subroutine that returns the i-th parameter name of the function
        procedure :: parameter_names_latex => LinearParametrizedParameterNamesLatex !< subroutine that returns the i-th parameter name of the function in latex format
        ! evaluation procedures:
        procedure :: value                 => LinearParametrizedValue               !< function that returns the value of the function
        procedure :: first_derivative      => LinearParametrizedFirstDerivative     !< function that returns the first derivative of the function
        procedure :: integral              => LinearParametrizedIntegral            !< function that returns the strange integral that we need for w_DE

    end type linear_parametrization

contains

    ! ---------------------------------------------------------------------------------------------
    ! Implementation of the linear function.
    ! ---------------------------------------------------------------------------------------------

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that initializes the linear parametrization
    subroutine LinearParametrizedInitialize( self, name, latexname )

        implicit none

        class(linear_parametrization)   :: self         !< the base class
        character(*), intent(in)        :: name         !< the name of the function
        character(*), intent(in)        :: latexname    !< the latex name of the function

        ! store the name of the function:
        self%name             = TRIM( name )
        ! store the latex name of the function:
        self%name_latex       = TRIM( latexname )
        ! initialize the number of parameters:
        self%parameter_number = 1

    end subroutine LinearParametrizedInitialize

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that reads a Ini file looking for the parameters of the function.
    subroutine LinearParametrizedInitFromFile( self, Ini )

        implicit none

        class(linear_parametrization)   :: self   !< the base class
        type(TIniFile)                  :: Ini    !< Input ini file

        self%linear_value = Ini_Read_Double_File( Ini, TRIM(self%name)//'_0', 0._dl )

    end subroutine LinearParametrizedInitFromFile

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that initializes the function parameters based on the values found in an input array.
    subroutine LinearParametrizedInit( self, array )

        implicit none

        class(linear_parametrization)                           :: self   !< the base class.
        real(dl), dimension(self%parameter_number), intent(in)  :: array  !< input array with the values of the parameters.

        self%linear_value = array(1)

    end subroutine LinearParametrizedInit

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that prints to screen the informations about the function.
    subroutine LinearParametrizedFeedback( self )

        implicit none

        class(linear_parametrization)   :: self         !< the base class

        integer                             :: i
        real(dl)                            :: param_value
        character(len=EFT_names_max_length) :: param_name

        write(*,*)     'Linear function: ', self%name
        do i=1, self%parameter_number
            call self%parameter_names( i, param_name  )
            call self%parameter_value( i, param_value )
            write(*,*) param_name, '=', param_value
        end do

    end subroutine LinearParametrizedFeedback

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the i-th parameter name
    subroutine LinearParametrizedParameterNames( self, i, name )

        implicit none

        class(linear_parametrization)   :: self   !< the base class
        integer     , intent(in)        :: i      !< the index of the parameter
        character(*), intent(out)       :: name   !< the output name of the i-th parameter

        select case (i)
            case(1)
                name = TRIM(self%name)//'0'
            case default
                write(*,*) 'Illegal index for parameter_names.'
                write(*,*) 'Maximum value is:', self%parameter_number
                stop
        end select

    end subroutine LinearParametrizedParameterNames

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the latex version of the i-th parameter name
    subroutine LinearParametrizedParameterNamesLatex( self, i, latexname )

        implicit none

        class(linear_parametrization)   :: self        !< the base class
        integer     , intent(in)        :: i           !< The index of the parameter
        character(*), intent(out)       :: latexname   !< the output latex name of the i-th parameter

        select case (i)
            case(1)
                latexname = TRIM(self%name_latex)//'_0'
            case default
                write(*,*) 'Illegal index for parameter_names.'
                write(*,*) 'Maximum value is:', self%parameter_number
                stop
        end select

    end subroutine LinearParametrizedParameterNamesLatex

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the value of the linear function.
    function LinearParametrizedValue( self, a )

        implicit none

        class(linear_parametrization)   :: self  !< the base class
        real(dl), intent(in)            :: a     !< the input scale factor
        real(dl) :: LinearParametrizedValue      !< the output value

        LinearParametrizedValue = self%linear_value*a

    end function LinearParametrizedValue

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the value of the linear function.
    function LinearParametrizedFirstDerivative( self, a )

        implicit none

        class(linear_parametrization)   :: self          !< the base class
        real(dl), intent(in)            :: a             !< the input scale factor
        real(dl) :: LinearParametrizedFirstDerivative    !< the output value

        LinearParametrizedFirstDerivative = self%linear_value

    end function LinearParametrizedFirstDerivative

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the integral of the linear function, as defined in the notes
    function LinearParametrizedIntegral( self, a )

        implicit none

        class(linear_parametrization)   :: self     !< the base class
        real(dl), intent(in)            :: a        !< the scale factor at which the integral is wanted
        real(dl) :: LinearParametrizedIntegral      !< the output value

        LinearParametrizedIntegral = Exp(-3._dl*(a-1._dl)*self%linear_value)/a

    end function LinearParametrizedIntegral

    ! ---------------------------------------------------------------------------------------------

end module EFTCAMB_linear_parametrizations

!----------------------------------------------------------------------------------------
