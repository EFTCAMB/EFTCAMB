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

!> @file 3p1_constant_parametrizations.f90
!! This file contains the definition of the constant parametrization, inheriting from
!! parametrized_function.


!----------------------------------------------------------------------------------------
!> This module contains the definition of the constant parametrization, inheriting from
!! parametrized_function.

!> @author Bin Hu, Marco Raveri

module EFTCAMB_constant_parametrization

    use precision
    use EFTDef
    use EFTCAMB_abstract_parametrizations

    implicit none

    ! ---------------------------------------------------------------------------------------------
    !> Type containing the constant function parametrization. Inherits from parametrized_function.
    !! Notice that the derivatives are not overridden since they are zero identically.
    type, extends ( parametrized_function ) :: constant_parametrization

        real(dl) :: constant_value

    contains

        ! utility functions:
        procedure :: feedback              => ConstantParametrizedFeedback            !< subroutine that prints to screen the informations about the function.
        procedure :: parameter_names       => ConstantParametrizedParameterNames      !< subroutine that returns the i-th parameter name of the function
        procedure :: parameter_names_latex => ConstantParametrizedParameterNamesLatex !< subroutine that returns the i-th parameter name of the function in latex format
        procedure :: parameter_value       => ConstantParametrizedParameterValues     !< subroutine that returns the value of the function i-th parameter.
        ! initialization:
        procedure :: init                  => ConstantParametrizedInitialize          !< subroutine that initializes the constant parametrization
        procedure :: init_from_file        => ConstantParametrizedInitFromFile        !< subroutine that reads a Ini file looking for the parameters of the function.
        procedure :: init_parameters       => ConstantParametrizedInit                !< subroutine that initializes the function parameters based on the values found in an input array.
        ! evaluation procedures:
        procedure :: value                 => ConstantParametrizedValue               !< function that returns the value of the function
        procedure :: integral              => ConstantParametrizedIntegral            !< function that returns the strange integral that we need for w_DE

    end type constant_parametrization

contains

    ! ---------------------------------------------------------------------------------------------
    ! Implementation of the constant function.
    ! ---------------------------------------------------------------------------------------------

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that prints to screen the informations about the function.
    subroutine ConstantParametrizedFeedback( self )

        implicit none

        class(constant_parametrization)     :: self       !< the base class

        integer                             :: i
        real(dl)                            :: param_value
        character(len=EFT_names_max_length) :: param_name

        write(*,*)     'Constant function: ', self%name
        do i=1, self%parameter_number
            call self%parameter_names( i, param_name  )
            call self%parameter_value( i, param_value )
            write(*,*) param_name, '=', param_value
        end do

    end subroutine ConstantParametrizedFeedback

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the i-th parameter name
    subroutine ConstantParametrizedParameterNames( self, i, name )

        implicit none

        class(constant_parametrization) :: self   !< the base class
        integer     , intent(in)        :: i      !< the index of the parameter
        character(*), intent(out)       :: name   !< the output name of the i-th parameter

        select case (i)
            case(1)
                name = TRIM(self%name)//'_0'
            case default
                write(*,*) 'Illegal index for parameter_names.'
                write(*,*) 'Maximum value is:', self%parameter_number
                stop
        end select

    end subroutine ConstantParametrizedParameterNames

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the latex version of the i-th parameter name
    subroutine ConstantParametrizedParameterNamesLatex( self, i, latexname )

        implicit none

        class(constant_parametrization) :: self        !< the base class
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

    end subroutine ConstantParametrizedParameterNamesLatex


    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the value of the function i-th parameter.
    subroutine ConstantParametrizedParameterValues( self, i, value )

        implicit none

        class(constant_parametrization) :: self        !< the base class
        integer     , intent(in)        :: i           !< The index of the parameter
        real(dl)    , intent(out)       :: value       !< the output value of the i-th parameter

        select case (i)
            case(1)
                value = self%constant_value
            case default
                write(*,*) 'Illegal index for parameter_names.'
                write(*,*) 'Maximum value is:', self%parameter_number
                stop
        end select

    end subroutine ConstantParametrizedParameterValues

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that initializes the constant parametrization
    subroutine ConstantParametrizedInitialize( self, name, latexname )

        implicit none

        class(constant_parametrization) :: self       !< the base class
        character(*), intent(in)        :: name       !< the name of the function
        character(*), intent(in)        :: latexname  !< the latex name of the function

        ! store the name of the function:
        self%name             = TRIM( name )
        ! store the latex name of the function:
        self%name_latex       = TRIM( latexname )
        ! initialize the number of parameters:
        self%parameter_number = 1

    end subroutine ConstantParametrizedInitialize

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that reads a Ini file looking for the parameters of the function.
    subroutine ConstantParametrizedInitFromFile( self, Ini )

        implicit none

        class(constant_parametrization) :: self   !< the base class
        type(TIniFile)                  :: Ini    !< Input ini file

        self%constant_value = Ini_Read_Double_File( Ini, TRIM(self%name)//'_0', 0._dl )

    end subroutine ConstantParametrizedInitFromFile


    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that initializes the function parameters based on the values found in an input array.
    subroutine ConstantParametrizedInit( self, array )

        implicit none

        class(constant_parametrization)                         :: self   !< the base class.
        real(dl), dimension(self%parameter_number), intent(in)  :: array  !< input array with the values of the parameters.

        self%constant_value = array(1)

    end subroutine ConstantParametrizedInit

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the value of the constant function.
    function ConstantParametrizedValue( self, a )

        implicit none

        class(constant_parametrization) :: self  !< the base class
        real(dl), intent(in)            :: a     !< the input scale factor
        real(dl) :: ConstantParametrizedValue    !< the output value

        ConstantParametrizedValue = self%constant_value

    end function ConstantParametrizedValue

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the integral of the constant function, as defined in the notes
    function ConstantParametrizedIntegral( self, a )

        implicit none

        class(constant_parametrization) :: self     !< the base class
        real(dl), intent(in)            :: a        !< the scale factor at which the integral is wanted
        real(dl) :: ConstantParametrizedIntegral    !< the output value

        ConstantParametrizedIntegral = a**(-1._dl-3._dl*self%constant_value)

    end function ConstantParametrizedIntegral

    ! ---------------------------------------------------------------------------------------------

end module EFTCAMB_constant_parametrization

!----------------------------------------------------------------------------------------
