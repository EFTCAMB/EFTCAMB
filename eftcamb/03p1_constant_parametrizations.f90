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

!> @file 03p1_constant_parametrizations.f90
!! This file contains the definition of the constant parametrization, inheriting from
!! parametrized_function_1D.


!----------------------------------------------------------------------------------------
!> This module contains the definition of the constant parametrization, inheriting from
!! parametrized_function_1D.

!> @author Bin Hu, Marco Raveri

module EFTCAMB_constant_parametrization_1D

    use precision
    use EFTDef
    use EFTCAMB_abstract_parametrizations_1D

    implicit none

    ! ---------------------------------------------------------------------------------------------
    !> Type containing the constant function parametrization. Inherits from parametrized_function_1D.
    !! Notice that the derivatives are not overridden since they are zero identically.
    type, extends ( parametrized_function_1D ) :: constant_parametrization_1D

        real(dl) :: constant_value

    contains

        ! initialization:
        procedure :: init                  => ConstantParametrized1DInitialize          !< subroutine that initializes the constant parametrization
        procedure :: init_from_file        => ConstantParametrized1DInitFromFile        !< subroutine that reads a Ini file looking for the parameters of the function.
        procedure :: init_parameters       => ConstantParametrized1DInit                !< subroutine that initializes the function parameters based on the values found in an input array.

        ! utility functions:
        procedure :: feedback              => ConstantParametrized1DFeedback            !< subroutine that prints to screen the informations about the function.
        procedure :: parameter_names       => ConstantParametrized1DParameterNames      !< subroutine that returns the i-th parameter name of the function
        procedure :: parameter_names_latex => ConstantParametrized1DParameterNamesLatex !< subroutine that returns the i-th parameter name of the function in latex format
        procedure :: parameter_value       => ConstantParametrized1DParameterValues     !< subroutine that returns the value of the function i-th parameter.

        ! evaluation procedures:
        procedure :: value                 => ConstantParametrized1DValue               !< function that returns the value of the function
        procedure :: integral              => ConstantParametrized1DIntegral            !< function that returns the strange integral that we need for w_DE

    end type constant_parametrization_1D

contains

    ! ---------------------------------------------------------------------------------------------
    ! Implementation of the constant function.
    ! ---------------------------------------------------------------------------------------------

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that initializes the constant parametrization
    subroutine ConstantParametrized1DInitialize( self, name, latexname )

        implicit none

        class(constant_parametrization_1D) :: self       !< the base class
        character(*), intent(in)           :: name       !< the name of the function
        character(*), intent(in)           :: latexname  !< the latex name of the function

        ! store the name of the function:
        self%name             = TRIM( name )
        ! store the latex name of the function:
        self%name_latex       = TRIM( latexname )
        ! initialize the number of parameters:
        self%parameter_number = 1

    end subroutine ConstantParametrized1DInitialize

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that reads a Ini file looking for the parameters of the function.
    subroutine ConstantParametrized1DInitFromFile( self, Ini )

        implicit none

        class(constant_parametrization_1D) :: self   !< the base class
        type(TIniFile)                     :: Ini    !< Input ini file

        character(len=EFT_names_max_length) :: param_name

        call self%parameter_names( 1, param_name )

        self%constant_value = Ini_Read_Double_File( Ini, TRIM(param_name), 0._dl )

    end subroutine ConstantParametrized1DInitFromFile

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that initializes the function parameters based on the values found in an input array.
    subroutine ConstantParametrized1DInit( self, array )

        implicit none

        class(constant_parametrization_1D)                      :: self   !< the base class.
        real(dl), dimension(self%parameter_number), intent(in)  :: array  !< input array with the values of the parameters.

        self%constant_value = array(1)

    end subroutine ConstantParametrized1DInit

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
    !> Subroutine that returns the i-th parameter name
    subroutine ConstantParametrized1DParameterNames( self, i, name )

        implicit none

        class(constant_parametrization_1D) :: self   !< the base class
        integer     , intent(in)           :: i      !< the index of the parameter
        character(*), intent(out)          :: name   !< the output name of the i-th parameter

        select case (i)
            case(1)
                name = TRIM(self%name)//'0'
            case default
                write(*,*) 'Illegal index for parameter_names.'
                write(*,*) 'Maximum value is:', self%parameter_number
                stop
        end select

    end subroutine ConstantParametrized1DParameterNames

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the latex version of the i-th parameter name
    subroutine ConstantParametrized1DParameterNamesLatex( self, i, latexname )

        implicit none

        class(constant_parametrization_1D) :: self        !< the base class
        integer     , intent(in)           :: i           !< The index of the parameter
        character(*), intent(out)          :: latexname   !< the output latex name of the i-th parameter

        select case (i)
            case(1)
                latexname = TRIM(self%name_latex)//'_0'
            case default
                write(*,*) 'Illegal index for parameter_names.'
                write(*,*) 'Maximum value is:', self%parameter_number
                stop
        end select

    end subroutine ConstantParametrized1DParameterNamesLatex


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
                write(*,*) 'Illegal index for parameter_names.'
                write(*,*) 'Maximum value is:', self%parameter_number
                stop
        end select

    end subroutine ConstantParametrized1DParameterValues

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the value of the constant function.
    function ConstantParametrized1DValue( self, x )

        implicit none

        class(constant_parametrization_1D) :: self  !< the base class
        real(dl), intent(in)               :: x     !< the input scale factor
        real(dl) :: ConstantParametrized1DValue     !< the output value

        ConstantParametrized1DValue = self%constant_value

    end function ConstantParametrized1DValue

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the integral of the constant function, as defined in the notes
    function ConstantParametrized1DIntegral( self, x )

        implicit none

        class(constant_parametrization_1D) :: self     !< the base class
        real(dl), intent(in)               :: x        !< the scale factor at which the integral is wanted
        real(dl) :: ConstantParametrized1DIntegral     !< the output value

        ConstantParametrized1DIntegral = x**(-1._dl-3._dl*self%constant_value)

    end function ConstantParametrized1DIntegral

    ! ---------------------------------------------------------------------------------------------

end module EFTCAMB_constant_parametrization_1D

!----------------------------------------------------------------------------------------
