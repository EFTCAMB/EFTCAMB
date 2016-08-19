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

!> @file 03p2_linear_parametrizations.f90
!! This file contains the definition of the linear parametrization, inheriting from
!! parametrized_function_1D.


!----------------------------------------------------------------------------------------
!> This module contains the definition of the linear parametrization, inheriting from
!! parametrized_function_1D.

!> @author Bin Hu, Marco Raveri

module EFTCAMB_linear_parametrizations_1D

    use precision
    use EFTDef
    use AMLutils
    use EFTCAMB_abstract_parametrizations_1D

    implicit none

    ! ---------------------------------------------------------------------------------------------
    !> Type containing the linear function parametrization. Inherits from parametrized_function_1D.
    !! Notice that the derivatives above the first are not overridden since they are zero identically.
    type, extends ( parametrized_function_1D ) :: linear_parametrization_1D

        real(dl) :: linear_value

    contains

        ! initialization:
        procedure :: init                  => LinearParametrized1DInitialize          !< subroutine that initializes the constant parametrization
        procedure :: init_from_file        => LinearParametrized1DInitFromFile        !< subroutine that reads a Ini file looking for the parameters of the function.
        procedure :: init_parameters       => LinearParametrized1DInit                !< subroutine that initializes the function parameters based on the values found in an input array.
        ! utility functions:
        procedure :: feedback              => LinearParametrized1DFeedback            !< subroutine that prints to screen the informations about the function.
        procedure :: parameter_names       => LinearParametrized1DParameterNames      !< subroutine that returns the i-th parameter name of the function
        procedure :: parameter_names_latex => LinearParametrized1DParameterNamesLatex !< subroutine that returns the i-th parameter name of the function in latex format
        procedure :: parameter_value       => LinearParametrized1DParameterValues     !< subroutine that returns the value of the function i-th parameter.

        ! evaluation procedures:
        procedure :: value                 => LinearParametrized1DValue               !< function that returns the value of the function
        procedure :: first_derivative      => LinearParametrized1DFirstDerivative     !< function that returns the first derivative of the function
        procedure :: integral              => LinearParametrized1DIntegral            !< function that returns the strange integral that we need for w_DE

    end type linear_parametrization_1D

contains

    ! ---------------------------------------------------------------------------------------------
    ! Implementation of the linear function.
    ! ---------------------------------------------------------------------------------------------

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that initializes the linear parametrization
    subroutine LinearParametrized1DInitialize( self, name, latexname )

        implicit none

        class(linear_parametrization_1D) :: self         !< the base class
        character(*), intent(in)         :: name         !< the name of the function
        character(*), intent(in)         :: latexname    !< the latex name of the function

        ! store the name of the function:
        self%name             = TRIM( name )
        ! store the latex name of the function:
        self%name_latex       = TRIM( latexname )
        ! initialize the number of parameters:
        self%parameter_number = 1

    end subroutine LinearParametrized1DInitialize

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that reads a Ini file looking for the parameters of the function.
    subroutine LinearParametrized1DInitFromFile( self, Ini )

        implicit none

        class(linear_parametrization_1D)  :: self   !< the base class
        type(TIniFile)                    :: Ini    !< Input ini file

        character(len=EFT_names_max_length) :: param_name

        call self%parameter_names( 1, param_name )

        self%linear_value = Ini_Read_Double_File( Ini, TRIM(param_name), 0._dl )

    end subroutine LinearParametrized1DInitFromFile

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that initializes the function parameters based on the values found in an input array.
    subroutine LinearParametrized1DInit( self, array )

        implicit none

        class(linear_parametrization_1D)                        :: self   !< the base class.
        real(dl), dimension(self%parameter_number), intent(in)  :: array  !< input array with the values of the parameters.

        self%linear_value = array(1)

    end subroutine LinearParametrized1DInit

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that prints to screen the informations about the function.
    subroutine LinearParametrized1DFeedback( self )

        implicit none

        class(linear_parametrization_1D)    :: self         !< the base class

        integer                             :: i
        real(dl)                            :: param_value
        character(len=EFT_names_max_length) :: param_name

        write(*,*)     'Linear function: ', self%name
        do i=1, self%parameter_number
            call self%parameter_names( i, param_name  )
            call self%parameter_value( i, param_value )
            write(*,*) param_name, '=', param_value
        end do

    end subroutine LinearParametrized1DFeedback

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the i-th parameter name
    subroutine LinearParametrized1DParameterNames( self, i, name )

        implicit none

        class(linear_parametrization_1D)  :: self   !< the base class
        integer     , intent(in)          :: i      !< the index of the parameter
        character(*), intent(out)         :: name   !< the output name of the i-th parameter

        select case (i)
            case(1)
                name = TRIM(self%name)//'0'
            case default
                write(*,*) 'Illegal index for parameter_names.'
                write(*,*) 'Maximum value is:', self%parameter_number
                call MpiStop('EFTCAMB error')
        end select

    end subroutine LinearParametrized1DParameterNames

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the latex version of the i-th parameter name
    subroutine LinearParametrized1DParameterNamesLatex( self, i, latexname )

        implicit none

        class(linear_parametrization_1D) :: self        !< the base class
        integer     , intent(in)         :: i           !< The index of the parameter
        character(*), intent(out)        :: latexname   !< the output latex name of the i-th parameter

        select case (i)
            case(1)
                latexname = TRIM(self%name_latex)//'_0'
            case default
                write(*,*) 'Illegal index for parameter_names.'
                write(*,*) 'Maximum value is:', self%parameter_number
                call MpiStop('EFTCAMB error')
        end select

    end subroutine LinearParametrized1DParameterNamesLatex

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the value of the function i-th parameter.
    subroutine LinearParametrized1DParameterValues( self, i, value )

        implicit none

        class(linear_parametrization_1D) :: self        !< the base class
        integer     , intent(in)         :: i           !< The index of the parameter
        real(dl)    , intent(out)        :: value       !< the output value of the i-th parameter

        select case (i)
            case(1)
                value = self%linear_value
            case default
                write(*,*) 'Illegal index for parameter_names.'
                write(*,*) 'Maximum value is:', self%parameter_number
                call MpiStop('EFTCAMB error')
        end select

    end subroutine LinearParametrized1DParameterValues

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the value of the linear function in the scale factor.
    function LinearParametrized1DValue( self, x )

        implicit none

        class(linear_parametrization_1D) :: self  !< the base class
        real(dl), intent(in)             :: x     !< the input scale factor
        real(dl) :: LinearParametrized1DValue     !< the output value

        LinearParametrized1DValue = self%linear_value*x

    end function LinearParametrized1DValue

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the value of the first derivative, wrt scale factor, of the linear function.
    function LinearParametrized1DFirstDerivative( self, x )

        implicit none

        class(linear_parametrization_1D) :: self          !< the base class
        real(dl), intent(in)             :: x             !< the input scale factor
        real(dl) :: LinearParametrized1DFirstDerivative   !< the output value

        LinearParametrized1DFirstDerivative = self%linear_value

    end function LinearParametrized1DFirstDerivative

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the integral of the linear function, as defined in the notes
    function LinearParametrized1DIntegral( self, x )

        implicit none

        class(linear_parametrization_1D) :: self     !< the base class
        real(dl), intent(in)             :: x        !< the scale factor at which the integral is wanted
        real(dl) :: LinearParametrized1DIntegral     !< the output value

        LinearParametrized1DIntegral = Exp(-3._dl*(x-1._dl)*self%linear_value)/x

    end function LinearParametrized1DIntegral

    ! ---------------------------------------------------------------------------------------------

end module EFTCAMB_linear_parametrizations_1D

!----------------------------------------------------------------------------------------
