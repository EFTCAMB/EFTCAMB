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

!> @file 04p2_linear_parametrizations_1D.f90
!! This file contains the definition of the linear parametrization, inheriting from
!! parametrized_function_1D.


!----------------------------------------------------------------------------------------
!> This module contains the definition of the linear parametrization, inheriting from
!! parametrized_function_1D.

!> @author Bin Hu, Marco Raveri, Simone Peirone

module EFTCAMB_linear_parametrizations_1D

    use precision
    use AMLutils
    use EFT_def
    use EFTCAMB_cache
    use EFTCAMB_abstract_parametrizations_1D

    implicit none

    private

    public linear_parametrization_1D

    ! ---------------------------------------------------------------------------------------------
    !> Type containing the linear function parametrization. Inherits from parametrized_function_1D.
    type, extends ( parametrized_function_1D ) :: linear_parametrization_1D

        real(dl) :: linear_value

    contains

        ! utility functions:
        procedure :: set_param_number      => LinearParametrized1DSetParamNumber      !< subroutine that sets the number of parameters of the linear parametrized function.
        procedure :: init_parameters       => LinearParametrized1DInitParams          !< subroutine that initializes the function parameters based on the values found in an input array.
        procedure :: parameter_value       => LinearParametrized1DParameterValues     !< subroutine that returns the value of the function i-th parameter.
        procedure :: feedback              => LinearParametrized1DFeedback            !< subroutine that prints to screen the informations about the function.

        ! evaluation procedures:
        procedure :: value                 => LinearParametrized1DValue               !< function that returns the value of the linear function.
        procedure :: first_derivative      => LinearParametrized1DFirstDerivative     !< function that returns the first derivative of the linear function.
        procedure :: second_derivative     => LinearParametrized1DSecondDerivative    !< function that returns the second derivative of the linear function.
        procedure :: third_derivative      => LinearParametrized1DThirdDerivative     !< function that returns the third derivative of the linear function.
        procedure :: integral              => LinearParametrized1DIntegral            !< function that returns the strange integral that we need for w_DE.

    end type linear_parametrization_1D

contains

    ! ---------------------------------------------------------------------------------------------
    ! Implementation of the linear function.
    ! ---------------------------------------------------------------------------------------------

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that sets the number of parameters of the linear parametrized function.
    subroutine LinearParametrized1DSetParamNumber( self )

        implicit none

        class(linear_parametrization_1D) :: self       !< the base class

        ! initialize the number of parameters:
        self%parameter_number = 1

    end subroutine LinearParametrized1DSetParamNumber

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that initializes the function parameters based on the values found in an input array.
    subroutine LinearParametrized1DInitParams( self, array )

        implicit none

        class(linear_parametrization_1D)                        :: self   !< the base class.
        real(dl), dimension(self%parameter_number), intent(in)  :: array  !< input array with the values of the parameters.

        self%linear_value = array(1)

    end subroutine LinearParametrized1DInitParams

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
            write(*,'(a23,a,F12.6)') param_name, '=', param_value
        end do

    end subroutine LinearParametrized1DFeedback

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the value of the linear function in the scale factor.
    function LinearParametrized1DValue( self, x, eft_cache )

        implicit none

        class(linear_parametrization_1D)                   :: self      !< the base class
        real(dl), intent(in)                               :: x         !< the input scale factor
        type(EFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: LinearParametrized1DValue                           !< the output value

        LinearParametrized1DValue = self%linear_value*x

    end function LinearParametrized1DValue

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the value of the first derivative, wrt scale factor, of the linear function.
    function LinearParametrized1DFirstDerivative( self, x, eft_cache )

        implicit none

        class(linear_parametrization_1D)                   :: self      !< the base class
        real(dl), intent(in)                               :: x         !< the input scale factor
        type(EFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: LinearParametrized1DFirstDerivative                 !< the output value

        LinearParametrized1DFirstDerivative = self%linear_value

    end function LinearParametrized1DFirstDerivative

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the second derivative of the linear function.
    function LinearParametrized1DSecondDerivative( self, x, eft_cache )

        implicit none

        class(linear_parametrization_1D)                   :: self      !< the base class
        real(dl), intent(in)                               :: x         !< the input scale factor
        type(EFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: LinearParametrized1DSecondDerivative                !< the output value

        LinearParametrized1DSecondDerivative = 0._dl

    end function LinearParametrized1DSecondDerivative

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the third derivative of the linear function.
    function LinearParametrized1DThirdDerivative( self, x, eft_cache )

        implicit none

        class(linear_parametrization_1D)                   :: self      !< the base class
        real(dl), intent(in)                               :: x         !< the input scale factor
        type(EFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: LinearParametrized1DThirdDerivative                 !< the output value

        LinearParametrized1DThirdDerivative = 0._dl

    end function LinearParametrized1DThirdDerivative

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the integral of the linear function, as defined in the notes.
    function LinearParametrized1DIntegral( self, x, eft_cache )

        implicit none

        class(linear_parametrization_1D)                   :: self      !< the base class
        real(dl), intent(in)                               :: x         !< the input scale factor
        type(EFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: LinearParametrized1DIntegral                        !< the output value

        LinearParametrized1DIntegral = Exp(-3._dl*(x-1._dl)*self%linear_value)/x

    end function LinearParametrized1DIntegral

    ! ---------------------------------------------------------------------------------------------

end module EFTCAMB_linear_parametrizations_1D

!----------------------------------------------------------------------------------------
