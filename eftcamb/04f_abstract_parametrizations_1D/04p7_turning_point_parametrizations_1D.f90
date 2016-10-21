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

!> @file 04p7_turning_point_parametrizations_1D.f90
!! This file contains the definition of the turning point parametrization,
!! inheriting from parametrized_function_1D.


!----------------------------------------------------------------------------------------
!> This module contains the definition of the turning point parametrization,
!! inheriting from parametrized_function_1D.

!> @author Bin Hu, Marco Raveri, Simone Peirone

module EFTCAMB_turning_point_parametrizations_1D

    use precision
    use AMLutils
    use EFT_def
    use EFTCAMB_cache
    use EFTCAMB_abstract_parametrizations_1D

    implicit none

    private

    public turning_point_parametrization_1D

    ! ---------------------------------------------------------------------------------------------
    !> Type containing the JBP function parametrization. Inherits from parametrized_function_1D.
    type, extends ( parametrized_function_1D ) :: turning_point_parametrization_1D

        real(dl) :: w0
        real(dl) :: wa
        real(dl) :: wat

    contains

        ! utility functions:
        procedure :: set_param_number      => TurningPointParametrized1DSetParamNumber      !< subroutine that sets the number of parameters of the turning point parametrized function.
        procedure :: init_parameters       => TurningPointParametrized1DInitParams          !< subroutine that initializes the function parameters based on the values found in an input array.
        procedure :: parameter_value       => TurningPointParametrized1DParameterValues     !< subroutine that returns the value of the function i-th parameter.
        procedure :: feedback              => TurningPointParametrized1DFeedback            !< subroutine that prints to screen the informations about the function.

        ! evaluation procedures:
        procedure :: value                 => TurningPointParametrized1DValue               !< function that returns the value of the turning point function.
        procedure :: first_derivative      => TurningPointParametrized1DFirstDerivative     !< function that returns the first derivative of the turning point function.
        procedure :: second_derivative     => TurningPointParametrized1DSecondDerivative    !< function that returns the second derivative of the turning point function.
        procedure :: third_derivative      => TurningPointParametrized1DThirdDerivative     !< function that returns the third derivative of the turning point function.
        procedure :: integral              => TurningPointParametrized1DIntegral            !< function that returns the strange integral that we need for w_DE.

    end type turning_point_parametrization_1D

contains

    ! ---------------------------------------------------------------------------------------------
    ! Implementation of the turning point function.
    ! ---------------------------------------------------------------------------------------------

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that sets the number of parameters of the turning point parametrized function.
    subroutine TurningPointParametrized1DSetParamNumber( self )

        implicit none

        class(turning_point_parametrization_1D) :: self       !< the base class

        ! initialize the number of parameters:
        self%parameter_number = 3

    end subroutine TurningPointParametrized1DSetParamNumber

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that initializes the function parameters based on the values found in an input array.
    subroutine TurningPointParametrized1DInitParams( self, array )

        implicit none

        class(turning_point_parametrization_1D)                 :: self   !< the base class.
        real(dl), dimension(self%parameter_number), intent(in)  :: array  !< input array with the values of the parameters.

        self%w0 = array(1)
        self%wa = array(2)
        self%wat= array(3)

    end subroutine TurningPointParametrized1DInitParams

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the value of the function i-th parameter.
    subroutine TurningPointParametrized1DParameterValues( self, i, value )

        implicit none

        class(turning_point_parametrization_1D):: self        !< the base class
        integer     , intent(in)               :: i           !< The index of the parameter
        real(dl)    , intent(out)              :: value       !< the output value of the i-th parameter

        select case (i)
            case(1)
                value = self%w0
            case(2)
                value = self%wa
            case(3)
                value = self%wat
            case default
                write(*,*) 'Illegal index for parameter_names.'
                write(*,*) 'Maximum value is:', self%parameter_number
                call MpiStop('EFTCAMB error')
        end select

    end subroutine TurningPointParametrized1DParameterValues

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that prints to screen the informations about the function.
    subroutine TurningPointParametrized1DFeedback( self )

        implicit none

        class(turning_point_parametrization_1D) :: self         !< the base class

        integer                                 :: i
        real(dl)                                :: param_value
        character(len=EFT_names_max_length)     :: param_name

        write(*,*)     'Turning Point parametrization: ', self%name
        do i=1, self%parameter_number
            call self%parameter_names( i, param_name  )
            call self%parameter_value( i, param_value )
            write(*,'(a23,a,F12.6)') param_name, '=', param_value
        end do

    end subroutine TurningPointParametrized1DFeedback

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the value of the function in the scale factor.
    function TurningPointParametrized1DValue( self, x, eft_cache )

        implicit none

        class(turning_point_parametrization_1D)            :: self      !< the base class
        real(dl), intent(in)                               :: x         !< the input scale factor
        type(EFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: TurningPointParametrized1DValue                     !< the output value

        TurningPointParametrized1DValue = self%w0 +self%wa*(self%wat - x)**2._dl

    end function TurningPointParametrized1DValue

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the value of the first derivative, wrt scale factor, of the function.
    function TurningPointParametrized1DFirstDerivative( self, x, eft_cache )

        implicit none

        class(turning_point_parametrization_1D)            :: self      !< the base class
        real(dl), intent(in)                               :: x         !< the input scale factor
        type(EFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: TurningPointParametrized1DFirstDerivative           !< the output value

        TurningPointParametrized1DFirstDerivative = 2._dl*self%wa*(x -self%wat)

    end function TurningPointParametrized1DFirstDerivative

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the second derivative of the function.
    function TurningPointParametrized1DSecondDerivative( self, x, eft_cache )

        implicit none

        class(turning_point_parametrization_1D)            :: self      !< the base class
        real(dl), intent(in)                               :: x         !< the input scale factor
        type(EFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: TurningPointParametrized1DSecondDerivative          !< the output value

        TurningPointParametrized1DSecondDerivative = 2._dl*self%wa

    end function TurningPointParametrized1DSecondDerivative

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the third derivative of the function.
    function TurningPointParametrized1DThirdDerivative( self, x, eft_cache )

        implicit none

        class(turning_point_parametrization_1D)            :: self      !< the base class
        real(dl), intent(in)                               :: x         !< the input scale factor
        type(EFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: TurningPointParametrized1DThirdDerivative           !< the output value

        TurningPointParametrized1DThirdDerivative = 0._dl

    end function TurningPointParametrized1DThirdDerivative

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the integral of the function, as defined in the notes.
    function TurningPointParametrized1DIntegral( self, x, eft_cache )

        implicit none

        class(turning_point_parametrization_1D)            :: self      !< the base class
        real(dl), intent(in)                               :: x         !< the input scale factor
        type(EFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: TurningPointParametrized1DIntegral                  !< the output value

        TurningPointParametrized1DIntegral = x**(2._dl -3._dl*(1._dl +self%w0 +self%wa*self%wat**2._dl))*Exp(-1.5_dl*(x -1._dl)*self%wa*(1._dl +x -4._dl*self%wat))

    end function TurningPointParametrized1DIntegral

    ! ---------------------------------------------------------------------------------------------

end module EFTCAMB_turning_point_parametrizations_1D

!----------------------------------------------------------------------------------------
