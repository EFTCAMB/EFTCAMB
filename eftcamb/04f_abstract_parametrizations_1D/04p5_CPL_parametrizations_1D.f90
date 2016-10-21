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

!> @file 04p5_CPL_parametrizations_1D.f90
!! This file contains the definition of the CPL parametrization,
!! inheriting from parametrized_function_1D.


!----------------------------------------------------------------------------------------
!> This module contains the definition of the CPL parametrization,
!! inheriting from parametrized_function_1D.

!> @author Bin Hu, Marco Raveri, Simone Peirone

module EFTCAMB_CPL_parametrizations_1D

    use precision
    use AMLutils
    use EFT_def
    use EFTCAMB_cache
    use EFTCAMB_abstract_parametrizations_1D

    implicit none

    private

    public CPL_parametrization_1D

    ! ---------------------------------------------------------------------------------------------
    !> Type containing the CPL function parametrization. Inherits from parametrized_function_1D.
    type, extends ( parametrized_function_1D ) :: CPL_parametrization_1D

        real(dl) :: w0
        real(dl) :: wa

    contains

        ! utility functions:
        procedure :: set_param_number      => CPLParametrized1DSetParamNumber      !< subroutine that sets the number of parameters of the CPL parametrized function.
        procedure :: init_parameters       => CPLParametrized1DInitParams          !< subroutine that initializes the function parameters based on the values found in an input array.
        procedure :: parameter_value       => CPLParametrized1DParameterValues     !< subroutine that returns the value of the function i-th parameter.
        procedure :: feedback              => CPLParametrized1DFeedback            !< subroutine that prints to screen the informations about the function.

        ! evaluation procedures:
        procedure :: value                 => CPLParametrized1DValue               !< function that returns the value of the CPL function.
        procedure :: first_derivative      => CPLParametrized1DFirstDerivative     !< function that returns the first derivative of the CPL function.
        procedure :: second_derivative     => CPLParametrized1DSecondDerivative    !< function that returns the second derivative of the CPL function.
        procedure :: third_derivative      => CPLParametrized1DThirdDerivative     !< function that returns the third derivative of the CPL function.
        procedure :: integral              => CPLParametrized1DIntegral            !< function that returns the strange integral that we need for w_DE.

    end type CPL_parametrization_1D

contains

    ! ---------------------------------------------------------------------------------------------
    ! Implementation of the CPL function.
    ! ---------------------------------------------------------------------------------------------

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that sets the number of parameters of the CPL parametrized function.
    subroutine CPLParametrized1DSetParamNumber( self )

        implicit none

        class(CPL_parametrization_1D) :: self       !< the base class

        ! initialize the number of parameters:
        self%parameter_number = 2

    end subroutine CPLParametrized1DSetParamNumber

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that initializes the function parameters based on the values found in an input array.
    subroutine CPLParametrized1DInitParams( self, array )

        implicit none

        class(CPL_parametrization_1D)                           :: self   !< the base class.
        real(dl), dimension(self%parameter_number), intent(in)  :: array  !< input array with the values of the parameters.

        self%w0 = array(1)
        self%wa = array(2)

    end subroutine CPLParametrized1DInitParams

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the value of the function i-th parameter.
    subroutine CPLParametrized1DParameterValues( self, i, value )

        implicit none

        class(CPL_parametrization_1D)        :: self        !< the base class
        integer     , intent(in)             :: i           !< The index of the parameter
        real(dl)    , intent(out)            :: value       !< the output value of the i-th parameter

        select case (i)
            case(1)
                value = self%w0
            case(2)
                value = self%wa
            case default
                write(*,*) 'Illegal index for parameter_names.'
                write(*,*) 'Maximum value is:', self%parameter_number
                call MpiStop('EFTCAMB error')
        end select

    end subroutine CPLParametrized1DParameterValues

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that prints to screen the informations about the function.
    subroutine CPLParametrized1DFeedback( self )

        implicit none

        class(CPL_parametrization_1D)         :: self         !< the base class

        integer                               :: i
        real(dl)                              :: param_value
        character(len=EFT_names_max_length)   :: param_name

        write(*,*)     'CPL parametrization: ', self%name
        do i=1, self%parameter_number
            call self%parameter_names( i, param_name  )
            call self%parameter_value( i, param_value )
            write(*,'(a23,a,F12.6)') param_name, '=', param_value
        end do

    end subroutine CPLParametrized1DFeedback

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the value of the function in the scale factor.
    function CPLParametrized1DValue( self, x, eft_cache )

        implicit none

        class(CPL_parametrization_1D)                      :: self      !< the base class
        real(dl), intent(in)                               :: x         !< the input scale factor
        type(EFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: CPLParametrized1DValue                              !< the output value

        CPLParametrized1DValue = self%w0 +(1._dl - x)*self%wa
    end function CPLParametrized1DValue

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the value of the first derivative, wrt scale factor, of the function.
    function CPLParametrized1DFirstDerivative( self, x, eft_cache )

        implicit none

        class(CPL_parametrization_1D)                      :: self      !< the base class
        real(dl), intent(in)                               :: x         !< the input scale factor
        type(EFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: CPLParametrized1DFirstDerivative                    !< the output value

        CPLParametrized1DFirstDerivative = -self%wa

    end function CPLParametrized1DFirstDerivative

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the second derivative of the function.
    function CPLParametrized1DSecondDerivative( self, x, eft_cache )

        implicit none

        class(CPL_parametrization_1D)                      :: self      !< the base class
        real(dl), intent(in)                               :: x         !< the input scale factor
        type(EFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: CPLParametrized1DSecondDerivative                   !< the output value

        CPLParametrized1DSecondDerivative = 0._dl

    end function CPLParametrized1DSecondDerivative

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the third derivative of the function.
    function CPLParametrized1DThirdDerivative( self, x, eft_cache )

        implicit none

        class(CPL_parametrization_1D)                      :: self      !< the base class
        real(dl), intent(in)                               :: x         !< the input scale factor
        type(EFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: CPLParametrized1DThirdDerivative                    !< the output value

        CPLParametrized1DThirdDerivative = 0._dl

    end function CPLParametrized1DThirdDerivative

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the integral of the function, as defined in the notes.
    function CPLParametrized1DIntegral( self, x, eft_cache )

        implicit none

        class(CPL_parametrization_1D)                      :: self      !< the base class
        real(dl), intent(in)                               :: x         !< the input scale factor
        type(EFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: CPLParametrized1DIntegral                           !< the output value

        CPLParametrized1DIntegral = x**(2._dl -3._dl*(1._dl +self%w0 +self%wa))*Exp(3._dl*(x-1._dl)*self%wa)

    end function CPLParametrized1DIntegral

    ! ---------------------------------------------------------------------------------------------

end module EFTCAMB_CPL_parametrizations_1D

!----------------------------------------------------------------------------------------
