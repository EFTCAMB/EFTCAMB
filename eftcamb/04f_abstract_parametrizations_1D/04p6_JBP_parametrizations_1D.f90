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

!> @file 04p6_JBP_parametrizations_1D.f90
!! This file contains the definition of the generalized Jassal-Bagla-Padmanabhan (JBP) parametrization,
!! inheriting from parametrized_function_1D.


!----------------------------------------------------------------------------------------
!> This module contains the definition of the generalized Jassal-Bagla-Padmanabhan (JBP) parametrization,
!! inheriting from parametrized_function_1D.

!> @author Bin Hu, Marco Raveri, Simone Peirone

module EFTCAMB_JBP_parametrizations_1D

    use precision
    use AMLutils
    use EFT_def
    use EFTCAMB_cache
    use EFTCAMB_abstract_parametrizations_1D

    implicit none

    private

    public JBP_parametrization_1D

    ! ---------------------------------------------------------------------------------------------
    !> Type containing the JBP function parametrization. Inherits from parametrized_function_1D.
    type, extends ( parametrized_function_1D ) :: JBP_parametrization_1D

        real(dl) :: w0
        real(dl) :: wa
        real(dl) :: wn

    contains

        ! utility functions:
        procedure :: set_param_number      => JBPParametrized1DSetParamNumber      !< subroutine that sets the number of parameters of the JBP parametrized function.
        procedure :: init_parameters       => JBPParametrized1DInitParams          !< subroutine that initializes the function parameters based on the values found in an input array.
        procedure :: parameter_value       => JBPParametrized1DParameterValues     !< subroutine that returns the value of the function i-th parameter.
        procedure :: feedback              => JBPParametrized1DFeedback            !< subroutine that prints to screen the informations about the function.

        ! evaluation procedures:
        procedure :: value                 => JBPParametrized1DValue               !< function that returns the value of the JBP function.
        procedure :: first_derivative      => JBPParametrized1DFirstDerivative     !< function that returns the first derivative of the JBP function.
        procedure :: second_derivative     => JBPParametrized1DSecondDerivative    !< function that returns the second derivative of the JBP function.
        procedure :: third_derivative      => JBPParametrized1DThirdDerivative     !< function that returns the third derivative of the JBP function.
        procedure :: integral              => JBPParametrized1DIntegral            !< function that returns the strange integral that we need for w_DE.

    end type JBP_parametrization_1D

contains

    ! ---------------------------------------------------------------------------------------------
    ! Implementation of the JBP function.
    ! ---------------------------------------------------------------------------------------------

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that sets the number of parameters of the JBP parametrized function.
    subroutine JBPParametrized1DSetParamNumber( self )

        implicit none

        class(JBP_parametrization_1D) :: self       !< the base class

        ! initialize the number of parameters:
        self%parameter_number = 3

    end subroutine JBPParametrized1DSetParamNumber

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that initializes the function parameters based on the values found in an input array.
    subroutine JBPParametrized1DInitParams( self, array )

        implicit none

        class(JBP_parametrization_1D)                           :: self   !< the base class.
        real(dl), dimension(self%parameter_number), intent(in)  :: array  !< input array with the values of the parameters.

        self%w0 = array(1)
        self%wa = array(2)
        self%wn = array(3)

    end subroutine JBPParametrized1DInitParams

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the value of the function i-th parameter.
    subroutine JBPParametrized1DParameterValues( self, i, value )

        implicit none

        class(JBP_parametrization_1D)        :: self        !< the base class
        integer     , intent(in)             :: i           !< The index of the parameter
        real(dl)    , intent(out)            :: value       !< the output value of the i-th parameter

        select case (i)
            case(1)
                value = self%w0
            case(2)
                value = self%wa
            case(3)
                value = self%wn
            case default
                write(*,*) 'Illegal index for parameter_names.'
                write(*,*) 'Maximum value is:', self%parameter_number
                call MpiStop('EFTCAMB error')
        end select

    end subroutine JBPParametrized1DParameterValues

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that prints to screen the informations about the function.
    subroutine JBPParametrized1DFeedback( self )

        implicit none

        class(JBP_parametrization_1D)         :: self         !< the base class

        integer                               :: i
        real(dl)                              :: param_value
        character(len=EFT_names_max_length)   :: param_name

        write(*,*)     'Generalized JBP parametrization: ', self%name
        do i=1, self%parameter_number
            call self%parameter_names( i, param_name  )
            call self%parameter_value( i, param_value )
            write(*,'(a23,a,F12.6)') param_name, '=', param_value
        end do

    end subroutine JBPParametrized1DFeedback

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the value of the function in the scale factor.
    function JBPParametrized1DValue( self, x, eft_cache )

        implicit none

        class(JBP_parametrization_1D)                      :: self      !< the base class
        real(dl), intent(in)                               :: x         !< the input scale factor
        type(EFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: JBPParametrized1DValue                              !< the output value

        JBPParametrized1DValue = self%w0 +(1._dl - x)*self%wa*x**(self%wn -1._dl)
    end function JBPParametrized1DValue

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the value of the first derivative, wrt scale factor, of the function.
    function JBPParametrized1DFirstDerivative( self, x, eft_cache )

        implicit none

        class(JBP_parametrization_1D)                      :: self      !< the base class
        real(dl), intent(in)                               :: x         !< the input scale factor
        type(EFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: JBPParametrized1DFirstDerivative                    !< the output value

        JBPParametrized1DFirstDerivative = x**(-2.0_dl +self%wn)*(-1._dl +self%wn -x*self%wn)*self%wa

    end function JBPParametrized1DFirstDerivative

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the second derivative of the function.
    function JBPParametrized1DSecondDerivative( self, x, eft_cache )

        implicit none

        class(JBP_parametrization_1D)                      :: self      !< the base class
        real(dl), intent(in)                               :: x         !< the input scale factor
        type(EFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: JBPParametrized1DSecondDerivative                   !< the output value

        JBPParametrized1DSecondDerivative = -x**(-3._dl +self%wn)*(-1.0_dl +self%wn)*(2._dl +self%wn*(-1._dl +x))*self%wa

    end function JBPParametrized1DSecondDerivative

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the third derivative of the function.
    function JBPParametrized1DThirdDerivative( self, x, eft_cache )

        implicit none

        class(JBP_parametrization_1D)                      :: self      !< the base class
        real(dl), intent(in)                               :: x         !< the input scale factor
        type(EFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: JBPParametrized1DThirdDerivative                    !< the output value

        JBPParametrized1DThirdDerivative = -x**(-4._dl +self%wn)*(-2.0_dl +self%wn)*(-1._dl +self%wn)*(3._dl +self%wn*(-1._dl +x))*self%wa

    end function JBPParametrized1DThirdDerivative

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the integral of the function, as defined in the notes.
    function JBPParametrized1DIntegral( self, x, eft_cache )

        implicit none

        class(JBP_parametrization_1D)                      :: self      !< the base class
        real(dl), intent(in)                               :: x         !< the input scale factor
        type(EFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: JBPParametrized1DIntegral                           !< the output value

        JBPParametrized1DIntegral = x**(-1._dl -3._dl*self%w0)*Exp(3._dl*(x +x**self%wn*(x*(self%wn -1._dl) -self%wn))*self%wa/(x*self%wn*(self%wn -1._dl)))

    end function JBPParametrized1DIntegral

    ! ---------------------------------------------------------------------------------------------

end module EFTCAMB_JBP_parametrizations_1D

!----------------------------------------------------------------------------------------
