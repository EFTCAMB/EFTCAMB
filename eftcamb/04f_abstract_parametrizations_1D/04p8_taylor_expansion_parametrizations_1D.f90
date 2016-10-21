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

!> @file 04p8_taylor_expansion_parametrizations_1D.f90
!! This file contains the definition of the Taylor expansion parametrization, around a=0,
!! up to third order, inheriting from parametrized_function_1D.


!----------------------------------------------------------------------------------------
!> This module contains the definition of the Taylor expansion parametrization, around a=0,
!! up to third order, inheriting from parametrized_function_1D.

!> @author Bin Hu, Marco Raveri, Simone Peirone

module EFTCAMB_taylor_parametrizations_1D

    use precision
    use AMLutils
    use EFT_def
    use EFTCAMB_cache
    use EFTCAMB_abstract_parametrizations_1D

    implicit none

    private

    public taylor_parametrization_1D

    ! ---------------------------------------------------------------------------------------------
    !> Type containing the Taylor expansion parametrization. Inherits from parametrized_function_1D.
    type, extends ( parametrized_function_1D ) :: taylor_parametrization_1D

        real(dl) :: w0
        real(dl) :: wa
        real(dl) :: w2
        real(dl) :: w3

    contains

        ! utility functions:
        procedure :: set_param_number      => TaylorParametrized1DSetParamNumber      !< subroutine that sets the number of parameters of the Taylor expansion parametrized function.
        procedure :: init_parameters       => TaylorParametrized1DInitParams          !< subroutine that initializes the function parameters based on the values found in an input array.
        procedure :: parameter_value       => TaylorParametrized1DParameterValues     !< subroutine that returns the value of the function i-th parameter.
        procedure :: feedback              => TaylorParametrized1DFeedback            !< subroutine that prints to screen the informations about the function.

        ! evaluation procedures:
        procedure :: value                 => TaylorParametrized1DValue               !< function that returns the value of the Taylor expansion.
        procedure :: first_derivative      => TaylorParametrized1DFirstDerivative     !< function that returns the first derivative of the Taylor expansion.
        procedure :: second_derivative     => TaylorParametrized1DSecondDerivative    !< function that returns the second derivative of the Taylor expansion.
        procedure :: third_derivative      => TaylorParametrized1DThirdDerivative     !< function that returns the third derivative of the Taylor expansion.
        procedure :: integral              => TaylorParametrized1DIntegral            !< function that returns the strange integral that we need for w_DE.

    end type taylor_parametrization_1D

contains

    ! ---------------------------------------------------------------------------------------------
    ! Implementation of the Taylor expansion parametrization.
    ! ---------------------------------------------------------------------------------------------

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that sets the number of parameters of the Taylor expansion parametrized function.
    subroutine TaylorParametrized1DSetParamNumber( self )

        implicit none

        class(taylor_parametrization_1D) :: self       !< the base class

        ! initialize the number of parameters:
        self%parameter_number = 4

    end subroutine TaylorParametrized1DSetParamNumber

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that initializes the function parameters based on the values found in an input array.
    subroutine TaylorParametrized1DInitParams( self, array )

        implicit none

        class(taylor_parametrization_1D)                        :: self   !< the base class.
        real(dl), dimension(self%parameter_number), intent(in)  :: array  !< input array with the values of the parameters.

        self%w0 = array(1)
        self%wa = array(2)
        self%w2 = array(3)
        self%w3 = array(4)

    end subroutine TaylorParametrized1DInitParams

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the value of the function i-th parameter.
    subroutine TaylorParametrized1DParameterValues( self, i, value )

        implicit none

        class(taylor_parametrization_1D)       :: self        !< the base class
        integer     , intent(in)               :: i           !< The index of the parameter
        real(dl)    , intent(out)              :: value       !< the output value of the i-th parameter

        select case (i)
            case(1)
                value = self%w0
            case(2)
                value = self%wa
            case(3)
                value = self%w2
            case(4)
                value = self%w3
            case default
                write(*,*) 'Illegal index for parameter_names.'
                write(*,*) 'Maximum value is:', self%parameter_number
                call MpiStop('EFTCAMB error')
        end select

    end subroutine TaylorParametrized1DParameterValues

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that prints to screen the informations about the function.
    subroutine TaylorParametrized1DFeedback( self )

        implicit none

        class(taylor_parametrization_1D)        :: self         !< the base class

        integer                                 :: i
        real(dl)                                :: param_value
        character(len=EFT_names_max_length)     :: param_name

        write(*,*)     'Taylor expansion parametrization: ', self%name
        do i=1, self%parameter_number
            call self%parameter_names( i, param_name  )
            call self%parameter_value( i, param_value )
            write(*,'(a23,a,F12.6)') param_name, '=', param_value
        end do

    end subroutine TaylorParametrized1DFeedback

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the value of the function in the scale factor.
    function TaylorParametrized1DValue( self, x, eft_cache )

        implicit none

        class(taylor_parametrization_1D)                   :: self      !< the base class
        real(dl), intent(in)                               :: x         !< the input scale factor
        type(EFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: TaylorParametrized1DValue                           !< the output value

        TaylorParametrized1DValue = self%w0 +self%wa*x +0.5_dl*self%w2*x**2 +self%w3/6._dl*x**3

    end function TaylorParametrized1DValue

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the value of the first derivative, wrt scale factor, of the function.
    function TaylorParametrized1DFirstDerivative( self, x, eft_cache )

        implicit none

        class(taylor_parametrization_1D)                   :: self      !< the base class
        real(dl), intent(in)                               :: x         !< the input scale factor
        type(EFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: TaylorParametrized1DFirstDerivative                 !< the output value

        TaylorParametrized1DFirstDerivative = self%wa +self%w2*x +0.5_dl*self%w3*x**2


    end function TaylorParametrized1DFirstDerivative

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the second derivative of the function.
    function TaylorParametrized1DSecondDerivative( self, x, eft_cache )

        implicit none

        class(taylor_parametrization_1D)                   :: self      !< the base class
        real(dl), intent(in)                               :: x         !< the input scale factor
        type(EFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: TaylorParametrized1DSecondDerivative                !< the output value

        TaylorParametrized1DSecondDerivative = self%w2 +self%w3*x

    end function TaylorParametrized1DSecondDerivative

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the third derivative of the function.
    function TaylorParametrized1DThirdDerivative( self, x, eft_cache )

        implicit none

        class(taylor_parametrization_1D)                   :: self      !< the base class
        real(dl), intent(in)                               :: x         !< the input scale factor
        type(EFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: TaylorParametrized1DThirdDerivative                 !< the output value

        TaylorParametrized1DThirdDerivative = self%w3

    end function TaylorParametrized1DThirdDerivative

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the integral of the function, as defined in the notes.
    function TaylorParametrized1DIntegral( self, x, eft_cache )

        implicit none

        class(taylor_parametrization_1D)                   :: self      !< the base class
        real(dl), intent(in)                               :: x         !< the input scale factor
        type(EFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: TaylorParametrized1DIntegral                        !< the output value

        TaylorParametrized1DIntegral = x**(-1._dl -3._dl*self%w0)*Exp(-1._dl/12._dl*(x -1._dl)*(9._dl*(1._dl+x)*self%w2 +2._dl*(1._dl +x +x**2)*self%w3 +36._dl*self%wa))

    end function TaylorParametrized1DIntegral

    ! ---------------------------------------------------------------------------------------------

end module EFTCAMB_taylor_parametrizations_1D

!----------------------------------------------------------------------------------------
