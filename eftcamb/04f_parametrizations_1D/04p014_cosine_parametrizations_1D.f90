!----------------------------------------------------------------------------------------
!
! This file is part of EFTCAMB.
!
! Copyright (C) 2013-2019 by the EFTCAMB authors
!
! The EFTCAMB code is free software;
! You can use it, redistribute it, and/or modify it under the terms
! of the GNU General Public License as published by the Free Software Foundation;
! either version 3 of the License, or (at your option) any later version.
! The full text of the license can be found in the file eftcamb/LICENSE at
! the top level of the EFTCAMB distribution.
!
!----------------------------------------------------------------------------------------

!> @file 04p017_cosine_parametrizations_1D.f90
!! This file contains the definition of the cosine parametrization defined as:
!! f(x) = f0*( Cos(alpha*x+beta) +gamma)
!! Inherits from parametrized_function_1D.

!----------------------------------------------------------------------------------------
!> This module contains the definition of the cosine parametrization defined as:
!! f(x) = f0*( Cos(alpha*x+beta) +gamma)
!! Inherits from parametrized_function_1D.

!> @author Marco Raveri

module EFTCAMB_cosine_parametrizations_1D

    use precision
    use constants, only : const_pi, const_twopi, const_fourpi
    use AMLutils
    use IniFile
    use EFT_def
    use EFTCAMB_mixed_algorithms
    use EFTCAMB_rootfind
    use linear_interpolation_1D
    use EFTCAMB_cache
    use EFTCAMB_abstract_parametrizations_1D

    implicit none

    private

    public cosine_parametrization_1D

    ! ---------------------------------------------------------------------------------------------
    !> Type containing the Cosine expansion parametrization. Inherits from parametrized_function_1D.
    type, extends ( parametrized_function_1D ) :: cosine_parametrization_1D

        real(dl)                           :: f0           !< amplitude of the cosine
        real(dl)                           :: alpha        !< frequency of the cosine
        real(dl)                           :: beta         !< phase shift of the cosine
        real(dl)                           :: gamma        !< constant offset

    contains

        ! utility functions:
        procedure :: set_param_number      => Cosine1DSetParamNumber         !< subroutine that sets the number of parameters of the Cosine function.
        procedure :: init_parameters       => Cosine1DInitParams             !< subroutine that initializes the function parameters based on the values found in an input array.
        procedure :: parameter_value       => Cosine1DParameterValues        !< subroutine that returns the value of the function i-th parameter.
        procedure :: feedback              => Cosine1DFeedback               !< subroutine that prints to screen the informations about the function.

        ! evaluation procedures:
        procedure :: value                 => Cosine1DValue               !< function that returns the value of the Cosine.
        procedure :: first_derivative      => Cosine1DFirstDerivative     !< function that returns the first derivative of the Cosine.
        procedure :: second_derivative     => Cosine1DSecondDerivative    !< function that returns the second derivative of the Cosine.
        procedure :: third_derivative      => Cosine1DThirdDerivative     !< function that returns the third derivative of the Cosine.
        procedure :: integral              => Cosine1DIntegral            !< function that returns the strange integral that we need for w_DE.

    end type cosine_parametrization_1D

contains

    ! ---------------------------------------------------------------------------------------------
    ! Implementation of the Cosine parametrization.
    ! ---------------------------------------------------------------------------------------------

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that sets the number of parameters of the Cosine expansion parametrization.
    subroutine Cosine1DSetParamNumber( self )

        implicit none

        class(cosine_parametrization_1D)  :: self       !< the base class

        ! initialize the number of parameters:
        self%parameter_number = 4

    end subroutine Cosine1DSetParamNumber

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that initializes the function parameters based on the values found in an input array.
    subroutine Cosine1DInitParams( self, array )

        implicit none

        class(cosine_parametrization_1D)                        :: self   !< the base class.
        real(dl), dimension(self%parameter_number), intent(in)  :: array  !< input array with the values of the parameters.

        integer  :: ind

        ! write parameters in:
        self%f0      = array(1)
        self%alpha   = array(2)
        self%beta    = array(3)
        self%gamma   = array(4)

    end subroutine Cosine1DInitParams

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the value of the function i-th parameter.
    subroutine Cosine1DParameterValues( self, i, value )

        implicit none

        class(cosine_parametrization_1D)    :: self        !< the base class
        integer     , intent(in)            :: i           !< The index of the parameter
        real(dl)    , intent(out)           :: value       !< the output value of the i-th parameter

        if ( i == 1 ) then
            value = self%f0
        else if ( i == 2 ) then
            value = self%alpha
        else if ( i == 3 ) then
            value = self%beta
        else if ( i == 4 ) then
            value = self%gamma
        else
            write(*,*) 'Illegal index for parameter_names.'
            write(*,*) 'Maximum value is:', self%parameter_number
            call MpiStop('EFTCAMB error')
        end if

    end subroutine Cosine1DParameterValues

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that prints to screen the informations about the function.
    subroutine Cosine1DFeedback( self, print_params )

        implicit none

        class(cosine_parametrization_1D)        :: self         !< the base class
        logical, optional                       :: print_params !< optional flag that decised whether to print numerical values
                                                                !! of the parameters.

        integer                                 :: i
        real(dl)                                :: param_value
        character(len=EFT_names_max_length)     :: param_name
        logical                                 :: print_params_temp

        if ( present(print_params) ) then
            print_params_temp = print_params
        else
            print_params_temp = .True.
        end if

        write(*,'(a,a)')      'Cosine parametrization: ', self%name

        if ( print_params_temp ) then
            do i=1, self%parameter_number
                call self%parameter_names( i, param_name  )
                call self%parameter_value( i, param_value )
                write(*,'(a23,a,F12.6)') param_name, '=', param_value
            end do
        end if

    end subroutine Cosine1DFeedback

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the value of the function in the scale factor.
    function Cosine1DValue( self, x, eft_cache )

        implicit none

        class(cosine_parametrization_1D)                   :: self      !< the base class
        real(dl), intent(in)                               :: x         !< the input scale factor
        type(EFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: Cosine1DValue                                      !< the output value

        Cosine1DValue = self%f0*( Cos( self%alpha*x +self%beta  ) +self%gamma )

    end function Cosine1DValue

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the value of the first derivative, wrt scale factor, of the function.
    function Cosine1DFirstDerivative( self, x, eft_cache )

        implicit none

        class(cosine_parametrization_1D)                   :: self      !< the base class
        real(dl), intent(in)                               :: x         !< the input scale factor
        type(EFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: Cosine1DFirstDerivative                            !< the output value

        Cosine1DFirstDerivative = -self%alpha*self%f0*Sin( self%beta + self%alpha*x )

    end function Cosine1DFirstDerivative

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the second derivative of the function.
    function Cosine1DSecondDerivative( self, x, eft_cache )

        implicit none

        class(cosine_parametrization_1D)                   :: self      !< the base class
        real(dl), intent(in)                               :: x         !< the input scale factor
        type(EFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: Cosine1DSecondDerivative                           !< the output value

        Cosine1DSecondDerivative = -self%alpha**2*self%f0*Cos( self%beta + self%alpha*x )

    end function Cosine1DSecondDerivative

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the third derivative of the function.
    function Cosine1DThirdDerivative( self, x, eft_cache )

        implicit none

        class(cosine_parametrization_1D)                   :: self      !< the base class
        real(dl), intent(in)                               :: x         !< the input scale factor
        type(EFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: Cosine1DThirdDerivative                            !< the output value

        Cosine1DThirdDerivative = self%alpha**3*self%f0*Sin( self%beta + self%alpha*x )

    end function Cosine1DThirdDerivative

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the integral of the function, as defined in the notes.
    function Cosine1DIntegral( self, x, eft_cache )

        implicit none

        class(cosine_parametrization_1D)                   :: self      !< the base class
        real(dl), intent(in)                               :: x         !< the input scale factor
        type(EFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: Cosine1DIntegral                                   !< the output value

        Cosine1DIntegral = 0._dl
        !< No analytic solution >!
        write(*,*) 'Cosine1DIntegral is not implemented.'
        write(*,*) 'Calculations cannot proceed.'
        call MpiStop('EFTCAMB error')

    end function Cosine1DIntegral

    ! ---------------------------------------------------------------------------------------------

end module EFTCAMB_cosine_parametrizations_1D

!----------------------------------------------------------------------------------------
