!----------------------------------------------------------------------------------------
!
! This file is part of EFTCAMB.
!
! Copyright (C) 2013-2023 by the EFTCAMB authors
!
! The EFTCAMB code is free software;
! You can use it, redistribute it, and/or modify it under the terms
! of the GNU General Public License as published by the Free Software Foundation;
! either version 3 of the License, or (at your option) any later version.
! The full text of the license can be found in the file eftcamb/LICENSE at
! the top level of the EFTCAMB distribution.
!
!----------------------------------------------------------------------------------------

!> @file 04p015_axion_potential.f90
!! This file contains the definition of the cosine potential for Axions defined as:
!! f(x) = f0*( 1 -Cos(alpha*x) )^n
!! Inherits from parametrized_function_1D.

!----------------------------------------------------------------------------------------
!> This module contains the definition of the cosine potential for Axions defined as:
!! f(x) = f0*( 1 -Cos(alpha*x) )^n
!! Inherits from parametrized_function_1D.

!> @author Marco Raveri

module EFTCAMB_axion_parametrizations_1D

    use precision
    use MpiUtils
    use constants, only : const_pi, const_twopi, const_fourpi
    use EFT_def
    use EFTCAMB_mixed_algorithms
    use EFTCAMB_rootfind
    use linear_interpolation_1D
    use EFTCAMB_cache
    use EFTCAMB_abstract_parametrizations_1D

    implicit none

    private

    public axion_parametrization_1D

    ! ---------------------------------------------------------------------------------------------
    !> Type containing the Axion potential parametrization. Inherits from parametrized_function_1D.
    type, extends ( parametrized_function_1D ) :: axion_parametrization_1D

        real(dl)                           :: f0           !< amplitude of the potential
        real(dl)                           :: alpha        !< frequency of the cosine
        real(dl)                           :: n            !< exponent of the cosine

    contains

        ! utility functions:
        procedure :: set_param_number      => Axion1DSetParamNumber      !< subroutine that sets the number of parameters of the Axion function.
        procedure :: init_parameters       => Axion1DInitParams          !< subroutine that initializes the function parameters based on the values found in an input array.
        procedure :: parameter_value       => Axion1DParameterValues     !< subroutine that returns the value of the function i-th parameter.
        procedure :: feedback              => Axion1DFeedback            !< subroutine that prints to screen the informations about the function.

        ! evaluation procedures:
        procedure :: value                 => Axion1DValue               !< function that returns the value of the Axion.
        procedure :: first_derivative      => Axion1DFirstDerivative     !< function that returns the first derivative of the Axion.
        procedure :: second_derivative     => Axion1DSecondDerivative    !< function that returns the second derivative of the Axion.
        procedure :: third_derivative      => Axion1DThirdDerivative     !< function that returns the third derivative of the Axion.
        procedure :: fourth_derivative     => Axion1DFourthDerivative    !< function that returns the fourth derivative of the Axion.
        procedure :: integral              => Axion1DIntegral            !< function that returns the strange integral that we need for w_DE.

    end type axion_parametrization_1D

contains

    ! ---------------------------------------------------------------------------------------------
    ! Implementation of the Axion parametrization.
    ! ---------------------------------------------------------------------------------------------

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that sets the number of parameters of the Axion expansion parametrization.
    subroutine Axion1DSetParamNumber( self )

        implicit none

        class(axion_parametrization_1D)  :: self       !< the base class

        ! initialize the number of parameters:
        self%parameter_number = 3

    end subroutine Axion1DSetParamNumber

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that initializes the function parameters based on the values found in an input array.
    subroutine Axion1DInitParams( self, array )

        implicit none

        class(axion_parametrization_1D)                        :: self   !< the base class.
        real(dl), dimension(self%parameter_number), intent(in)  :: array  !< input array with the values of the parameters.

        integer  :: ind

        ! write parameters in:
        self%f0      = array(1)
        self%alpha   = array(2)
        self%n       = array(3)

    end subroutine Axion1DInitParams

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the value of the function i-th parameter.
    subroutine Axion1DParameterValues( self, i, value )

        implicit none

        class(axion_parametrization_1D)    :: self        !< the base class
        integer     , intent(in)            :: i           !< The index of the parameter
        real(dl)    , intent(out)           :: value       !< the output value of the i-th parameter

        if ( i == 1 ) then
            value = self%f0
        else if ( i == 2 ) then
            value = self%alpha
        else if ( i == 3 ) then
            value = self%n
        else
            write(*,*) 'Illegal index for parameter_names.'
            write(*,*) 'Maximum value is:', self%parameter_number
            call MpiStop('EFTCAMB error')
        end if

    end subroutine Axion1DParameterValues

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that prints to screen the informations about the function.
    subroutine Axion1DFeedback( self, print_params )

        implicit none

        class(axion_parametrization_1D)        :: self         !< the base class
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

        write(*,'(a,a)')      'Axion parametrization: ', self%name

        if ( print_params_temp ) then
            do i=1, self%parameter_number
                call self%parameter_names( i, param_name  )
                call self%parameter_value( i, param_value )
                write(*,'(a23,a,F12.6)') param_name, '=', param_value
            end do
        end if

    end subroutine Axion1DFeedback

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the value of the function in the scale factor.
    function Axion1DValue( self, x, eft_cache )

        implicit none

        class(axion_parametrization_1D)                   :: self      !< the base class
        real(dl), intent(in)                               :: x         !< the input scale factor
        type(TEFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: Axion1DValue                                      !< the output value

        Axion1DValue = self%f0*( 1._dl -Cos( self%alpha*x ) )**self%n

    end function Axion1DValue

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the value of the first derivative, wrt scale factor, of the function.
    function Axion1DFirstDerivative( self, x, eft_cache )

        implicit none

        class(axion_parametrization_1D)                   :: self      !< the base class
        real(dl), intent(in)                               :: x         !< the input scale factor
        type(TEFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: Axion1DFirstDerivative                            !< the output value

        Axion1DFirstDerivative = self%f0*self%n*self%alpha*(1._dl -Cos( self%alpha*x ))**(self%n-1._dl)*Sin(self%alpha*x)

    end function Axion1DFirstDerivative

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the second derivative of the function.
    function Axion1DSecondDerivative( self, x, eft_cache )

        implicit none

        class(axion_parametrization_1D)                   :: self      !< the base class
        real(dl), intent(in)                               :: x         !< the input scale factor
        type(TEFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: Axion1DSecondDerivative                           !< the output value

        Axion1DSecondDerivative = self%n*self%f0*self%alpha**2*( 1._dl-Cos(x*self%alpha))**(self%n-1._dl)*Cos(x*self%alpha) &
            & +(-1._dl+self%n)*self%n*self%f0*self%alpha**2*(1._dl-Cos(x*self%alpha))**(-2._dl+self%n)*Sin(x*self%alpha)**2

    end function Axion1DSecondDerivative

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the third derivative of the function.
    function Axion1DThirdDerivative( self, x, eft_cache )

        implicit none

        class(axion_parametrization_1D)                   :: self      !< the base class
        real(dl), intent(in)                               :: x         !< the input scale factor
        type(TEFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: Axion1DThirdDerivative                            !< the output value

        Axion1DThirdDerivative = -(self%n*self%f0*self%alpha**3*(1._dl-Cos(x*self%alpha))**(self%n-1._dl)*Sin(x*self%alpha)) &
            & +3._dl*(-1._dl+self%n)*self%n*self%f0*self%alpha**3*(1._dl-Cos(x*self%alpha))**(-2._dl+self%n)*Cos(x*self%alpha)*Sin(x*self%alpha) &
            & +(-2._dl+self%n)*(-1._dl+self%n)*self%n*self%f0*self%alpha**3*(1._dl-Cos(x*self%alpha))**(-3._dl+self%n)*Sin(x*self%alpha)**3

    end function Axion1DThirdDerivative

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the fourth derivative of the function.
    function Axion1DFourthDerivative( self, x, eft_cache )

        implicit none

        class(axion_parametrization_1D)                   :: self      !< the base class
        real(dl), intent(in)                               :: x         !< the input scale factor
        type(TEFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: Axion1DFourthDerivative                            !< the output value

        Axion1DFourthDerivative = 0._dl

        !<2._dl*self%alpha**4*self%n*(Sin(self%alpha*x/2._dl))**4*(1._dl-Cos*(self%alpha*x))**(-4._dl+self%n)&
        !<&*(self%n**3*Cos(2._dl*self%alpha*x)+2._dl*(2._dl*self%n*(-2._dl+self%n)*(-1._dl+self%n)-1._dl)*Cos(self%alpha*x)+(-2._dl+self%n)*(3._dl*self%n*(-2._dl+self%n)+2._dl))*self%f0

    end function Axion1DFourthDerivative

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the integral of the function, as defined in the notes.
    function Axion1DIntegral( self, x, eft_cache )

        implicit none

        class(axion_parametrization_1D)                   :: self      !< the base class
        real(dl), intent(in)                               :: x         !< the input scale factor
        type(TEFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: Axion1DIntegral                                   !< the output value

        Axion1DIntegral = 0._dl
        !< No analytic solution >!
        write(*,*) 'Axion1DIntegral is not implemented.'
        write(*,*) 'Calculations cannot proceed.'
        call MpiStop('EFTCAMB error')

    end function Axion1DIntegral

    ! ---------------------------------------------------------------------------------------------

end module EFTCAMB_axion_parametrizations_1D

!----------------------------------------------------------------------------------------
