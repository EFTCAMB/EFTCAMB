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

!> @file 04p017_step_parametrizations_1D.f90
!! This file contains the definition of the step parametrization in scale factor,
!! inheriting from parametrized_function_1D.

!----------------------------------------------------------------------------------------
!> This module contains the definition of the step parametrization,
!! inheriting from parametrized_function_1D.

!> @author Marco Raveri, Meng-Xiang Lin

module EFTCAMB_step_parametrizations_1D

    use precision
    use MpiUtils
    use EFT_def
    use EFTCAMB_cache
    use EFTCAMB_abstract_parametrizations_1D

    implicit none

    private

    public step_parametrization_1D

    ! ---------------------------------------------------------------------------------------------
    !> Type containing the step function parametrization. Inherits from parametrized_function_1D.
    type, extends ( parametrized_function_1D ) :: step_parametrization_1D

        real(dl) :: v1      ! The value of the function when a is -infinity (in cosmology we use a=0)
        real(dl) :: v2      ! The value of the function when a is infinity (in cosmology we use a=1)
        real(dl) :: aT      ! The time of transition
        real(dl) :: Delta   ! The width of transition

    contains

        ! utility functions:
        procedure :: set_param_number      => stepParametrized1DSetParamNumber      !< subroutine that sets the number of parameters of the step parametrized function.
        procedure :: init_parameters       => stepParametrized1DInitParams          !< subroutine that initializes the function parameters based on the values found in an input array.
        procedure :: parameter_value       => stepParametrized1DParameterValues     !< subroutine that returns the value of the function i-th parameter.
        procedure :: feedback              => stepParametrized1DFeedback            !< subroutine that prints to screen the informations about the function.
        procedure :: parameter_names       => stepParameterNames                    !< subroutine that returns the i-th parameter name of the function.
        procedure :: parameter_names_latex => stepParameterNamesLatex               !< subroutine that returns the i-th parameter name of the function in latex format.

        ! evaluation procedures:
        procedure :: value                 => stepParametrized1DValue               !< function that returns the value of the step function.
        procedure :: first_derivative      => stepParametrized1DFirstDerivative     !< function that returns the first derivative of the step function.
        procedure :: second_derivative     => stepParametrized1DSecondDerivative    !< function that returns the second derivative of the step function.
        procedure :: third_derivative      => stepParametrized1DThirdDerivative     !< function that returns the third derivative of the step function.
        procedure :: fourth_derivative     => stepParametrized1DFourthDerivative    !< function that returns the fourth derivative of the step function.
        procedure :: integral              => stepParametrized1DIntegral            !< function that returns the strange integral that we need for w_DE.

    end type step_parametrization_1D

contains

    ! ---------------------------------------------------------------------------------------------
    ! Implementation of the step function.
    ! ---------------------------------------------------------------------------------------------

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that sets the number of parameters of the step parametrized function.
    subroutine stepParametrized1DSetParamNumber( self )

        implicit none

        class(step_parametrization_1D) :: self       !< the base class

        ! initialize the number of parameters:
        self%parameter_number = 4

    end subroutine stepParametrized1DSetParamNumber

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that initializes the function parameters based on the values found in an input array.
    subroutine stepParametrized1DInitParams( self, array )

        implicit none

        class(step_parametrization_1D)                           :: self   !< the base class.
        real(dl), dimension(self%parameter_number), intent(in)  :: array  !< input array with the values of the parameters.

        self%v1    = array(1)
        self%v2    = array(2)
        self%aT    = array(3)
        self%Delta = array(4)

    end subroutine stepParametrized1DInitParams

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the value of the function i-th parameter.
    subroutine stepParametrized1DParameterValues( self, i, value )

        implicit none

        class(step_parametrization_1D)        :: self        !< the base class
        integer     , intent(in)             :: i           !< The index of the parameter
        real(dl)    , intent(out)            :: value       !< the output value of the i-th parameter

        select case (i)
            case(1)
                value = self%v1
            case(2)
                value = self%v2
            case(3)
                value = self%aT
            case(4)
                value = self%Delta
            case default
                write(*,*) 'Illegal index for parameter_names.'
                write(*,*) 'Maximum value is:', self%parameter_number
                call MpiStop('EFTCAMB error')
        end select

    end subroutine stepParametrized1DParameterValues

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that prints to screen the informations about the function.
    subroutine stepParametrized1DFeedback( self, print_params )

        implicit none

        class(step_parametrization_1D)   :: self         !< the base class
        logical, optional                :: print_params !< optional flag that decised whether to print numerical values
                                                         !! of the parameters.

        integer                               :: i
        real(dl)                              :: param_value
        character(len=EFT_names_max_length)   :: param_name
        logical                               :: print_params_temp

        if ( present(print_params) ) then
            print_params_temp = print_params
        else
            print_params_temp = .True.
        end if

        write(*,*)     'step parametrization: ', self%name
        if ( print_params_temp ) then
            do i=1, self%parameter_number
                call self%parameter_names( i, param_name  )
                call self%parameter_value( i, param_value )
                write(*,'(a23,a,F12.6)') param_name, '=', param_value
            end do
        end if

    end subroutine stepParametrized1DFeedback

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the i-th parameter name.
    subroutine stepParameterNames( self, i, name )

        implicit none

        class(step_parametrization_1D)     :: self   !< the base class
        integer     , intent(in)           :: i      !< the index of the parameter
        character(*), intent(out)          :: name   !< the output name of the i-th parameter

        ! check the input index:
        if ( i>self%parameter_number ) then
            write(*,*) 'In parametrized_function_1D:', self%name
            write(*,*) 'Illegal index for parameter_names.'
            write(*,*) 'Maximum value is:', self%parameter_number
            call MpiStop('EFTCAMB error')
        end if
        ! return the parameter name:
        if ( allocated(self%param_names) ) then
            name = self%param_names(i)%string
        else
            if ( i==1 ) then
                name = TRIM(self%name)//'_v1'
            else if ( i==2 ) then
                name = TRIM(self%name)//'_v2'
            else if ( i==3 ) then
                name = TRIM(self%name)//'_at'
            else if ( i==4 ) then
                name = TRIM(self%name)//'_delta'
            end if
        end if

    end subroutine stepParameterNames

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the latex version of the i-th parameter name.
    subroutine stepParameterNamesLatex( self, i, latexname )

        implicit none

        class(step_parametrization_1D)      :: self       !< the base class
        integer     , intent(in)            :: i          !< the index of the parameter
        character(*), intent(out)           :: latexname  !< the output latex name of the i-th parameter

        ! check the input index:
        if ( i>self%parameter_number ) then
            write(*,*) 'In parametrized_function_1D:', self%name
            write(*,*) 'Illegal index for parameter_names.'
            write(*,*) 'Maximum value is:', self%parameter_number
            call MpiStop('EFTCAMB error')
        end if
        ! return the parameter name:
        if ( allocated(self%param_names_latex) ) then
            latexname = self%param_names_latex(i)%string
        else
            if ( i==1 ) then
                latexname = TRIM(self%name_latex)//'_1'
            else if ( i==2 ) then
                latexname = TRIM(self%name_latex)//'_2'
            else if ( i==3 ) then
                latexname = TRIM(self%name_latex)//'_{a_T}'
            else if ( i==4 ) then
                latexname = TRIM(self%name_latex)//'_{\delta}'
            end if
        end if

    end subroutine stepParameterNamesLatex

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the value of the function in the scale factor.
    function stepParametrized1DValue( self, x, eft_cache )

        implicit none

        class(step_parametrization_1D)                     :: self      !< the base class
        real(dl), intent(in)                               :: x         !< the input scale factor
        type(TEFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: stepParametrized1DValue                              !< the output value

        stepParametrized1DValue = 0.5_dl*(self%v2-self%v1)*((x-self%aT)/self%Delta) /sqrt(1+((x-self%aT)/self%Delta)**2) + 0.5_dl*(self%v1+self%v2)

        if(isnan(stepParametrized1DValue)) then
            write(*,*)  'a:', x
            write(*,*)  'params:', self%v1, self%v2, self%aT, self%Delta
        end if

    end function stepParametrized1DValue

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the value of the first derivative, wrt scale factor, of the function.
    function stepParametrized1DFirstDerivative( self, x, eft_cache )

        implicit none

        class(step_parametrization_1D)                      :: self      !< the base class
        real(dl), intent(in)                               :: x         !< the input scale factor
        type(TEFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: stepParametrized1DFirstDerivative                    !< the output value

        stepParametrized1DFirstDerivative = 0.5_dl*(self%v2-self%v1)/self%Delta /(1+((x-self%aT)/self%Delta)**2)**1.5_dl

    end function stepParametrized1DFirstDerivative

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the second derivative of the function.
    function stepParametrized1DSecondDerivative( self, x, eft_cache )

        implicit none

        class(step_parametrization_1D)                      :: self      !< the base class
        real(dl), intent(in)                               :: x         !< the input scale factor
        type(TEFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: stepParametrized1DSecondDerivative                   !< the output value

        stepParametrized1DSecondDerivative = 1.5_dl*(x-self%aT)*(self%v1-self%v2)*self%Delta /sqrt(1+((x-self%aT)/self%Delta)**2) /((x-self%aT)**2+self%Delta**2)**2

    end function stepParametrized1DSecondDerivative

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the third derivative of the function.
    function stepParametrized1DThirdDerivative( self, x, eft_cache )

        implicit none

        class(step_parametrization_1D)                      :: self      !< the base class
        real(dl), intent(in)                               :: x         !< the input scale factor
        type(TEFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: stepParametrized1DThirdDerivative                    !< the output value

        stepParametrized1DThirdDerivative = 1.5_dl*(self%v1-self%v2)*self%Delta*(-4*(x-self%aT)**2+self%Delta**2) /sqrt(1+((x-self%aT)/self%Delta)**2) /((x-self%aT)**2+self%Delta**2)**3

    end function stepParametrized1DThirdDerivative

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the fourth derivative of the function.
    function stepParametrized1DFourthDerivative( self, x, eft_cache )

        implicit none

        class(step_parametrization_1D)                      :: self      !< the base class
        real(dl), intent(in)                               :: x         !< the input scale factor
        type(TEFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: stepParametrized1DFourthDerivative                    !< the output value

        stepParametrized1DFourthDerivative = (15._dl/2._dl)*(self%v2-self%v1)*(self%aT-x)*(4._dl*(self%aT-x)**2 - 3._dl*self%Delta**2)&
        &*(1._dl/(self%Delta**7 * ((self%aT-x)**2 + self%Delta**2)**(4.5_dl)))

    end function stepParametrized1DFourthDerivative

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the integral of the function, as defined in the notes.
    function stepParametrized1DIntegral( self, x, eft_cache )

        implicit none

        class(step_parametrization_1D)                      :: self      !< the base class
        real(dl), intent(in)                               :: x         !< the input scale factor
        type(TEFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: stepParametrized1DIntegral                           !< the output value
        real(dl) :: term1, term2, term3, term4

        if (x==0) then
            stepParametrized1DIntegral = 0._dl
        else
            term1 = -asinh(self%aT/self%Delta) - asinh((x-self%aT)/self%Delta)
            term2 = log((1-self%aT+sqrt((-1+self%aT)**2+self%Delta**2))*(self%aT+sqrt(self%aT**2+self%Delta**2))/self%Delta**2)
            term3 = self%aT*(log(((-1+self%aT)*self%aT+self%Delta**2+sqrt(((-1+self%aT)**2+self%Delta**2)*(self%aT**2+self%Delta**2)))/(-x*self%aT+self%aT**2+self%Delta**2+sqrt(((x-self%aT)**2+self%Delta**2)*(self%aT**2+self%Delta**2)))))/sqrt(self%aT**2+self%Delta**2)
            term4 = 1.5_dl*(self%v2-self%v1) * self%aT/sqrt(self%aT**2+self%Delta**2)

            stepParametrized1DIntegral = x**(-1._dl-1.5_dl*(self%v1+self%v2)+term4) *Exp( 1.5_dl*(self%v2-self%v1)*(term1+term2+term3) )
        end if

        if (isnan(stepParametrized1DIntegral)) then
            write(*,*)  'a:', x
            write(*,*)  'params:', self%v1, self%v2, self%aT, self%Delta
            write(*,*)  'term1:', term1
            write(*,*)  'term2:', term2
            write(*,*)  'term3:', term3
            write(*,*)  'term4:', term4
        end if

    end function stepParametrized1DIntegral

    ! ---------------------------------------------------------------------------------------------

end module EFTCAMB_step_parametrizations_1D

!----------------------------------------------------------------------------------------
