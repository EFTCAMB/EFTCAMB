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

!> @file 04p011_fourier_parametrizations_1D.f90
!! This file contains the definition of the fourier series parametrization.
!! Inherits from parametrized_function_1D.

!----------------------------------------------------------------------------------------
!> This module contains the definition of the fourier series parametrization.
!! Inherits from parametrized_function_1D.

!> @author Marco Raveri

module EFTCAMB_fourier_parametrizations_1D

    use precision
    use IniObjects
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

    public fourier_parametrization_1D

    ! ---------------------------------------------------------------------------------------------
    !> Type containing the Fourier expansion parametrization. Inherits from parametrized_function_1D.
    type, extends ( parametrized_function_1D ) :: fourier_parametrization_1D

        integer                               :: order_cos    !< order of the cosine expansion
        integer                               :: order_sin    !< order of the sine expansion
        real(dl), dimension(:), allocatable   :: coeff_cos    !< array with the coefficients of the cosine expansion
        real(dl), dimension(:), allocatable   :: coeff_sin    !< array with the coefficients of the sine expansion
        real(dl)                              :: a0           !< first order term (a constant)
        real(dl)                              :: phase        !< an overall phase shift
        logical                               :: min_a0       !< whether a0 is the minimum over a period of the series or the average over a period
        real(dl)                              :: offset       !< the offset value that makes a0 the minimum or the DC value

    contains

        ! utility functions:
        procedure :: set_param_number      => Fourier1DSetParamNumber         !< subroutine that sets the number of parameters of the Fourier function.
        procedure :: init_parameters       => Fourier1DInitParams             !< subroutine that initializes the function parameters based on the values found in an input array.
        procedure :: init_func_from_file   => Fourier1DInitFunctionFromFile   !< subroutine that reads a Ini file looking for initialization parameters for the function..
        procedure :: parameter_value       => Fourier1DParameterValues        !< subroutine that returns the value of the function i-th parameter.
        procedure :: feedback              => Fourier1DFeedback               !< subroutine that prints to screen the informations about the function.
        procedure :: parameter_names       => Fourier1DParameterNames         !< subroutine that returns the i-th parameter name of the function.
        procedure :: parameter_names_latex => Fourier1DParameterNamesLatex    !< subroutine that returns the i-th parameter name of the function in latex format.
        procedure :: minimum               => FourierMinimum                  !< finds the minimum value of the series

        ! evaluation procedures:
        procedure :: value                 => Fourier1DValue               !< function that returns the value of the Fourier.
        procedure :: first_derivative      => Fourier1DFirstDerivative     !< function that returns the first derivative of the Fourier.
        procedure :: second_derivative     => Fourier1DSecondDerivative    !< function that returns the second derivative of the Fourier.
        procedure :: third_derivative      => Fourier1DThirdDerivative     !< function that returns the third derivative of the Fourier.
        procedure :: fourth_derivative     => Fourier1DFourthDerivative    !< function that returns the fourth derivative of the Fourier.
        procedure :: integral              => Fourier1DIntegral            !< function that returns the strange integral that we need for w_DE.

    end type fourier_parametrization_1D

contains

    ! ---------------------------------------------------------------------------------------------
    ! Implementation of the Fourier parametrization.
    ! ---------------------------------------------------------------------------------------------

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that sets the number of parameters of the Fourier expansion parametrization.
    subroutine Fourier1DSetParamNumber( self )

        implicit none

        class(fourier_parametrization_1D) :: self       !< the base class

        ! initialize the number of parameters:
        self%parameter_number = 2 ! phase and offset
        self%parameter_number = self%parameter_number +self%order_cos
        self%parameter_number = self%parameter_number +self%order_sin

    end subroutine Fourier1DSetParamNumber

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that initializes the function parameters based on the values found in an input array.
    subroutine Fourier1DInitParams( self, array )

        implicit none

        class(fourier_parametrization_1D)                       :: self   !< the base class.
        real(dl), dimension(self%parameter_number), intent(in)  :: array  !< input array with the values of the parameters.

        integer  :: ind

        ! write parameters in:
        self%a0        = array(1)
        self%phase     = array(2)
        self%coeff_cos = array( 2+1:2+self%order_cos )
        self%coeff_sin = array( 2+self%order_cos+1:2+self%order_cos+self%order_sin )
        ! find the minimum if needed and initialize offset so that a0 is the minimum value:
        self%offset = 0._dl
        if ( self%min_a0 ) then
            self%a0        = 0._dl
            self%offset    = -self%minimum()
            self%a0        = array(1)
        end if

    end subroutine Fourier1DInitParams

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that reads a Ini file looking for initialization parameters for the function..
    subroutine Fourier1DInitFunctionFromFile( self, Ini, eft_error )

        implicit none

        class(fourier_parametrization_1D)  :: self      !< the base class
        type(TIniFile)                     :: Ini       !< Input ini file
        integer                            :: eft_error !< error code: 0 all fine, 1 initialization failed

        ! read order of the expansion:
        self%order_cos = Ini%Read_Int( TRIM(self%name)//'_Fourier_order_cos', 0 )
        self%order_sin = Ini%Read_Int( TRIM(self%name)//'_Fourier_order_sin', 0 )
        self%min_a0    = Ini%Read_Logical( TRIM(self%name)//'_min_a0', .False.  )
        ! allocate:
        if ( allocated(self%coeff_cos) ) deallocate( self%coeff_cos )
        allocate( self%coeff_cos(self%order_cos) )
        if ( allocated(self%coeff_sin) ) deallocate( self%coeff_sin )
        allocate( self%coeff_sin(self%order_sin) )
        ! ensure parameter number remains consistent:
        call self%set_param_number()

    end subroutine Fourier1DInitFunctionFromFile

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the value of the function i-th parameter.
    subroutine Fourier1DParameterValues( self, i, value )

        implicit none

        class(fourier_parametrization_1D)   :: self        !< the base class
        integer     , intent(in)            :: i           !< The index of the parameter
        real(dl)    , intent(out)           :: value       !< the output value of the i-th parameter

        if ( i == 1 ) then
            value = self%a0
        else if ( i == 2 ) then
            value = self%phase
        else if ( i>2 .and. i<3+self%order_cos ) then
            value = self%coeff_cos( i-2 )
        else if ( i>2+self%order_cos .and. i<=self%parameter_number ) then
            value = self%coeff_sin( i-2-self%order_cos )
        else
            write(*,*) 'Illegal index for parameter_names.'
            write(*,*) 'Maximum value is:', self%parameter_number
            call MpiStop('EFTCAMB error')
        end if

    end subroutine Fourier1DParameterValues

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that prints to screen the informations about the function.
    subroutine Fourier1DFeedback( self, print_params )

        implicit none

        class(fourier_parametrization_1D)       :: self         !< the base class
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

        write(*,'(a,a)')      'Fourier parametrization: ', self%name
        write(*,'(a,a,a,I3)') ' order of the cos expansion (',TRIM(self%name)//'_Fourier_order_cos','): ', self%order_cos
        write(*,'(a,a,a,I3)') ' order of the sin expansion (',TRIM(self%name)//'_Fourier_order_sin','): ', self%order_sin
        if ( self%min_a0 ) then
            write(*,'(a)') ' The parameter '//TRIM(self%name)//'_a0 is the function minimum'
        end if

        if ( print_params_temp ) then
            if ( self%min_a0 ) then
                write(*,'(a23,a,F12.6)') '   Offset              ', '=', self%offset
            end if
            do i=1, self%parameter_number
                call self%parameter_names( i, param_name  )
                call self%parameter_value( i, param_value )
                write(*,'(a23,a,F12.6)') param_name, '=', param_value
            end do
        end if

    end subroutine Fourier1DFeedback

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the i-th parameter name.
    subroutine Fourier1DParameterNames( self, i, name )

        implicit none

        class(fourier_parametrization_1D)          :: self   !< the base class
        integer     , intent(in)                   :: i      !< the index of the parameter
        character(*), intent(out)                  :: name   !< the output name of the i-th parameter

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
                name = TRIM(self%name)//'_a0'
            else if ( i==2 ) then
                name = TRIM(self%name)//'_phase'
            else if ( i>2 .and. i<3+self%order_cos ) then
                name = TRIM(self%name)//'_a'//integer_to_string(i-2)
            else if ( i>2+self%order_cos .and. i<=self%parameter_number ) then
                name = TRIM(self%name)//'_b'//integer_to_string(i-2-self%order_cos)
            end if
        end if

    end subroutine Fourier1DParameterNames

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the latex version of the i-th parameter name.
    subroutine Fourier1DParameterNamesLatex( self, i, latexname )

        implicit none

        class(fourier_parametrization_1D)          :: self       !< the base class
        integer     , intent(in)                   :: i          !< the index of the parameter
        character(*), intent(out)                  :: latexname  !< the output latex name of the i-th parameter

        ! check the input index:
        if ( i>self%parameter_number ) then
            write(*,*) 'In parametrized_function_1D:', self%name
            write(*,*) 'Illegal index for parameter_names.'
            write(*,*) 'Maximum value is:', self%parameter_number
            call MpiStop('EFTCAMB error')
        end if
        ! return the parameter name:
        if ( allocated(self%param_names) ) then
            latexname = self%param_names(i)%string
        else
            if ( i==1 ) then
                latexname = TRIM(self%name)//'_a_0'
            else if ( i==2 ) then
                latexname = TRIM(self%name)//'_phase'
            else if ( i>2 .and. i<3+self%order_cos ) then
                latexname = TRIM(self%name)//'_a_'//integer_to_string(i-2)
            else if ( i>2+self%order_cos .and. i<=self%parameter_number ) then
                latexname = TRIM(self%name)//'_b_'//integer_to_string(i-2-self%order_cos)
            end if
        end if

    end subroutine Fourier1DParameterNamesLatex

    ! ---------------------------------------------------------------------------------------------
    !> Function that computes the minimum of the Fourier series
    function FourierMinimum( self )

        implicit none

        class(fourier_parametrization_1D) :: self  !< the base class

        real(dl) :: FourierMinimum
        real(dl) :: min_sample, max_sample, global_minimum, initial_value, final_value
        real(dl) :: func_zero, new_min
        integer  :: num_points, i
        logical  :: success

        ! 1) decide how fine to presample the value grid over a period:
        min_sample = self%phase
        max_sample = 1._dl +self%phase
        num_points = 4*2*max( self%order_cos, self%order_sin ) +1 ! try to be well above the Shannon rate

        ! 2) search for the roots of the first derivative and look for the global minimum:
        global_minimum = 1.d16 ! just initialize to a big number
        do i=1, num_points-1
            initial_value = min_sample +REAL(i-1)/REAL(num_points-1)*( max_sample -min_sample )
            final_value   = min_sample +REAL(i)/REAL(num_points-1)*( max_sample -min_sample )
            func_zero     = zbrent( helper_min, initial_value, final_value, 1.d-10, 0._dl, success)
            if ( success ) then
                new_min = self%value( func_zero )
                if ( new_min < global_minimum ) then
                    global_minimum = new_min
                end if
            end if
        end do

        ! 3) return the value:
        FourierMinimum = global_minimum

    contains

        ! small interface function to feed the root finder:
        function helper_min(x)
            implicit none
            real(dl) :: x, helper_min
            helper_min = self%first_derivative( x )
        end function helper_min

    end function FourierMinimum

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the value of the function in the scale factor.
    function Fourier1DValue( self, x, eft_cache )

        implicit none

        class(fourier_parametrization_1D)                  :: self      !< the base class
        real(dl), intent(in)                               :: x         !< the input scale factor
        type(TEFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: Fourier1DValue                                      !< the output value

        real(dl) :: x_temp
        integer  :: ind

        ! phase:
        x_temp = x-self%phase
        ! offset:
        Fourier1DValue = self%a0 +self%offset
        ! cosine sum:
        do ind=1, self%order_cos
            Fourier1DValue = Fourier1DValue +self%coeff_cos(ind)*Cos( const_twopi*REAL(ind)*x_temp )
        end do
        ! sine sum:
        do ind=1, self%order_sin
            Fourier1DValue = Fourier1DValue +self%coeff_sin(ind)*Sin( const_twopi*REAL(ind)*x_temp )
        end do

    end function Fourier1DValue

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the value of the first derivative, wrt scale factor, of the function.
    function Fourier1DFirstDerivative( self, x, eft_cache )

        implicit none

        class(fourier_parametrization_1D)                  :: self      !< the base class
        real(dl), intent(in)                               :: x         !< the input scale factor
        type(TEFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: Fourier1DFirstDerivative                            !< the output value

        real(dl) :: x_temp
        integer  :: ind

        ! phase:
        x_temp = x-self%phase
        Fourier1DFirstDerivative = 0._dl
        ! cosine sum:
        do ind=1, self%order_cos
            Fourier1DFirstDerivative = Fourier1DFirstDerivative -const_twopi*REAL(ind)*self%coeff_cos(ind)*Sin( const_twopi*REAL(ind)*x_temp )
        end do
        ! sine sum:
        do ind=1, self%order_sin
            Fourier1DFirstDerivative = Fourier1DFirstDerivative +const_twopi*REAL(ind)*self%coeff_sin(ind)*Cos( const_twopi*REAL(ind)*x_temp )
        end do

    end function Fourier1DFirstDerivative

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the second derivative of the function.
    function Fourier1DSecondDerivative( self, x, eft_cache )

        implicit none

        class(fourier_parametrization_1D)                  :: self      !< the base class
        real(dl), intent(in)                               :: x         !< the input scale factor
        type(TEFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: Fourier1DSecondDerivative                           !< the output value

        real(dl) :: x_temp
        integer  :: ind

        ! phase:
        x_temp = x-self%phase
        Fourier1DSecondDerivative = 0._dl
        ! cosine sum:
        do ind=1, self%order_cos
            Fourier1DSecondDerivative = Fourier1DSecondDerivative -(const_twopi*REAL(ind))**2*self%coeff_cos(ind)*Cos( const_twopi*REAL(ind)*x_temp )
        end do
        ! sine sum:
        do ind=1, self%order_sin
            Fourier1DSecondDerivative = Fourier1DSecondDerivative -(const_twopi*REAL(ind))**2*self%coeff_sin(ind)*Sin( const_twopi*REAL(ind)*x_temp )
        end do

    end function Fourier1DSecondDerivative

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the third derivative of the function.
    function Fourier1DThirdDerivative( self, x, eft_cache )

        implicit none

        class(fourier_parametrization_1D)                  :: self      !< the base class
        real(dl), intent(in)                               :: x         !< the input scale factor
        type(TEFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: Fourier1DThirdDerivative                            !< the output value

        real(dl) :: x_temp
        integer  :: ind

        ! phase:
        x_temp = x-self%phase
        Fourier1DThirdDerivative = 0._dl
        ! cosine sum:
        do ind=1, self%order_cos
            Fourier1DThirdDerivative = Fourier1DThirdDerivative +(const_twopi*REAL(ind))**3*self%coeff_cos(ind)*Sin( const_twopi*REAL(ind)*x_temp )
        end do
        ! sine sum:
        do ind=1, self%order_sin
            Fourier1DThirdDerivative = Fourier1DThirdDerivative -(const_twopi*REAL(ind))**3*self%coeff_sin(ind)*Cos( const_twopi*REAL(ind)*x_temp )
        end do

    end function Fourier1DThirdDerivative

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the fourth derivative of the function.
    function Fourier1DFourthDerivative( self, x, eft_cache )

        implicit none

        class(fourier_parametrization_1D)                  :: self      !< the base class
        real(dl), intent(in)                               :: x         !< the input scale factor
        type(TEFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: Fourier1DFourthDerivative                            !< the output value

        real(dl) :: x_temp
        integer  :: ind

        ! phase:
        x_temp = x-self%phase
        Fourier1DFourthDerivative = 0._dl
        ! cosine sum:
        do ind=1, self%order_cos
            Fourier1DFourthDerivative = Fourier1DFourthDerivative +(const_twopi*REAL(ind))**4*self%coeff_cos(ind)*Cos( const_twopi*REAL(ind)*x_temp )
        end do
        ! sine sum:
        do ind=1, self%order_sin
            Fourier1DFourthDerivative = Fourier1DFourthDerivative +(const_twopi*REAL(ind))**4*self%coeff_sin(ind)*Sin( const_twopi*REAL(ind)*x_temp )
        end do

    end function Fourier1DFourthDerivative

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the integral of the function, as defined in the notes.
    function Fourier1DIntegral( self, x, eft_cache )

        implicit none

        class(fourier_parametrization_1D)                  :: self      !< the base class
        real(dl), intent(in)                               :: x         !< the input scale factor
        type(TEFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: Fourier1DIntegral                                   !< the output value

        Fourier1DIntegral = 0._dl
        !< No analytic solution >!
        write(*,*) 'Fourier1DIntegral is not implemented.'
        write(*,*) 'Calculations cannot proceed.'
        call MpiStop('EFTCAMB error')

    end function Fourier1DIntegral

    ! ---------------------------------------------------------------------------------------------

end module EFTCAMB_fourier_parametrizations_1D

!----------------------------------------------------------------------------------------
