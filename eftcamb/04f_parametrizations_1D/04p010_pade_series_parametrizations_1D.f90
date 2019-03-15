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

!> @file 04p010_pade_series_parametrizations_1D.f90
!! This file contains the definition of the Pade series parametrization, in the scale
!! factor, around an arbitrary expansion point. For best performances Horner’s algorithm is used.
!! Inherits from parametrized_function_1D.

!----------------------------------------------------------------------------------------
!> This module contains the definition of the Pade series parametrization, in the scale
!! factor, around an arbitrary expansion point. For best performances Horner’s algorithm is used.
!! Inherits from parametrized_function_1D.

!> @author Marco Raveri, Simone Peirone

module EFTCAMB_padeseries_parametrizations_1D

    use precision
    use AMLutils
    use IniFile
    use EFT_def
    use EFTCAMB_mixed_algorithms
    use EFTCAMB_cache
    use EFTCAMB_abstract_parametrizations_1D
    use EFTCAMB_taylorseries_parametrizations_1D

    implicit none

    private

    public padeseries_parametrization_1D

    ! ---------------------------------------------------------------------------------------------
    !> Type containing the Pade expansion parametrization. Inherits from parametrized_function_1D.
    type, extends ( parametrized_function_1D ) :: padeseries_parametrization_1D

        integer                               :: order_up    !< order of numerator polinomial in the Pade expansion
        integer                               :: order_down  !< order of denominator polinomial in the Pade expansion
        real(dl)                              :: a0          !< expansion point
        real(dl), dimension(:), allocatable   :: coeff_up    !< array with the coefficients of the numerator polinomial in the Pade expansion
        real(dl), dimension(:), allocatable   :: coeff_down  !< array with the coefficients of the denominator polinomial in the Pade expansion

        type(taylorseries_parametrization_1D) :: taylor_up   !< auxiliary Taylor expansion defining the numerator
        type(taylorseries_parametrization_1D) :: taylor_down !< auxiliary Taylor expansion defining the denominator

    contains

        ! utility functions:
        procedure :: set_param_number      => PadeSeriesParametrized1DSetParamNumber         !< subroutine that sets the number of parameters of the Pade expansion parametrized function.
        procedure :: init_parameters       => PadeSeriesParametrized1DInitParams             !< subroutine that initializes the function parameters based on the values found in an input array.
        procedure :: init_func_from_file   => PadeSeriesParametrized1DInitFunctionFromFile   !< subroutine that reads a Ini file looking for initialization parameters for the function..
        procedure :: parameter_value       => PadeSeriesParametrized1DParameterValues        !< subroutine that returns the value of the function i-th parameter.
        procedure :: feedback              => PadeSeriesParametrized1DFeedback               !< subroutine that prints to screen the informations about the function.
        procedure :: parameter_names       => PadeSeriesParametrized1DParameterNames         !< subroutine that returns the i-th parameter name of the function.
        procedure :: parameter_names_latex => PadeSeriesParametrized1DParameterNamesLatex    !< subroutine that returns the i-th parameter name of the function in latex format.

        ! evaluation procedures:
        procedure :: value                 => PadeSeriesParametrized1DValue               !< function that returns the value of the Pade expansion.
        procedure :: first_derivative      => PadeSeriesParametrized1DFirstDerivative     !< function that returns the first derivative of the Pade expansion.
        procedure :: second_derivative     => PadeSeriesParametrized1DSecondDerivative    !< function that returns the second derivative of the Pade expansion.
        procedure :: third_derivative      => PadeSeriesParametrized1DThirdDerivative     !< function that returns the third derivative of the Pade expansion.
        procedure :: integral              => PadeSeriesParametrized1DIntegral            !< function that returns the strange integral that we need for w_DE.

    end type padeseries_parametrization_1D

contains

    ! ---------------------------------------------------------------------------------------------
    ! Implementation of the Pade series parametrization.
    ! ---------------------------------------------------------------------------------------------

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that sets the number of parameters of the Pade expansion parametrized function.
    subroutine PadeSeriesParametrized1DSetParamNumber( self )

        implicit none

        class(Padeseries_parametrization_1D) :: self       !< the base class

        ! initialize the number of parameters:
        self%parameter_number = self%order_up +self%order_down +2

    end subroutine PadeSeriesParametrized1DSetParamNumber

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that initializes the function parameters based on the values found in an input array.
    subroutine PadeSeriesParametrized1DInitParams( self, array )

        implicit none

        class(Padeseries_parametrization_1D)                    :: self   !< the base class.
        real(dl), dimension(self%parameter_number), intent(in)  :: array  !< input array with the values of the parameters.

        integer  :: ind

        ! read in parameters:
        self%a0             = array(1)
        self%coeff_up(:)    = array(2:self%order_up+2)
        self%coeff_down(1)  = 1._dl
        self%coeff_down(2:) = array(self%order_up+3:)
        ! pass to the two Taylor series manually:
        self%taylor_up%a0      = self%a0
        self%taylor_down%a0    = self%a0
        self%taylor_up%coeff   = self%coeff_up
        self%taylor_down%coeff = self%coeff_down

    end subroutine PadeSeriesParametrized1DInitParams

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that reads a Ini file looking for initialization parameters for the function..
    subroutine PadeSeriesParametrized1DInitFunctionFromFile( self, Ini )

        implicit none

        class(Padeseries_parametrization_1D)      :: self   !< the base class
        type(TIniFile)                            :: Ini    !< Input ini file

        ! read order of the Pade expansion:
        self%order_up   = Ini_Read_Int_File( Ini, TRIM(self%name)//'_Pade_order_N', 0 )
        self%order_down = Ini_Read_Int_File( Ini, TRIM(self%name)//'_Pade_order_D', 0 )
        ! allocate:
        if ( allocated(self%coeff_up) ) deallocate(self%coeff_up)
        allocate( self%coeff_up( self%order_up+1 ) )
        self%coeff_up = 0._dl
        if ( allocated(self%coeff_down) ) deallocate(self%coeff_down)
        allocate( self%coeff_down( self%order_down+1 ) )
        self%coeff_down = 0._dl
        ! allocate the auxiliary Taylor series:
        self%taylor_up%order = self%order_up
        if ( allocated(self%taylor_up%coeff) ) deallocate(self%taylor_up%coeff)
        allocate( self%taylor_up%coeff( self%order_up+1 ) )
        self%taylor_up%coeff = 0._dl
        self%taylor_down%order = self%order_down
        if ( allocated(self%taylor_down%coeff) ) deallocate(self%taylor_down%coeff)
        allocate( self%taylor_down%coeff( self%order_down+1 ) )
        self%taylor_down%coeff = 0._dl
        ! ensure parameter number remains consistent:
        call self%set_param_number()

    end subroutine PadeSeriesParametrized1DInitFunctionFromFile

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the value of the function i-th parameter.
    subroutine PadeSeriesParametrized1DParameterValues( self, i, value )

        implicit none

        class(Padeseries_parametrization_1D)   :: self        !< the base class
        integer     , intent(in)               :: i           !< The index of the parameter
        real(dl)    , intent(out)              :: value       !< the output value of the i-th parameter

        if ( i==1 ) then
            value = self%a0
        else if ( i>1 .and. i<=self%order_up+2 ) then
            value = self%coeff_up( i-1 )
        else if ( i>self%order_up+2 .and. i<=self%parameter_number ) then
            value = self%coeff_down( i-self%order_up-1 )
        else
            write(*,*) 'Illegal index for parameter_names.'
            write(*,*) 'Maximum value is:', self%parameter_number
            call MpiStop('EFTCAMB error')
        end if

    end subroutine PadeSeriesParametrized1DParameterValues

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that prints to screen the informations about the function.
    subroutine PadeSeriesParametrized1DFeedback( self, print_params )

        implicit none

        class(Padeseries_parametrization_1D)    :: self         !< the base class
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

        write(*,'(a,a)')      'Pade series parametrization: ', self%name
        write(*,'(a,a,a,I3)') ' numerator order (',TRIM(self%name)//'_Pade_order_N','): ', self%order_up
        write(*,'(a,a,a,I3)') ' denominator order (',TRIM(self%name)//'_Pade_order_D','): ', self%order_down
        if ( print_params_temp ) then
            do i=1, self%parameter_number
                call self%parameter_names( i, param_name  )
                call self%parameter_value( i, param_value )
                write(*,'(a23,a,F12.6)') param_name, '=', param_value
            end do
        end if

    end subroutine PadeSeriesParametrized1DFeedback

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the i-th parameter name.
    subroutine PadeSeriesParametrized1DParameterNames( self, i, name )

        implicit none

        class(Padeseries_parametrization_1D)   :: self   !< the base class
        integer     , intent(in)               :: i      !< the index of the parameter
        character(*), intent(out)              :: name   !< the output name of the i-th parameter

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
                name = TRIM(self%name)//'a0'
            else if ( i>1 .and. i<=self%order_up+2 ) then
                name = TRIM(self%name)//'N'//integer_to_string(i-2)
            else if ( i>self%order_up+2 .and. i<=self%parameter_number ) then
                name = TRIM(self%name)//'D'//integer_to_string(i-self%order_up-2)
            end if
        end if

    end subroutine PadeSeriesParametrized1DParameterNames

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the latex version of the i-th parameter name.
    subroutine PadeSeriesParametrized1DParameterNamesLatex( self, i, latexname )

        implicit none

        class(Padeseries_parametrization_1D)   :: self       !< the base class
        integer     , intent(in)               :: i          !< the index of the parameter
        character(*), intent(out)              :: latexname  !< the output latex name of the i-th parameter

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
                latexname = TRIM(self%name)//'_'//'a_0'
            else if ( i>1 .and. i<=self%order_up+2 ) then
                latexname = TRIM(self%name)//'_N_'//integer_to_string(i-2)
            else if ( i>self%order_up+2 .and. i<=self%parameter_number ) then
                latexname = TRIM(self%name)//'_D_'//integer_to_string(i-self%order_up-2)
            end if
        end if

    end subroutine PadeSeriesParametrized1DParameterNamesLatex

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the value of the function in the scale factor.
    function PadeSeriesParametrized1DValue( self, x, eft_cache )

        implicit none

        class(Padeseries_parametrization_1D)               :: self      !< the base class
        real(dl), intent(in)                               :: x         !< the input scale factor
        type(EFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: PadeSeriesParametrized1DValue                     !< the output value

        real(dl) :: x_temp, up, down
        integer  :: ind

        up   = self%taylor_up%value( x, eft_cache )
        down = self%taylor_down%value( x, eft_cache )

        PadeSeriesParametrized1DValue = up/down

    end function PadeSeriesParametrized1DValue

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the value of the first derivative, wrt scale factor, of the function.
    function PadeSeriesParametrized1DFirstDerivative( self, x, eft_cache )

        implicit none

        class(Padeseries_parametrization_1D)               :: self      !< the base class
        real(dl), intent(in)                               :: x         !< the input scale factor
        type(EFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: PadeSeriesParametrized1DFirstDerivative           !< the output value

        real(dl) :: x_temp, up, dup, down, ddown
        integer  :: ind

        up    = self%taylor_up%value( x, eft_cache )
        dup   = self%taylor_up%first_derivative( x, eft_cache )
        down  = self%taylor_down%value( x, eft_cache )
        ddown = self%taylor_down%first_derivative( x, eft_cache )

        PadeSeriesParametrized1DFirstDerivative = ( dup*down - up*ddown )/(down)**2

    end function PadeSeriesParametrized1DFirstDerivative

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the second derivative of the function.
    function PadeSeriesParametrized1DSecondDerivative( self, x, eft_cache )

        implicit none

        class(Padeseries_parametrization_1D)               :: self      !< the base class
        real(dl), intent(in)                               :: x         !< the input scale factor
        type(EFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: PadeSeriesParametrized1DSecondDerivative                !< the output value

        real(dl) :: x_temp, up, dup, ddup, down, ddown, dddown
        integer  :: ind

        up     = self%taylor_up%value( x, eft_cache )
        dup    = self%taylor_up%first_derivative( x, eft_cache )
        ddup   = self%taylor_up%second_derivative( x, eft_cache )
        down   = self%taylor_down%value( x, eft_cache )
        ddown  = self%taylor_down%first_derivative( x, eft_cache )
        dddown = self%taylor_down%second_derivative( x, eft_cache )

        PadeSeriesParametrized1DSecondDerivative = (ddup*down**2 - 2._dl*ddown*down*dup + 2._dl*ddown**2*up - dddown*down*up)/down**3

    end function PadeSeriesParametrized1DSecondDerivative

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the third derivative of the function.
    function PadeSeriesParametrized1DThirdDerivative( self, x, eft_cache )

        implicit none

        class(Padeseries_parametrization_1D)               :: self      !< the base class
        real(dl), intent(in)                               :: x         !< the input scale factor
        type(EFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: PadeSeriesParametrized1DThirdDerivative           !< the output value

        real(dl) :: x_temp, up, dup, ddup, dddup, down, ddown, dddown, ddddown
        integer  :: ind

        up      = self%taylor_up%value( x, eft_cache )
        dup     = self%taylor_up%first_derivative( x, eft_cache )
        ddup    = self%taylor_up%second_derivative( x, eft_cache )
        dddup   = self%taylor_up%third_derivative( x, eft_cache )
        down    = self%taylor_down%value( x, eft_cache )
        ddown   = self%taylor_down%first_derivative( x, eft_cache )
        dddown  = self%taylor_down%second_derivative( x, eft_cache )
        ddddown = self%taylor_down%third_derivative( x, eft_cache )

        PadeSeriesParametrized1DThirdDerivative = (-3._dl*ddown*ddup*down**2 + dddup*down**3 + 6._dl*ddown**2*down*dup - 3._dl*dddown*down**2*dup - 6._dl*ddown**3*up + 6._dl*dddown*ddown*down*up - ddddown*down**2*up)/down**4

    end function PadeSeriesParametrized1DThirdDerivative

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the integral of the function, as defined in the notes.
    function PadeSeriesParametrized1DIntegral( self, x, eft_cache )

        implicit none

        class(Padeseries_parametrization_1D)               :: self      !< the base class
        real(dl), intent(in)                               :: x         !< the input scale factor
        type(EFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: PadeSeriesParametrized1DIntegral                  !< the output value

        PadeSeriesParametrized1DIntegral = 0._dl
        !< No analytic solution >!
        write(*,*) 'PadeSeriesParametrized1DIntegral is not implemented.'
        write(*,*) 'Calculations cannot proceed.'
        call MpiStop('EFTCAMB error')

    end function PadeSeriesParametrized1DIntegral

    ! ---------------------------------------------------------------------------------------------

end module EFTCAMB_padeseries_parametrizations_1D

!----------------------------------------------------------------------------------------
