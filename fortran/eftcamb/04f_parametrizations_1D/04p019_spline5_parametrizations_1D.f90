!----------------------------------------------------------------------------------------
!
! This file is part of EFTCAMB.
!
! Copyright (C) 2013-2020 by the EFTCAMB authors
!
! The EFTCAMB code is free software;
! You can use it, redistribute it, and/or modify it under the terms
! of the GNU General Public License as published by the Free Software Foundation;
! either version 3 of the License, or (at your option) any later version.
! The full text of the license can be found in the file eftcamb/LICENSE at
! the top level of the EFTCAMB distribution.
!
!----------------------------------------------------------------------------------------

!> @file 04p019_spline5_parametrizations_1D.f90
!! This file contains the definition of a splined 1D function, up to polynomial order 5.
!! Inherits from parametrized_function_1D.

!----------------------------------------------------------------------------------------
!> This module contains the definition of a splined 1D function, up to polynomial order 5.
!! Inherits from parametrized_function_1D.

!> @author Gen Ye

module EFTCAMB_spline5_parametrizations_1D

    use precision
    use IniObjects
    use MpiUtils
    use EFT_splines, only: spline3pars, poly3, dpoly3, d2poly3
    use splines5
    use EFT_def
    use EFTCAMB_mixed_algorithms
    use linear_interpolation_1D
    use EFTCAMB_cache
    use EFTCAMB_abstract_parametrizations_1D

    implicit none

    private

    public spline5_parametrization_1D

    ! ---------------------------------------------------------------------------------------------
    !> Type containing the spline parametrization, up to polynomial order 5. Inherits from parametrized_function_1D.
    type, extends ( parametrized_function_1D ) :: spline5_parametrization_1D

        logical                               :: logx_flag      !< If true we do interpolation in log(x) but return derivative wrt x.
        logical                               :: logy_flag      !< If true we do interpolation in log(y) but return derivative wrt y.
        integer                               :: spl_order      !< order of spline e.g: spl_order = 3 for cubic spline
        integer                               :: num_pixels     !< number of reconstruction points
        real(dl)                              :: null_value     !< the value to report when function called outside the interpolation range
        real(dl), dimension(:), allocatable   :: x_coord        !< array with the value of the coordinates where the reconstruction is defined
        real(dl), dimension(:), allocatable   :: pixels         !< array with the function value in each pixel

        ! internal variables:
        real(dl)                              :: x_initial      !< initial point of interpolation
        real(dl)                              :: x_final        !< ond point of interpolation
        real(dl), dimension(:,:), allocatable :: bcoef          !< internal coefficients for the spline interpolation

        real(dl), allocatable, dimension(:)   :: yint           !< array containing the values of the function w DE integral \f$ yint_i= \exp\left(-3\int_1^{x_i} \frac{1+f(x)}{x} \, dx \right) \f$.

    contains

        ! utility functions:
        procedure :: set_param_number      => SplineQ1DSetParamNumber         !< subroutine that sets the number of parameters of the reconstruction function.
        procedure :: init_parameters       => SplineQ1DInitParams             !< subroutine that initializes the function parameters based on the values found in an input array.
        procedure :: init_from_file        => SplineQ1DInitFromFile           !< subroutine that reads a Ini file looking for the parameters of the function.
        procedure :: init_func_from_file   => SplineQ1DInitFunctionFromFile   !< subroutine that reads a Ini file looking for initialization parameters for the function..
        procedure :: parameter_value       => SplineQ1DParameterValues        !< subroutine that returns the value of the function i-th parameter.
        procedure :: feedback              => SplineQ1DFeedback               !< subroutine that prints to screen the informations about the function.
        procedure :: parameter_names       => SplineQ1DParameterNames         !< subroutine that returns the i-th parameter name of the function.
        procedure :: parameter_names_latex => SplineQ1DParameterNamesLatex    !< subroutine that returns the i-th parameter name of the function in latex format.

        ! evaluation procedures:
        procedure :: value                 => SplineQ1DValue               !< function that returns the value of the reconstruction.
        procedure :: first_derivative      => SplineQ1DFirstDerivative     !< function that returns the first derivative of the reconstruction.
        procedure :: second_derivative     => SplineQ1DSecondDerivative    !< function that returns the second derivative of the reconstruction.
        procedure :: third_derivative      => SplineQ1DThirdDerivative     !< function that returns the third derivative of the reconstruction.
        procedure :: fourth_derivative     => SplineQ1DFourthDerivative    !< function that returns the fourth derivative of the reconstruction.
        procedure :: integral              => SplineQ1DIntegral            !< function that returns the strange integral that we need for w_DE.

    end type spline5_parametrization_1D

contains

    ! ---------------------------------------------------------------------------------------------
    ! Implementation of the reconstructioin parametrization.
    ! ---------------------------------------------------------------------------------------------

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that sets the number of parameters of the splined function.
    subroutine SplineQ1DSetParamNumber( self )

        implicit none

        class(spline5_parametrization_1D) :: self       !< the base class

        ! initialize the number of parameters:
        self%parameter_number = 2*self%num_pixels

    end subroutine SplineQ1DSetParamNumber

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that initializes the function parameters based on the values found in an input array.
    subroutine SplineQ1DInitParams( self, array )

        implicit none

        class(spline5_parametrization_1D)                    :: self   !< the base class.
        real(dl), dimension(self%parameter_number), intent(in)  :: array  !< input array with the values of the parameters.

        integer  :: ind
        ! real(dl) :: x

        self%x_coord   = array(1:self%num_pixels)
        if (self%logx_flag) then
            do ind = 1, self%num_pixels
                if (self%x_coord(ind)<=0) call MpiStop(TRIM(self%name)//' error: spline x coordinate must be non-negative when '//TRIM(self%name)//'_Logx_Flag is True.')
                self%x_coord(ind) = Log(self%x_coord(ind))
            end do
        end if
        self%pixels    = array(self%num_pixels+1:)
        if (self%logy_flag) then
            do ind = 1, self%num_pixels
                if (self%pixels(ind)<=0) call MpiStop(TRIM(self%name)//' error: spline y values must be non-negative when '//TRIM(self%name)//'_Logy_Flag is True.')
                self%pixels(ind) = Log(self%pixels(ind))
            end do
        end if
        self%x_initial = self%x_coord(1)
        self%x_final   = self%x_coord(self%num_pixels)
        ! check the coordinates:
        do ind = 1, self%num_pixels-1
          if (self%x_coord(ind)>=self%x_coord(ind+1)) call MpiStop(TRIM(self%name)//' error: spline x coordinate must be increasing. Pass correct values for '//TRIM(self%name)//'xN')
        end do
        ! create the spline coefficients:
        if ( allocated(self%bcoef)         ) deallocate( self%bcoef )
        if ( allocated(self%yint)          ) deallocate( self%yint )
        allocate( self%bcoef(self%spl_order+2, self%num_pixels-1) )
        ! compute spline polynomial coefficients
        select case ( self%spl_order )
            case (3)
                call spline3pars( self%x_coord, self%pixels, [2, 2], [0._dl, 0._dl], self%bcoef )
            case (5)
                call spline5pars( self%x_coord, self%pixels, [2, 2], [0._dl, 0._dl, 0._dl, 0._dl], self%bcoef )
        end select

        ! write(*,*) "-------------------------"
        ! do ind=1, self%num_pixels
        !     x = exp(self%x_coord(ind))
        !     write(*,*) x, self%pixels(ind), self%first_derivative(x), self%second_derivative(x), self%third_derivative(x), self%fourth_derivative(x)
        ! end do
        ! write(*,*) "=========================="
        
    end subroutine SplineQ1DInitParams

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that reads a Ini file looking for the parameters of the function.
    subroutine SplineQ1DInitFromFile( self, Ini, eft_error )

        implicit none

        class(spline5_parametrization_1D)   :: self      !< the base class
        type(TIniFile)                         :: Ini       !< Input ini file
        integer                                :: eft_error !< error code: 0 all fine, 1 initialization failed

        character(len=EFT_names_max_length)          :: param_name
        real(dl), dimension( self%parameter_number ) :: parameters

        integer  :: i

        ! ensure that the number of parameters is properly associated:
        call self%set_param_number()

        ! read the parameters and store them in a vector:
        do i=1, self%parameter_number
            call self%parameter_names( i, param_name )
            if ( i<=self%num_pixels ) then
              parameters(i) = Ini%Read_Double( TRIM(param_name), 1._dl*real(i-1)/real(self%num_pixels-1) )
            else
              parameters(i) = Ini%Read_Double( TRIM(param_name), 0._dl )
            end if
        end do

        ! initialize the function parameters from the vector:
        call self%init_parameters( parameters )

    end subroutine SplineQ1DInitFromFile

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that reads a Ini file looking for initialization parameters for the function..
    subroutine SplineQ1DInitFunctionFromFile( self, Ini, eft_error )

        implicit none

        class(spline5_parametrization_1D)  :: self      !< the base class
        type(TIniFile)                        :: Ini       !< Input ini file
        integer                               :: eft_error !< error code: 0 all fine, 1 initialization failed

        ! read spline configuration
        self%logx_flag = Ini%Read_Logical( TRIM(self%name)//'_Logx_Flag', .false. )
        self%logy_flag = Ini%Read_Logical( TRIM(self%name)//'_Logy_Flag', .false. )
        ! read spline order
        self%spl_order = Ini%Read_Int( TRIM(self%name)//'_Spline_Order', 5 )
        if ( self%spl_order/=3 .and. self%spl_order/=5 ) then
            write(*,*) TRIM(self%name)//'_Spline_Order =', self%spl_order, 'not supported. Use 3 or 5.'
            eft_error = 1
            return
        end if
        ! read number of points:
        self%num_pixels   = Ini%Read_Int( TRIM(self%name)//'_Spline_Pixels', 6 )
        if ( self%num_pixels < self%spl_order+1 ) then
            write(*,*) 'Spline requires ', TRIM(self%name)//'_Spline_Pixels to be no smaller than ', TRIM(self%name)//'_Spline_Order + 1 = ', self%spl_order+1
            eft_error = 1
            return
        end if
        self%null_value   = Ini%Read_Double( TRIM(self%name)//'_null_value', 0._dl )
        ! allocate:
        if ( allocated(self%x_coord) ) deallocate( self%x_coord )
        allocate( self%x_coord(self%num_pixels) )
        if ( allocated(self%pixels) ) deallocate( self%pixels )
        allocate( self%pixels(self%num_pixels) )
        ! ensure parameter number remains consistent:
        call self%set_param_number()

    end subroutine SplineQ1DInitFunctionFromFile

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the value of the function i-th parameter.
    subroutine SplineQ1DParameterValues( self, i, value )

        implicit none

        class(spline5_parametrization_1D)   :: self        !< the base class
        integer     , intent(in)               :: i           !< The index of the parameter
        real(dl)    , intent(out)              :: value       !< the output value of the i-th parameter

        if ( i<=self%num_pixels ) then
            value = self%x_coord( i )
        else if ( i>self%num_pixels .and. i<=self%parameter_number ) then
            value = self%pixels( i-self%num_pixels )
        else
            write(*,*) 'Illegal index for parameter_names.'
            write(*,*) 'Maximum value is:', self%parameter_number
            call MpiStop('EFTCAMB error')
        end if

    end subroutine SplineQ1DParameterValues

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that prints to screen the informations about the function.
    subroutine SplineQ1DFeedback( self, print_params )

        implicit none

        class(spline5_parametrization_1D)     :: self         !< the base class
        logical, optional                        :: print_params !< optional flag that decised whether to print numerical values
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

        write(*,'(a,a)')      ' alternative spline parametrization: ', self%name
        write(*,'(a,a)')      ' whether interpolate in terms of log(x): ', self%logx_flag
        write(*,'(a,a)')      ' whether interpolate in terms of log(y): ', self%logy_flag
        write(*,'(a,a,a,I3)') ' number of pixels (',TRIM(self%name)//'_Spline_Pixels','): ', self%num_pixels
        write(*,'(a,a,a,I3)') ' spline order (',TRIM(self%name)//'_Spline_Order','): ', self%spl_order
        write(*,'(a,F12.6)')  ' null value: ', self%null_value
        if ( print_params_temp ) then
            do i=1, self%parameter_number
                call self%parameter_names( i, param_name  )
                call self%parameter_value( i, param_value )
                write(*,'(a23,a,F12.6)') param_name, '=', param_value
            end do
        end if

    end subroutine SplineQ1DFeedback

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the i-th parameter name.
    subroutine SplineQ1DParameterNames( self, i, name )

        implicit none

        class(spline5_parametrization_1D)   :: self   !< the base class
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
            if ( i<=self%num_pixels ) then
                name = TRIM(self%name)//'x'//integer_to_string(i)
            else if ( i>self%num_pixels .and. i<=self%parameter_number ) then
                name = TRIM(self%name)//'v'//integer_to_string(i-self%num_pixels)
            end if
        end if

    end subroutine SplineQ1DParameterNames

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the latex version of the i-th parameter name.
    subroutine SplineQ1DParameterNamesLatex( self, i, latexname )

        implicit none

        class(spline5_parametrization_1D)       :: self       !< the base class
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
        if ( allocated(self%param_names_latex) ) then
            latexname = self%param_names_latex(i)%string
        else
            if ( i<=self%num_pixels ) then
                latexname = TRIM(self%name)//'x_'//integer_to_string(i)
            else if ( i>self%num_pixels .and. i<=self%parameter_number ) then
                latexname = TRIM(self%name)//'v_'//integer_to_string(i-self%num_pixels)
            end if
        end if

    end subroutine SplineQ1DParameterNamesLatex

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the value of the function in the scale factor.
    function SplineQ1DValue( self, x, eft_cache )

        implicit none

        class(spline5_parametrization_1D)                :: self      !< the base class
        real(dl), intent(in)                                :: x         !< the input scale factor
        type(TEFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: SplineQ1DValue                                     !< the output value

        integer  :: ip
        real(dl) :: a, x_
        real(dl) :: y, f

        a  = x
        x_ = a
        if (self%logx_flag) then
            x_ = Log(a)
        end if
        ! Do not extrapolate:
        SplineQ1DValue = self%null_value
        if ( x_ < self%x_initial .or. x_ > self%x_final ) return
        ! if inside return the spline value:
        ip = iixmin(x_, self%x_coord, 0 )
        select case ( self%spl_order )
        case (3)
            y = poly3( x_, self%bcoef(:, ip) )
        case (5)
            y = poly5( x_, self%bcoef(:, ip) )
        end select
        f = y
        if (self%logy_flag) then
            f = exp(y)
        end if
        SplineQ1DValue = f

    end function SplineQ1DValue

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the value of the first derivative, wrt scale factor, of the function.
    function SplineQ1DFirstDerivative( self, x, eft_cache )

        implicit none

        class(spline5_parametrization_1D)                :: self      !< the base class
        real(dl), intent(in)                                :: x         !< the input scale factor
        type(TEFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: SplineQ1DFirstDerivative                           !< the output value

        integer  :: ip
        real(dl) :: a, x_
        real(dl) :: y, f, yx, fx
        
        a  = x
        x_ = a
        if (self%logx_flag) then
            x_ = Log(a)
        end if
        ! Do not extrapolate:
        SplineQ1DFirstDerivative = self%null_value
        if ( x_ < self%x_initial .or. x_ > self%x_final ) return
        ! if inside return the spline value:
        ip = iixmin(x_, self%x_coord, 0 )
        select case ( self%spl_order )
        case (3)
            if (self%logx_flag .or. self%logy_flag) then
                y = poly3( x_, self%bcoef(:, ip) )
            end if
            yx = dpoly3( x_, self%bcoef(:, ip) )
        case (5)
            if (self%logx_flag .or. self%logy_flag) then
                y = poly5( x_, self%bcoef(:, ip) )
            end if
            yx = dpoly5( x_, self%bcoef(:, ip) )
        end select
        ! convert interpolation function y(x) to f(a) 
        f  = y
        fx = yx
        if (self%logy_flag) then
            f  = exp(y)
            fx = yx*f
        end if
        SplineQ1DFirstDerivative = fx
        if (self%logx_flag) then
            SplineQ1DFirstDerivative = fx/a
        end if


    end function SplineQ1DFirstDerivative

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the second derivative of the function.
    function SplineQ1DSecondDerivative( self, x, eft_cache )

        implicit none

        class(spline5_parametrization_1D)                :: self      !< the base class
        real(dl), intent(in)                                :: x         !< the input scale factor
        type(TEFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: SplineQ1DSecondDerivative                          !< the output value

        integer  :: ip
        real(dl) :: a, x_
        real(dl) :: y, f, yx, fx, yxx, fxx
        
        a  = x
        x_ = a
        if (self%logx_flag) then
            x_ = Log(a)
        end if
        ! Do not extrapolate:
        SplineQ1DSecondDerivative = self%null_value
        if ( x_ < self%x_initial .or. x_ > self%x_final ) return
        ! if inside return the spline value:
        ip = iixmin(x_, self%x_coord, 0 )
        select case ( self%spl_order )
        case (3)
            if (self%logx_flag .or. self%logy_flag) then
                y  = poly3( x_, self%bcoef(:, ip) )
                yx = dpoly3( x_, self%bcoef(:, ip) )
            end if
            yxx = d2poly3( x_, self%bcoef(:, ip) )
        case (5)
            if (self%logx_flag .or. self%logy_flag) then
                y  = poly5( x_, self%bcoef(:, ip) )
                yx = dpoly5( x_, self%bcoef(:, ip) )
            end if
            yxx = d2poly5( x_, self%bcoef(:, ip) )
        end select
        ! convert interpolation function y(x) to f(a) 
        f   = y
        fx  = yx
        fxx = yxx
        if (self%logy_flag) then
            f   = exp(y)
            fx  = yx*f
            fxx = (yxx + yx**2)*f
        end if
        SplineQ1DSecondDerivative = fxx
        if (self%logx_flag) then
            SplineQ1DSecondDerivative = (fxx - fx)/a**2
        end if

    end function SplineQ1DSecondDerivative

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the third derivative of the function.
    function SplineQ1DThirdDerivative( self, x, eft_cache )

        implicit none

        class(spline5_parametrization_1D)                :: self      !< the base class
        real(dl), intent(in)                                :: x         !< the input scale factor
        type(TEFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: SplineQ1DThirdDerivative                           !< the output value

        integer  :: ip
        real(dl) :: a, x_
        real(dl) :: y, f, yx, fx, yxx, fxx, yxxx, fxxx

        a  = x
        x_ = a
        if (self%logx_flag) then
            x_ = Log(a)
        end if
        ! Do not extrapolate:
        SplineQ1DThirdDerivative = self%null_value
        if ( x_ < self%x_initial .or. x_ > self%x_final ) return
        ! if inside return the spline value:
        ip = iixmin(x_, self%x_coord, 0 )
        select case ( self%spl_order )
        case (3)
            write(*,'(a,a,a,I3)') 'Third derivative invalid for ',TRIM(self%name)//'_Spline_Order = ', self%spl_order
            write(*,*) 'Calculations cannot proceed.'
            call MpiStop('EFTCAMB error')
        case (5)
            if (self%logx_flag .or. self%logy_flag) then
                y   = poly5( x_, self%bcoef(:, ip) )
                yx  = dpoly5( x_, self%bcoef(:, ip) )
                yxx = d2poly5( x_, self%bcoef(:, ip) )
            end if
            yxxx = d3poly5( x_, self%bcoef(:, ip) )
        end select
        ! convert interpolation function y(x) to f(a) 
        f    = y
        fx   = yx
        fxx  = yxx
        fxxx = yxxx
        if (self%logy_flag) then
            f    = exp(y)
            fx   = yx*f
            fxx  = (yxx + yx**2)*f
            fxxx = (yxxx + 3._dl*yx*yxx + yx**3)*f
        end if
        SplineQ1DThirdDerivative = fxxx
        if (self%logx_flag) then
            SplineQ1DThirdDerivative = (fxxx - 3._dl*fxx + 2._dl*fx)/a**3
        end if

    end function SplineQ1DThirdDerivative

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the fourth derivative of the function.
    function SplineQ1DFourthDerivative( self, x, eft_cache )

        implicit none

        class(spline5_parametrization_1D)                :: self      !< the base class
        real(dl), intent(in)                                :: x         !< the input scale factor
        type(TEFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: SplineQ1DFourthDerivative                          !< the output value

        integer  :: ip
        real(dl) :: a, x_
        real(dl) :: y, f, yx, fx, yxx, fxx, yxxx, fxxx, yxxxx, fxxxx
        
        a  = x
        x_ = a
        if (self%logx_flag) then
            x_ = Log(a)
        end if
        ! Do not extrapolate:
        SplineQ1DFourthDerivative = self%null_value
        if ( x_ < self%x_initial .or. x_ > self%x_final ) return
        ! if inside return the spline value:
        ip = iixmin(x_, self%x_coord, 0 )
        select case ( self%spl_order )
        case (3)
            write(*,'(a,a,a,I3)') 'Fourth derivative invalid for ',TRIM(self%name)//'_Spline_Order = ', self%spl_order
            write(*,*) 'Calculations cannot proceed.'
            call MpiStop('EFTCAMB error')
        case (5)
            if (self%logx_flag .or. self%logy_flag) then
                y    = poly5( x_, self%bcoef(:, ip) )
                yx   = dpoly5( x_, self%bcoef(:, ip) )
                yxx  = d2poly5( x_, self%bcoef(:, ip) )
                yxxx = d3poly5( x_, self%bcoef(:, ip) )
            end if
            yxxxx = d4poly5( x_, self%bcoef(:, ip) )
        end select
        ! convert interpolation function y(x) to f(a) 
        f     = y
        fx    = yx
        fxx   = yxx
        fxxx  = yxxx
        fxxxx = yxxxx
        if (self%logy_flag) then
            f     = exp(y)
            fx    = yx*f
            fxx   = (yxx + yx**2)*f
            fxxx  = (yxxx + 3._dl*yx*yxx + yx**3)*f
            fxxxx = (yxxxx + 6._dl*yx**2*yxx + 3._dl*yxx**2 + 4._dl*yx*yxxx + yx**4)*f
        end if
        SplineQ1DFourthDerivative = fxxxx
        if (self%logx_flag) then
            SplineQ1DFourthDerivative = (fxxxx - 6._dl*fxxx + 11._dl*fxx - 6._dl*fx)/a**4
        end if

    end function SplineQ1DFourthDerivative

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the integral of the function, as defined in the notes. Only meaningful for w_DE(a), by Gen
    function SplineQ1DIntegral( self, x, eft_cache )

        implicit none

        class(spline5_parametrization_1D)                :: self      !< the base class
        real(dl), intent(in)                                :: x         !< the input scale factor
        type(TEFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: SplineQ1DIntegral                                  !< the output value

        real(dl) :: tmp, a, x_, x0, dx
        integer  :: ip, j

        SplineQ1DIntegral = 0._dl

        ! initialize integration vector when first called
        if ( .not. allocated(self%yint) ) then
            call SplineQ1DInitIntegration( self )
        end if

        a  = x
        x_ = a
        if (self%logx_flag) then
            x_ = Log(a)
        end if
        ! Do not extrapolate:
        if ( x_ < self%x_initial ) then
            SplineQ1DIntegral = self%yint(1)
        elseif ( x_ > self%x_final ) then
            SplineQ1DIntegral = 1._dl
        else
            ! simply use linear interpolation
            ip = iixmin(x_, self%x_coord, 0 )
            x0 = self%x_coord(ip)
            dx = x_ - x0
            tmp = 0._dl
            do j=1, size( self%bcoef(:, ip) )
                tmp = tmp + self%bcoef(j,ip)*dx**j/real(j)
            end do
            SplineQ1DIntegral = self%yint(ip)*(x_/x0)**3*exp(tmp)
        end if

    end function SplineQ1DIntegral

    ! ---------------------------------------------------------------------------------------------
    !> Private subroutine to initialize the integration vecotr. \f$ yint_i= \exp\left(-3\int_1^{x_i} \frac{1+f(x)}{x} \, dx \right) \f$. Only meaningful when we are parameterizing w_DE(a). by Gen
    subroutine SplineQ1DInitIntegration( self )

        implicit none

        class(spline5_parametrization_1D)                    :: self      !< the base class

        integer                              :: i, j
        real(dl)                             :: x1, x2, dx, tmp

        if ( self%logy_flag ) then
            write(*,*) 'Spline integral is not implemented for ', TRIM(self%name)//'_Logy_Flag = True.'
            write(*,*) 'Calculations cannot proceed.'
            call MpiStop('EFTCAMB error')
        end if
        if ( .not. self%logx_flag ) then
            write(*,*) 'Spline integral is not implemented for ', TRIM(self%name)//'_Logx_Flag = False.'
            write(*,*) 'Calculations cannot proceed.'
            call MpiStop('EFTCAMB error')
        end if
        
        allocate( self%yint(self%num_pixels) )
        
        tmp = 0._dl
        self%yint(self%num_pixels) = 1._dl
        do i=self%num_pixels-1, 1, -1
            x1 = self%x_coord(i)
            x2 = self%x_coord(i+1)
            dx = x2 - x1
            do j=1, size( self%bcoef(:, i) )
                tmp = tmp + self%bcoef(j,i)*dx**j/real(j)
            end do
            self%yint(i) = self%yint(i+1)*(x1/x2)**3*exp(-tmp)
        end do
    
    end subroutine SplineQ1DInitIntegration

end module EFTCAMB_spline5_parametrizations_1D

!----------------------------------------------------------------------------------------