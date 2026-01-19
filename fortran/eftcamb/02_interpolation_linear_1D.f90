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

!> @file 02_interpolation_linear_1D.f90
!! This file contains the class that can be used for 1D linearly interpolated functions
!! on a general grid.


!----------------------------------------------------------------------------------------
!> This module contains the class that can be used for 1D linearly interpolated functions
!! on a general grid.

!> @author Bin Hu, Marco Raveri

module linear_interpolation_1D

    use precision
    use EFT_splines
    use EFTCAMB_mixed_algorithms
    use MpiUtils

    implicit none

    private

    public linear_interpolate_function_1D

    !----------------------------------------------------------------------------------------
    !> This is the type that can be used for the 1D linear interpolation.
    type :: linear_interpolate_function_1D

        ! options:
        integer                             :: num_points     !< number of points of the interpolating function.

        ! parameters:
        real(dl)                            :: x_initial      !< first value of x.
        real(dl)                            :: x_final        !< last value of x.
        real(dl)                            :: null_value     !< value that is returned if a point outside the range of interpolation is requested.
        logical                             :: has_null_value !< wether to use the null value outside interpolation range. If no null value is passed to initialization assume that function is constant outside interpolation range.

        ! arrays with the values:
        real(dl), allocatable, dimension(:) :: x              !< array containing the values of x.
        real(dl), allocatable, dimension(:) :: y              !< array containing the values of the function \f$ y_i=f(x_i) \f$.
        real(dl), allocatable, dimension(:) :: yp             !< array containing the values of the function derivative \f$ yp_i= \frac{d f(x_i)}{dx} \f$.
        real(dl), allocatable, dimension(:) :: ypp            !< array containing the values of the function second derivative \f$ ypp_i= \frac{d^2 f(x_i)}{dx^2} \f$.
        real(dl), allocatable, dimension(:) :: yppp           !< array containing the values of the function third derivative \f$ yppp_i= \frac{d^3 f(x_i)}{dx^3} \f$.
        real(dl), allocatable, dimension(:) :: ypppp          !< array containing the values of the function fourth derivative \f$ ypppp_i=\frac{d^3 f(x_i)}{dx^4} \f$.
        real(dl), allocatable, dimension(:) :: yint           !< array containing the values of the function w DE integral \f$ yint_i= \exp\left(-3\int_1^{x_i} \frac{1+f(x)}{x} \, dx \right) \f$.

    contains

        procedure :: initialize             => LinearIntepolateFunction1DInit               !< subroutine that initialize the interpolating function.
        procedure :: precompute             => LinearIntepolateFunction1DPrecompute         !< subroutine that does precomputations for the interpolation. Usefull when calling values and derivatives.
        procedure :: value                  => LinearIntepolateFunction1DValue              !< function that gives the value of the function at a given coordinate x.
        procedure :: first_derivative       => LinearIntepolateFunction1DFirstDerivative    !< function that gives the value of the function first derivative at a given coordinate x.
        procedure :: second_derivative      => LinearIntepolateFunction1DSecondDerivative   !< function that gives the value of the function second derivative at a given coordinate x.
        procedure :: third_derivative       => LinearIntepolateFunction1DThirdDerivative    !< function that gives the value of the function third derivative at a given coordinate x.
        procedure :: fourth_derivative      => LinearIntepolateFunction1DFourthDerivative   !< function that gives the value of the function fourth derivative at a given coordinate x.
        procedure :: integral               => LinearIntepolateFunction1DIntegral           !< function that gives the value of the interpolated w DE integral at a given coordinate x.
        procedure :: initialize_derivatives => LinearIntepolateFunction1DInitDerivatives    !< subroutine that initializes the derivatives if the derivatives vectors are not initialized. The derivative are inferred from the function itself.
        procedure :: initialize_integration => LinearIntepolateFunction1DInitIntegration    !< subroutine that initializes the integration if the integration vector is not initialized. The integration is inferred from the function itself.

    end type linear_interpolate_function_1D

    ! ---------------------------------------------------------------------------------------------

contains

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that initialize the interpolating function.
    subroutine LinearIntepolateFunction1DInit( self, x_values, null_value )

        implicit none

        class(linear_interpolate_function_1D) :: self         !< the base class

        real(dl), intent(in), dimension(:)    :: x_values     !< first value of x. Optional.
        real(dl), intent(in), optional        :: null_value   !< value that is returned if a point outside the range of interpolation is requested. Optional. If not passed to the constructor it is assumed to be zero.

        integer  :: i

        ! initialize the null value:
        if ( present(null_value) ) then
            self%has_null_value = .true.
            self%null_value     = null_value
        else
            self%has_null_value = .false.
            self%null_value     = 0._dl
        end if

        ! allocate the x vector:
        if ( allocated(self%x) ) deallocate( self%x )

        self%num_points = SIZE(x_values)

        allocate( self%x, source=x_values )

        ! check the ordering of the x array:
        do i=1, self%num_points-1
            if (self%x(i) > self%x(i+1)) then
                write(*,*) 'ERROR: the x vector is not ordered'
                call MpiStop('EFTCAMB error')
            end if
        end do

        ! store initial and final times:
        self%x_initial  = self%x(1)
        self%x_final    = self%x(self%num_points)

        ! allocate the other vectors:
        if ( allocated(self%y)    ) deallocate( self%y    )
        if ( allocated(self%yp)   ) deallocate( self%yp   )
        if ( allocated(self%ypp)  ) deallocate( self%ypp  )
        if ( allocated(self%yppp) ) deallocate( self%yppp )
        if ( allocated(self%ypppp) ) deallocate( self%ypppp )
        if ( allocated(self%yint) ) deallocate( self%yint )

        allocate( self%y   ( self%num_points ) )
        allocate( self%yp  ( self%num_points ) )
        allocate( self%ypp ( self%num_points ) )
        allocate( self%yppp( self%num_points ) )
        allocate( self%ypppp( self%num_points ) )
        allocate( self%yint( self%num_points ) )

    end subroutine LinearIntepolateFunction1DInit

    ! ---------------------------------------------------------------------------------------------
    !> Function that computes the main interpolation index.
    !! Usefull when calling value and derivatives in the same place.
    subroutine LinearIntepolateFunction1DPrecompute( self, x, ind, mu )

        implicit none

        class(linear_interpolate_function_1D)  :: self       !< the base class
        real(dl), intent(in)                   :: x          !< the value of x at which the function is required
        integer , intent(out)                  :: ind        !< the main interpolation index
        real(dl), intent(out)                  :: mu         !< the interpolation coefficient

        real(dl) :: x1, x2

        ! compute the interpolation index:
        call hunt( self%x, self%num_points, x, ind)
        ! store the x values:
        x1  = self%x(ind)
        x2  = self%x(ind+1)
        ! compute the linear interpolation coefficient:
        mu  = (x-x1)/(x2-x1)

    end subroutine LinearIntepolateFunction1DPrecompute

    ! ---------------------------------------------------------------------------------------------
    !> Function that gives the value of the function at a given coordinate x.
    function LinearIntepolateFunction1DValue( self, x, index, coeff )

        implicit none

        class(linear_interpolate_function_1D)          :: self       !< the base class
        real(dl), intent(in)                           :: x          !< the value of x at which the function is required
        integer , intent(in), optional                 :: index      !< optional precomputed value of the interpolation index
        real(dl), intent(in), optional                 :: coeff      !< optional precomputed value of the interpolation coefficient
        real(dl) :: LinearIntepolateFunction1DValue                  !< the output value of the function
        integer  :: ind
        real(dl) :: x1, x2, y1, y2, mu

        ! initialize to null value:
        LinearIntepolateFunction1DValue = self%null_value
        if ( self%has_null_value ) then
            ! if outside the interpolation range return the null value:
            if ( x <= self%x_initial .or. x >= self%x_final ) return
        else
            ! if below the interpolation range return the first value:
            if ( x <= self%x_initial ) then
                LinearIntepolateFunction1DValue = self%y(1)
                return
            end if
            ! if above the interpolation range return the first value:
            if ( x >= self%x_final   ) then
                LinearIntepolateFunction1DValue = self%y(self%num_points)
                return
            end if
        end if

        ! return the index of the point:
        if ( present(index) ) then
            ind = index
        else
            call hunt( self%x, self%num_points, x, ind)
        end if

        ! get the interpolation coefficient:
        if ( present(coeff) ) then
            mu = coeff
        else
            ! store the x values:
            x1  = self%x(ind)
            x2  = self%x(ind+1)
            ! compute the linear interpolation coefficient:
            mu  = (x-x1)/(x2-x1)
        end if

        ! store the y values:
        y1  = self%y(ind)
        y2  = self%y(ind+1)

        ! compute the linear interpolation:
        LinearIntepolateFunction1DValue = y1*( 1._dl -mu ) +y2*mu

    end function LinearIntepolateFunction1DValue

    ! ---------------------------------------------------------------------------------------------
    !> Function that gives the value of the function first derivative at a given coordinate x.
    function LinearIntepolateFunction1DFirstDerivative( self, x, index, coeff )

        implicit none

        class(linear_interpolate_function_1D)             :: self        !< the base class
        real(dl), intent(in)                              :: x           !< the value of x at which the function derivative is required
        integer , intent(in), optional                    :: index       !< optional precomputed value of the interpolation index
        real(dl), intent(in), optional                    :: coeff       !< optional precomputed value of the interpolation coefficient
        real(dl) :: LinearIntepolateFunction1DFirstDerivative  !< the output value of the function

        integer  :: ind
        real(dl) :: x1, x2, y1, y2, mu

        ! initialize to null value:
        LinearIntepolateFunction1DFirstDerivative = self%null_value
        if ( self%has_null_value ) then
            ! if outside the interpolation range return the null value:
            if ( x <= self%x_initial .or. x >= self%x_final ) return
        else
            ! if below the interpolation range return the first value:
            if ( x <= self%x_initial ) then
                LinearIntepolateFunction1DFirstDerivative = self%yp(1)
                return
            end if
            ! if above the interpolation range return the first value:
            if ( x >= self%x_final   ) then
                LinearIntepolateFunction1DFirstDerivative = self%yp(self%num_points)
                return
            end if
        end if

        ! return the index of the point:
        if ( present(index) ) then
            ind = index
        else
            call hunt( self%x, self%num_points, x, ind)
        end if

        ! get the interpolation coefficient:
        if ( present(coeff) ) then
            mu = coeff
        else
            ! store the x values:
            x1  = self%x(ind)
            x2  = self%x(ind+1)
            ! compute the linear interpolation coefficient:
            mu  = (x-x1)/(x2-x1)
        end if

        ! store the y values:
        y1  = self%yp(ind)
        y2  = self%yp(ind+1)

        ! compute the linear interpolation:
        LinearIntepolateFunction1DFirstDerivative = y1*( 1._dl -mu ) +y2*mu

    end function LinearIntepolateFunction1DFirstDerivative

    ! ---------------------------------------------------------------------------------------------
    !> Function that gives the value of the function second derivative at a given coordinate x.
    function LinearIntepolateFunction1DSecondDerivative( self, x, index, coeff )

        implicit none

        class(linear_interpolate_function_1D)             :: self        !< the base class
        real(dl), intent(in)                              :: x           !< the value of x at which the function derivative is required
        integer , intent(in), optional                    :: index       !< optional precomputed value of the interpolation index
        real(dl), intent(in), optional                    :: coeff       !< optional precomputed value of the interpolation coefficient
        real(dl) :: LinearIntepolateFunction1DSecondDerivative !< the output value of the function

        integer  :: ind
        real(dl) :: x1, x2, y1, y2, mu

        ! initialize to null value:
        LinearIntepolateFunction1DSecondDerivative = self%null_value
        if ( self%has_null_value ) then
            ! if outside the interpolation range return the null value:
            if ( x <= self%x_initial .or. x >= self%x_final ) return
        else
            ! if below the interpolation range return the first value:
            if ( x <= self%x_initial ) then
                LinearIntepolateFunction1DSecondDerivative = self%ypp(1)
                return
            end if
            ! if above the interpolation range return the first value:
            if ( x >= self%x_final   ) then
                LinearIntepolateFunction1DSecondDerivative = self%ypp(self%num_points)
                return
            end if
        end if

        ! return the index of the point:
        if ( present(index) ) then
            ind = index
        else
            call hunt( self%x, self%num_points, x, ind)
        end if

        ! get the interpolation coefficient:
        if ( present(coeff) ) then
            mu = coeff
        else
            ! store the x values:
            x1  = self%x(ind)
            x2  = self%x(ind+1)
            ! compute the linear interpolation coefficient:
            mu  = (x-x1)/(x2-x1)
        end if

        ! store the y values:
        y1  = self%ypp(ind)
        y2  = self%ypp(ind+1)

        ! compute the linear interpolation:
        LinearIntepolateFunction1DSecondDerivative = y1*( 1._dl -mu ) +y2*mu

    end function LinearIntepolateFunction1DSecondDerivative

    ! ---------------------------------------------------------------------------------------------
    !> Function that gives the value of the function third derivative at a given coordinate x.
    function LinearIntepolateFunction1DThirdDerivative( self, x, index, coeff )

        implicit none

        class(linear_interpolate_function_1D)             :: self        !< the base class
        real(dl), intent(in)                              :: x           !< the value of x at which the function derivative is required
        integer , intent(in), optional                    :: index       !< optional precomputed value of the interpolation index
        real(dl), intent(in), optional                    :: coeff       !< optional precomputed value of the interpolation coefficient
        real(dl) :: LinearIntepolateFunction1DThirdDerivative  !< the output value of the function

        integer  :: ind
        real(dl) :: x1, x2, y1, y2, mu

        ! initialize to null value:
        LinearIntepolateFunction1DThirdDerivative = self%null_value
        if ( self%has_null_value ) then
            ! if outside the interpolation range return the null value:
            if ( x <= self%x_initial .or. x >= self%x_final ) return
        else
            ! if below the interpolation range return the first value:
            if ( x <= self%x_initial ) then
                LinearIntepolateFunction1DThirdDerivative = self%yppp(1)
                return
            end if
            ! if above the interpolation range return the first value:
            if ( x >= self%x_final   ) then
                LinearIntepolateFunction1DThirdDerivative = self%yppp(self%num_points)
                return
            end if
        end if

        ! return the index of the point:
        if ( present(index) ) then
            ind = index
        else
            call hunt( self%x, self%num_points, x, ind)
        end if

        ! get the interpolation coefficient:
        if ( present(coeff) ) then
            mu = coeff
        else
            ! store the x values:
            x1  = self%x(ind)
            x2  = self%x(ind+1)
            ! compute the linear interpolation coefficient:
            mu  = (x-x1)/(x2-x1)
        end if

        ! store the y values:
        y1  = self%yppp(ind)
        y2  = self%yppp(ind+1)

        ! compute the linear interpolation:
        LinearIntepolateFunction1DThirdDerivative = y1*( 1._dl -mu ) +y2*mu

    end function LinearIntepolateFunction1DThirdDerivative

    ! ---------------------------------------------------------------------------------------------
    !> Function that gives the value of the function third derivative at a given coordinate x.
    function LinearIntepolateFunction1DFourthDerivative( self, x, index, coeff )

        implicit none

        class(linear_interpolate_function_1D)             :: self        !< the base class
        real(dl), intent(in)                              :: x           !< the value of x at which the function derivative is required
        integer , intent(in), optional                    :: index       !< optional precomputed value of the interpolation index
        real(dl), intent(in), optional                    :: coeff       !< optional precomputed value of the interpolation coefficient
        real(dl) :: LinearIntepolateFunction1DFourthDerivative  !< the output value of the function

        integer  :: ind
        real(dl) :: x1, x2, y1, y2, mu

        ! initialize to null value:
        LinearIntepolateFunction1DFourthDerivative = self%null_value
        if ( self%has_null_value ) then
           ! if outside the interpolation range return the null value:
           if ( x <= self%x_initial .or. x >= self%x_final ) return
        else
           ! if below the interpolation range return the first value:
           if ( x <= self%x_initial ) then
               LinearIntepolateFunction1DFourthDerivative = self%ypppp(1)
               return
           end if
           ! if above the interpolation range return the first value:
           if ( x >= self%x_final   ) then
               LinearIntepolateFunction1DFourthDerivative = self%ypppp(self%num_points)
               return
           end if
        end if

        ! return the index of the point:
        if ( present(index) ) then
           ind = index
        else
           call hunt( self%x, self%num_points, x, ind)
        end if

        ! get the interpolation coefficient:
        if ( present(coeff) ) then
           mu = coeff
        else
           ! store the x values:
           x1  = self%x(ind)
           x2  = self%x(ind+1)
           ! compute the linear interpolation coefficient:
           mu  = (x-x1)/(x2-x1)
        end if

        ! store the y values:
        y1  = self%ypppp(ind)
        y2  = self%ypppp(ind+1)

        ! compute the linear interpolation:
        LinearIntepolateFunction1DFourthDerivative = y1*( 1._dl -mu ) +y2*mu
        ! LinearIntepolateFunction1DFourthDerivative = 0._dl

    end function LinearIntepolateFunction1DFourthDerivative

    ! ---------------------------------------------------------------------------------------------
    !> Function that gives the value of the interpolated w DE integral at a given coordinate x.
    function LinearIntepolateFunction1DIntegral( self, x, index, coeff )

        implicit none

        class(linear_interpolate_function_1D)             :: self        !< the base class
        real(dl), intent(in)                              :: x           !< the value of x at which the function derivative is required
        integer , intent(in), optional                    :: index       !< optional precomputed value of the interpolation index
        real(dl), intent(in), optional                    :: coeff       !< optional precomputed value of the interpolation coefficient
        real(dl) :: LinearIntepolateFunction1DIntegral         !< the output value of the function

        integer  :: ind
        real(dl) :: x1, x2, y1, y2, mu

        ! initialize to null value:
        LinearIntepolateFunction1DIntegral = self%null_value
        if ( self%has_null_value ) then
            ! if outside the interpolation range return the null value:
            if ( x <= self%x_initial .or. x >= self%x_final ) return
        else
            ! if below the interpolation range return the first value:
            if ( x <= self%x_initial ) then
                LinearIntepolateFunction1DIntegral = self%yint(1)
                return
            end if
            ! if above the interpolation range return the first value:
            if ( x >= self%x_final   ) then
                LinearIntepolateFunction1DIntegral = self%yint(self%num_points)
                return
            end if
        end if

        ! return the index of the point:
        if ( present(index) ) then
            ind = index
        else
            call hunt( self%x, self%num_points, x, ind)
        end if

        ! get the interpolation coefficient:
        if ( present(coeff) ) then
            mu = coeff
        else
            ! store the x values:
            x1  = self%x(ind)
            x2  = self%x(ind+1)
            ! compute the linear interpolation coefficient:
            mu  = (x-x1)/(x2-x1)
        end if

        ! store the y values:
        y1  = self%yint(ind)
        y2  = self%yint(ind+1)

        ! compute the linear interpolation:
        LinearIntepolateFunction1DIntegral = y1*( 1._dl -mu ) +y2*mu

    end function LinearIntepolateFunction1DIntegral

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that initializes the derivatives if the derivatives vectors are not initialized.
    !! The derivative are inferred from the function itself.
    subroutine LinearIntepolateFunction1DInitDerivatives( self, jacobian )

        implicit none

        class(linear_interpolate_function_1D)             :: self        !< the base class
        real(dl), dimension(self%num_points), optional    :: jacobian    !< Jacobian of the transformation to use if we want the derivative wrt to another variable

        real(dl), dimension(self%num_points) :: dummy, dummy2
        integer  :: i

        ! compute the first derivative:
        call spline3ders( self%x, self%y, self%x, dummy, self%yp, self%ypp )
        ! compute the second derivative:
        call spline3ders( self%x, self%yp, self%x, dummy, self%ypp, self%yppp )
        ! compute the third derivative:
        call spline3ders( self%x, self%ypp, self%x, dummy, self%yppp, dummy2 )
        ! compute the fourth derivative:
        call spline3ders( self%x, self%yppp, self%x, dummy, self%ypppp, dummy2 )

        ! apply the Jacobian:
        if ( present(jacobian) ) then
            self%yp = jacobian*self%yp
            self%ypp = jacobian*self%ypp
            self%yppp = jacobian*self%yppp
            self%ypppp = jacobian*self%ypppp
        end if

        ! set the derivatives to zero at the boundary for continuity:
        self%yp(1)   = 0._dl
        self%ypp(1)  = 0._dl
        self%yppp(1) = 0._dl
        self%ypppp(1) = 0._dl
        self%yp(self%num_points)   = 0._dl
        self%ypp(self%num_points)  = 0._dl
        self%yppp(self%num_points) = 0._dl
        self%ypppp(self%num_points) = 0._dl

    end subroutine LinearIntepolateFunction1DInitDerivatives

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that initializes the integration vector of \f$ yint_i= \exp\left(-3\int_1^{x_i} \frac{1+f(x)}{x} \, dx \right) \f$. Only meaningful if we are parameterizing w_DE(a)
    subroutine LinearIntepolateFunction1DInitIntegration( self )

        implicit none

        class(linear_interpolate_function_1D)  :: self        !< the base class

        integer                              :: i
        real(dl)                             :: b

        self%yint(self%num_points) = 1._dl
        do i=self%num_points-1, 1, -1
            b = (self%x(i)*self%y(i+1) - self%x(i+1)*self%y(i))/(self%x(i) - self%x(i+1))
            self%yint(i) = self%yint(i+1)*(self%x(i)/self%x(i+1))**(-3._dl*(1._dl+b))*exp(-3._dl*(self%y(i)-self%y(i+1)))
        end do


    end subroutine LinearIntepolateFunction1DInitIntegration

end module linear_interpolation_1D

!----------------------------------------------------------------------------------------
