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

!> @file 02_equispaced_interpolation_linear_1D.f90
!! This file contains the class that can be used for 1D linearly interpolated functions
!! on an equispaced grid.


!----------------------------------------------------------------------------------------
!> This module contains the class that can be used for 1D linearly interpolated functions
!! on an equispaced grid.

!> @author Bin Hu, Marco Raveri

module equispaced_linear_interpolation_1D

    use precision
    use AMLutils
    use EFTCAMB_mixed_algorithms

    implicit none

    private

    public equispaced_linear_interpolate_function_1D

    !----------------------------------------------------------------------------------------
    !> This is the type that can be used for the 1D linear interpolation.
    type :: equispaced_linear_interpolate_function_1D

        ! options:
        integer                             :: num_points   !< number of points of the interpolating function.

        ! parameters:
        real(dl)                            :: x_initial      !< first value of x.
        real(dl)                            :: x_final        !< last value of x.
        real(dl)                            :: null_value     !< value that is returned if a point outside the range of interpolation is requested.
        logical                             :: has_null_value !< wether to use the null value outside interpolation range. If no null value is passed to initialization assume that function is constant outside interpolation range.
        real(dl)                            :: grid_width     !< the width of the interpolating grid.

        ! arrays with the values:
        real(dl), allocatable, dimension(:) :: x              !< array containing the values of x.
        real(dl), allocatable, dimension(:) :: y              !< array containing the values of the function \f$ y_i=f(x_i) \f$.
        real(dl), allocatable, dimension(:) :: yp             !< array containing the values of the function derivative \f$ yp_i= \frac{d f(x_i)}{dx} \f$.
        real(dl), allocatable, dimension(:) :: ypp            !< array containing the values of the function second derivative \f$ ypp_i= \frac{d^2 f(x_i)}{dx^2} \f$.
        real(dl), allocatable, dimension(:) :: yppp           !< array containing the values of the function third derivative \f$ yppp_i= \frac{d^3 f(x_i)}{dx^3} \f$.

    contains

        procedure :: initialize             => EquispacedLinearIntepolateFunction1DInit               !< subroutine that initialize the interpolating function.
        procedure :: value                  => EquispacedLinearIntepolateFunction1DValue              !< function that gives the value of the function at a given coordinate x.
        procedure :: first_derivative       => EquispacedLinearIntepolateFunction1DFirstDerivative    !< function that gives the value of the function first derivative at a given coordinate x.
        procedure :: second_derivative      => EquispacedLinearIntepolateFunction1DSecondDerivative   !< function that gives the value of the function second derivative at a given coordinate x.
        procedure :: third_derivative       => EquispacedLinearIntepolateFunction1DThirdDerivative    !< function that gives the value of the function third derivative at a given coordinate x.
        procedure :: initialize_derivatives => EquispacedLinearIntepolateFunction1DInitDerivatives    !< subroutine that initializes the derivatives if the derivatives vectors are not initialized. The derivative are inferred from the function itself.

    end type equispaced_linear_interpolate_function_1D

    ! ---------------------------------------------------------------------------------------------

contains

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that initialize the interpolating function.
    subroutine EquispacedLinearIntepolateFunction1DInit( self, num_points, x_initial, x_final, null_value )

        implicit none

        class(equispaced_linear_interpolate_function_1D) :: self         !< the base class

        integer , intent(in)                             :: num_points   !< number of points in the interpolation. Optional. If not passed it is inferred from x.
        real(dl), intent(in)                             :: x_initial    !< first value of x. Optional.
        real(dl), intent(in)                             :: x_final      !< last value of x. Optional.
        real(dl), intent(in), optional                   :: null_value   !< value that is returned if a point outside the range of interpolation is requested. Optional. If not passed to the constructor it is assumed to be zero.

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

        self%num_points = num_points

        allocate( self%x( self%num_points ) )

        ! fill the x vector:
        self%x_initial  = x_initial
        self%x_final    = x_final
        ! fill in the equispaced x array:
        do i=1, self%num_points
            self%x(i)   = self%x_initial + REAL(i-1)/REAL(self%num_points-1)*( self%x_final-self%x_initial )
        end do

        ! compute the width of the interpolation:
        self%grid_width = ( self%x_final -self%x_initial )/REAL( self%num_points -1 )

        ! allocate the other vectors:
        if ( allocated(self%y)    ) deallocate( self%y    )
        if ( allocated(self%yp)   ) deallocate( self%yp   )
        if ( allocated(self%ypp)  ) deallocate( self%ypp  )
        if ( allocated(self%yppp) ) deallocate( self%yppp )

        allocate( self%y( self%num_points )    )
        allocate( self%yp( self%num_points )   )
        allocate( self%ypp( self%num_points )  )
        allocate( self%yppp( self%num_points ) )

    end subroutine EquispacedLinearIntepolateFunction1DInit

    ! ---------------------------------------------------------------------------------------------
    !> Function that gives the value of the function at a given coordinate x.
    function EquispacedLinearIntepolateFunction1DValue( self, x )

        implicit none

        class(equispaced_linear_interpolate_function_1D)  :: self       !< the base class
        real(dl), intent(in)                              :: x          !< the value of x at which the function is required
        real(dl) :: EquispacedLinearIntepolateFunction1DValue           !< the output value of the function

        integer  :: ind
        real(dl) :: x1, x2, y1, y2, mu

        ! initialize to null value:
        EquispacedLinearIntepolateFunction1DValue = self%null_value
        if ( self%has_null_value ) then
            ! if outside the interpolation range return the null value:
            if ( x .lt. self%x_initial .or. x .gt. self%x_final ) return
        else
            ! if below the interpolation range return the first value:
            if ( x .lt. self%x_initial ) then
                EquispacedLinearIntepolateFunction1DValue = self%y(1)
                return
            end if
            ! if above the interpolation range return the first value:
            if ( x .gt. self%x_final   ) then
                EquispacedLinearIntepolateFunction1DValue = self%y(self%num_points)
                return
            end if
        end if

        ! return the index of the point:
        ind = int( ( x-self%x_initial)/self%grid_width ) +1

        ! store the values:
        x1  = self%x(ind)
        x2  = self%x(ind+1)
        y1  = self%y(ind)
        y2  = self%y(ind+1)
        ! compute the linear interpolation:
        mu  = (x-x1)/(x2-x1)
        EquispacedLinearIntepolateFunction1DValue = y1*( 1._dl -mu ) +y2*mu

    end function EquispacedLinearIntepolateFunction1DValue

    ! ---------------------------------------------------------------------------------------------
    !> Function that gives the value of the function first derivative at a given coordinate x.
    function EquispacedLinearIntepolateFunction1DFirstDerivative( self, x )

        implicit none

        class(equispaced_linear_interpolate_function_1D)  :: self        !< the base class
        real(dl), intent(in)                              :: x           !< the value of x at which the function derivative is required
        real(dl) :: EquispacedLinearIntepolateFunction1DFirstDerivative  !< the output value of the function

        integer  :: ind
        real(dl) :: x1, x2, y1, y2, mu

        ! initialize to null value:
        EquispacedLinearIntepolateFunction1DFirstDerivative = self%null_value
        if ( self%has_null_value ) then
            ! if outside the interpolation range return the null value:
            if ( x .lt. self%x_initial .or. x .gt. self%x_final ) return
        else
            ! if below the interpolation range return the first value:
            if ( x .lt. self%x_initial ) then
                EquispacedLinearIntepolateFunction1DFirstDerivative = self%yp(1)
                return
            end if
            ! if above the interpolation range return the first value:
            if ( x .gt. self%x_final   ) then
                EquispacedLinearIntepolateFunction1DFirstDerivative = self%yp(self%num_points)
                return
            end if
        end if

        ! return the index of the point:
        ind = int( ( x-self%x_initial)/self%grid_width ) +1

        ! store the values:
        x1  = self%x(ind)
        x2  = self%x(ind+1)
        y1  = self%yp(ind)
        y2  = self%yp(ind+1)
        ! compute the linear interpolation:
        mu  = (x-x1)/(x2-x1)
        EquispacedLinearIntepolateFunction1DFirstDerivative = y1*( 1._dl -mu ) +y2*mu

    end function EquispacedLinearIntepolateFunction1DFirstDerivative

    ! ---------------------------------------------------------------------------------------------
    !> Function that gives the value of the function second derivative at a given coordinate x.
    function EquispacedLinearIntepolateFunction1DSecondDerivative( self, x )

        implicit none

        class(equispaced_linear_interpolate_function_1D)  :: self        !< the base class
        real(dl), intent(in)                              :: x           !< the value of x at which the function derivative is required
        real(dl) :: EquispacedLinearIntepolateFunction1DSecondDerivative !< the output value of the function

        integer  :: ind
        real(dl) :: x1, x2, y1, y2, mu

        ! initialize to null value:
        EquispacedLinearIntepolateFunction1DSecondDerivative = self%null_value
        if ( self%has_null_value ) then
            ! if outside the interpolation range return the null value:
            if ( x .lt. self%x_initial .or. x .gt. self%x_final ) return
        else
            ! if below the interpolation range return the first value:
            if ( x .lt. self%x_initial ) then
                EquispacedLinearIntepolateFunction1DSecondDerivative = self%ypp(1)
                return
            end if
            ! if above the interpolation range return the first value:
            if ( x .gt. self%x_final   ) then
                EquispacedLinearIntepolateFunction1DSecondDerivative = self%ypp(self%num_points)
                return
            end if
        end if

        ! return the index of the point:
        ind = int( ( x-self%x_initial)/self%grid_width ) +1

        ! store the values:
        x1  = self%x(ind)
        x2  = self%x(ind+1)
        y1  = self%ypp(ind)
        y2  = self%ypp(ind+1)

        ! compute the linear interpolation:
        mu  = (x-x1)/(x2-x1)
        EquispacedLinearIntepolateFunction1DSecondDerivative = y1*( 1._dl -mu ) +y2*mu

    end function EquispacedLinearIntepolateFunction1DSecondDerivative

    ! ---------------------------------------------------------------------------------------------
    !> Function that gives the value of the function third derivative at a given coordinate x.
    function EquispacedLinearIntepolateFunction1DThirdDerivative( self, x )

        implicit none

        class(equispaced_linear_interpolate_function_1D)  :: self        !< the base class
        real(dl), intent(in)                              :: x           !< the value of x at which the function derivative is required
        real(dl) :: EquispacedLinearIntepolateFunction1DThirdDerivative  !< the output value of the function

        integer  :: ind
        real(dl) :: x1, x2, y1, y2, mu

        ! initialize to null value:
        EquispacedLinearIntepolateFunction1DThirdDerivative = self%null_value
        if ( self%has_null_value ) then
            ! if outside the interpolation range return the null value:
            if ( x .lt. self%x_initial .or. x .gt. self%x_final ) return
        else
            ! if below the interpolation range return the first value:
            if ( x .lt. self%x_initial ) then
                EquispacedLinearIntepolateFunction1DThirdDerivative = self%yppp(1)
                return
            end if
            ! if above the interpolation range return the first value:
            if ( x .gt. self%x_final   ) then
                EquispacedLinearIntepolateFunction1DThirdDerivative = self%yppp(self%num_points)
                return
            end if
        end if

        ! return the index of the point:
        ind = int( ( x-self%x_initial)/self%grid_width ) +1

        ! store the values:
        x1  = self%x(ind)
        x2  = self%x(ind+1)
        y1  = self%yppp(ind)
        y2  = self%yppp(ind+1)
        ! compute the linear interpolation:
        mu  = (x-x1)/(x2-x1)
        EquispacedLinearIntepolateFunction1DThirdDerivative = y1*( 1._dl -mu ) +y2*mu

    end function EquispacedLinearIntepolateFunction1DThirdDerivative

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that initializes the derivatives if the derivatives vectors are not initialized.
    !! The derivative are inferred from the function itself.
    subroutine EquispacedLinearIntepolateFunction1DInitDerivatives( self )

        implicit none

        class(equispaced_linear_interpolate_function_1D)  :: self        !< the base class

        write(*,*) 'IW'
        stop

    end subroutine EquispacedLinearIntepolateFunction1DInitDerivatives

    ! ---------------------------------------------------------------------------------------------

end module equispaced_linear_interpolation_1D

!----------------------------------------------------------------------------------------
