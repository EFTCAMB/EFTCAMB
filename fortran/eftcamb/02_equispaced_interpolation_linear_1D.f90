! ----------------------------------------------------------------------------------------
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
! ----------------------------------------------------------------------------------------

!> @file 02_equispaced_interpolation_linear_1D.f90
!! This file contains the class that can be used for 1D linearly interpolated functions
!! on an equispaced grid.


! ----------------------------------------------------------------------------------------
!> This module contains the class that can be used for 1D linearly interpolated functions
!! on an equispaced grid.

!> @author Bin Hu, Marco Raveri, Minsu Park

    module equispaced_linear_interpolation_1D

    use precision
    use EFTCAMB_mixed_algorithms
    use MpiUtils
    use EFT_splines

    implicit none

    private

    public equispaced_linear_interpolate_function_1D

    ! ----------------------------------------------------------------------------------------
    !> This is the type that can be used for the 1D linear interpolation.
    type :: equispaced_linear_interpolate_function_1D

        ! options:
        integer                             :: num_points     !< number of points of the interpolating function.

        ! parameters:
        real(dl)                            :: x_initial      !< first value of x.
        real(dl)                            :: x_final        !< last value of x.
        !< value that is returned if a point outside the range of interpolation is requested.
        real(dl)                            :: null_value
        !< wether to use the null value outside interpolation range. If no null value is passed to initialization assume that function is constant outside interpolation range.
        logical                             :: has_null_value
        real(dl)                            :: grid_width     !< the width of the interpolating grid.

        ! arrays with the values:
        real(dl), allocatable, dimension(:) :: x              !< array containing the values of x.
        !< array containing the values of the function \f$ y_i=f(x_i) \f$.
        real(dl), allocatable, dimension(:) :: y
        !< array containing the values of the function derivative \f$ yp_i= \frac{d f(x_i)}{dx} \f$.
        real(dl), allocatable, dimension(:) :: yp
        !< array containing the values of the function second derivative \f$ ypp_i= \frac{d^2 f(x_i)}{dx^2} \f$.
        real(dl), allocatable, dimension(:) :: ypp
        !< array containing the values of the function third derivative \f$ yppp_i= \frac{d^3 f(x_i)}{dx^3} \f$.
        real(dl), allocatable, dimension(:) :: yppp
        !< array containing the values of the function fourth derivative \f$ypppp_i= \frac{d^4 f(x_i)}{dx^4} \f$.
        real(dl), allocatable, dimension(:) :: ypppp
        !< array containing the values of the function w DE integral \f$ yint_i= \exp\left(-3\int_1^{x_i} \frac{1+f(x)}{x} \, dx \right) \f$.
        real(dl), allocatable, dimension(:) :: yint

    contains

    !< subroutine that initialize the interpolating function.
    procedure :: initialize => EquispacedLinearIntepolateFunction1DInit
    !< subroutine that does precomputations for the interpolation. Usefull when calling values and derivatives.
    procedure :: precompute => EquispacedLinearIntepolateFunction1DPrecompute
    !< function that gives the value of the function at a given coordinate x.
    procedure :: value => EquispacedLinearIntepolateFunction1DValue
    !< function that gives the value of the function first derivative at a given coordinate x.
    procedure :: first_derivative => EquispacedLinearIntepolateFunction1DFirstDerivative
    !< function that gives the value of the function second derivative at a given coordinate x.
    procedure :: second_derivative => EquispacedLinearIntepolateFunction1DSecondDerivative
    !< function that gives the value of the function third derivative at a given coordinate x.
    procedure :: third_derivative => EquispacedLinearIntepolateFunction1DThirdDerivative
    !< function that gives the value of the function fourth derivative at a given coordinate x.
    procedure :: fourth_derivative => EquispacedLinearIntepolateFunction1DFourthDerivative
    !< function that gives the value of the interpolated w DE integral at a given coordinate x.
    procedure :: integral => EquispacedLinearIntepolateFunction1DIntegral
    !< subroutine that initializes the derivatives if the derivatives vectors are not initialized. The derivative are inferred from the function itself.
    procedure :: initialize_derivatives => EquispacedLinearIntepolateFunction1DInitDerivatives
    !< subroutine that initializes the integration vector. Jacobian is not supported now. by Gen
    procedure :: initialize_integration => EquispacedLinearIntepolateFunction1DInitIntegration

    end type equispaced_linear_interpolate_function_1D

    ! ---------------------------------------------------------------------------------------------

    contains

    ! ----------------------------------------------------------------------------------------
    !> Equispaced cubic spline derivative helpers, from the CAMB splines module
    !! (CAMB fortran/subroutines.f90, removed upstream).
    subroutine splder(y, dy, n, g)
    ! Splder fits a cubic spline to y and returns the first derivatives at
    ! the grid points in dy. Dy is equivalent to a 4th-order Pade
    ! difference formula for dy/di.
    integer, intent(in) :: n
    real(dl), intent(in) :: y(n), g(n)
    real(dl), intent(out) :: dy(n)
    integer :: n1, i
    real(dl), allocatable, dimension(:) :: f

    allocate(f(n))
    n1 = n - 1
    ! Quartic fit to dy/di at boundaries, assuming d3y/di3=0.
    f(1) = (-10._dl*y(1) + 15._dl*y(2) - 6._dl*y(3) + y(4))/6._dl
    f(n) = (10._dl*y(n) - 15._dl*y(n1) + 6._dl*y(n - 2) - y(n - 3))/6._dl
    ! Solve the tridiagonal system
    ! dy(i-1)+4*dy(i)+dy(i+1)=3*(y(i+1)-y(i-1)), i=2,3,...,n1,
    ! with dy(1)=f(1), dy(n)=f(n).
    do i = 2, n1
        f(i) = g(i)*(3._dl*(y(i + 1) - y(i - 1)) - f(i - 1))
    end do
    dy(n) = f(n)
    do i = n1, 1, -1
        dy(i) = f(i) - g(i)*dy(i + 1)
    end do

    end subroutine splder

    subroutine splini(g, n)
    ! Splini must be called before splder to initialize array g.
    integer, intent(in) :: n
    real(dl), intent(out) :: g(n)
    integer :: i

    g(1) = 0._dl
    do i = 2, n
        g(i) = 1/(4._dl - g(i - 1))
    end do

    end subroutine splini

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that initialize the interpolating function.
    subroutine EquispacedLinearIntepolateFunction1DInit(self, num_points, x_initial, x_final, null_value)

    implicit none

    class(equispaced_linear_interpolate_function_1D) :: self         !< the base class

    !< number of points in the interpolation. Optional. If not passed it is inferred from x.
    integer, intent(in)            :: num_points
    real(dl), intent(in)           :: x_initial    !< first value of x. Optional.
    real(dl), intent(in)           :: x_final      !< last value of x. Optional.
    !< value that is returned if a point outside the range of interpolation is requested. Optional. If not passed to the constructor it is assumed to be zero.
    real(dl), intent(in), optional :: null_value

    integer :: i

    if (x_initial > x_final) then
        write(*, *) "x_initial must be smaller than x_final"
        call MpiStop('EFTCAMB error')
    end if

    ! initialize the null value:
    if (present(null_value)) then
        self%has_null_value = .true.
        self%null_value = null_value
    else
        self%has_null_value = .false.
        self%null_value = 0._dl
    end if

    ! allocate the x vector:
    if (allocated(self%x)) deallocate(self%x)

    self%num_points = num_points

    allocate(self%x(self%num_points))

    ! store initial and final times:
    self%x_initial = x_initial
    self%x_final = x_final

    ! compute the width of the interpolation:
    self%grid_width = (self%x_final - self%x_initial)/real(self%num_points - 1)

    ! fill in the equispaced x array:
    do i = 1, self%num_points
        self%x(i) = self%x_initial + real(i - 1)*self%grid_width
    end do

    ! allocate the other vectors:
    if (allocated(self%y)) deallocate(self%y)
    if (allocated(self%yp)) deallocate(self%yp)
    if (allocated(self%ypp)) deallocate(self%ypp)
    if (allocated(self%yppp)) deallocate(self%yppp)
    if (allocated(self%ypppp)) deallocate(self%ypppp)
    if (allocated(self%yint)) deallocate(self%yint)

    allocate(self%y(self%num_points))
    allocate(self%yp(self%num_points))
    allocate(self%ypp(self%num_points))
    allocate(self%yppp(self%num_points))
    allocate(self%ypppp(self%num_points))
    allocate(self%yint(self%num_points))

    end subroutine EquispacedLinearIntepolateFunction1DInit

    ! ---------------------------------------------------------------------------------------------
    !> Function that computes the main interpolation index.
    !! Usefull when calling value and derivatives in the same place.
    subroutine EquispacedLinearIntepolateFunction1DPrecompute(self, x, ind, mu)

    implicit none

    class(equispaced_linear_interpolate_function_1D) :: self       !< the base class
    real(dl), intent(in)                             :: x          !< the value of x at which the function is required
    integer, intent(out)                             :: ind        !< the main interpolation index
    real(dl), intent(out)                            :: mu         !< the interpolation coefficient

    real(dl) :: x1, x2

    ! check input:
    if (x <= self%x_initial .or. x >= self%x_final) then
        ind = self%num_points
        x1 = 0._dl
        x2 = 0._dl
        mu = 0._dl
        return
    end if
    ! compute the interpolation index:
    ind = int((x - self%x_initial)/self%grid_width) + 1
    ! store the x values:
    x1 = self%x(ind)
    x2 = self%x(ind + 1)
    ! compute the linear interpolation coefficient:
    mu = (x - x1)/(x2 - x1)

    end subroutine EquispacedLinearIntepolateFunction1DPrecompute

    ! ---------------------------------------------------------------------------------------------
    !> Function that gives the value of the function at a given coordinate x.
    function EquispacedLinearIntepolateFunction1DValue(self, x, index, coeff)

    implicit none

    class(equispaced_linear_interpolate_function_1D) :: self       !< the base class
    real(dl), intent(in)                             :: x          !< the value of x at which the function is required
    !< optional precomputed value of the interpolation index
    integer, intent(in), optional                    :: index
    !< optional precomputed value of the interpolation coefficient
    real(dl), intent(in), optional                   :: coeff
    real(dl) :: EquispacedLinearIntepolateFunction1DValue           !< the output value of the function
    integer  :: ind
    real(dl) :: x1, x2, y1, y2, mu

    ! initialize to null value:
    EquispacedLinearIntepolateFunction1DValue = self%null_value
    if (self%has_null_value) then
        ! if outside the interpolation range return the null value:
        if (x <= self%x_initial .or. x >= self%x_final) return
    else
        ! if below the interpolation range return the first value:
        if (x <= self%x_initial) then
            EquispacedLinearIntepolateFunction1DValue = self%y(1)
            return
        end if
        ! if above the interpolation range return the first value:
        if (x >= self%x_final) then
            EquispacedLinearIntepolateFunction1DValue = self%y(self%num_points)
            return
        end if
    end if

    ! return the index of the point:
    if (present(index)) then
        ind = index
    else
        ind = int((x - self%x_initial)/self%grid_width) + 1
    end if

    ! get the interpolation coefficient:
    if (present(coeff)) then
        mu = coeff
    else
        ! store the x values:
        x1 = self%x(ind)
        x2 = self%x(ind + 1)
        ! compute the linear interpolation coefficient:
        mu = (x - x1)/(x2 - x1)
    end if

    ! store the y values:
    y1 = self%y(ind)
    y2 = self%y(ind + 1)

    ! compute the linear interpolation:
    EquispacedLinearIntepolateFunction1DValue = y1*(1._dl - mu) + y2*mu

    end function EquispacedLinearIntepolateFunction1DValue

    ! ---------------------------------------------------------------------------------------------
    !> Function that gives the value of the function first derivative at a given coordinate x.
    function EquispacedLinearIntepolateFunction1DFirstDerivative(self, x, index, coeff)

    implicit none

    class(equispaced_linear_interpolate_function_1D) :: self        !< the base class
    !< the value of x at which the function derivative is required
    real(dl), intent(in)                             :: x
    !< optional precomputed value of the interpolation index
    integer, intent(in), optional                    :: index
    !< optional precomputed value of the interpolation coefficient
    real(dl), intent(in), optional                   :: coeff
    real(dl) :: EquispacedLinearIntepolateFunction1DFirstDerivative  !< the output value of the function

    integer  :: ind
    real(dl) :: x1, x2, y1, y2, mu

    ! initialize to null value:
    EquispacedLinearIntepolateFunction1DFirstDerivative = self%null_value
    if (self%has_null_value) then
        ! if outside the interpolation range return the null value:
        if (x <= self%x_initial .or. x >= self%x_final) return
    else
        ! if below the interpolation range return the first value:
        if (x <= self%x_initial) then
            EquispacedLinearIntepolateFunction1DFirstDerivative = self%yp(1)
            return
        end if
        ! if above the interpolation range return the first value:
        if (x >= self%x_final) then
            EquispacedLinearIntepolateFunction1DFirstDerivative = self%yp(self%num_points)
            return
        end if
    end if

    ! return the index of the point:
    if (present(index)) then
        ind = index
    else
        ind = int((x - self%x_initial)/self%grid_width) + 1
    end if

    ! get the interpolation coefficient:
    if (present(coeff)) then
        mu = coeff
    else
        ! store the x values:
        x1 = self%x(ind)
        x2 = self%x(ind + 1)
        ! compute the linear interpolation coefficient:
        mu = (x - x1)/(x2 - x1)
    end if

    ! store the y values:
    y1 = self%yp(ind)
    y2 = self%yp(ind + 1)

    ! compute the linear interpolation:
    EquispacedLinearIntepolateFunction1DFirstDerivative = y1*(1._dl - mu) + y2*mu

    end function EquispacedLinearIntepolateFunction1DFirstDerivative

    ! ---------------------------------------------------------------------------------------------
    !> Function that gives the value of the function second derivative at a given coordinate x.
    function EquispacedLinearIntepolateFunction1DSecondDerivative(self, x, index, coeff)

    implicit none

    class(equispaced_linear_interpolate_function_1D) :: self        !< the base class
    !< the value of x at which the function derivative is required
    real(dl), intent(in)                             :: x
    !< optional precomputed value of the interpolation index
    integer, intent(in), optional                    :: index
    !< optional precomputed value of the interpolation coefficient
    real(dl), intent(in), optional                   :: coeff
    real(dl) :: EquispacedLinearIntepolateFunction1DSecondDerivative !< the output value of the function

    integer  :: ind
    real(dl) :: x1, x2, y1, y2, mu

    ! initialize to null value:
    EquispacedLinearIntepolateFunction1DSecondDerivative = self%null_value
    if (self%has_null_value) then
        ! if outside the interpolation range return the null value:
        if (x <= self%x_initial .or. x >= self%x_final) return
    else
        ! if below the interpolation range return the first value:
        if (x <= self%x_initial) then
            EquispacedLinearIntepolateFunction1DSecondDerivative = self%ypp(1)
            return
        end if
        ! if above the interpolation range return the first value:
        if (x >= self%x_final) then
            EquispacedLinearIntepolateFunction1DSecondDerivative = self%ypp(self%num_points)
            return
        end if
    end if

    ! return the index of the point:
    if (present(index)) then
        ind = index
    else
        ind = int((x - self%x_initial)/self%grid_width) + 1
    end if

    ! get the interpolation coefficient:
    if (present(coeff)) then
        mu = coeff
    else
        ! store the x values:
        x1 = self%x(ind)
        x2 = self%x(ind + 1)
        ! compute the linear interpolation coefficient:
        mu = (x - x1)/(x2 - x1)
    end if

    ! store the y values:
    y1 = self%ypp(ind)
    y2 = self%ypp(ind + 1)

    ! compute the linear interpolation:
    EquispacedLinearIntepolateFunction1DSecondDerivative = y1*(1._dl - mu) + y2*mu

    end function EquispacedLinearIntepolateFunction1DSecondDerivative

    ! ---------------------------------------------------------------------------------------------
    !> Function that gives the value of the function third derivative at a given coordinate x.
    function EquispacedLinearIntepolateFunction1DThirdDerivative(self, x, index, coeff)

    implicit none

    class(equispaced_linear_interpolate_function_1D) :: self        !< the base class
    !< the value of x at which the function derivative is required
    real(dl), intent(in)                             :: x
    !< optional precomputed value of the interpolation index
    integer, intent(in), optional                    :: index
    !< optional precomputed value of the interpolation coefficient
    real(dl), intent(in), optional                   :: coeff
    real(dl) :: EquispacedLinearIntepolateFunction1DThirdDerivative  !< the output value of the function

    integer  :: ind
    real(dl) :: x1, x2, y1, y2, mu

    ! initialize to null value:
    EquispacedLinearIntepolateFunction1DThirdDerivative = self%null_value
    if (self%has_null_value) then
        ! if outside the interpolation range return the null value:
        if (x <= self%x_initial .or. x >= self%x_final) return
    else
        ! if below the interpolation range return the first value:
        if (x <= self%x_initial) then
            EquispacedLinearIntepolateFunction1DThirdDerivative = self%yppp(1)
            return
        end if
        ! if above the interpolation range return the first value:
        if (x >= self%x_final) then
            EquispacedLinearIntepolateFunction1DThirdDerivative = self%yppp(self%num_points)
            return
        end if
    end if

    ! return the index of the point:
    if (present(index)) then
        ind = index
    else
        ind = int((x - self%x_initial)/self%grid_width) + 1
    end if

    ! get the interpolation coefficient:
    if (present(coeff)) then
        mu = coeff
    else
        ! store the x values:
        x1 = self%x(ind)
        x2 = self%x(ind + 1)
        ! compute the linear interpolation coefficient:
        mu = (x - x1)/(x2 - x1)
    end if

    ! store the y values:
    y1 = self%yppp(ind)
    y2 = self%yppp(ind + 1)

    ! compute the linear interpolation:
    EquispacedLinearIntepolateFunction1DThirdDerivative = y1*(1._dl - mu) + y2*mu

    end function EquispacedLinearIntepolateFunction1DThirdDerivative

    ! ---------------------------------------------------------------------------------------------
    !> Function that gives the value of the function fourth derivative at a given coordinate x.
    function EquispacedLinearIntepolateFunction1DFourthDerivative(self, x, index, coeff)

    implicit none

    class(equispaced_linear_interpolate_function_1D) :: self        !< the base class
    !< the value of x at which the function derivative is required
    real(dl), intent(in)                             :: x
    !< optional precomputed value of the interpolation index
    integer, intent(in), optional                    :: index
    !< optional precomputed value of the interpolation coefficient
    real(dl), intent(in), optional                   :: coeff
    real(dl) :: EquispacedLinearIntepolateFunction1DFourthDerivative  !< the output value of the function

    integer  :: ind
    real(dl) :: x1, x2, y1, y2, mu

    ! initialize to null value:
    EquispacedLinearIntepolateFunction1DFourthDerivative = self%null_value
    if (self%has_null_value) then
        ! if outside the interpolation range return the null value:
        if (x <= self%x_initial .or. x >= self%x_final) return
    else
        ! if below the interpolation range return the first value:
        if (x <= self%x_initial) then
            EquispacedLinearIntepolateFunction1DFourthDerivative = self%ypppp(1)
            return
        end if
        ! if above the interpolation range return the first value:
        if (x >= self%x_final) then
            EquispacedLinearIntepolateFunction1DFourthDerivative = self%ypppp(self%num_points)
            return
        end if
    end if

    ! return the index of the point:
    if (present(index)) then
        ind = index
    else
        ind = int((x - self%x_initial)/self%grid_width) + 1
    end if

    ! get the interpolation coefficient:
    if (present(coeff)) then
        mu = coeff
    else
        ! store the x values:
        x1 = self%x(ind)
        x2 = self%x(ind + 1)
        ! compute the linear interpolation coefficient:
        mu = (x - x1)/(x2 - x1)
    end if

    ! store the y values:
    y1 = self%ypppp(ind)
    y2 = self%ypppp(ind + 1)

    ! compute the linear interpolation:
    EquispacedLinearIntepolateFunction1DFourthDerivative = y1*(1._dl - mu) + y2*mu

    end function EquispacedLinearIntepolateFunction1DFourthDerivative

    ! ---------------------------------------------------------------------------------------------
    !> Function that gives the value of the interpolated w DE integral at a given coordinate x.
    function EquispacedLinearIntepolateFunction1DIntegral(self, x, index, coeff)

    implicit none

    class(equispaced_linear_interpolate_function_1D) :: self        !< the base class
    !< the value of x at which the function derivative is required
    real(dl), intent(in)                             :: x
    !< optional precomputed value of the interpolation index
    integer, intent(in), optional                    :: index
    !< optional precomputed value of the interpolation coefficient
    real(dl), intent(in), optional                   :: coeff
    real(dl) :: EquispacedLinearIntepolateFunction1DIntegral         !< the output value of the function

    integer  :: ind
    real(dl) :: x1, x2, y1, y2, mu

    ! initialize to null value:
    EquispacedLinearIntepolateFunction1DIntegral = self%null_value
    if (self%has_null_value) then
        ! if outside the interpolation range return the null value:
        if (x <= self%x_initial .or. x >= self%x_final) return
    else
        ! if below the interpolation range return the first value:
        if (x <= self%x_initial) then
            EquispacedLinearIntepolateFunction1DIntegral = self%yint(1)
            return
        end if
        ! if above the interpolation range return the first value:
        if (x >= self%x_final) then
            EquispacedLinearIntepolateFunction1DIntegral = self%yint(self%num_points)
            return
        end if
    end if

    ! return the index of the point:
    if (present(index)) then
        ind = index
    else
        ind = int((x - self%x_initial)/self%grid_width) + 1
    end if

    ! get the interpolation coefficient:
    if (present(coeff)) then
        mu = coeff
    else
        ! store the x values:
        x1 = self%x(ind)
        x2 = self%x(ind + 1)
        ! compute the linear interpolation coefficient:
        mu = (x - x1)/(x2 - x1)
    end if

    ! store the y values:
    y1 = self%yint(ind)
    y2 = self%yint(ind + 1)

    ! compute the linear interpolation:
    EquispacedLinearIntepolateFunction1DIntegral = y1*(1._dl - mu) + y2*mu

    end function EquispacedLinearIntepolateFunction1DIntegral

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that initializes the derivatives if the derivatives vectors are not initialized.
    !! The derivative are inferred from the function itself.
    subroutine EquispacedLinearIntepolateFunction1DInitDerivatives(self, jacobian)

    implicit none

    class(equispaced_linear_interpolate_function_1D) :: self        !< the base class
    !< Jacobian of the transformation to use if we want the derivative wrt to another variable
    real(dl), dimension(self%num_points), optional   :: jacobian

    real(dl), dimension(self%num_points) :: spline_workspace

    ! initialize the calculation:
    call splini(spline_workspace, self%num_points)
    ! compute the numerical first derivative:
    call splder(self%y, self%yp, self%num_points, spline_workspace)
    self%yp = self%yp/self%grid_width
    if (present(jacobian)) then
        self%yp = jacobian*self%yp
    end if
    ! compute the numerical second derivative:
    call splder(self%yp, self%ypp, self%num_points, spline_workspace)
    self%ypp = self%ypp/self%grid_width
    if (present(jacobian)) then
        self%ypp = jacobian*self%ypp
    end if
    ! compute the numerical third derivative:
    call splder(self%ypp, self%yppp, self%num_points, spline_workspace)
    self%yppp = self%yppp/self%grid_width
    if (present(jacobian)) then
        self%yppp = jacobian*self%yppp
    end if

    ! compute the numerical fourth derivative:
    call splder(self%yppp, self%ypppp, self%num_points, spline_workspace)
    self%ypppp = self%ypppp/self%grid_width
    if (present(jacobian)) then
        self%ypppp = jacobian*self%ypppp
    end if

    end subroutine EquispacedLinearIntepolateFunction1DInitDerivatives

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that initializes the integration vector of \f$ yint_i= \exp\left(-3\int_1^{x_i} \frac{1+f(x)}{x} \, dx \right) \f$. Only meaningful if we are parameterizing w_DE(a)
    subroutine EquispacedLinearIntepolateFunction1DInitIntegration(self)

    implicit none

    class(equispaced_linear_interpolate_function_1D) :: self        !< the base class

    integer  :: i
    real(dl) :: b

    self%yint(self%num_points) = 1._dl
    do i = self%num_points - 1, 1, -1
        b = (self%x(i)*self%y(i + 1) - self%x(i + 1)*self%y(i))/(self%x(i) - self%x(i + 1))
        self%yint(i) = &
            self%yint(i + 1)*(self%x(i)/self%x(i + 1))**(-3._dl*(1._dl + b))*exp(-3._dl*(self%y(i) - self%y(i + 1)))
    end do

    end subroutine EquispacedLinearIntepolateFunction1DInitIntegration

    ! ---------------------------------------------------------------------------------------------

    end module equispaced_linear_interpolation_1D

! ----------------------------------------------------------------------------------------
