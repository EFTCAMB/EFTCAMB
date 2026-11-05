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

!> @file 04p016_spline_parametrizations_1D.f90
!! This file contains the definition of a 1D function given by a natural cubic spline of its nodes, with
!! approximate derivatives (spline_rec).
!! Inherits from spline3_parametrization_1D.

! ----------------------------------------------------------------------------------------
!> This module contains the definition of a 1D function given by a natural cubic spline of its nodes, with
!! approximate derivatives (spline_rec).
!! The value, the w_DE integral and the parameter handling are those of spline3_parametrization_1D.
!! The derivatives are not the derivatives of the value spline: each one is a separate natural cubic spline
!! through nodal derivatives, computed by repeated cubic spline differentiation of the nodal values and set to
!! zero at the first and last node. They are zero outside the nodes, where the function is the constant null_value.

!> @author Marco Raveri

    module EFTCAMB_spline_parametrizations_1D

    use precision
    use EFT_splines
    use EFT_def
    use EFTCAMB_cache
    use EFTCAMB_spline3_parametrizations_1D

    implicit none

    private

    public spline_rec_parametrization_1D

    ! ---------------------------------------------------------------------------------------------
    !> Type containing the approximate cubic spline parametrization. Inherits from spline3_parametrization_1D.
    type, extends (spline3_parametrization_1D) :: spline_rec_parametrization_1D

        ! internal variables:
        !< internal coefficients of the splines of the nodal derivatives
        real(dl), dimension(:, :), allocatable :: cyp, cypp, cyppp, cypppp

    contains

    ! utility functions:
    !< subroutine that initializes the function parameters based on the values found in an input array.
    procedure :: init_parameters => SplineRec1DInitParams
    !< subroutine that prints to screen the informations about the function.
    procedure :: feedback => SplineRec1DFeedback

    ! evaluation procedures:
    !< function that returns the first derivative of the reconstruction.
    procedure :: first_derivative => SplineRec1DFirstDerivative
    !< function that returns the second derivative of the reconstruction.
    procedure :: second_derivative => SplineRec1DSecondDerivative
    !< function that returns the third derivative of the reconstruction.
    procedure :: third_derivative => SplineRec1DThirdDerivative
    !< function that returns the fourth derivative of the reconstruction.
    procedure :: fourth_derivative => SplineRec1DFourthDerivative

    end type spline_rec_parametrization_1D

    contains

    ! ---------------------------------------------------------------------------------------------
    ! Implementation of the reconstructioin parametrization.
    ! ---------------------------------------------------------------------------------------------

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that initializes the function parameters based on the values found in an input array.
    subroutine SplineRec1DInitParams(self, array)

    implicit none

    class(spline_rec_parametrization_1D)                   :: self   !< the base class.
    real(dl), dimension(self%parameter_number), intent(in) :: array  !< input array with the values of the parameters.

    integer :: ind
    logical :: valid

    ! nodes, value spline and node checks:
    call self%spline3_parametrization_1D%init_parameters(array)
    valid = .true.
    do ind = 1, self%num_pixels - 1
        if (self%x_coord(ind) >= self%x_coord(ind + 1)) valid = .false.
    end do
    ! create the splines of the nodal derivatives:
    if (allocated(self%cyp)) deallocate(self%cyp)
    if (allocated(self%cypp)) deallocate(self%cypp)
    if (allocated(self%cyppp)) deallocate(self%cyppp)
    if (allocated(self%cypppp)) deallocate(self%cypppp)
    allocate(self%cyp   (5, self%num_pixels - 1))
    allocate(self%cypp  (5, self%num_pixels - 1))
    allocate(self%cyppp (5, self%num_pixels - 1))
    allocate(self%cypppp(5, self%num_pixels - 1))
    if (valid) then
        ! (the error for invalid nodes has been raised by spline3, the coefficients are left at zero)
        call self%interpolation%initialize_derivatives()
        call spline3pars(self%x_coord, self%interpolation%yp, [2, 2], [0._dl, 0._dl], self%cyp)
        call spline3pars(self%x_coord, self%interpolation%ypp, [2, 2], [0._dl, 0._dl], self%cypp)
        call spline3pars(self%x_coord, self%interpolation%yppp, [2, 2], [0._dl, 0._dl], self%cyppp)
        call spline3pars(self%x_coord, self%interpolation%ypppp, [2, 2], [0._dl, 0._dl], self%cypppp)
    else
        self%cyp = 0._dl
        self%cypp = 0._dl
        self%cyppp = 0._dl
        self%cypppp = 0._dl
    end if

    end subroutine SplineRec1DInitParams

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that prints to screen the informations about the function.
    subroutine SplineRec1DFeedback(self, print_params)

    implicit none

    class(spline_rec_parametrization_1D) :: self         !< the base class
    logical, optional                    :: print_params !< optional flag that decised whether to print numerical
    !! values of the parameters.

    integer                             :: i
    real(dl)                            :: param_value
    character(len=EFT_names_max_length) :: param_name
    logical                             :: print_params_temp

    if (present(print_params)) then
        print_params_temp = print_params
    else
        print_params_temp = .true.
    end if

    write(*, '(a,a)')      'spline_rec parametrization: ', self%name
    write(*, '(a,a,a,I3)') ' number of pixels (', trim(self%name) // '_Spline_Pixels', '): ', self%num_pixels
    write(*, '(a,F12.6)')  ' null value: ', self%null_value
    if (print_params_temp) then
        do i = 1, self%parameter_number
            call self%parameter_names(i, param_name)
            call self%parameter_value(i, param_value)
            write(*, '(a23,a,F12.6)') param_name, '=', param_value
        end do
    end if

    end subroutine SplineRec1DFeedback

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the value of the first derivative, wrt scale factor, of the function.
    function SplineRec1DFirstDerivative(self, x, eft_cache)

    implicit none

    class(spline_rec_parametrization_1D)                :: self      !< the base class
    real(dl), intent(in)                                :: x         !< the input scale factor
    type(TEFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
    real(dl) :: SplineRec1DFirstDerivative                           !< the output value

    integer :: ip

    ! outside the nodes the function is the constant null_value:
    SplineRec1DFirstDerivative = 0._dl
    if (x < self%interpolation%x_initial .or. x > self%interpolation%x_final) return
    ! if inside return the spline of the nodal derivative:
    ip = iixmin(x, self%x_coord, 0)
    SplineRec1DFirstDerivative = poly3 (x, self%cyp(:, ip))

    end function SplineRec1DFirstDerivative

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the second derivative of the function.
    function SplineRec1DSecondDerivative(self, x, eft_cache)

    implicit none

    class(spline_rec_parametrization_1D)                :: self      !< the base class
    real(dl), intent(in)                                :: x         !< the input scale factor
    type(TEFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
    real(dl) :: SplineRec1DSecondDerivative                          !< the output value

    integer :: ip

    ! outside the nodes the function is the constant null_value:
    SplineRec1DSecondDerivative = 0._dl
    if (x < self%interpolation%x_initial .or. x > self%interpolation%x_final) return
    ! if inside return the spline of the nodal derivative:
    ip = iixmin(x, self%x_coord, 0)
    SplineRec1DSecondDerivative = poly3 (x, self%cypp(:, ip))

    end function SplineRec1DSecondDerivative

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the third derivative of the function.
    function SplineRec1DThirdDerivative(self, x, eft_cache)

    implicit none

    class(spline_rec_parametrization_1D)                :: self      !< the base class
    real(dl), intent(in)                                :: x         !< the input scale factor
    type(TEFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
    real(dl) :: SplineRec1DThirdDerivative                           !< the output value

    integer :: ip

    ! outside the nodes the function is the constant null_value:
    SplineRec1DThirdDerivative = 0._dl
    if (x < self%interpolation%x_initial .or. x > self%interpolation%x_final) return
    ! if inside return the spline of the nodal derivative:
    ip = iixmin(x, self%x_coord, 0)
    SplineRec1DThirdDerivative = poly3 (x, self%cyppp(:, ip))

    end function SplineRec1DThirdDerivative

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the fourth derivative of the function.
    function SplineRec1DFourthDerivative(self, x, eft_cache)

    implicit none

    class(spline_rec_parametrization_1D)                :: self      !< the base class
    real(dl), intent(in)                                :: x         !< the input scale factor
    type(TEFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
    real(dl) :: SplineRec1DFourthDerivative                          !< the output value

    integer :: ip

    ! outside the nodes the function is the constant null_value:
    SplineRec1DFourthDerivative = 0._dl
    if (x < self%interpolation%x_initial .or. x > self%interpolation%x_final) return
    ! if inside return the spline of the nodal derivative:
    ip = iixmin(x, self%x_coord, 0)
    SplineRec1DFourthDerivative = poly3 (x, self%cypppp(:, ip))

    end function SplineRec1DFourthDerivative

    ! ---------------------------------------------------------------------------------------------

    end module EFTCAMB_spline_parametrizations_1D

! ----------------------------------------------------------------------------------------
