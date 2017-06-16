!----------------------------------------------------------------------------------------
!
! This file is part of EFTCAMB.
!
! Copyright (C) 2013-2017 by the EFTCAMB authors
!
! The EFTCAMB code is free software;
! You can use it, redistribute it, and/or modify it under the terms
! of the GNU General Public License as published by the Free Software Foundation;
! either version 3 of the License, or (at your option) any later version.
! The full text of the license can be found in the file eftcamb/LICENSE at
! the top level of the EFTCAMB distribution.
!
!----------------------------------------------------------------------------------------

program tester

    use precision
    use equispaced_linear_interpolation_1D
    use AMLutils

    implicit none

    integer , parameter :: num_points = 1000                           !< Number of points
    real(dl), parameter :: x_initial  = 0.0_dl                       !< initial point
    real(dl), parameter :: x_final    = 6.24_dl                       !< final point
    type(equispaced_linear_interpolate_function_1D) :: test_function

    integer :: i
    real(dl) :: x, y
    real(dl), dimension(num_points) :: spline_internal

    call test_function%initialize( num_points, x_initial, x_final )


    do i=1, num_points
        x = test_function%x(i)
        y = dummy_function(x)

        test_function%y(i) = y

    end do

    call test_function%initialize_derivatives()

    ! output to file:
    call CreateTxtFile( 'test.dat'  ,666 )
    do i=1, num_points
        write(666,'(18e18.10)') test_function%x(i), test_function%y(i), test_function%yp(i), test_function%ypp(i), test_function%yppp(i)
    end do
    close(666)

contains

    function dummy_function(x)
        implicit none
        real(dl), intent(in) :: x
        real(dl)  :: dummy_function
        dummy_function = sin(x)
    end function dummy_function

end program tester
