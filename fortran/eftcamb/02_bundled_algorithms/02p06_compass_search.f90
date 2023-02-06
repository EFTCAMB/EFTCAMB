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

!> @file 02p06_compass_search.f90
!! This file contains the relevant code for the double precision compass search algorithm.
!! This code was developed by many authors that retain the copyright for the following code.
!! This source file was modified to use it with the EFTCAMB code.


!----------------------------------------------------------------------------------------
!> This module contains the relevant code for the double precision compass search algorithm.

!> @author Marco Raveri

module EFTCAMB_compass_search

    implicit none

contains

    !----------------------------------------------------------------------------------------

    subroutine compass_search ( function_handle, m, x0, delta_tol, delta_init, &
        k_max, x, fx, k )

        !*****************************************************************************80
        !
        !! COMPASS_SEARCH carries out a direct search minimization algorithm.
        !
        !  Licensing:
        !
        !    This code is distributed under the GNU LGPL license.
        !
        !  Modified:
        !
        !    05 January 2012
        !
        !  Author:
        !
        !    John Burkardt
        !
        !  Reference:
        !
        !    Tamara Kolda, Robert Michael Lewis, Virginia Torczon,
        !    Optimization by Direct Search: New Perspectives on Some Classical
        !    and Modern Methods,
        !    SIAM Review,
        !    Volume 45, Number 3, 2003, pages 385-482.
        !
        !  Parameters:
        !
        !    Input, external real ( kind = 8 ) FUNCTION_HANDLE, the name of
        !    a FORTRAN90 function which evaluates the function to be minimized, of the
        !    form FUNCTION FUNCTION_HANDLE ( M, X ).
        !
        !    Input, integer ( kind = 4 ) M, the number of variables.
        !
        !    Input, real ( kind = 8 ) X0(M), a starting estimate for the minimizer.
        !
        !    Input, real ( kind = 8 ) DELTA_TOL, the smallest step size that is allowed.
        !
        !    Input, real ( kind = 8 ) DELTA_INIT, the starting stepsize.
        !
        !    Input, integer ( kind = 4 ) K_MAX, the maximum number of steps allowed.
        !
        !    Output, real ( kind = 8 ) X(M), the estimated minimizer.
        !
        !    Output, real ( kind = 8 ) FX, the function value at X.
        !
        !    Output, integer ( kind = 4 ) K, the number of steps taken.
        !
        implicit none

        integer ( kind = 4 ) m

        logical decrease
        real ( kind = 8 ) delta
        real ( kind = 8 ) delta_init
        real ( kind = 8 ) delta_tol
        real ( kind = 8 ), external :: function_handle
        real ( kind = 8 ) fx
        real ( kind = 8 ) fxd
        integer ( kind = 4 ) i
        integer ( kind = 4 ) ii
        integer ( kind = 4 ) k
        integer ( kind = 4 ) k_max
        real ( kind = 8 ) s
        real ( kind = 8 ) x(m)
        real ( kind = 8 ) x0(m)
        real ( kind = 8 ) xd(m)

        k = 0
        x(1:m) = x0(1:m)
        fx = function_handle ( m, x )

        if ( delta_tol <= 0 ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'COMPASS_SEARCH - Fatal error!'
            write ( *, '(a)' ) '  DELTA_TOL <= 0.0.'
            write ( *, '(a,g14.6)' ) '  DELTA_TOL = ', delta_tol
            stop
        end if

        if ( delta_init <= delta_tol ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'COMPASS_SEARCH - Fatal error!'
            write ( *, '(a)' ) '  DELTA_INIT < DELTA_TOL.'
            write ( *, '(a,g14.6)' ) '  DELTA_INIT = ', delta_init
            write ( *, '(a,g14.6)' ) '  DELTA_TOL = ', delta_tol
            stop
        end if

        delta = delta_init

        do while ( k < k_max )

            k = k + 1
            !
            !  For each coordinate direction I, seek a lower function value
            !  by increasing or decreasing X(I) by DELTA.
            !
            decrease = .false.
            s = + 1.0D+00
            i = 1

            do ii = 1, 2 * m

                xd = x
                xd(i) = xd(i) + s * delta
                fxd = function_handle ( m, xd )
                !
                !  As soon as a decrease is noticed, accept the new point.
                !
                if ( fxd < fx ) then
                    x = xd
                    fx = fxd
                    decrease = .true.
                    exit
                end if

                s = - s
                if ( s == + 1.0D+00 ) then
                    i = i + 1
                end if

            end do
            !
            !  If no decrease occurred, reduce DELTA.
            !
            if ( .not. decrease ) then
                delta = delta / 2.0D+00
                if ( delta < delta_tol ) then
                    exit
                end if
            end if

        end do

        return
    end

    !----------------------------------------------------------------------------------------

    subroutine r8vec_print ( n, a, title )

        !*****************************************************************************80
        !
        !! R8VEC_PRINT prints an R8VEC.
        !
        !  Discussion:
        !
        !    An R8VEC is a vector of R8's.
        !
        !  Licensing:
        !
        !    This code is distributed under the GNU LGPL license.
        !
        !  Modified:
        !
        !    22 August 2000
        !
        !  Author:
        !
        !    John Burkardt
        !
        !  Parameters:
        !
        !    Input, integer ( kind = 4 ) N, the number of components of the vector.
        !
        !    Input, real ( kind = 8 ) A(N), the vector to be printed.
        !
        !    Input, character ( len = * ) TITLE, a title.
        !
        implicit none

        integer ( kind = 4 ) n

        real ( kind = 8 ) a(n)
        integer ( kind = 4 ) i
        character ( len = * ) title

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) trim ( title )
        write ( *, '(a)' ) ' '

        do i = 1, n
            write ( *, '(2x,i8,a,1x,g16.8)' ) i, ':', a(i)
        end do

        return
    end

    !----------------------------------------------------------------------------------------

end module EFTCAMB_compass_search

!----------------------------------------------------------------------------------------
