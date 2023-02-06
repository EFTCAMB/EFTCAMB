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

!> @file 10_EFTCAMB_background_output.f90
!! This file contains the code that prints to file the background evolution of observables.


!----------------------------------------------------------------------------------------
!> This module contains the code that prints to file the background evolution of observables.

!> @author Marco Raveri

module EFTCAMB_background_output

    use precision, only: dl
    use constants, only: c
    use FileUtils
    use results

    implicit none

    private

    public output_background

contains

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that saves to file the background evolution. Assumes that CAMB was already called
    !! to perform all initialization.
    subroutine output_background( outroot, state )

        implicit none

        character(LEN=*), intent(in)  :: outroot  !< the output file root
        Type(CAMBdata)                :: state    !< camb data state

        ! hard coded parameters:
        real(dl), parameter :: minimum_scale_factor_1 = 1.e-5  !< minimum scale factor at which the code computes background observables log spaced
        real(dl), parameter :: minimum_scale_factor_2 = 0.01   !< minimum scale factor at which the code computes background observables linear spaced
        real(dl), parameter :: maximum_scale_factor   = 1.0    !< maximum scale factor at which the code computes background observables
        integer , parameter :: number_points          = 20     !< number of points, evenly spaced in logarithmic space in a from minimum_scale_factor_1 to minimum_scale_factor_2 and linearly in a to maximum_scale_factor

        integer  :: ind, unit
        logical  :: is_open
        real(dl) :: rs_drag
        real(dl), dimension(number_points) :: a, z, tau, Hz, r, D_L, D_A, DVoRs, DMoRs
        type(TTextFile) :: out_file

        ! open the output file:
        call out_file%CreateFile( TRIM(outroot)//'background.dat' )
        ! compute the a spacing log and then linear:
        do ind=1, number_points/2
            a(ind) = exp( log(minimum_scale_factor_1) +REAL(ind-1)/REAL(number_points/2-1)*( log(minimum_scale_factor_2)-log(minimum_scale_factor_1) ) )
        end do
        do ind=1, number_points/2
            a(number_points/2+ind) = minimum_scale_factor_2 +REAL(ind)/REAL(number_points/2)*( maximum_scale_factor - minimum_scale_factor_2 )
        end do

        ! compute redshift:
        z = 1._dl/a -1._dl
        ! compute tau:
        call state%TimeOfzArr( tau, z, number_points )
        ! compute H(z):
        call state%HofzArr( Hz, z, number_points )
        Hz = Hz*c/1000._dl
        ! compute radial distance:
        call state%ComovingRadialDistanceArr( r, z(number_points:1:-1), number_points, 1e-4_dl )
        r = r(number_points:1:-1)
        ! angular diameter distance:
        call state%AngularDiameterDistanceArr( D_A, z(number_points:1:-1), number_points )
        D_A = D_A(number_points:1:-1)
        ! luminosity distance:
        D_L = D_A*(1+z)**2
        ! get sound horizon:
        rs_drag = state%ThermoDerivedParams( derived_rdrag )
        ! get the BAO DV:
        do ind=1, number_points
            DVoRs(ind) = BAO_D_v_from_DA_H( z(ind), D_A(ind), Hz(ind)*1000._dl/c )
        end do
        ! get DM:
        DMoRs   = (1+z)*D_A

        ! save to file:
        write(out_file%unit,'(a)') '# a z tau r Hz DL DA DV DM HzRs DAoRs DVoRs DMoRs'
        do ind=1, number_points
            write(out_file%unit,'(200ES20.10E3)') a(ind), z(ind), tau(ind), r(ind), Hz(ind), D_L(ind), D_A(ind), DVoRs(ind), DMoRs(ind), Hz(ind)*rs_drag, D_A(ind)/rs_drag, DVoRs(ind)/rs_drag, DMoRs(ind)/rs_drag
        end do

        call out_file%Close()

    end subroutine output_background

    !----------------------------------------------------------------------------------------

end module EFTCAMB_background_output

!----------------------------------------------------------------------------------------
