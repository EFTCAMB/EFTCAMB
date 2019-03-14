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

!> @file 02_utilities.f90
!! This file contains various generic algorithms that are useful to EFTCAMB.


!----------------------------------------------------------------------------------------
!> This module contains various generic algorithms that are useful to EFTCAMB.

!> @author Bin Hu, Marco Raveri

module EFTCAMB_mixed_algorithms

    use precision

    implicit none

    private

    public hunt, integer_to_string, string, dfridr

    !----------------------------------------------------------------------------------------
    !> This is a utility type that allows to handle strings arrays. It is used because
    !! of partial coverage of these F2003 features that would prevent compiling with gfortran
    type :: string
        character(len=:), allocatable :: string
    end type string

contains

    ! ---------------------------------------------------------------------------------------------
    !> Hunting algorithm: This is used to efficiently search an ordered table by means of a
    !> hunting and a bisection algorithm.
    subroutine hunt(xx,n,x,jlo)

        implicit none

        integer  :: n      !< the length of the table
        real(dl) :: xx(n)  !< the table to be searched
        real(dl) :: x      !< the requested value
        integer  :: jlo    !< the index of the closest (from below) entry of the table

        integer inc,jhi,jm
        logical ascnd
        ascnd=xx(n).ge.xx(1)
        if(jlo.le.0.or.jlo.gt.n)then
            jlo=0
            jhi=n+1
            goto 3
        endif
        inc=1
        if(x.ge.xx(jlo).eqv.ascnd)then
1           jhi=jlo+inc
            if(jhi.gt.n)then
                jhi=n+1
            else if(x.ge.xx(jhi).eqv.ascnd)then
                jlo=jhi
                inc=inc+inc
                goto 1
            endif
        else
            jhi=jlo
2           jlo=jhi-inc
            if(jlo.lt.1)then
                jlo=0
            else if(x.lt.xx(jlo).eqv.ascnd)then
                jhi=jlo
                inc=inc+inc
                goto 2
            endif
        endif
3       if(jhi-jlo.eq.1)then
            if(x.eq.xx(n))jlo=n-1
            if(x.eq.xx(1))jlo=1
            return
        endif
        jm=(jhi+jlo)/2
        if(x.ge.xx(jm).eqv.ascnd)then
            jlo=jm
        else
            jhi=jm
        endif
        goto 3

    end subroutine hunt

    ! ---------------------------------------------------------------------------------------------
    !> This function converts an integer to a string. Usefull for numbered files output.
    function integer_to_string( number )

        implicit none

        integer, intent(in) :: number               !< Input integer number
        character(10)       :: integer_to_string    !< Output string with the number

        write( integer_to_string, '(i10)' ) number

        integer_to_string = TRIM(ADJUSTL( integer_to_string ))

    end function integer_to_string

    !----------------------------------------------------------------------------------------
    !> Algorithm to compute the numerical derivative as in the Numerical Recepies.
    function dfridr( func, x, h, err )

        implicit none

        integer,parameter :: ntab = 100
        real(dl) dfridr,err,h,x,func
        external func

        real(dl), parameter :: CON=1.4_dl       ! decrease of the stepsize.
        real(dl), parameter :: CON2=CON*CON
        real(dl), parameter :: BIG=1.d+30
        real(dl), parameter :: SAFE=2._dl

        integer i,j
        real(dl) errt, fac, hh
        real(dl), dimension(ntab,ntab) :: a

        if (h.eq.0._dl) h = 1.d-8

        hh=h
        a(1,1)=(func(x+hh)-func(x-hh))/(2.0_dl*hh)
        err=BIG

        do 12 i=2,NTAB
            hh=hh/CON
            a(1,i)=(func(x+hh)-func(x-hh))/(2.0_dl*hh)
            fac=CON2
            do 11 j=2,i
                a(j,i)=(a(j-1,i)*fac-a(j-1,i-1))/(fac-1._dl)
                fac=CON2*fac
                errt=max(abs(a(j,i)-a(j-1,i)),abs(a(j,i)-a(j-1,i-1)))
                if (errt.le.err) then
                    err=errt
                    dfridr=a(j,i)
                endif
11          continue
            if(abs(a(i,i)-a(i-1,i-1)).ge.SAFE*err)return
12      continue
        return
    end function dfridr

    !----------------------------------------------------------------------------------------

end module EFTCAMB_mixed_algorithms

!----------------------------------------------------------------------------------------
