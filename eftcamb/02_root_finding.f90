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

!> @file 02_root_finding.f90
!! This file contains the EFTCAMB root finding algorithms.


!----------------------------------------------------------------------------------------
!> This module contains the definitions of all the EFTCAMB root finding algorithms.

!> @author Bin Hu, Marco Raveri

module EFTCAMB_rootfind

    use precision

    implicit none

    private

    public zbrac, zbrent

contains

    ! ---------------------------------------------------------------------------------------------
    !> Bracketing subroutine: This subroutine does a outward search for the smallest intervall
    !! containing a root of the equation func=funcZero
    subroutine zbrac( func, x1, x2, succes, funcZero )

        implicit none

        real(dl) func                         !< function to find the root of
        real(dl), intent(inout)  :: x1        !< in input lower bound of the interval in which to find the root,
                                              !< in output lower bound of the bracketing interval
        real(dl), intent(inout)  :: x2        !< in input upper bound of the interval in which to find the root,
                                              !< in output upper bound of the bracketing interval
        logical , intent(out)    :: succes    !< true if the algorithm succeeds false if not
        real(dl), intent(in)     :: funcZero  !< the value desired func = funcZero

        integer , parameter      :: ntry = 1000     !< number of subsequent enlargement of the intervall
        real(dl), parameter      :: FACTOR = 5._dl  !< factor that is used to expand the intervall

        real(dl) delta, temp, f1,f2
        integer  j
        external func

        ! initial check:
        if ( x1.eq.x2 ) stop 'you have to guess an initial range in zbrac'
        ! get an ordered interval:
        if (x2<x1) then
            temp=x2
            x2=x1
            x1=temp
        end if
        ! compute the function:
        f1 = func(x1) -funcZero
        f2 = func(x2) -funcZero

        succes=.true.
        do j=1, ntry
            ! check wether the interval is bracketing:
            if ( f1*f2 .lt. 0._dl ) return
            ! compute the interval and enlarge it:
            delta=ABS(x2-x1)
            x1 = x1 -FACTOR*delta
            x2 = x2 +FACTOR*delta
            ! recompute the functions:
            f1 = func(x1) -funcZero
            f2 = func(x2) -funcZero
        end do
        succes=.false.

        return

    end subroutine zbrac

    ! ---------------------------------------------------------------------------------------------
    !> Brent root finding algorithm:  This is used to solve numerically the equation
    !! func=funcZero by means of the Brent method. Notice that the initial interval has to bracket
    !! the root of the function.
    function zbrent(func,x1,x2,tol,funcZero,succes)

        implicit none

        real(dl) :: func      !< function to find the root of
        real(dl) :: x1        !< in input lower bound of the interval in which to find the root. Has to be braketing.
        real(dl) :: x2        !< in input upper bound of the interval in which to find the root. Has to be braketing.
        real(dl) :: tol       !< numerical absolute tollerance. Notice that internally the algorithm mitigates this with
                              !< absoute tollerance. Values below 1.d-15 are useless.
        real(dl) :: funcZero  !< the value desired func = funcZero
        logical  :: succes    !< true if the algorithm succeeds false if not

        real(dl) :: zbrent    !< returns the value of x at which the function satisfies func = funcZero up to numerical accuracy

        external func

        integer ,parameter :: ITMAX = 2000   !< Max Number of iterations.
        real(dl),parameter :: EPS   = 3.d-15 !< Relative tollerance.
        integer iter
        real(dl) a,b,c,d,e,fa,fb,fc,p,q,r,s,tol1,xm

        succes = .true.
        a=x1
        b=x2
        fa=func(a)-funcZero
        fb=func(b)-funcZero
        if((fa.gt.0..and.fb.gt.0.).or.(fa.lt.0..and.fb.lt.0.)) then
            succes=.false.
            return
        end if
        c=b
        fc=fb
        do 11 iter=1,ITMAX
            if((fb.gt.0..and.fc.gt.0.).or.(fb.lt.0..and.fc.lt.0.))then
                c=a
                fc=fa
                d=b-a
                e=d
            endif
            if(abs(fc).lt.abs(fb)) then
                a=b
                b=c
                c=a
                fa=fb
                fb=fc
                fc=fa
            endif
            tol1=2.*EPS*abs(b)+0.5*tol
            xm=.5*(c-b)
            if(abs(xm).le.tol1 .or. fb.eq.0.)then
                zbrent=b
                return
            endif
            if(abs(e).ge.tol1 .and. abs(fa).gt.abs(fb)) then
                s=fb/fa
                if(a.eq.c) then
                    p=2.*xm*s
                    q=1.-s
                else
                    q=fa/fc
                    r=fb/fc
                    p=s*(2.*xm*q*(q-r)-(b-a)*(r-1.))
                    q=(q-1.)*(r-1.)*(s-1.)
                endif
                if(p.gt.0.) q=-q
                p=abs(p)
                if(2.*p .lt. min(3.*xm*q-abs(tol1*q),abs(e*q))) then
                    e=d
                    d=p/q
                else
                    d=xm
                    e=d
                endif
            else
                d=xm
                e=d
            endif
            a=b
            fa=fb
            if(abs(d) .gt. tol1) then
                b=b+d
            else
                b=b+sign(tol1,xm)
            endif
            fb=func(b)-funcZero
11      continue
        succes = .false.
        zbrent=b
        return

    end function zbrent

end module EFTCAMB_rootfind

!----------------------------------------------------------------------------------------
