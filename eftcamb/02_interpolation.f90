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

!> @file 02_interpolation.f90
!! This file contains the EFTCAMB interpolation algorithms.


!----------------------------------------------------------------------------------------
!> This module contains the definitions of all the EFTCAMB interpolation algorithms.

!> @author Bin Hu, Marco Raveri

module EFTCAMB_interpolation

    use precision

    implicit none

    private

    public Polint

contains

    ! ---------------------------------------------------------------------------------------------
    !> Neville interpolator: computes polynomial interpolation of a set of points
    subroutine Polint(n,xa,ya,xpl,ypl,dypl)

        implicit none

        integer , intent(in) :: n     !< number of points in the table
        real(dl), intent(in) :: xa(n) !< first coordinate of the points to be interpolated
        real(dl), intent(in) :: ya(n) !< second coordinate of the points to be interpolated y=f(x)
        real(dl) :: xpl               !< requested value of x
        real(dl) :: ypl               !< output value of the interpolated function at xpl
        real(dl) :: dypl              !< output error estimate on the value of the interpolated function at xpl

        integer  :: i,m,ns
        real(dl) :: den,dif,dift,ho,hp,wpl,cc(n),d(n)

        ns  = 1
        dif = abs( xpl-xa(1) )
        do i=1,n
            dift=abs(xpl-xa(i))
            if (dift<dif) then
                ns=i
                dif=dift
            endif
            cc(i)=ya(i)
            d(i)=ya(i)
        end do
        ypl = ya(ns)
        ns  = ns -1
        do m=1, n-1
            do i=1, n-m
                ho  = xa(i)-xpl
                hp  = xa(i+m)-xpl
                wpl = cc(i+1)-d(i)
                den = ho-hp
                if (den==0.) then
                    write(*,*) 'failure in polint'
                    stop
                end if
                den   = wpl/den
                d(i)  = hp*den
                cc(i) = ho*den
            end do
            if (2*ns<n-m) then
                dypl=cc(ns+1)
            else
                dypl=d(ns)
                ns=ns-1
            endif
            ypl=ypl+dypl
        end do

        return

    end subroutine Polint

end module EFTCAMB_interpolation

!----------------------------------------------------------------------------------------
