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

!> @file 02_utilities.f90
!! This file contains various generic algorithms that are useful to EFTCAMB.


!----------------------------------------------------------------------------------------
!> This module contains various generic algorithms that are useful to EFTCAMB.

!> @author Bin Hu, Marco Raveri

module EFTCAMB_mixed_algorithms

    use precision

    implicit none

    private

    public hunt

contains


    ! ---------------------------------------------------------------------------------------------
    !> Hunting algorithm: This is used to efficiently search an ordered table by means of a
    !> hunting and a bisection algorithm.
    subroutine hunt(xx,n,x,jlo)

        implicit none

        integer  n     !< the length of the table
        real(dl) xx(n) !< the table to be searched
        real(dl) x     !< the requested value
        real(dl) jlo   !< the index of the closest (from below) entry of the table

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

end module EFTCAMB_mixed_algorithms

!----------------------------------------------------------------------------------------
