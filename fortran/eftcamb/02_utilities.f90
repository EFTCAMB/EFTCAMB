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

!> @file 02_utilities.f90
!! This file contains various generic algorithms that are useful to EFTCAMB.


!----------------------------------------------------------------------------------------
!> This module contains various generic algorithms that are useful to EFTCAMB.

!> @author Bin Hu, Marco Raveri

module EFTCAMB_mixed_algorithms

    use precision

    implicit none

    private

    public hunt, integer_to_string, string, dfridr, double_NaN, log1p

    !----------------------------------------------------------------------------------------
    !> This is a utility type that allows to handle strings arrays. It is used because
    !! of partial coverage of these F2003 features that would prevent compiling with gfortran
    type :: string
        character(len=:), allocatable :: string
    end type string

    ! real number that gets represented as a Nan:
    real(dl), parameter :: double_NaN = TRANSFER((/ real(Z'00000000'), real(Z'7FF80000') /),1.0_8) !< there are situations where we need to initialize a variable to Nan...

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
        dfridr=0._dl
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
    !> Algorithm to compute log(1+x)
    function log1p( x0 )
        !    
        !   Double precision function program
        !   to compute a special exponential function,
        !       log1p(x) = log(1+x),
        !   accurately by minimax rational approximation
        !
        !   "dlog1p" keeps 16 digit accuracy for arbitrary argument
        !   while runs 20 % faster than dlog(x) itself
        !
        !   Author: Fukushima, T. <Toshio.Fukushima@nao.ac.jp>
        !   Date: 2019/10/30
        !
        implicit none
        
        real(dl), intent(in) :: x0
        real(dl) x1, x2, x4, log1p
        !
        real(dl), parameter  :: a1=0.99999999999999999405_dl
        real(dl), parameter  :: a2=2.45235728562912886048_dl
        real(dl), parameter  :: a3=2.17053627298972253249_dl
        real(dl), parameter  :: a4=0.83928994566440838378_dl
        real(dl), parameter  :: a5=0.13520496594993836479_dl
        real(dl), parameter  :: a6=0.00682631751459270270_dl
        real(dl), parameter  :: a7=0.00002291289324181940_dl
        real(dl), parameter  :: b1=2.95235728562912599232_dl
        real(dl), parameter  :: b2=3.31338158247117791600_dl
        real(dl), parameter  :: b3=1.76186164168333482938_dl
        real(dl), parameter  :: b4=0.44976458082070468584_dl
        real(dl), parameter  :: b5=0.04896199808811261680_dl
        real(dl), parameter  :: b6=0.00157389087429218809_dl
        !
        !   bit loss occurs when -1/2 <= x <= 1
        !
        if(x0.lt.-0.5_dl.or.x0.gt.1._dl) then
            log1p=log(1._dl+x0)
        !
        !   1/2 <= x < 0: log1p(x)=-log1p(x/(1+x))
        !
        elseif(x0.lt.0._dl) then
            x1=-x0/(1._dl+x0)
            x2=x1*x1; x4=x2*x2
            log1p=-x1*(((a1+x1*a2)+x2*(a3+x1*a4))+x4*((a5+x1*a6)+x2*a7)) &
                /(((1._dl+x1*b1)+x2*(b2+x1*b3))+x4*((b4+x1*b5)+x2*b6))
        !
        !   0 <= x <= 1
        !
        else
            x1=x0
            x2=x1*x1; x4=x2*x2
            log1p=x1*(((a1+x1*a2)+x2*(a3+x1*a4))+x4*((a5+x1*a6)+x2*a7)) &
                /(((1._dl+x1*b1)+x2*(b2+x1*b3))+x4*((b4+x1*b5)+x2*b6))
        endif
    end function log1p

    !----------------------------------------------------------------------------------------

end module EFTCAMB_mixed_algorithms

!----------------------------------------------------------------------------------------
