! EFTCAMB MOD START

! -------------------------------------------------------------------------------------------------

!   EFTCAMB numerical subroutines.

!   Includes several algorithms that are used by the EFTCAMB code.
!   They are written in a general form and can be used at different purposes.

! -------------------------------------------------------------------------------------------------

!   1) Brent root finding algorithm.
!      This is used to solve numerically the equation func=funcZero by means of the Brent method.

function zbrent(func,x1,x2,tol,funcZero,succes)
    use precision
    implicit none

    real(dl) :: tol
    real(dl) :: x1,x2
    real(dl) :: funcZero
    logical  :: succes
    real(dl) :: func

    real(dl) zbrent
    external func
    real(dl),parameter :: EPS = 3.d-15
    integer ,parameter :: ITMAX = 2000  !Max Number of iterations.
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
11  continue
    succes = .false.
    zbrent=b
    return

end function zbrent

! -------------------------------------------------------------------------------------------------

!   2) Hunting algorithm.
!      This is used to efficiently search an ordered table by means of a hunting and
!      a bisection algorithm.

subroutine hunt(xx,n,x,jlo)
    ! xx  = the table to be searched (in);
    ! n   = the length of the table  (in);
    ! x   = the requested value      (in);
    ! jlo = the closest (from below) entry of the table (out).
    use precision
    implicit none
    integer jlo,n
    real(dl) x,xx(n)
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
1       jhi=jlo+inc
        if(jhi.gt.n)then
            jhi=n+1
        else if(x.ge.xx(jhi).eqv.ascnd)then
            jlo=jhi
            inc=inc+inc
            goto 1
        endif
    else
        jhi=jlo
2       jlo=jhi-inc
        if(jlo.lt.1)then
            jlo=0
        else if(x.lt.xx(jlo).eqv.ascnd)then
            jhi=jlo
            inc=inc+inc
            goto 2
        endif
    endif
3   if(jhi-jlo.eq.1)then
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

! -------------------------------------------------------------------------------------------------

!   3) Bracketing subroutine.
!      This subroutine does a outward search for the smallest intervall containing a root
!      of the equation func=funcZero

subroutine zbrac(func,x1,x2,succes,funcZero)
    ! func     = function to find the root of                           (in);
    ! x1       = lower bound of the interval in which to find the root  (in);
    ! x2       = upper bound of the interval in which to find the root.
    !            The intervall must bracket the root                    (in);
    ! success  = true if the algorithm succeeds false if not            (inout);
    ! funcZero = the value desired func = funcZero                      (in).
    use precision
    implicit none

    real(dl) func, funcZero
    real(dl) x1,x2
    logical  succes

    integer, parameter :: NTRY = 1000
    real(dl), parameter :: FACTOR = 5._dl

    real(dl) delta, temp
    real(dl) f1,f2
    integer j
    external func

    if (x1.eq.x2) stop 'you have to guess an initial range in zbrac'
    if (x2<x1) then
        temp=x2
        x2=x1
        x1=temp
    end if
    f1=func(x1)-funcZero
    f2=func(x2)-funcZero
    succes=.true.
    do 11 j=1,NTRY
        if(f1*f2.lt.0.)return
        delta=ABS(x2-x1)
        x1=x1-FACTOR*delta
        x2=x2+FACTOR*delta

        f1=func(x1)-funcZero
        f2=func(x2)-funcZero
11  continue
    succes=.false.
    return

end subroutine zbrac

! -------------------------------------------------------------------------------------------------

!   4) Neville interpolator.
!      This is used to interpolate the EFT functions once the designer/mapping
!      code has sampled them

subroutine Polint(n,xa,ya,xpl,ypl,dypl)
    ! n    = number of points in the table                          (in)
    ! xa   = first coordinate of the points to be interpolated      (in)
    ! ya   = second coordinate of the points to be interpolated     (in)
    ! xpl  = requested value of x                                   (in)
    ! ypl  = value of the interpolated function at xpl              (out)
    ! dypl =
    use precision
    implicit none

    integer , intent(in) :: n               ! Length of the table of points
    real(dl), intent(in) :: xa(n),ya(n)     ! The tables of points to be interpolated as y=f(x)
    real(dl) :: dypl,xpl,ypl

    integer  :: i,m,ns
    real(dl) :: den,dif,dift,ho,hp,wpl,cc(n),d(n)

    ns=1
    dif=abs(xpl-xa(1))
    do i=1,n
        dift=abs(xpl-xa(i))
        if (dift<dif) then
            ns=i
            dif=dift
        endif
        cc(i)=ya(i)
        d(i)=ya(i)
    end do
    ypl=ya(ns)
    ns=ns-1
    do m=1,n-1
        do i=1,n-m
            ho=xa(i)-xpl
            hp=xa(i+m)-xpl
            wpl=cc(i+1)-d(i)
            den=ho-hp
            if (den==0.) then
                write(*,*) 'failure in polint'
                stop
            end if
            den=wpl/den
            d(i)=hp*den
            cc(i)=ho*den
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

! -------------------------------------------------------------------------------------------------

!   5) Fourth order Runge-Kutta.
!      This is a very simple algorithm that is used to solve the designer equation.
!

subroutine EFT_rk4(n, y, dydx, x, h, yout, deriv)
    ! n     = dimensionality of the problem;
    ! y     = 'position' at t=x;
    ! dydx  = 'velocity' at t=x;
    ! x     = initial time;
    ! h     = time step;
    ! yout  = 'position' at t=x+h computed using fourth order Runge-Kutta;
    ! deriv = name of the subroutine that computes dydx.
    use precision
    implicit none

    real(dl) :: x, h
    integer  :: n
    real(dl), dimension(n) :: y, dydx, yout

    interface
        subroutine deriv(n, x, y, dydx)
        use precision
            implicit none
            real(dl) :: x
            integer  :: n
            real(dl), dimension(n) :: y
            real(dl), dimension(n) :: dydx
        end subroutine deriv
    end interface

    real(dl), dimension(n) :: yt, dyt,dym
    real(dl) :: hh,h6,xh
    integer :: i

    hh=h*0.5_dl
    h6=h/6._dl

    xh=x+hh
    yt=y+hh*dydx

    call deriv(n, xh, yt, dyt)

    yt=y+hh*dyt

    call deriv(n, xh, yt, dym)

    yt  =y+h*dym
    dym =dyt+dym

    call deriv(n, x+h, yt, dyt)

    yout=y+h6*(dydx+dyt+2.0*dym)

    return
end subroutine EFT_rk4

! EFTCAMB MOD END
