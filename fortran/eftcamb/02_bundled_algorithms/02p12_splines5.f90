module splines5

    ! This module is based on cubic spline code written by John E. Pask, LLNL and at https://github.com/certik/fortran-utils

    ! Raised to quintic spline by Gen Ye 
    
    ! Splines are fully specified by the interpolation points, except that
    ! at the ends, we have the freedom to prescribe the third and fourth derivatives.
    ! If we know a derivative at an end (exactly), then best is to impose that.
    ! Otherwise, it is better to use the "consistent" end conditions: the third and fourth
    ! derivatives are determined such that it is smooth (i.e: interpolating a polynomial with base {x^n} using the first 6 points counting from the end and then use its derivatives at boundary).
    !
    ! High level API: spline5, spline5ders.
    ! Low level API: the rest of public soubroutines.
    !
    ! Use the high level API to obtain quintic spline fit with consistent boundary
    ! conditions and optionally the derivatives. Use the low level API if more fine
    ! grained control is needed.
    

    use precision, only: dl
    !use lapack, only: dgesv, dgbsv
    implicit none
    private
    public spline5pars, spline5valder, iix, iixmin, iixun, iixexp, poly5, dpoly5, &
        d2poly5, d3poly5, d4poly5, spline5, spline5ders

contains

    function spline5(x, y, xnew) result(ynew)
        ! Takes the function values 'y' on the grid 'x' and returns new values 'ynew'
        ! at the given grid 'xnew' using cubic splines interpolation with such
        ! boundary conditions so that the 2nd derivative is consistent with the
        ! interpolating cubic.
        real(dl), intent(in) :: x(:), y(:), xnew(:)
        real(dl) :: ynew(size(xnew))
        real(dl) :: c(0:6, size(x)-1)
        integer :: i, ip
        ! get spline parameters: 2nd derivs at ends determined by cubic interpolation
        call spline5pars(x, y, [2, 2], [0._dl, 0._dl, 0._dl, 0._dl], c)

        ip = 0
        do i = 1, size(xnew)
            ip = iixmin(xnew(i), x, ip)
            ynew(i) = poly5(xnew(i), c(:, ip))
        end do
    end function

    subroutine spline5ders(x, y, xnew, ynew, dynew, d2ynew, d3ynew, d4ynew)
        ! Just like 'spline', but also calculate 1st and 2nd derivatives
        real(dl), intent(in) :: x(:), y(:), xnew(:)
        real(dl), intent(out), optional :: ynew(:), dynew(:), d2ynew(:), d3ynew(:), d4ynew(:)
        real(dl) :: c(0:6, size(x)-1)
        integer :: i, ip
        call spline5pars(x, y, [2, 2], [0._dl, 0._dl, 0._dl, 0._dl], c)

        ip = 0
        do i = 1, size(xnew)
            ip = iixmin(xnew(i), x, ip)
            if (present(  ynew))   ynew(i) =   poly5(xnew(i), c(:, ip))
            if (present( dynew))  dynew(i) =  dpoly5(xnew(i), c(:, ip))
            if (present(d2ynew)) d2ynew(i) = d2poly5(xnew(i), c(:, ip))
            if (present(d3ynew)) d3ynew(i) = d3poly5(xnew(i), c(:, ip))
            if (present(d4ynew)) d4ynew(i) = d4poly5(xnew(i), c(:, ip))
        end do
    end subroutine

    subroutine spline5pars(xi,yi,bctype,bcval,c)
        ! Returns parameters c defining quintic spline interpolating x-y data xi, yi, with
        ! boundary conditions specified by bcytpe, bcvals
        real(dl), intent(in):: xi(:)        ! x values of data
        real(dl), intent(in):: yi(:)        ! y values of data
        integer, intent(in):: bctype(2)     ! type of boundary condition at each end:
           ! bctype(1) = type at left end, bctype(2) = type at right end.
           ! 1 = specified 1st and 2nd derivatives, 2 = 1st and 2nd derivatives consistent with interpolating quintic.
        real(dl), intent(in):: bcval(4)     ! boundary condition values at each end:
           ! bcval(1:2) = 3rd and 4th derivative values at left end, bcval(3:4) = values at right end
        real(dl), intent(out):: c(0:,:)     ! parameters defining spline: c(i,j) = ith parameter of jth
           ! spline polynomial, p_j = sum_{i=1}^6 c(i,j) (x-c(0,j))^(i-1), j = 1..n-1, n = # of data pts.
           ! dimensions: c(0:6,1:n-1)
        real(dl) As(4*size(c,2),4*size(c,2))   ! spline eq. matrix
        real(dl) bs(4*size(c,2))               ! spline eq. rhs vector
        real(dl) cs(4*size(c,2))               ! spline eq. solution vector
        real(dl) hi(size(c,2))                 ! spline intervals
        real(dl) dyi(size(c,2))                ! dyi(i)=yi(i+1)-yi(i)
        real(dl) Ae(6,6)                       ! end-quintic eq. matrix
        real(dl) be(6)                         ! end-quintic eq. rhs vector
        real(dl) ce(6)                         ! end-quintic eq. solution vector
        real(dl) xe(6),ye(6)                   ! x,y values at ends
        real(dl) d1p1,d2p1,d1pn,d2pn           ! 1st and 2nd derivatives at ends
        real(dl) x0                            ! expansion center
        real(dl) c1,c2,c3,c4,c5,c6             ! expansion coefficients
        integer n                              ! number of data points
        integer i,j,i4
        ! lapack variables
        integer ipiv(6),ipiv2(4*size(c,2))
        real(dl) bemat(6,1),bmat(4*size(c,2),1)
        integer info

        ! check input parameters
        if (bctype(1) < 1 .or. bctype(1) > 2) stop "spline5pars error: bctype /= 1 or 2."
        if (bctype(2) < 1 .or. bctype(2) > 2) stop "spline5pars error: bctype /= 1 or 2."
        if (size(c,1) /= 7) stop "spline5pars error: size(c,1) /= 7."
        if (size(c,2) /= size(xi)-1) stop "spline5pars error: size(c,2) /= size(xi)-1."
        if (size(xi) /= size(yi)) stop "spline5pars error: size(xi) /= size(yi)"

        ! To get rid of compiler warnings:
        d1p1 = 0
        d1pn = 0
        d2p1 = 0
        d2pn = 0

        ! initializations
        n=size(xi)
        do i=1,n-1
            hi(i)=xi(i+1)-xi(i)
            dyi(i)=yi(i+1)-yi(i)
        end do

        ! compute interpolating-cubic 3rd and 4th derivs at ends, if required
           ! left end
        if(bctype(1)==2) then
            if (n < 6) stop "spline5pars error: n < 6"
            xe=xi(1:6)
            ye=yi(1:6)
            x0=xe(1) ! center at end
            do i=1,6
                do j=1,6
                    Ae(i,j) = (xe(i)-x0)**(j-1)
                end do
            end do
            Ae(:,1) = 1    ! set 0^0 = 1
            be=ye; bemat(:,1)=be
            call dgesv(6, 1, Ae, 6, ipiv, bemat, 6, info)
            if (info /= 0) stop "spline5pars error: dgesv error."
            ce=bemat(:,1)
            d1p1=ce(2)
            d2p1=2._dl*ce(3)
        end if
           ! right end
        if(bctype(2)==2) then
            if (n < 6) stop "spline5pars error: n < 6"
            xe=xi(n-5:n)
            ye=yi(n-5:n)
            x0=xe(6) ! center at end
            do i=1,6
                do j=1,6
                    Ae(i,j) = (xe(i)-x0)**(j-1)
                end do
            end do
            Ae(:,1) = 1    ! set 0^0 = 1
            be=ye; bemat(:,1)=be
            call dgesv(6, 1, Ae, 6, ipiv, bemat, 6, info)
            if (info /= 0) stop "spline5pars error: dgesv error."
            ce=bemat(:,1)
            d1pn=ce(2)
            d2pn=2._dl*ce(3)
        end if

        ! set 1st and 2nd derivs at ends
        if(bctype(1)==1) then
            d1p1=bcval(1)
            d2p1=bcval(2)
        end if
        if(bctype(2)==1) then
            d1pn=bcval(3)
            d2pn=bcval(4)
        end if
        !write(*,*) d1p1,d2p1,d1pn,d2pn

        ! construct spline equations -- LAPACK GE matrix form
        ! basis: 
        !    * H1(t) = 1 - 10*t^3 + 15*t^4 - 6*t^5
        !    * H2(t) = t - 6*t^3 + 8*t^4 - 3*t^5
        !    * H3(t) = t^2/2 - 3*t^3/2 + 3*t^4/2 - t^5/2
        !    * H4(t) = H3(1-t)
        !    * H5(t) = -H2(1-t)
        !    * H6(t) = H1(1-t)
        ! where t=(x - x_i)/h_i is defined on interval [x_i,x_{i+1}] of length h_i = x_{i+1}-x_i
        
        As=0._dl
            ! boundary condition
        As(1,1)=1._dl/hi(1)
        bs(1)=d1p1
        As(2,2)=1._dl/hi(1)**2
        bs(2)=d2p1
        As(3,4*(n-1)-1)=1._dl/hi(n-1)**2
        bs(3)=d2pn
        As(4,4*(n-1))=1._dl/hi(n-1)
        bs(4)=d1pn
           ! internal knot conditions
        do i=2,n-1
            i4=4*i
            As(i4,i4)     = 168._dl/hi(i)**4
            As(i4,i4-1)   = -24._dl/hi(i)**4
            As(i4,i4-2)   = 36._dl/hi(i)**4
            As(i4,i4-3)   = 192._dl/hi(i)**4
            As(i4,i4-4)   = 192._dl/hi(i-1)**4
            As(i4,i4-5)   = -36._dl/hi(i-1)**4
            As(i4,i4-6)   = 24._dl/hi(i-1)**4
            As(i4,i4-7)   = 168._dl/hi(i-1)**4
            bs(i4)        = 360._dl*(dyi(i-1)/hi(i-1)**4 + dyi(i)/hi(i)**4)
            As(i4-1,i4)   = -24._dl/hi(i)**3
            As(i4-1,i4-1) = 3._dl/hi(i)**3
            As(i4-1,i4-2) = -9._dl/hi(i)**3
            As(i4-1,i4-3) = -36._dl/hi(i)**3
            As(i4-1,i4-4) = 36._dl/hi(i-1)**3
            As(i4-1,i4-5) = -9._dl/hi(i-1)**3
            As(i4-1,i4-6) = 3._dl/hi(i-1)**3
            As(i4-1,i4-7) = 24._dl/hi(i-1)**3
            bs(i4-1)      = 60._dl*(dyi(i-1)/hi(i-1)**3 - dyi(i)/hi(i)**3)
            As(i4-2,i4-2) = 1._dl/hi(i)**2
            As(i4-2,i4-5) = -1._dl/hi(i-1)**2
            bs(i4-2)      = 0._dl
            As(i4-3,i4-3) = 1._dl/hi(i)
            As(i4-3,i4-4) = -1._dl/hi(i-1)
            bs(i4-3)      = 0._dl
        end do

        ! solve spline equations -- full matrix
        bmat(:,1)=bs
        call dgesv(4*(n-1), 1, As, 4*(n-1), ipiv2, bmat, 4*(n-1), info)
        if (info /= 0) stop "spline5pars error: dgesv error."
        cs=bmat(:,1)
        !write(*,*) cs(1:6)

        ! transform to (x-x0)^(i-1) basis and return
        do i=1,n-1
            i4 = 4*i
            ! coefficients in spline basis:
            c1=yi(i)
            c6=yi(i+1)
            c2=cs(i4-3)
            c3=cs(i4-2)
            c4=cs(i4-1)
            c5=cs(i4)
            ! coefficients in (x-x0)^(i-1) basis
            c(0,i)=xi(i)
            c(1,i)=c1
            c(2,i)=c2/hi(i)
            c(3,i)=c3/hi(i)**2/2._dl
            c(4,i)=(-10._dl*c1 - 6._dl*c2 - (3._dl*c3)/2._dl + c4/2._dl - 4._dl*c5 + 10._dl*c6)/hi(i)**3
            c(5,i)=(15._dl*c1 + 8._dl*c2 + (3._dl*c3)/2._dl - c4 + 7._dl*c5 - 15._dl*c6)/hi(i)**4
            c(6,i)=(-6._dl*c1 - 3._dl*c2 - c3/2._dl + c4/2._dl - 3._dl*c5 + 6._dl*c6)/hi(i)**5
        end do
    end subroutine

    !--------------------------------------------------------------------------------------------------!

    subroutine spline5valder(x,xi,c,val,der)
        ! Returns value and 1st derivative of spline defined by knots xi and parameters c
        ! returned by spline5pars
        real(dl), intent(in):: x            ! point at which to evaluate spline
        real(dl), intent(in):: xi(:)        ! spline knots (x values of data)
        real(dl), intent(in):: c(0:,:)      ! spline parameters: c(i,j) = ith parameter of jth
           ! spline polynomial, p_j = sum_{i=1}^4 c(i,j) (x-c(0,j))^(i-1), j = 1..n-1, n = # of data pts.
           ! dimensions: c(0:4,1:n-1)
        real(dl), intent(out):: val         ! value of spline at x
        real(dl), intent(out):: der         ! 1st derivative of spline at x
        integer n                           ! number of knots
        integer i1

        ! initialize, check input parameters
        n=size(xi)
        if (size(c,1) /= 7) stop "spline5 error: size(c,1) /= 7."
        if (size(c,2) /= size(xi)-1) stop "spline3 error: size(c,2) /= size(xi)-1."
        ! find interval containing x
        i1=iix(x,xi)
        ! return value and derivative
        val=poly5(x,c(:,i1))
        der=dpoly5(x,c(:,i1))
    end subroutine

    !--------------------------------------------------------------------------------------------------!

    integer function iix(x, xi) result(i1)
        ! Returns index i of interval [xi(i),xi(i+1)] containing x in mesh xi,
        ! with intervals indexed by left-most points.
        ! N.B.: x outside [x1,xn] are indexed to nearest end.
        ! Uses bisection, except if "x" lies in the first or second elements (which is
        ! often the case)
        real(dl), intent(in) :: x            ! target value
        real(dl), intent(in) :: xi(:)        ! mesh, xi(i) < xi(i+1)
        integer n                            ! number of mesh points
        integer i2, ic

        n = size(xi)
        i1 = 1
        if (n < 2) then
            stop "error in iix: n < 2"
        elseif (n == 2) then
            i1 = 1
        elseif (n == 3) then
            if (x <= xi(2)) then ! first element
                i1 = 1
            else
                i1 = 2
            end if
        elseif (x <= xi(1)) then ! left end
            i1 = 1
        elseif (x <= xi(2)) then ! first element
            i1 = 1
        elseif (x <= xi(3)) then ! second element
            i1 = 2
        elseif (x >= xi(n)) then  ! right end
            i1 = n-1
        else
            ! bisection: xi(i1) <= x < xi(i2)
            i1 = 3; i2 = n
            do
                if (i2 - i1 == 1) exit
                ic = i1 + (i2 - i1)/2
                if (x >= xi(ic)) then
                    i1 = ic
                else
                    i2 = ic
                endif
            end do
        end if
    end function

    integer function iixmin(x, xi, i_min) result(ip)
        ! Just like iix, but assumes that x >= xi(i_min)
        real(dl), intent(in) :: x, xi(:)
        integer, intent(in) :: i_min
        if (i_min >= 1 .and. i_min <= size(xi)-1) then
            ip = iix(x, xi(i_min:)) + i_min - 1
        else
            ip = iix(x, xi)
        end if
    end function

    !--------------------------------------------------------------------------------------------------!

    function iixun(x,n,x1,xn)
        ! Returns index i of interval [x(i),x(i+1)] containing x in uniform mesh defined by
        !   x(i) = x1 + (i-1)/(n-1)*(xn-x1), i = 1 .. n,
        ! with intervals indexed by left-most points.
        ! N.B.: x outside [x1,xn] are indexed to nearest end.
        integer iixun                       ! index i of interval [x(i),x(i+1)] containing x
        real(dl), intent(in):: x            ! target value
        integer, intent(in):: n             ! number of mesh points
        real(dl), intent(in):: x1           ! initial point of mesh
        real(dl), intent(in):: xn           ! final point of mesh
        integer i

        ! compute index
        i=int((x-x1)/(xn-x1)*(n-1))+1
        ! reset if ouside 1..n
        if (i<1) i=1
        if (i>n-1) i=n-1
        iixun=i
    end function

    !--------------------------------------------------------------------------------------------------!

    function iixexp(x,n,x1,alpha,beta)
        ! Returns index i of interval [x(i),x(i+1)] containing x in exponential mesh defined by
        !   x(i) = x1 + alpha [ exp(beta(i-1)) - 1 ], i = 1 .. n,
        ! where alpha = (x(n) - x(1))/[ exp(beta(n-1)) - 1 ],
        ! beta = log(r)/(n-2), r = (x(n)-x(n-1))/(x(2)-x(1)) = ratio of last to first interval,
        ! and intervals indexed by left-most points.
        ! N.B.: x outside [x1,xn] are indexed to nearest end.
        integer iixexp                      ! index i of interval [x(i),x(i+1)] containing x
        real(dl), intent(in):: x            ! target value
        integer , intent(in):: n             ! number of mesh points
        real(dl), intent(in):: x1           ! initial point of mesh
        real(dl), intent(in):: alpha        ! mesh parameter:
        !   x(i) = x1 + alpha [ exp(beta(i-1)) - 1 ], i = 1 .. n,
        ! where alpha = (x(n) - x(1))/[ exp(beta(n-1)) - 1 ],
        ! beta = log(r)/(n-2), r = (x(n)-x(n-1))/(x(2)-x(1)) = ratio of last to first interval,
        real(dl), intent(in):: beta         ! mesh parameter
        integer i

        ! compute index
        i=int(log((x-x1)/alpha + 1)/beta) + 1
        ! reset if outside 1..n
        if (i<1) i=1
        if (i>n-1) i=n-1
        iixexp=i
    end function

    !--------------------------------------------------------------------------------------------------!

    function poly5(x,c)
        ! returns value of polynomial \sum_{i=1}^6 c(i) (x-c(0))^(i-1)
        real(dl) poly5
        real(dl), intent(in):: x      ! point at which to evaluate polynomial
        real(dl), intent(in):: c(0:)  ! coefficients: poly = \sum_{i=1}^6 c(i) (x-c(0))^(i-1)
        real(dl) dx
        dx=x-c(0)
        poly5=c(1)+c(2)*dx+c(3)*dx**2+c(4)*dx**3+c(5)*dx**4+c(6)*dx**6
    end function

    !--------------------------------------------------------------------------------------------------!

    function dpoly5(x,c)
        ! returns 1st derivative of polynomial \sum_{i=1}^6 c(i) (x-c(0))^(i-1)
        real(dl) dpoly5
        real(dl), intent(in):: x      ! point at which to evaluate polynomial
        real(dl), intent(in):: c(0:)  ! coefficients: poly = \sum_{i=1}^4 c(i) (x-c(0))^(i-1)
        real(dl) dx
        dx=x-c(0)
        dpoly5=c(2)+2._dl*c(3)*dx+3._dl*c(4)*dx**2+4._dl*c(5)*dx**3+5._dl*c(6)*dx**4
    end function

    !--------------------------------------------------------------------------------------------------!

    function d2poly5(x,c)
        ! returns 2nd derivative of polynomial \sum_{i=1}^6 c(i) (x-c(0))^(i-1)
        real(dl) d2poly5
        real(dl), intent(in):: x      ! point at which to evaluate polynomial
        real(dl), intent(in):: c(0:)  ! coefficients: poly = \sum_{i=1}^6 c(i) (x-c(0))^(i-1)
        real(dl) dx
        dx=x-c(0)
        d2poly5=2._dl*c(3)+6._dl*c(4)*dx+12._dl*c(5)*dx**2+20._dl*c(6)*dx**3
    end function

    !--------------------------------------------------------------------------------------------------!

    function d3poly5(x,c)
        ! returns 3rd derivative of polynomial \sum_{i=1}^6 c(i) (x-c(0))^(i-1)
        real(dl) d3poly5
        real(dl), intent(in):: x      ! point at which to evaluate polynomial
        real(dl), intent(in):: c(0:)  ! coefficients: poly = \sum_{i=1}^6 c(i) (x-c(0))^(i-1)
        real(dl) dx
        dx=x-c(0)
        d3poly5=6._dl*c(4)+24._dl*c(5)*dx+60._dl*c(6)*dx**2
    end function

    !--------------------------------------------------------------------------------------------------!

    function d4poly5(x,c)
        ! returns 4th derivative of polynomial \sum_{i=1}^6 c(i) (x-c(0))^(i-1)
        real(dl) d4poly5
        real(dl), intent(in):: x      ! point at which to evaluate polynomial
        real(dl), intent(in):: c(0:)  ! coefficients: poly = \sum_{i=1}^6 c(i) (x-c(0))^(i-1)
        real(dl) dx
        dx=x-c(0)
        d4poly5=24._dl*c(5)+120._dl*c(6)*dx
    end function

end module
