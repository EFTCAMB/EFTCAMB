program interp_xfinal_unit
    use precision
    use equispaced_linear_interpolation_1D
    implicit none
    type(equispaced_linear_interpolate_function_1D) :: f
    integer :: ind, i, k, n
    real(dl) :: mu, x, v, xi
    n = 6000
    xi = log(1.d-10)
    call f%initialize(n, xi, 0._dl)
    do i = 1, n
        f%y(i) = real(i, dl)
    end do
    ! find x < x_final with formula index == N (no Fortran call that reads out of bounds yet)
    do k = 1, 2000
        x = -real(k, dl)*1.d-18
        ind = int((x - f%x_initial)/f%grid_width) + 1
        if (ind >= n) then
            write(*,'(a,es12.4,a,i6)') 'formula: x = ', x, ' (< x_final) gives ind = ', ind
            exit
        end if
    end do
    x = -1.d-16
    ind = int((x - f%x_initial)/f%grid_width) + 1
    write(*,'(a,es12.4,a,i6)') 'formula: x = ', x, ' gives ind = ', ind
    write(*,*) 'calling precompute at x=', x
    call f%precompute(x, ind, mu)
    write(*,*) 'precompute returned ind, mu =', ind, mu
    v = f%value(x)
    write(*,*) 'value =', v
end program
