! Verifier unit test (spline5-order3-errorstop): drive spline5_parametrization_1D directly, linked against
! the frozen A0 camblib.so, through the same init path EFTCAMB uses (init_func_from_file + init_from_file).
! usage: spline5_o3_unit <ini> <mode>   mode = 0 (value/d1/d2 checks + candidate d3), 3 (call third), 4 (call fourth)
program spline5_o3_unit
    use precision
    use IniObjects
    use EFTCAMB_spline5_parametrizations_1D
    implicit none
    type(spline5_parametrization_1D) :: f
    type(TIniFile) :: Ini
    character(len=512) :: ininame, smode
    integer :: eft_error, mode, i, ip
    real(dl) :: a, h, fd1, fd2, fd3, d3cand, maxe1, maxe2, maxe3, r

    call get_command_argument(1, ininame)
    call get_command_argument(2, smode)
    read(smode, *) mode
    call Ini%Open(trim(ininame))
    call f%set_name('EFTOmega', '\Omega')
    eft_error = 0
    call f%init_func_from_file(Ini, eft_error)
    write(*, '(a,i3,a,i3,a,i3)') 'init_func_from_file: eft_error=', eft_error, '  spl_order=', f%spl_order, &
        '  num_pixels=', f%num_pixels
    call f%init_from_file(Ini, eft_error)
    write(*, '(a,i3,a,2i3)') 'init_from_file: eft_error=', eft_error, '  shape(bcoef)=', shape(f%bcoef)
    select case (mode)
    case (0)
        h = 1.e-4_dl
        maxe1 = 0; maxe2 = 0; maxe3 = 0
        do i = 1, 40
            a = 0.013_dl + 0.024_dl*(i - 1)   ! interior points, never within 2h of a node
            if (minval(abs(a - f%x_coord)) < 3*h) cycle
            fd1 = (f%value(a + h) - f%value(a - h))/(2*h)
            fd2 = (f%first_derivative(a + h) - f%first_derivative(a - h))/(2*h)
            fd3 = (f%second_derivative(a + h) - f%second_derivative(a - h))/(2*h)
            ! candidate fix from the finding: cubic third derivative = 6*bcoef(5,ip)
            ip = count(f%x_coord <= a)
            d3cand = 6._dl*f%bcoef(5, ip)
            maxe1 = max(maxe1, abs(fd1 - f%first_derivative(a)))
            maxe2 = max(maxe2, abs(fd2 - f%second_derivative(a)))
            maxe3 = max(maxe3, abs(fd3 - d3cand))
            if (mod(i, 8) == 1) write(*, '(a,f7.4,a,es12.4,a,es12.4,a,es12.4,a,es12.4)') 'a=', a, ' f=', f%value(a), &
                ' f1=', f%first_derivative(a), ' f2=', f%second_derivative(a), ' 6*c4=', d3cand
        end do
        write(*, '(a,es10.2,a,es10.2,a,es10.2)') 'max|FD-d1|=', maxe1, '  max|FD-d2|=', maxe2, &
            '  max|FD(d2)-6*bcoef(5,ip)|=', maxe3
    case (3)
        write(*, *) 'calling third_derivative(0.5) ...'
        flush(6)
        r = f%third_derivative(0.5_dl)
        write(*, '(a,es14.6)') 'third_derivative returned ', r
    case (4)
        write(*, *) 'calling fourth_derivative(0.5) ...'
        flush(6)
        r = f%fourth_derivative(0.5_dl)
        write(*, '(a,es14.6)') 'fourth_derivative returned ', r
    end select
end program spline5_o3_unit
