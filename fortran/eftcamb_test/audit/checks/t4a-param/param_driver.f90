! Track t4a-param: driver that instantiates EFTCAMB 1D/2D parametrizations through the
! abstract interfaces and evaluates every method on a list of points.
!
! usage: param_driver <case_file>
!
! case file (list directed, one item per line):
!   kind            'alloc' (1D via allocate_parametrized_1D_function), 'type' (1D by type name), '2d'
!   key             allocator integer flag (alloc) or type name (type / 2d)
!   fname           function name
!   flabel          function latex label
!   inipath         ini file path read by init_func_from_file / init_from_file
!   do_d3 do_d4 do_int do_feedback_noarg   logicals: evaluate third / fourth derivative / integral,
!                   call feedback() without the optional argument
!   npts
!   npts lines of x (1D) or x y (2D)
!
! Output on stdout, one record per line, prefixed by a tag.
    program param_driver

    use precision
    use IniObjects
    use EFT_def
    use EFTCAMB_abstract_parametrizations_1D
    use EFTCAMB_abstract_parametrizations_2D
    use EFTCAMB_parametrizations_1D
    use EFTCAMB_step_parametrizations_1D
    use EFTCAMB_power_law_sum_parametrizations_1D
    use EFTCAMB_neutral_parametrization_2D
    use EFTCAMB_constant_parametrization_2D
    use EFTCAMB_bilinear_parametrizations_2D

    implicit none

    character(len=512) :: case_file, kind, key, fname, flabel, inipath
    logical :: do_d3, do_d4, do_int, do_fb_noarg
    integer :: npts, i, flag, eft_error, u, ios
    real(dl), allocatable :: xs(:), ys(:), params(:)
    class(parametrized_function_1D), allocatable :: f1
    class(parametrized_function_2D), allocatable :: f2
    type(TIniFile) :: Ini
    character(len=200) :: long_name, long_latex
    character(len=EFT_names_max_length) :: short_name
    real(dl) :: pv, v, d1, d2, d3, d4, vint, before, after
    logical :: stab

    call get_command_argument(1, case_file)
    open(newunit=u, file=trim(case_file), status='old', action='read')
    read(u, '(a)') kind
    read(u, '(a)') key
    read(u, '(a)') fname
    read(u, '(a)') flabel
    read(u, '(a)') inipath
    read(u, *) do_d3, do_d4, do_int, do_fb_noarg
    read(u, *) npts
    allocate(xs(npts), ys(npts))
    ys = 0._dl
    do i = 1, npts
        if (trim(kind) == '2d') then
            read(u, *) xs(i), ys(i)
        else
            read(u, *) xs(i)
        end if
    end do
    close(u)

    eft_error = 0
    call Ini%Open(trim(inipath))

    if (trim(kind) /= '2d') then
        ! ------------------------------------------------------------------ 1D
        if (trim(kind) == 'alloc') then
            read(key, *) flag
            call allocate_parametrized_1D_function(f1, flag, trim(fname), trim(flabel), eft_error, 1)
            write(*, '(a,1x,i0)') 'ALLOC_ERROR', eft_error
            if (eft_error /= 0) stop
        else
            select case (trim(key))
            case ('zero');      allocate(zero_parametrization_1D :: f1)
            case ('wlcdm');     allocate(wDE_LCDM_parametrization_1D :: f1)
            case ('constant');  allocate(constant_parametrization_1D :: f1)
            case ('cpl');       allocate(CPL_parametrization_1D :: f1)
            case ('jbp');       allocate(JBP_parametrization_1D :: f1)
            case ('turning');   allocate(turning_point_parametrization_1D :: f1)
            case ('taylor');    allocate(taylor_parametrization_1D :: f1)
            case ('dexp');      allocate(double_exponential_parametrization_1D :: f1)
            case ('cosine');    allocate(cosine_parametrization_1D :: f1)
            case ('axion');     allocate(axion_parametrization_1D :: f1)
            case ('powsum');    allocate(power_law_sum_parametrization_1D :: f1)
            case ('step');      allocate(step_parametrization_1D :: f1)
            case ('exp2');      allocate(exponential_parametrization_2_1D :: f1)
            case default
                write(*, '(a)') 'UNKNOWN_TYPE ' // trim(key)
                stop
            end select
            call f1%set_name(trim(fname), trim(flabel))
        end if

        write(*, '(a,1x,a)') 'TYPE', trim(type_name_1d(f1))

        call f1%init_func_from_file(Ini, eft_error)
        write(*, '(a,1x,i0)') 'INITFUNC_ERROR', eft_error
        if (eft_error /= 0) stop
        call f1%init_from_file(Ini, eft_error)
        write(*, '(a,1x,i0)') 'INIT_ERROR', eft_error

        write(*, '(a,1x,i0,1x,i0)') 'NPARAM', f1%parameter_number, f1%param_number()
        allocate(params(f1%parameter_number))
        do i = 1, f1%parameter_number
            long_name = ''
            long_latex = ''
            short_name = ''
            call f1%parameter_names(i, long_name)
            call f1%parameter_names(i, short_name)
            call f1%parameter_names_latex(i, long_latex)
            call f1%parameter_value(i, pv)
            params(i) = pv
            write(*, '(a,1x,i0,1x,es26.17e3,1x,a,a,a,a,a,a)') 'PARAM', i, pv, '|', trim(long_name), '|', &
                trim(short_name), '|', trim(long_latex)
        end do

        ! feedback:
        write(*, '(a)') 'FEEDBACK_BEGIN'
        call f1%feedback(.true.)
        call f1%feedback(.false.)
        if (do_fb_noarg) call f1%feedback()
        write(*, '(a)') 'FEEDBACK_END'

        do i = 1, npts
            v = f1%value(xs(i))
            d1 = f1%first_derivative(xs(i))
            d2 = f1%second_derivative(xs(i))
            d3 = 0._dl
            d4 = 0._dl
            vint = 0._dl
            if (do_d3) d3 = f1%third_derivative(xs(i))
            if (do_d4) d4 = f1%fourth_derivative(xs(i))
            if (do_int .and. xs(i) > 0._dl) vint = f1%integral(xs(i))
            stab = f1%compute_function_stability(xs(i))
            write(*, '(a,7(1x,es26.17e3),1x,l1)') 'D', xs(i), v, d1, d2, d3, d4, vint, stab
        end do
        write(*, '(a)') 'DONE_EVAL'
        ! set/get round trip through init_parameters (after evaluation: it may stop the program):
        before = f1%value(xs(1))
        call f1%init_parameters(params)
        after = f1%value(xs(1))
        write(*, '(a,1x,es26.17e3,1x,es26.17e3)') 'ROUNDTRIP', before, after
    else
        ! ------------------------------------------------------------------ 2D
        select case (trim(key))
        case ('zero');     allocate(zero_parametrization_2D :: f2)
        case ('constant'); allocate(constant_parametrization_2D :: f2)
        case ('bilinear'); allocate(bilinear_parametrization_2D :: f2)
        case default
            write(*, '(a)') 'UNKNOWN_TYPE ' // trim(key)
            stop
        end select
        call f2%set_name(trim(fname), trim(flabel))
        write(*, '(a,1x,a)') 'TYPE', trim(key)
        call f2%init_from_file(Ini, eft_error)
        write(*, '(a,1x,i0)') 'INIT_ERROR', eft_error
        write(*, '(a,1x,i0,1x,i0)') 'NPARAM', f2%parameter_number, f2%param_number()
        allocate(params(f2%parameter_number))
        do i = 1, f2%parameter_number
            long_name = ''
            long_latex = ''
            short_name = ''
            call f2%parameter_names(i, long_name)
            call f2%parameter_names(i, short_name)
            call f2%parameter_names_latex(i, long_latex)
            call f2%parameter_value(i, pv)
            params(i) = pv
            write(*, '(a,1x,i0,1x,es26.17e3,1x,a,a,a,a,a,a)') 'PARAM', i, pv, '|', trim(long_name), '|', &
                trim(short_name), '|', trim(long_latex)
        end do
        write(*, '(a)') 'FEEDBACK_BEGIN'
        call f2%feedback(.true.)
        call f2%feedback(.false.)
        if (do_fb_noarg) call f2%feedback()
        write(*, '(a)') 'FEEDBACK_END'
        before = f2%value(xs(1), ys(1))
        call f2%init_parameters(params)
        after = f2%value(xs(1), ys(1))
        write(*, '(a,1x,es26.17e3,1x,es26.17e3)') 'ROUNDTRIP', before, after
        do i = 1, npts
            write(*, '(a,8(1x,es26.17e3))') 'D2', xs(i), ys(i), f2%value(xs(i), ys(i)), &
                f2%first_derivative_x(xs(i), ys(i)), f2%first_derivative_y(xs(i), ys(i)), &
                f2%second_derivative_x(xs(i), ys(i)), f2%second_derivative_y(xs(i), ys(i)), &
                f2%second_derivative_xy(xs(i), ys(i))
        end do
        write(*, '(a)') 'DONE_EVAL'
    end if
    write(*, '(a)') 'DONE'

    contains

    function type_name_1d(f) result(tn)
    class(parametrized_function_1D) :: f
    character(len=64) :: tn
    select type (f)
    type is (zero_parametrization_1D);              tn = 'zero_parametrization_1D'
    type is (wDE_LCDM_parametrization_1D);          tn = 'wDE_LCDM_parametrization_1D'
    type is (constant_parametrization_1D);          tn = 'constant_parametrization_1D'
    type is (linear_parametrization_1D);            tn = 'linear_parametrization_1D'
    type is (power_law_parametrization_1D);         tn = 'power_law_parametrization_1D'
    type is (exponential_parametrization_1D);       tn = 'exponential_parametrization_1D'
    type is (CPL_parametrization_1D);               tn = 'CPL_parametrization_1D'
    type is (JBP_parametrization_1D);               tn = 'JBP_parametrization_1D'
    type is (turning_point_parametrization_1D);     tn = 'turning_point_parametrization_1D'
    type is (taylor_parametrization_1D);            tn = 'taylor_parametrization_1D'
    type is (taylorseries_parametrization_1D);      tn = 'taylorseries_parametrization_1D'
    type is (padeseries_parametrization_1D);        tn = 'padeseries_parametrization_1D'
    type is (fourier_parametrization_1D);           tn = 'fourier_parametrization_1D'
    type is (exponential_parametrization_2_1D);     tn = 'exponential_parametrization_2_1D'
    type is (double_exponential_parametrization_1D); tn = 'double_exponential_parametrization_1D'
    type is (cosine_parametrization_1D);            tn = 'cosine_parametrization_1D'
    type is (axion_parametrization_1D);             tn = 'axion_parametrization_1D'
    type is (spline_parametrization_1D);            tn = 'spline_parametrization_1D'
    type is (step_parametrization_1D);              tn = 'step_parametrization_1D'
    type is (steplog_parametrization_1D);           tn = 'steplog_parametrization_1D'
    type is (power_law_sum_parametrization_1D);     tn = 'power_law_sum_parametrization_1D'
    type is (spline5_parametrization_1D);           tn = 'spline5_parametrization_1D'
    class default;                                  tn = 'UNKNOWN'
    end select
    end function type_name_1d

    end program param_driver
