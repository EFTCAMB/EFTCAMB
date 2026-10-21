! Verifier (finding stability-cache-cdotdot-uninitialised): independent pure-Fortran check, linked against a frozen
! release camblib.so (no ctypes, no special -finit flags in the library).
!
! usage: cdotdot_stack <test> <fill> [ini]
!   test = init      : fill a TEFTCAMB_timestep_cache with <fill>, call %initialize(), list the fields %is_nan flags
!   test = setparams : read ini, force stability check on, dirty 1 MB of stack with <fill>, call CAMBdata%SetParams
!   test = direct    : read ini, SetParams with skip_stability, dirty stack with <fill>, call EFTCAMB_Stability_Check
!   fill = zero | nan | one
module cdd_dirty
    implicit none
contains
    subroutine dirty(val)
        real(8), intent(in) :: val
        ! fixed-size local; -fopenmp implies -frecursive so this lives on the stack
        real(8), volatile :: buf(131072)
        integer :: i
        do i = 1, size(buf)
            buf(i) = val
        end do
    end subroutine dirty
end module cdd_dirty

program cdotdot_stack
    use, intrinsic :: ieee_arithmetic
    use, intrinsic :: iso_c_binding, only: c_loc, c_f_pointer
    use precision
    use config
    use CAMB
    use results
    use EFTCAMB_cache
    use EFTCAMB_stability
    use cdd_dirty
    implicit none
    character(len=1024) :: ini, test, fillstr
    type(CAMBparams) :: P
    type(CAMBdata), target :: Res
    type(TEFTCAMB_timestep_cache), target :: tc
    real(dl), pointer :: raw(:)
    real(dl) :: fill, kmax
    logical :: ok, have, success
    integer :: err, nd, i

    call get_command_argument(1, test)
    call get_command_argument(2, fillstr)
    call get_command_argument(3, ini)
    select case (trim(fillstr))
    case ('nan')
        fill = ieee_value(1._dl, ieee_quiet_nan)
    case ('one')
        fill = 1._dl
    case default
        fill = 0._dl
    end select

    select case (trim(test))
    case ('init')
        nd = storage_size(tc)/storage_size(1._dl)
        call c_f_pointer(c_loc(tc), raw, [nd])   ! raw aliases the bytes of tc (the name transfer is shadowed by CAMB)
        raw = fill
        call tc%initialize()
        write(*, '(a,i0,a,i0)') 'RESULT init ndoubles=', nd, ' n_nonzero_after_initialize=', count(raw /= 0._dl)
        write(*, '(a,*(1x,i0))') 'RESULT init nonzero_indices=', pack([(i, i=1, nd)], raw /= 0._dl)
        write(*, '(a,2(1x,i0))') 'RESULT init index_of_EFTcdotdot_EFTGamma1PP=', &
            int((loc(tc%EFTcdotdot) - loc(tc))/8) + 1, int((loc(tc%EFTGamma1PP) - loc(tc))/8) + 1
        call tc%is_nan(have, 3)
        write(*, '(a,L2)') 'RESULT init HaveNan=', have
    case ('setparams', 'direct')
        ok = CAMB_ReadParamFile(P, trim(ini), len_trim(ini))
        if (.not. ok) stop 'read failed'
        if (trim(test) == 'setparams') then
            P%EFTCAMB%EFTCAMB_skip_stability = .false.
            call dirty(fill)
            call Res%SetParams(P, err)
            write(*, '(a,i0,2a)') 'RESULT setparams err=', err, ' msg=', trim(global_error_message)
        else
            P%EFTCAMB%EFTCAMB_skip_stability = .true.
            call Res%SetParams(P, err)
            if (err /= 0) stop 'SetParams failed'
            kmax = 10._dl
            call dirty(fill)
            call EFTCAMB_Stability_Check(success, Res%CP%EFTCAMB, Res%CP%eft_par_cache, &
                Res%CP%EFTCAMB%EFTCAMB_stability_time, 1._dl, kmax)
            write(*, '(a,L2)') 'RESULT direct stable=', success
        end if
    end select
end program cdotdot_stack
