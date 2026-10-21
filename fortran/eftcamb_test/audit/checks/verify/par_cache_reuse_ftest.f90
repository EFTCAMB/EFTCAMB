! Independent pure-Fortran check of finding eft-par-cache-realloc-abort (no Python wrapper involved).
! Usage: par_cache_reuse_ftest <ini> <mode>
!   mode = reuse          : second SetParams with P2 = State%CP (what results.Params.copy() does)
!   mode = reuse_dealloc  : same, but grhormass/nu_masses of the copy deallocated before the call
!   mode = user           : second SetParams with the original user parameters P
program par_cache_reuse_ftest
    use precision
    use constants, only: c
    use model
    use results
    use camb
    implicit none
    type(CAMBparams) :: P, P2
    type(CAMBdata) :: S
    character(len=1024) :: ini, mode
    integer :: err
    logical :: ok

    call get_command_argument(1, ini)
    call get_command_argument(2, mode)
    call CAMB_SetDefParams(P)
    ok = CAMB_ReadParamFile(P, trim(ini), len_trim(ini))
    if (.not. ok) stop 'ini read failed'
    write(*, '(a,i3)') 'EFTflag = ', P%EFTCAMB%EFTflag
    write(*, '(a,l2)') 'user P: allocated(eft_par_cache%grhormass) = ', allocated(P%eft_par_cache%grhormass)
    err = 0
    call S%SetParams(P, error=err)
    write(*, '(a,i3,a,l2)') 'first SetParams: error = ', err, '; allocated(State%CP%eft_par_cache%grhormass) = ', &
        allocated(S%CP%eft_par_cache%grhormass)
    write(*, '(a,f20.12)') 'first SetParams: H(z=1) [km/s/Mpc] = ', S%Hofz(1._dl)*c/1000._dl
    P2 = S%CP
    write(*, '(a,l2)') 'copy P2 = State%CP: allocated(eft_par_cache%grhormass) = ', allocated(P2%eft_par_cache%grhormass)
    select case (trim(mode))
    case ('reuse_dealloc')
        if (allocated(P2%eft_par_cache%grhormass)) deallocate(P2%eft_par_cache%grhormass)
        if (allocated(P2%eft_par_cache%nu_masses)) deallocate(P2%eft_par_cache%nu_masses)
        write(*, '(a)') 'deallocated grhormass/nu_masses of the copy'
    case ('user')
        P2 = P
        write(*, '(a)') 'second call uses the original user P'
    end select
    flush(6)
    call S%SetParams(P2, error=err)
    write(*, '(a,i3)') 'second SetParams: error = ', err
    write(*, '(a,f20.12)') 'second SetParams: H(z=1) [km/s/Mpc] = ', S%Hofz(1._dl)*c/1000._dl
    write(*, '(a)') 'DONE'
end program par_cache_reuse_ftest
