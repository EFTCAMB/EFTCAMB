! Prints storage size of TEFTCAMB_timestep_cache and offsets of the fields checked by the uninit-cache verifier,
! to compare with the python ctypes mapping in camb/eftcamb.py.
program uninit_cache_layout
    use EFTCAMB_cache
    implicit none
    type(TEFTCAMB_timestep_cache), target :: c
    integer(8) :: base
    base = loc(c)
    write(*, '(a,i0)') 'size ', storage_size(c)/8
    write(*, '(a,i0)') 'EFTcdotdot ', loc(c%EFTcdotdot) - base
    write(*, '(a,i0)') 'EFTGamma1PP ', loc(c%EFTGamma1PP) - base
    write(*, '(a,i0)') 'V_scf ', loc(c%V_scf) - base
    write(*, '(a,i0)') 'EFTcdot ', loc(c%EFTcdot) - base
end program uninit_cache_layout
