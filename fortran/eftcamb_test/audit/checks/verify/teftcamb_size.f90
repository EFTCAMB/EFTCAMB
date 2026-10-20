! Prints the byte size of TEFTCAMB and the offset of its components, to compare with the python ctypes mapping.
program teftcamb_size
    use EFTCAMB_main
    use iso_c_binding
    implicit none
    type(TEFTCAMB), target :: x
    integer(c_intptr_t) :: base
    base = loc(x)
    write(*, '(a,i0)') 'storage_size(TEFTCAMB)/8 = ', storage_size(x)/8
    write(*, '(a,i0)') 'offset EFTCAMB_write_background = ', loc(x%EFTCAMB_write_background) - base
    write(*, '(a,i0)') 'offset outroot (data)           = ', loc(x%outroot) - base
    write(*, '(a,i0)') 'offset read_param_keys          = ', loc(x%read_param_keys) - base
end program teftcamb_size
