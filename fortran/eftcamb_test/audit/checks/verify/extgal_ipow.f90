! verifier: value of the integer expression 10**(-1) used as the bisection tolerance in
! 008p7_Extended_Galileon.f90:479, and the iteration count of that loop condition.
program extgal_ipow
    implicit none
    integer, parameter :: dl = kind(1.d0)
    real(dl) :: f
    integer :: counter
    print *, '10**(-1) =', 10**(-1)
    f = 1.e-3_dl
    counter = 0
    do while (sqrt(f**2) > 10**(-1) .and. counter < 50**1)
        counter = counter + 1
    end do
    print *, 'iterations with |f| = 1e-3 (tolerance would stop at once if it were 0.1):', counter
end program extgal_ipow
