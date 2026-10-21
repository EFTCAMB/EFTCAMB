! Track 10a: accuracy of the massive-neutrino pressure derivatives used by EFTCAMB.
! Compares ThermalNuBackground%pidot/pidotdot/pidotdotdot (with adotoa=1, Hdot=Hdotdot=0, i.e. derivatives with
! respect to log(a m)) and %drho with high-order finite differences of %rho_P on a log grid of a*m.
! Build against a tree's camblib.so and module directory (-DB1API for the pre-rebase API that needs init).
program nu_deriv_test
    use constants
    use MassiveNu
    implicit none
    integer, parameter :: n = 600
    integer :: i
    real(dl) :: lam, am, h, r0, p0, d1n, d2n, d3n, d1a, d2a, d3a, drn, dra
    real(dl) :: rp(-3:3), pp(-3:3)
    integer :: j
#ifdef B1API
    call ThermalNuBackground%init()
#endif
    h = 2.e-3_dl
    write(*, '(a)') '# am P dP1_ana dP1_num dP2_ana dP2_num dP3_ana dP3_num drho_ana drho_num'
    do i = 0, n
        lam = log(1.e-3_dl) + (log(3.e3_dl) - log(1.e-3_dl))*i/real(n, dl)
        am = exp(lam)
        do j = -3, 3
            call ThermalNuBackground%rho_P(exp(lam + j*h), rp(j), pp(j))
        end do
        r0 = rp(0); p0 = pp(0)
        ! 6th-order central differences in log(am)
        d1n = (-pp(-3) + 9*pp(-2) - 45*pp(-1) + 45*pp(1) - 9*pp(2) + pp(3))/(60*h)
        d2n = (2*pp(-3) - 27*pp(-2) + 270*pp(-1) - 490*pp(0) + 270*pp(1) - 27*pp(2) + 2*pp(3))/(180*h**2)
        d3n = (pp(-3) - 8*pp(-2) + 13*pp(-1) - 13*pp(1) + 8*pp(2) - pp(3))/(8*h**3)
        drn = (-rp(-3) + 9*rp(-2) - 45*rp(-1) + 45*rp(1) - 9*rp(2) + rp(3))/(60*h)
        d1a = ThermalNuBackground%pidot(am, 1._dl, p0)
        d2a = ThermalNuBackground%pidotdot(am, 1._dl, 0._dl, p0, d1a)
        d3a = ThermalNuBackground%pidotdotdot(am, 1._dl, 0._dl, 0._dl, p0, d1a, d2a)
        dra = ThermalNuBackground%drho(am, 1._dl)
        write(*, '(11es24.15)') am, p0, d1a, d1n, d2a, d2n, d3a, d3n, dra, drn, r0
    end do
end program nu_deriv_test
