! Audit check t8b-integration: compare the analytic massive-neutrino pressure derivatives
! (EFTCAMB MOD pidot/pidotdot/pidotdotdot in massive_neutrinos.f90) with finite differences
! of ThermalNuBackground%rho_P in log(a m), over the full am range including branch boundaries.
program nu_pressure_derivs
    use precision
    use MassiveNu
    implicit none
    integer :: i, n
    real(dl) :: am, lnam, h, rho, p(-3:3), d1fd, d2fd, d3fd, d1, d2, d3, adotoa, hdot, hdotdot
    real(dl) :: maxe1, maxe2, maxe3, e1, e2, e3, am_w1, am_w2, am_w3
    real(dl) :: amb(2), jl(3), jr(3)
    integer :: j
    n = 4000
    maxe1 = 0; maxe2 = 0; maxe3 = 0
    adotoa = 1._dl; hdot = 0._dl; hdotdot = 0._dl
    h = 1.e-3_dl
    do i = 0, n
        lnam = log(1.e-3_dl) + (log(1.e4_dl) - log(1.e-3_dl))*i/real(n, dl)
        am = exp(lnam)
        ! skip points whose stencil straddles a branch boundary
        if (abs(lnam - log(0.42_dl)) < 4*h .or. abs(lnam - log(70._dl)) < 4*h) cycle
        do j = -3, 3
            call ThermalNuBackground%rho_P(exp(lnam + j*h), rho, p(j))
        end do
        d1fd = (p(1) - p(-1))/(2*h)
        d2fd = (p(1) - 2*p(0) + p(-1))/h**2
        d3fd = (p(2) - 2*p(1) + 2*p(-1) - p(-2))/(2*h**3)
        ! with adotoa=1, Hdot=Hdotdot=0, pidot = dP/dln(am), pidotdot = d2P/dln(am)^2, ...
        d1 = ThermalNuBackground%pidot(am, adotoa, p(0))
        d2 = ThermalNuBackground%pidotdot(am, adotoa, hdot, p(0), d1)
        d3 = ThermalNuBackground%pidotdotdot(am, adotoa, hdot, hdotdot, p(0), d1, d2)
        e1 = abs(d1 - d1fd)/max(abs(d1fd), 1e-6_dl*abs(p(0)))
        e2 = abs(d2 - d2fd)/max(abs(d2fd), 1e-4_dl*abs(p(0)))
        e3 = abs(d3 - d3fd)/max(abs(d3fd), 1e-3_dl*abs(p(0)))
        if (e1 > maxe1) then; maxe1 = e1; am_w1 = am; end if
        if (e2 > maxe2) then; maxe2 = e2; am_w2 = am; end if
        if (e3 > maxe3) then; maxe3 = e3; am_w3 = am; end if
    end do
    write(*, '(a,es10.3,a,es10.3)') 'max rel err dP/dlnam      = ', maxe1, ' at am = ', am_w1
    write(*, '(a,es10.3,a,es10.3)') 'max rel err d2P/dlnam2    = ', maxe2, ' at am = ', am_w2
    write(*, '(a,es10.3,a,es10.3)') 'max rel err d3P/dlnam3    = ', maxe3, ' at am = ', am_w3
    ! chain rule check: with nonzero Hdot, Hdotdot, compare against numerical tau-derivatives
    ! P(tau) with ln a(tau) = adotoa*tau + Hdot*tau^2/2 + Hdotdot*tau^3/6 (so d(adotoa)/dtau = Hdot ...)
    adotoa = 0.7_dl; hdot = -0.3_dl; hdotdot = 0.2_dl
    am = 3._dl
    h = 1.e-3_dl
    do j = -3, 3
        call ThermalNuBackground%rho_P(am*exp(adotoa*j*h + hdot*(j*h)**2/2 + hdotdot*(j*h)**3/6), rho, p(j))
    end do
    d1fd = (p(1) - p(-1))/(2*h)
    d2fd = (p(1) - 2*p(0) + p(-1))/h**2
    d3fd = (p(2) - 2*p(1) + 2*p(-1) - p(-2))/(2*h**3)
    d1 = ThermalNuBackground%pidot(am, adotoa, p(0))
    d2 = ThermalNuBackground%pidotdot(am, adotoa, hdot, p(0), d1)
    d3 = ThermalNuBackground%pidotdotdot(am, adotoa, hdot, hdotdot, p(0), d1, d2)
    write(*, '(a,3es12.4)') 'tau-chain rule analytic   = ', d1, d2, d3
    write(*, '(a,3es12.4)') 'tau-chain rule finite diff= ', d1fd, d2fd, d3fd
    ! branch-boundary jumps of the derivatives (value just left vs just right)
    amb = [0.42_dl, 70._dl]
    adotoa = 1._dl; hdot = 0._dl; hdotdot = 0._dl
    do i = 1, 2
        am = amb(i)*(1 - 1e-9_dl)
        call ThermalNuBackground%rho_P(am, rho, p(0))
        jl(1) = ThermalNuBackground%pidot(am, adotoa, p(0))
        jl(2) = ThermalNuBackground%pidotdot(am, adotoa, hdot, p(0), jl(1))
        jl(3) = ThermalNuBackground%pidotdotdot(am, adotoa, hdot, hdotdot, p(0), jl(1), jl(2))
        am = amb(i)*(1 + 1e-9_dl)
        call ThermalNuBackground%rho_P(am, rho, p(1))
        jr(1) = ThermalNuBackground%pidot(am, adotoa, p(1))
        jr(2) = ThermalNuBackground%pidotdot(am, adotoa, hdot, p(1), jr(1))
        jr(3) = ThermalNuBackground%pidotdotdot(am, adotoa, hdot, hdotdot, p(1), jr(1), jr(2))
        write(*, '(a,f6.2,a,es10.3,a,3es11.3)') 'boundary am=', amb(i), ' P jump(rel)=', abs(p(1) - p(0))/p(0), &
            ' dP1,dP2,dP3 rel jumps=', abs(jr - jl)/max(abs(jl), 1e-30_dl)
    end do
end program nu_pressure_derivs
