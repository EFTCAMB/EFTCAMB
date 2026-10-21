! Check ThermalNuBackground%pidot / pidotdot (EFTCAMB additions used by the EFT blocks of equations.f90)
! against finite differences of the pressure returned by ThermalNuBackground%rho_P, and drho against rho.
! Background: a(tau) with H = adotoa assumed constant over the stencil in log a; we test the log(am) derivatives:
!   pidot(am, H=1)            = dP/dln(am)
!   pidotdot(am, H=1, Hdot=0) = d2P/dln(am)^2
program nu_pidot_test
    use precision
    use MassiveNu
    implicit none
    real(dl) :: am, h, rp, pp, rm, pm, r0, p0, d1fd, d2fd, d1, d2, maxe1, maxe2, maxe3, drfd, dr
    real(dl) :: amw
    integer :: i
    maxe1 = 0; maxe2 = 0; maxe3 = 0
    h = 1.d-3
    do i = 0, 400
        am = 10._dl**(-3._dl + 7._dl*i/400._dl)
        call ThermalNuBackground%rho_P(am*exp(h), rp, pp)
        call ThermalNuBackground%rho_P(am*exp(-h), rm, pm)
        call ThermalNuBackground%rho_P(am, r0, p0)
        d1fd = (pp - pm)/(2*h)
        d2fd = (pp - 2*p0 + pm)/h**2
        drfd = (rp - rm)/(2*h)
        d1 = ThermalNuBackground%pidot(am, 1._dl, p0)
        d2 = ThermalNuBackground%pidotdot(am, 1._dl, 0._dl, p0, d1)
        dr = ThermalNuBackground%drho(am, 1._dl)
        if (abs(d1 - d1fd)/max(abs(p0), 1.d-30) > maxe1) then
            maxe1 = abs(d1 - d1fd)/max(abs(p0), 1.d-30); amw = am
        end if
        maxe2 = max(maxe2, abs(d2 - d2fd)/max(abs(p0), 1.d-30))
        maxe3 = max(maxe3, abs(dr - drfd)/max(abs(r0), 1.d-30))
        if (mod(i, 50) == 0) write(*, '(a,es10.3,a,es12.4,a,es12.4,a,es12.4,a,es12.4)') 'am=', am, ' P=', p0, &
            ' dP/dlnam: analytic=', d1, ' fd=', d1fd, '  rel.err(d2)=', abs(d2 - d2fd)/abs(p0)
    end do
    write(*, '(a,es10.3,a,es10.3)') 'max |pidot - fd|/P    = ', maxe1, ' at am=', amw
    write(*, '(a,es10.3)') 'max |pidotdot - fd|/P = ', maxe2
    write(*, '(a,es10.3)') 'max |drho - fd|/rho   = ', maxe3
end program nu_pidot_test
