"""Apply the t4a-param candidate fixes to a COPY of an audit tree (never to the repo).

usage: python make_fix_tree.py TREE_COPY [--only name,name]
Each patch is an exact (old, new) replacement and must match once.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

P1 = "fortran/eftcamb/04f_parametrizations_1D/"

PATCHES: dict[str, list[tuple[str, str, str]]] = {
    # Pade fourth derivative: the 4 u' g''' term lost its u' factor
    "pade_d4": [(
        P1 + "04p010_pade_series_1D.f90",
        "        4._dl*(-ddddown/down**2 - 6._dl*ddown**3/down**4 + 6._dl*ddown*dddown/down**3) + &",
        "        4._dl*dup*(-ddddown/down**2 - 6._dl*ddown**3/down**4 + 6._dl*ddown*dddown/down**3) + &",
    )],
    # quintic spline value: last term must be c(6)*dx**5
    "poly5": [(
        "fortran/eftcamb/02_bundled_algorithms/02p12_splines5.f90",
        "    poly5 = c(1) + c(2)*dx + c(3)*dx**2 + c(4)*dx**3 + c(5)*dx**4 + c(6)*dx**6",
        "    poly5 = c(1) + c(2)*dx + c(3)*dx**2 + c(4)*dx**3 + c(5)*dx**4 + c(6)*dx**5",
    )],
    # steplog integral: no division by log(x) (exact a=1 gave 1^Inf = 1 for one factor only)
    "steplog_int": [(
        P1 + "04p018_steplog_parametrizations_1D.f90",
        """        term1 = x**(-1._dl - 1.5_dl*(self%v1 + self%v2) - &
            1.5_dl*(self%v2 - self%v1)*(log(x) - 2*log(self%aT) + log(self%aT)**2/log(x))/ &
            sqrt(self%Delta**2 + (log(x) - log(self%aT))**2))
        term2 = exp(1.5_dl*(self%v1 - self%v2)* &
            (self%Delta**2/sqrt(self%Delta**2 + (log(x) - log(self%aT))**2) - sqrt(self%Delta**2 + log(self%aT)**2)))
        steplogParametrized1DIntegral = term1*term2""",
        """        term1 = x**(-1._dl - 1.5_dl*(self%v1 + self%v2))
        term2 = exp(-1.5_dl*(self%v2 - self%v1)* &
            (sqrt(self%Delta**2 + (log(x) - log(self%aT))**2) - sqrt(self%Delta**2 + log(self%aT)**2)))
        steplogParametrized1DIntegral = term1*term2""",
    )],
    # steplog fourth derivative: chain rule in N = ln x
    "steplog_d4": [(
        P1 + "04p018_steplog_parametrizations_1D.f90",
        """        steplogParametrized1DFourthDerivative = (1._dl/2._dl)*(self%v2 - self%v1)* &
            (-(3._dl*log(x/self%aT))/(self%Delta**2*x**4*(1._dl - log(x/self%aT)/self%Delta**2)**(1.5_dl)) - &
            6._dl/(x**4*sqrt(1._dl - log(x/self%aT)/self%Delta**2)) + &
            11._dl/(self%Delta**2*x**4*(1._dl - log(x/self%aT)/self%Delta**2)**(1.5_dl)) + &
            (105._dl*log(x/self%aT))/(16._dl*self%Delta**8*x**4*(1._dl - log(x/self%aT)/self%Delta**2)**(4.5_dl)) - &
            (45._dl*log(x/self%aT))/(4._dl*self%Delta**6*x**4*(1._dl - log(x/self%aT)/self%Delta**2)**(3.5_dl)) + &
            15._dl/(2._dl*self%Delta**6*x**4*(1._dl - log(x/self%aT)/self%Delta**2)**(3.5_dl)) + &
            (33._dl*log(x/self%aT))/(4._dl*self%Delta**4*x**4*(1._dl - log(x/self%aT)/self%Delta**2)**(2.5_dl)) - &
            27._dl/(2._dl*self%Delta**4*x**4*(1._dl - log(x/self%aT)/self%Delta**2)**(2.5_dl)))/self%Delta""",
        """        term1 = (log(x) - log(self%aT))/self%Delta
        steplogParametrized1DFourthDerivative = 0.5_dl*(self%v2 - self%v1)/x**4*( &
            15._dl*term1*(3._dl - 4._dl*term1**2)/(1._dl + term1**2)**4.5_dl/self%Delta**4 &
            - 6._dl*(12._dl*term1**2 - 3._dl)/(1._dl + term1**2)**3.5_dl/self%Delta**3 &
            - 11._dl*3._dl*term1/(1._dl + term1**2)**2.5_dl/self%Delta**2 &
            - 6._dl/(1._dl + term1**2)**1.5_dl/self%Delta)""",
    )],
    # step fourth derivative: Delta^2 in the numerator, not Delta^-7
    "step_d4": [(
        P1 + "04p017_step_parametrizations_1D.f90",
        "        (1._dl/(self%Delta**7*((self%aT - x)**2 + self%Delta**2)**(4.5_dl)))",
        "        (self%Delta**2/((self%aT - x)**2 + self%Delta**2)**(4.5_dl))",
    )],
    # cubic spline: evaluate the value spline at the end nodes and differentiate the value spline itself
    "spline": [
        (P1 + "04p016_spline_parametrizations_1D.f90",
         """    Spline1DValue = self%null_value
    if (x <= self%interpolation%x_initial .or. x >= self%interpolation%x_final) return""",
         """    Spline1DValue = self%null_value
    if (x < self%interpolation%x_initial .or. x > self%interpolation%x_final) return"""),
        (P1 + "04p016_spline_parametrizations_1D.f90",
         """    Spline1DFirstDerivative = self%null_value
    if (x <= self%interpolation%x_initial .or. x >= self%interpolation%x_final) return
    ! if inside return the cubic spline value:
    ip = iixmin(x, self%x_coord, 0)
    Spline1DFirstDerivative = poly3 (x, self%cyp(:, ip))""",
         """    Spline1DFirstDerivative = 0._dl
    if (x < self%interpolation%x_initial .or. x > self%interpolation%x_final) return
    ! if inside return the cubic spline value:
    ip = iixmin(x, self%x_coord, 0)
    Spline1DFirstDerivative = dpoly3 (x, self%cy(:, ip))"""),
        (P1 + "04p016_spline_parametrizations_1D.f90",
         """    Spline1DSecondDerivative = self%null_value
    if (x <= self%interpolation%x_initial .or. x >= self%interpolation%x_final) return
    ! if inside return the cubic spline value:
    ip = iixmin(x, self%x_coord, 0)
    Spline1DSecondDerivative = poly3 (x, self%cypp(:, ip))""",
         """    Spline1DSecondDerivative = 0._dl
    if (x < self%interpolation%x_initial .or. x > self%interpolation%x_final) return
    ! if inside return the cubic spline value:
    ip = iixmin(x, self%x_coord, 0)
    Spline1DSecondDerivative = d2poly3 (x, self%cy(:, ip))"""),
        (P1 + "04p016_spline_parametrizations_1D.f90",
         """    Spline1DThirdDerivative = self%null_value
    if (x <= self%interpolation%x_initial .or. x >= self%interpolation%x_final) return
    ! if inside return the cubic spline value:
    ip = iixmin(x, self%x_coord, 0)
    Spline1DThirdDerivative = poly3 (x, self%cyppp(:, ip))""",
         """    Spline1DThirdDerivative = 0._dl
    if (x < self%interpolation%x_initial .or. x > self%interpolation%x_final) return
    ! if inside return the cubic spline value:
    ip = iixmin(x, self%x_coord, 0)
    Spline1DThirdDerivative = 6._dl*self%cy(5, ip)"""),
        (P1 + "04p016_spline_parametrizations_1D.f90",
         """    Spline1DFourthDerivative = self%null_value
    if (x <= self%interpolation%x_initial .or. x >= self%interpolation%x_final) return
    ! if inside return the cubic spline value:
    ip = iixmin(x, self%x_coord, 0)
    Spline1DFourthDerivative = poly3 (x, self%cypppp(:, ip))""",
         """    Spline1DFourthDerivative = 0._dl
    if (x < self%interpolation%x_initial .or. x > self%interpolation%x_final) return
    ! a cubic spline has zero fourth derivative inside each interval:
    ip = iixmin(x, self%x_coord, 0)
    Spline1DFourthDerivative = 0._dl"""),
    ],
    # spline integrals: EFTCAMB convention (a^2 rho_DE/rho_DE0) computed by Gauss-Legendre from value()
    "spline_int": [
        (P1 + "04p016_spline_parametrizations_1D.f90", '    real(dl) :: x_temp\n    integer  :: ind, ip\n\n    Spline1DIntegral = 0._dl\n    if (x <= self%interpolation%x_initial .or. x >= self%interpolation%x_final) return\n    ! if inside return the cubic spline value:\n    ip = iixmin(x, self%x_coord, 0)\n    Spline1DIntegral = poly3 (x, self%cyint(:, ip)) !< the output value\n', '    real(dl), parameter :: gx(4) = [0.1834346424956498_dl, 0.5255324099163290_dl, 0.7966664774136267_dl, &\n        0.9602898564975363_dl]\n    real(dl), parameter :: gw(4) = [0.3626837833783620_dl, 0.3137066458778873_dl, 0.2223810344533745_dl, &\n        0.1012285362903763_dl]\n    real(dl) :: bounds(size(self%x_coord) + 2), n_lo, n_hi, nb_i, c, h, s\n    integer  :: nb, i, j, k, nsub\n\n    ! I(a) = a^2 exp(3 int_{ln a}^0 (1 + w(N)) dN): composite 8-point Gauss-Legendre in N = ln a, split at the nodes\n    Spline1DIntegral = 0._dl\n    if (x <= 0._dl) return\n    n_lo = min(log(x), 0._dl)\n    n_hi = max(log(x), 0._dl)\n    nb = 1\n    bounds(1) = n_lo\n    do i = 1, size(self%x_coord)\n        nb_i = -huge(1._dl)\n        if (self%x_coord(i) > 0._dl) nb_i = log(self%x_coord(i))\n        if (nb_i > n_lo .and. nb_i < n_hi) then\n            nb = nb + 1\n            bounds(nb) = nb_i\n        end if\n    end do\n    nb = nb + 1\n    bounds(nb) = n_hi\n    s = 0._dl\n    do i = 1, nb - 1\n        nsub = max(1, ceiling((bounds(i + 1) - bounds(i))/0.25_dl))\n        h = (bounds(i + 1) - bounds(i))/nsub\n        do k = 1, nsub\n            c = bounds(i) + (k - 0.5_dl)*h\n            do j = 1, 4\n                s = s + gw(j)*0.5_dl*h*(2._dl + self%value(exp(c + 0.5_dl*h*gx(j))) + &\n                    self%value(exp(c - 0.5_dl*h*gx(j))))\n            end do\n        end do\n    end do\n    if (log(x) > 0._dl) s = -s\n    Spline1DIntegral = x**2*exp(3._dl*s)\n'),
        (P1 + "04p019_spline5_parametrizations_1D.f90", '    real(dl) :: tmp, a, x_, x0, dx\n    integer  :: ip, j\n\n    SplineQ1DIntegral = 0._dl\n\n    ! initialize integration vector when first called\n    if (.not. allocated(self%yint)) then\n        call SplineQ1DInitIntegration(self)\n    end if\n\n    a = x\n    x_ = a\n    if (self%logx_flag) then\n        x_ = log(a)\n    end if\n    ! Do not extrapolate:\n    if (x_ < self%x_initial) then\n        SplineQ1DIntegral = self%yint(1)\n    else if (x_ > self%x_final) then\n        SplineQ1DIntegral = 1._dl\n    else\n        ! simply use linear interpolation\n        ip = iixmin(x_, self%x_coord, 0)\n        x0 = self%x_coord(ip)\n        dx = x_ - x0\n        tmp = 0._dl\n        do j = 1, size(self%bcoef(:, ip))\n            tmp = tmp + self%bcoef(j, ip)*dx**j/real(j)\n        end do\n        SplineQ1DIntegral = self%yint(ip)*(x_/x0)**3*exp(tmp)\n    end if\n', '    real(dl), parameter :: gx(4) = [0.1834346424956498_dl, 0.5255324099163290_dl, 0.7966664774136267_dl, &\n        0.9602898564975363_dl]\n    real(dl), parameter :: gw(4) = [0.3626837833783620_dl, 0.3137066458778873_dl, 0.2223810344533745_dl, &\n        0.1012285362903763_dl]\n    real(dl) :: bounds(size(self%x_coord) + 2), n_lo, n_hi, nb_i, c, h, s\n    integer  :: nb, i, j, k, nsub\n\n    ! I(a) = a^2 exp(3 int_{ln a}^0 (1 + w(N)) dN): composite 8-point Gauss-Legendre in N = ln a, split at the nodes\n    SplineQ1DIntegral = 0._dl\n    if (x <= 0._dl) return\n    n_lo = min(log(x), 0._dl)\n    n_hi = max(log(x), 0._dl)\n    nb = 1\n    bounds(1) = n_lo\n    do i = 1, size(self%x_coord)\n        nb_i = -huge(1._dl)\n        if (self%logx_flag) then\n            nb_i = self%x_coord(i)\n        else if (self%x_coord(i) > 0._dl) then\n            nb_i = log(self%x_coord(i))\n        end if\n        if (nb_i > n_lo .and. nb_i < n_hi) then\n            nb = nb + 1\n            bounds(nb) = nb_i\n        end if\n    end do\n    nb = nb + 1\n    bounds(nb) = n_hi\n    s = 0._dl\n    do i = 1, nb - 1\n        nsub = max(1, ceiling((bounds(i + 1) - bounds(i))/0.25_dl))\n        h = (bounds(i + 1) - bounds(i))/nsub\n        do k = 1, nsub\n            c = bounds(i) + (k - 0.5_dl)*h\n            do j = 1, 4\n                s = s + gw(j)*0.5_dl*h*(2._dl + self%value(exp(c + 0.5_dl*h*gx(j))) + &\n                    self%value(exp(c - 0.5_dl*h*gx(j))))\n            end do\n        end do\n    end do\n    if (log(x) > 0._dl) s = -s\n    SplineQ1DIntegral = x**2*exp(3._dl*s)\n'),
    ],
    # quintic spline: derivatives of the constant null_value outside the nodes are 0, not null_value
    "spline5_outside": [
        (P1 + "04p019_spline5_parametrizations_1D.f90", "    SplineQ1DFirstDerivative = self%null_value\n",
         "    SplineQ1DFirstDerivative = 0._dl\n"),
        (P1 + "04p019_spline5_parametrizations_1D.f90", "    SplineQ1DSecondDerivative = self%null_value\n",
         "    SplineQ1DSecondDerivative = 0._dl\n"),
        (P1 + "04p019_spline5_parametrizations_1D.f90", "    SplineQ1DThirdDerivative = self%null_value\n",
         "    SplineQ1DThirdDerivative = 0._dl\n"),
        (P1 + "04p019_spline5_parametrizations_1D.f90", "    SplineQ1DFourthDerivative = self%null_value\n",
         "    SplineQ1DFourthDerivative = 0._dl\n"),
    ],
    # quintic-spline type with the accepted Spline_Order = 3: cubic third derivative, zero fourth derivative
    "spline5_o3": [
        (P1 + "04p019_spline5_parametrizations_1D.f90",
         """        write(*, '(a,a,a,I3)') 'Third derivative invalid for ', trim(self%name) // '_Spline_Order = ', self%spl_order
        write(*, *) 'Calculations cannot proceed.'
        call MpiStop('EFTCAMB error')
""",
         """        if (self%logx_flag .or. self%logy_flag) then
            y = poly3(x_, self%bcoef(:, ip))
            yx = dpoly3(x_, self%bcoef(:, ip))
            yxx = d2poly3(x_, self%bcoef(:, ip))
        end if
        yxxx = 6._dl*self%bcoef(5, ip)
"""),
        (P1 + "04p019_spline5_parametrizations_1D.f90",
         """        write(*, '(a,a,a,I3)') 'Fourth derivative invalid for ', trim(self%name) // '_Spline_Order = ', self%spl_order
        write(*, *) 'Calculations cannot proceed.'
        call MpiStop('EFTCAMB error')
""",
         """        if (self%logx_flag .or. self%logy_flag) then
            y = poly3(x_, self%bcoef(:, ip))
            yx = dpoly3(x_, self%bcoef(:, ip))
            yxx = d2poly3(x_, self%bcoef(:, ip))
            yxxx = 6._dl*self%bcoef(5, ip)
        end if
        yxxxx = 0._dl
"""),
    ],
    # 2D abstract feedback assigned to an absent optional argument
    "feedback2d": [
        ("fortran/eftcamb/05_abstract_parametrizations_2D.f90",
         """    character(len=EFT_names_max_length) :: param_name

    if (.not. present(print_params)) print_params = .true.

    if (self%parameter_number > 0) then
        write(*, *)     'Parametrized function 2D: ', self%name
        if (print_params) then""",
         """    character(len=EFT_names_max_length) :: param_name
    logical                             :: print_params_temp

    print_params_temp = .true.
    if (present(print_params)) print_params_temp = print_params

    if (self%parameter_number > 0) then
        write(*, *)     'Parametrized function 2D: ', self%name
        if (print_params_temp) then"""),
    ],
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("tree")
    ap.add_argument("--only", default="")
    args = ap.parse_args()
    tree = Path(args.tree)
    if "eftcamb_update/EFTCAMB" in str(tree.resolve()) and "scratchpad" not in str(tree.resolve()):
        sys.exit("refusing to patch the repository checkout")
    sel = [s for s in args.only.split(",") if s] or list(PATCHES)
    for name in sel:
        for rel, old, new in PATCHES[name]:
            p = tree / rel
            s = p.read_text()
            n = s.count(old)
            if n != 1:
                sys.exit(f"{name}: pattern found {n} times in {rel}")
            p.write_text(s.replace(old, new))
        print("applied", name)


if __name__ == "__main__":
    main()
