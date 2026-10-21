"""Independent check of steplog (04p018) fourth derivative and integral:
mpmath transcription of the Fortran formulas vs mpmath numerical differentiation / quadrature."""
import mpmath as mp

mp.mp.dps = 40
v1, v2, aT, D = mp.mpf(-1), mp.mpf("-0.8"), mp.mpf("0.3"), mp.mpf("0.5")
lg = mp.log


def val(x):
    L = (lg(x) - lg(aT)) / D
    return 0.5 * (v2 - v1) * L / mp.sqrt(1 + L**2) + 0.5 * (v1 + v2)


def d4_fortran(x):  # verbatim transcription of lines 340-349 (real arithmetic; negative base -> nan)
    l = lg(x / aT)
    b = 1 - l / D**2
    if b < 0:
        return mp.nan
    x4 = x**4
    r = (-(3 * l) / (D**2 * x4 * b**1.5) - 6 / (x4 * mp.sqrt(b)) + 11 / (D**2 * x4 * b**1.5)
         + (105 * l) / (16 * D**8 * x4 * b**4.5) - (45 * l) / (4 * D**6 * x4 * b**3.5) + 15 / (2 * D**6 * x4 * b**3.5)
         + (33 * l) / (4 * D**4 * x4 * b**2.5) - 27 / (2 * D**4 * x4 * b**2.5))
    return 0.5 * (v2 - v1) * r / D


def integ_fortran(x, prec_double=True):
    lx, la = lg(x), lg(aT)
    s = mp.sqrt(D**2 + (lx - la) ** 2)
    if lx == 0:
        expo = mp.inf  # log(aT)**2/log(1) -> +Inf in IEEE; sign of prefactor
        # exponent = -1-1.5(v1+v2) - 1.5(v2-v1)*(0-2la+la^2/0)/s ; la^2/0 = +Inf, times -1.5(v2-v1)<0 -> -Inf
        t1 = mp.mpf(1)  # 1**(-Inf) = 1 in IEEE
    else:
        t1 = x ** (-1 - 1.5 * (v1 + v2) - 1.5 * (v2 - v1) * (lx - 2 * la + la**2 / lx) / s)
    t2 = mp.exp(1.5 * (v1 - v2) * (D**2 / s - mp.sqrt(D**2 + la**2)))
    return t1 * t2


def integ_exact(x):  # a^2 rho_DE(a)/rho_DE0 = a^2 exp(-3 int_1^a (1+w)/a' da')
    return x**2 * mp.exp(-3 * mp.quad(lambda t: (1 + val(t)) / t, [1, x]))


print(f"{'a':>9} | {'d4 fortran':>14} {'d4 numeric':>14} | {'I fortran':>14} {'I quad':>14}")
for a in ["0.05", "0.2", "0.3", "0.385", "0.39", "0.5", "0.9", "0.999999", "1"]:
    x = mp.mpf(a)
    d4n = mp.diff(val, x, 4)
    print(f"{a:>9} | {mp.nstr(d4_fortran(x), 8):>14} {mp.nstr(d4n, 8):>14} | "
          f"{mp.nstr(integ_fortran(x), 10):>14} {mp.nstr(integ_exact(x), 10):>14}")
print("NaN threshold a > aT*exp(Delta^2) =", mp.nstr(aT * mp.exp(D**2), 6))
