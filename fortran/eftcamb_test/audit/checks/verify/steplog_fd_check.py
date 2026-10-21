"""Verifier: reimplement 04p018 steplog d3/d4/integral formulas in mpmath and compare with high-precision
finite differences of the value function / direct quadrature of the integral definition."""
import mpmath as mp

mp.mp.dps = 50
v1, v2, aT, D = mp.mpf(-1), mp.mpf("-0.8"), mp.mpf("0.3"), mp.mpf("0.5")
lg = mp.log


def val(x):
    L = (lg(x) - lg(aT)) / D
    return 0.5 * (v2 - v1) * L / mp.sqrt(1 + L**2) + 0.5 * (v1 + v2)


def d3_code(x):
    t1 = (9 * D**2 + 2 * lg(x)**3 + lg(x)**2 * (9 - 6 * lg(aT)) - 4 * lg(aT) * (3 + D**2) + 9 * lg(aT)**2
          - 2 * lg(aT)**3 + 2 * lg(x) * (6 + 2 * D**2 + 3 * lg(aT) * (-3 + lg(aT))))
    return (0.5 * (v2 - v1) * D * (D**2 * (-3 + 2 * D**2) + (lg(x) - lg(aT)) * t1) / x**3
            / mp.sqrt(1 + ((lg(x) - lg(aT)) / D)**2) / ((lg(x) - lg(aT))**2 + D**2)**3)


def d4_code(x):
    l = lg(x / aT)
    b = 1 - l / D**2
    try:
        r = (-(3 * l) / (D**2 * x**4 * b**1.5) - 6 / (x**4 * mp.sqrt(b)) + 11 / (D**2 * x**4 * b**1.5)
             + 105 * l / (16 * D**8 * x**4 * b**4.5) - 45 * l / (4 * D**6 * x**4 * b**3.5)
             + 15 / (2 * D**6 * x**4 * b**3.5) + 33 * l / (4 * D**4 * x**4 * b**2.5)
             - 27 / (2 * D**4 * x**4 * b**2.5))
        r = 0.5 * (v2 - v1) * r / D
        return r if mp.im(r) == 0 else mp.nan
    except Exception:
        return mp.nan


def integ_exact(x):
    return x**2 * mp.exp(3 * mp.quad(lambda n: 1 + val(mp.exp(n)), [lg(x), 0]))


def integ_code_float(x):
    import math
    v1f, v2f, aTf, Df = -1.0, -0.8, 0.3, 0.5
    lx = math.log(x)
    with __import__("numpy").errstate(all="ignore"):
        import numpy as np
        lxn = np.float64(lx)
        ratio = np.float64(math.log(aTf)**2) / lxn
    e = -1 - 1.5 * (v1f + v2f) - 1.5 * (v2f - v1f) * (lx - 2 * math.log(aTf) + ratio) / math.sqrt(Df**2 + (lx - math.log(aTf))**2)
    t1 = np.float64(x) ** np.float64(e)
    t2 = math.exp(1.5 * (v1f - v2f) * (Df**2 / math.sqrt(Df**2 + (lx - math.log(aTf))**2) - math.sqrt(Df**2 + math.log(aTf)**2)))
    return float(t1 * t2)


print(f"{'a':>8} {'d3 code':>14} {'d3 FD':>14} {'d4 code':>14} {'d4 FD':>14} {'I code':>12} {'I exact':>12}")
for x in ["0.05", "0.2", "0.3", "0.5", "0.9", "1.0"]:
    x = mp.mpf(x)
    fd3 = mp.diff(val, x, 3)
    fd4 = mp.diff(val, x, 4)
    print(f"{float(x):8.3f} {float(d3_code(x)):14.6g} {float(fd3):14.6g} {float(d4_code(x)) if d4_code(x)==d4_code(x) else float('nan'):14.6g}"
          f" {float(fd4):14.6g} {integ_code_float(float(x)):12.8f} {float(integ_exact(x)):12.8f}")
