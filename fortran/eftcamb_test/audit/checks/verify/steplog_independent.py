"""Independent check of steplog d4 / integral: numpy replica of the Fortran expressions vs mpmath
high-precision numerical differentiation / quadrature of the documented w(a)."""
import numpy as np
import mpmath as mp

mp.mp.dps = 40
v1, v2, aT, D = -1.0, -0.8, 0.3, 0.5


def w(a):
    L = (mp.log(a) - mp.log(aT)) / D
    return 0.5 * (v2 - v1) * L / mp.sqrt(1 + L**2) + 0.5 * (v1 + v2)


def d4_code(x):
    with np.errstate(all="ignore"):
        l = np.log(x / aT); b = 1.0 - l / D**2
        return 0.5 * (v2 - v1) * (-(3 * l) / (D**2 * x**4 * b**1.5) - 6 / (x**4 * np.sqrt(b)) + 11 / (D**2 * x**4 * b**1.5)
                                 + 105 * l / (16 * D**8 * x**4 * b**4.5) - 45 * l / (4 * D**6 * x**4 * b**3.5)
                                 + 15 / (2 * D**6 * x**4 * b**3.5) + 33 * l / (4 * D**4 * x**4 * b**2.5)
                                 - 27 / (2 * D**4 * x**4 * b**2.5)) / D


def integral_code(x):
    with np.errstate(all="ignore"):
        lx, la = np.log(x), np.log(aT)
        t1 = x ** (-1 - 1.5 * (v1 + v2) - 1.5 * (v2 - v1) * (lx - 2 * la + la**2 / lx) / np.sqrt(D**2 + (lx - la) ** 2))
        t2 = np.exp(1.5 * (v1 - v2) * (D**2 / np.sqrt(D**2 + (lx - la) ** 2) - np.sqrt(D**2 + la**2)))
        return t1 * t2


def integral_exact(x):
    return x**2 * mp.exp(3 * mp.quad(lambda n: 1 + w(mp.exp(n)), [mp.log(x), 0]))


for x in [0.05, 0.2, 0.3, 0.385, 0.4, 0.5, 0.9, 1.0]:
    e4 = mp.diff(w, mp.mpf(x), 4)
    print(f"a={x:6.3f} d4 code={d4_code(np.float64(x)):14.6g} mpdiff={float(e4):14.6g} | "
          f"I code={integral_code(np.float64(x)):.12f} exact={float(integral_exact(mp.mpf(x))):.12f}")
for x in [1 - 1e-6, 1 - 1e-12, 1 - 1e-15, 1.0]:
    print(f"x=1-{1-x:.0e}: I code={integral_code(np.float64(x)):.12f}")
