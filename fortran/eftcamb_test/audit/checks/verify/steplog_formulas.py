"""Independent check of steplog d4 / integral formulas (proposed fix vs sympy/quad)."""
import math
import sympy as sp
from scipy.integrate import quad

v1, v2, aT, de = -1.0, -0.8, 0.3, 0.5
a = sp.symbols("a", positive=True)
L = (sp.log(a) - math.log(aT)) / de
w = (v2 - v1) / 2 * L / sp.sqrt(1 + L**2) + (v1 + v2) / 2
d4 = sp.lambdify(a, sp.diff(w, a, 4))
wf = sp.lambdify(a, w)


def d4_fix(x):
    L = math.log(x / aT) / de
    q = 1 + L * L
    return 0.5 * (v2 - v1) / x**4 * (15 * L * (3 - 4 * L * L) / q**4.5 / de**4 - 6 * (12 * L * L - 3) / q**3.5 / de**3
                                     - 33 * L / q**2.5 / de**2 - 6 / q**1.5 / de)


def I_fix(x):
    t1 = x ** (-1 - 1.5 * (v1 + v2))
    t2 = math.exp(-1.5 * (v2 - v1) * (math.sqrt(de**2 + (math.log(x) - math.log(aT)) ** 2) - math.sqrt(de**2 + math.log(aT) ** 2)))
    return t1 * t2


def I_code(x):
    lx, la = math.log(x), math.log(aT)
    try:
        e = -1 - 1.5 * (v1 + v2) - 1.5 * (v2 - v1) * (lx - 2 * la + la**2 / lx) / math.sqrt(de**2 + (lx - la) ** 2)
    except ZeroDivisionError:
        e = -math.inf
    t1 = x**e
    t2 = math.exp(1.5 * (v1 - v2) * (de**2 / math.sqrt(de**2 + (lx - la) ** 2) - math.sqrt(de**2 + la**2)))
    return t1 * t2


def I_ex(x):
    return x * x * math.exp(3 * quad(lambda n: 1 + wf(math.exp(n)), math.log(x), 0, epsabs=1e-13, epsrel=1e-13)[0])


for x in [0.01, 0.05, 0.2, 0.3, 0.5, 0.9, 0.999999, 1.0, 2.0]:
    print(f"{x:9.6f} d4 sympy {d4(x):14.6g} fix {d4_fix(x):14.6g} | I quad {I_ex(x):.12f} fix {I_fix(x):.12f} code(py) {I_code(x):.12f}")
