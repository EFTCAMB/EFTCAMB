"""Check root-cause claim for steplog d4: (a) code formula with base (1 - L/D^2) (as in A0) -> reproduces A0 numbers;
(b) same formula with base 1 + (L/D)^2 -> is it then correct?; (c) proposed fix formula vs mpmath.diff."""
import mpmath as mp
mp.mp.dps = 40


def code_d4(x, v1, v2, at, D, base):
    L = mp.log(x / at)
    B = base(L, D)
    return mp.mpf(1) / 2 * (v2 - v1) * (-(3 * L) / (D**2 * x**4 * B**1.5) - 6 / (x**4 * mp.sqrt(B))
                                        + 11 / (D**2 * x**4 * B**1.5) + (105 * L) / (16 * D**8 * x**4 * B**4.5)
                                        - (45 * L) / (4 * D**6 * x**4 * B**3.5) + 15 / (2 * D**6 * x**4 * B**3.5)
                                        + (33 * L) / (4 * D**4 * x**4 * B**2.5) - 27 / (2 * D**4 * x**4 * B**2.5)) / D


def fix_d4(x, v1, v2, at, D):
    t = (mp.log(x) - mp.log(at)) / D
    return mp.mpf(1) / 2 * (v2 - v1) / x**4 * (15 * t * (3 - 4 * t**2) / (1 + t**2)**4.5 / D**4
                                              - 6 * (12 * t**2 - 3) / (1 + t**2)**3.5 / D**3
                                              - 33 * t / (1 + t**2)**2.5 / D**2 - 6 / (1 + t**2)**1.5 / D)


for (v1, v2, at, D) in [(-1.0, -0.8, 0.3, 0.5), (-0.9, -1.1, 0.7, 1.0), (0.0, 1.5, 0.5, 0.5)]:
    v1, v2, at, D = map(lambda s: mp.mpf(str(s)), (v1, v2, at, D))

    def w(a):
        u = (mp.log(a) - mp.log(at)) / D
        return (v2 - v1) / 2 * u / mp.sqrt(1 + u * u) + (v1 + v2) / 2
    print(f"params {float(v1)}, {float(v2)}, {float(at)}, {float(D)}")
    for x in ["0.05", "0.2", "0.3", "0.5", "0.9", "1.0"]:
        x = mp.mpf(x)
        ex = mp.diff(w, x, 4)
        a = code_d4(x, v1, v2, at, D, lambda L, D: 1 - L / D**2)
        b = code_d4(x, v1, v2, at, D, lambda L, D: 1 + (L / D)**2)
        c = fix_d4(x, v1, v2, at, D)
        fa = mp.nstr(a, 8) if mp.im(a) == 0 else "complex(NaN)"
        print(f"  a={float(x):4.2f} exact {mp.nstr(ex, 10):>14} | code-as-is {fa:>14} | code w/ base 1+(L/D)^2 "
              f"{mp.nstr(b, 8):>14} | proposed fix {mp.nstr(c, 10):>14} relerr {float(abs(c - ex) / abs(ex)):.1e}")
