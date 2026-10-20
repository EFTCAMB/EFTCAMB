"""Independent reference for the quintic spline: scipy make_interp_spline(k=5) with the same end conditions that
spline5pars uses (1st and 2nd derivative of the quintic interpolating the 6 end points). Compared with the
parametrization output of the A0 and poly5-fix builds (via verify_spline5_param_eval.evaluate).
usage: python verify_spline5_scipy_ref.py
"""

import numpy as np
from scipy.interpolate import make_interp_spline

from verify_spline5_param_eval import evaluate


def end_derivs(x, y, left=True):
    xe, ye = (x[:6], y[:6]) if left else (x[-6:], y[-6:])
    x0 = xe[0] if left else xe[-1]
    c = np.linalg.solve(np.vander(xe - x0, 6, increasing=True), ye)
    return c[1], 2 * c[2]


def check(name, nodes, vals, logx, xs):
    u = np.log(nodes) if logx else np.asarray(nodes, float)
    y = np.asarray(vals, float)
    d1l, d2l = end_derivs(u, y, True)
    d1r, d2r = end_derivs(u, y, False)
    s = make_interp_spline(u, y, k=5, bc_type=([(1, d1l), (2, d2l)], [(1, d1r), (2, d2r)]))
    uq = np.log(xs) if logx else np.asarray(xs)
    ref_v = s(uq)
    ref_d = s(uq, 1) / (np.asarray(xs) if logx else 1.0)
    scale = np.max(np.abs(y))
    print(f"{name}: scipy quintic spline (same end conditions) vs EFTCAMB spline5 parametrization")
    for t in ("A0", "FIX"):
        r = np.array(evaluate(t, nodes, vals, logx, xs))
        inside = (uq >= u[0]) & (uq <= u[-1])
        ev = np.max(np.abs(r[inside, 1] - ref_v[inside])) / scale
        ed = np.max(np.abs(r[inside, 2] - ref_d[inside])) / np.max(np.abs(ref_d[inside]))
        print(f"   {t}: max|value - ref|/max|y| = {ev:.3e}   max|d1 - ref|/max|d1| = {ed:.3e}")


def main():
    n15 = [10.0 ** (-3 + 3 * i / 14) for i in range(15)]
    xs = [10 ** (-3 + 3 * k / 280) for k in range(1, 280)]
    check("log15_a2 (Logx=T)", n15, [a**2 for a in n15], True, xs)
    n9 = [1e-3, 0.01, 0.05, 0.1, 0.25, 0.4, 0.6, 0.8, 1.0]
    v9 = [0.2 * np.sin(4 * a) + 0.1 for a in n9]
    xs9 = list(np.exp(np.linspace(np.log(1e-3), 0, 400))[1:-1])
    check("sin9 (Logx=T)", n9, v9, True, xs9)
    check("sin9 (Logx=F)", n9, v9, False, list(np.linspace(1e-3, 1, 400)[1:-1]))


if __name__ == "__main__":
    main()
