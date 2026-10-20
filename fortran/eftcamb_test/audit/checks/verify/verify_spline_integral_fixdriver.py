"""Verifier (spline-integral-convention): same Fortran driver as verify_spline_integral_formula.py, linked against the
integral-only patched tree (patch_spline_integral_only.py), to check the patched integral() against the exact
EFTCAMB convention I(a) = a^2 exp(-3 int_1^a (1+w)/a' da') before using it in model runs.
usage: /Users/marco/Software/miniforge3/bin/python verify_spline_integral_fixdriver.py
"""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import verify_spline_integral_formula as F  # noqa: E402

F.BUILD = F.AUDIT / "work" / "verify" / "spline_integral" / "tree_intfix_build"
F.WORK = F.AUDIT / "work" / "verify" / "spline_integral" / "fixdriver"


def main():
    F.build()
    a = np.array([1e-8, 3e-6, 1e-4, 0.003, 0.05, 0.15, 0.35, 0.55, 0.85, 0.95, 1.0])
    w0, wa = -0.9, 0.2
    cases = {
        "cubic const w, nodes 1e-8,0.2..1": (9, [1e-8, 0.2, 0.4, 0.6, 0.8, 1.0], lambda t: -0.9 + 0 * t, ""),
        "cubic linear w, nodes 0.05..1 (null=w(x1))": (9, list(np.linspace(0.05, 1, 11)), None, ""),
        "quintic Logx const w, nodes 1e-8..1": (10, [10.0 ** e for e in range(-8, 1)], lambda t: -0.9 + 0 * t,
                                                 "EFTw_Logx_Flag = T\n"),
    }
    for name, (flag, xs, wf, extra) in cases.items():
        xs = np.array(xs)
        if wf is None:
            ws = w0 + wa * (1 - xs)
            # w(a) as the code defines it: spline inside, null_value = w(x1) outside
            wfun = lambda t, xs=xs: np.where(t < xs[0], w0 + wa * (1 - xs[0]), w0 + wa * (1 - t))  # noqa: E731
        else:
            ws = wf(xs)
            wfun = wf
        v, i_f = F.run(flag, xs, ws, extra, a=a)
        ex = F.exact_I(lambda t: float(wfun(np.array(t))), a)
        rel = np.abs(i_f - ex) / ex
        print(f"{name}: max rel |I_patched - I_exact| = {rel.max():.3e}  (I(1)={i_f[-1]:.15g}, "
              f"I(1e-8) patched {i_f[0]:.6e} exact {ex[0]:.6e})")


if __name__ == "__main__":
    main()
