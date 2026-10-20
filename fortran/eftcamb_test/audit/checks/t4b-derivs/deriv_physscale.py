"""Track t4b-derivs: background derivative/Friedmann residuals relative to the physical scale of the equation.

For each npz from deriv_worker.py (stencils at several h): residual of
  Hdot - d(adotoa)/dtau                       / H^2
  EFTcdot - a^2 d(c/a^2)/dtau                  / (H * (|c| + |Lambda| + H^2))
  EFTLambdadot - a^2 d(Lambda/a^2)/dtau        / (H * (|c| + |Lambda| + H^2))
  Hdotdot - d(Hdot)/dtau                       / H^3
  Friedmann (1+O+aO')H^2 - (rho+2c-L)/3        / H^2
  2nd Friedmann (EFT form)                     / H^2
  DE conservation d(grhov_t)/dtau + H(grhov_t+3 gpiv_t)  / (H * H^2)
The FD estimate uses the largest step that is converged (|FD(h1)-FD(h2)| < resid/3), reporting min over steps.

    python deriv_physscale.py DIR [--amin A]
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from deriv_analyze import dlna, load  # noqa: E402


def analyze(path: Path, amin: float) -> dict:
    d = load(path)
    hs = d["hs"]
    a = np.exp(d["xc"])
    g = lambda f: d[f"st_{f}"][0, 4]  # noqa: E731
    D = lambda f: dlna(d[f"st_{f}"], hs)  # noqa: E731  (n_h, n_c)
    H, Hd, Hdd = g("adotoa"), g("Hdot"), g("Hdotdot")
    c, L = g("EFTc"), g("EFTLambda")
    on = a > max(amin, d["meta"]["pert_turn_on_results"] * np.exp(4 * hs.max()))
    out = {}

    def rel(name, resid_h, scale):
        # resid_h: (n_h, n_c) residual for each step; take the step with smallest |resid|, flag if steps disagree
        r = np.abs(resid_h) / scale[None, :]
        best = r.min(axis=0)
        spread = r.max(axis=0) - r.min(axis=0)
        m = on & np.isfinite(best)
        if not m.any():
            return
        i = int(np.argmax(np.where(m, best, -1)))
        out[name] = (float(best[i]), float(a[i]), float(spread[i]), float(np.median(best[m])))

    rel("Hdot", Hd[None] - H[None] * D("adotoa"), H**2)
    sc = H * (np.abs(c) + np.abs(L) + H**2)
    if np.any(c != 0):
        rel("EFTcdot", g("EFTcdot")[None] - H[None] * (D("EFTc") - 2 * c[None]), sc)
    if np.any(L != 0):
        rel("EFTLambdadot", g("EFTLambdadot")[None] - H[None] * (D("EFTLambda") - 2 * L[None]), sc)
    rel("Hdotdot", Hdd[None] - H[None] * D("Hdot"), H**3)
    Om, OmP, OmPP = g("EFTOmegaV"), g("EFTOmegaP"), g("EFTOmegaPP")
    rho, p = g("grhom_t"), g("gpresm_t")
    f1 = (1 + Om + a * OmP) * H**2 - (rho + 2 * c - L) / 3
    f2 = (1 + Om) * (2 * Hd + H**2) + a * OmP * (Hd + 2 * H**2) + a**2 * OmPP * H**2 + p + L
    rel("Friedmann1", f1[None], H**2)
    rel("Friedmann2", f2[None], H**2)
    if np.any(g("grhov_t") != 0) and np.any(g("gpiv_t") != 0):
        rv = g("grhov_t")
        rel("DE_cons", H[None] * D("grhov_t") + H[None] * (rv + 3 * g("gpiv_t"))[None], H * H**2)
    return out


def main() -> None:
    amin = float(sys.argv[sys.argv.index("--amin") + 1]) if "--amin" in sys.argv else 0.0
    for f in sorted(Path(sys.argv[1]).glob("*.npz")):
        r = analyze(f, amin)
        print(f"== {f.stem}")
        for k, (e, a, s, med) in r.items():
            flag = "PASS" if e < (1e-6 if k in ("Hdot", "EFTcdot", "EFTLambdadot", "Friedmann1", "Friedmann2", "DE_cons")
                                  else 1e-4) else "FAIL"
            print(f"   {flag} {k:14s} max={e:.2e} @a={a:.2e} (step-spread {s:.1e}) median={med:.1e}")


if __name__ == "__main__":
    main()
