"""Track t4b-derivs: offline check of a corrected full-mapping Hdotdotdot formula (06p1 compute_H_derivs).

Hdot D = N with D = 1+O+aO'/2, N = -B H^2/2 - (p+L)/2, B = 1+O+2aO'+a^2O''  =>
Hdotdotdot = (N'' - 2 Hdotdot D' - Hdot D'')/D  (' = d/dtau).  Inputs not filled by get_eft_functions for
full-mapping models (EFTLambdadotdot, EFTOmegaPPPP) are replaced by finite differences of the stored lower
derivatives.  Compares (i) the stored Hdotdotdot, (ii) this formula, (iii) a finite difference of stored Hdotdot.

    python hdddot_formula_check.py NPZ [NPZ ...]     (deriv_worker.py output)
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from deriv_analyze import dlna, load  # noqa: E402

for f in sys.argv[1:]:
    d = load(Path(f))
    hs, a = d["hs"], np.exp(d["xc"])
    g = lambda k: d[f"st_{k}"][1, 4]  # noqa: E731  (centre, h index 1)
    H, Hd, Hdd = g("adotoa"), g("Hdot"), g("Hdotdot")
    O, O1, O2, O3 = g("EFTOmegaV"), g("EFTOmegaP"), g("EFTOmegaPP"), g("EFTOmegaPPP")
    O4 = dlna(d["st_EFTOmegaPPP"], hs)[1] / a
    L, Ld = g("EFTLambda"), g("EFTLambdadot")
    Ldd = H * (dlna(d["st_EFTLambdadot"], hs)[1] - 2 * Ld)  # a^2 d(Ldot/a^2)/dtau
    p, pd, pdd = g("gpresm_t"), g("gpresdotm_t"), g("gpresdotdotm_t")
    D = 1 + O + 0.5 * a * O1
    D1 = a * H * (1.5 * O1 + 0.5 * a * O2)
    D2 = a * (H**2 + Hd) * (1.5 * O1 + 0.5 * a * O2) + a**2 * H**2 * (2 * O2 + 0.5 * a * O3)
    B = 1 + O + 2 * a * O1 + a**2 * O2
    B1 = a * H * (3 * O1 + 4 * a * O2 + a**2 * O3)
    B2 = a * (H**2 + Hd) * (3 * O1 + 4 * a * O2 + a**2 * O3) + a**2 * H**2 * (7 * O2 + 6 * a * O3 + a**2 * O4)
    X2 = lambda x, xd, xdd: xdd + 4 * H * xd + (2 * Hd + 4 * H**2) * x  # noqa: E731
    P2 = X2(p, pd, pdd) + X2(L, Ld, Ldd)
    N2 = -0.5 * B2 * H**2 - 2 * B1 * H * Hd - B * (Hd**2 + H * Hdd) - 0.5 * P2
    Hddd = (N2 - 2 * Hdd * D1 - Hd * D2) / D
    fd = H * dlna(d["st_Hdotdot"], hs)[1]
    sc = np.maximum(np.abs(fd), H**4)
    m = a > 1e-3
    e_new = np.max(np.abs(Hddd - fd)[m] / sc[m])
    e_old = np.max(np.abs(g("Hdotdotdot") - fd)[m] / sc[m])
    print(f"{Path(f).stem:24s} a>1e-3: max|stored-FD|/max(|FD|,H^4) = {e_old:.2e}   max|corrected-FD|/... = {e_new:.2e}")
