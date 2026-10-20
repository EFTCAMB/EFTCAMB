"""Side check: 5_scg analytic EFTLambdadotdot (008p6) vs numerical a^2 d2(EFTLambda/a^2)/dtau^2 and vs
d(EFTLambdadot)/dtau - 2H EFTLambdadot, from get_eft_functions on A0.  Usage: python minktach_scg_ldd.py"""
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness  # noqa: E402

os.environ["OMP_NUM_THREADS"] = "1"
ROOT = "/private/tmp/claude-501/-Users-marco-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit"
sys.path.insert(0, f"{ROOT}/trees/A0")
PAR = "/Users/marco/IW/eftcamb_update/EFTCAMB/fortran/eftcamb_test/parameters"
with harness.slot(1):
    import numpy as np
    import camb
    cp = camb.read_ini(f"{PAR}/5_scg.ini", no_validate=True)
    cp.EFTCAMB.EFTCAMB_feedback_level = 0
    cp.EFTCAMB.EFTCAMB_write_background = False
    res = camb.get_background(cp)
    a = np.logspace(-3, 0, 6001)
    c = res.Params.EFTCAMB.get_eft_functions(res, a)[1]
    g = lambda f: np.array(c[f], dtype=float)
    H, Hd, L, Ld, Ldd = g("adotoa"), g("Hdot"), g("EFTLambda"), g("EFTLambdadot"), g("EFTLambdadotdot")
    N = np.log(a)
    lam = L / a**2
    LN = np.gradient(lam, N, edge_order=2)
    Ldd_v = a**2 * (Hd * LN + H**2 * np.gradient(LN, N, edge_order=2))
    Ldd_y = H * np.gradient(Ld, N, edge_order=2) - 2 * H * Ld
    print("convention max|a^2 H Lambda_N - EFTLambdadot|/max|EFTLambdadot| =",
          float(np.abs(a**2 * H * LN - Ld)[5:-5].max() / np.abs(Ld[5:-5]).max()))
    for q in (0.01, 0.1, 0.3, 0.5, 0.8, 0.99):
        i = int(np.argmin(abs(a - q)))
        print(f"a={a[i]:.3f}: EFTLambdadotdot code={Ldd[i]: .4e}  num(values)={Ldd_v[i]: .4e}  "
              f"num(EFTLambdadot)={Ldd_y[i]: .4e}  H*EFTLambdadot={H[i]*Ld[i]: .4e}")
