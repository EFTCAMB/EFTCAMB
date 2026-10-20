"""Verifier (independent, numeric): compare the Hdotdotdot field from get_eft_functions with
  (i)  d(Hdotdot)/dtau from a cubic spline of the Hdotdot field against the *tau field* on a dense a grid
       (different method from the repro: no a-derivative, no Richardson),
  (ii) the sympy derivative of the shipped Hdotdot (lambdified from verify_hdddot_sympy.py's construction),
       fed with the cache fields + spline values of EFTLambdadotdot and EFTOmegaPPPP,
  (iii) the shipped Hdotdotdot expression lambdified with Ldd = O4 = 0 (what the code sees on the python path),
       to confirm the field is exactly what that expression gives.
Usage: python verify_hdddot_numeric.py [--tree PATH] [--models m1,m2]
"""

import argparse
import os
import sys
import tempfile
from pathlib import Path

import numpy as np
from scipy.interpolate import CubicSpline

REPO = Path(__file__).resolve().parents[5]
ROOT = "/private/tmp/claude-501/-Users-marco-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit"
ap = argparse.ArgumentParser()
ap.add_argument("--tree", default=f"{ROOT}/trees/A0")
ap.add_argument("--models", default="5_Kmouflage_1,5_Kmouflage_2,5_quint_exponential,5_ADE_1,4_Des_F_R_LCDM_1")
args = ap.parse_args()
os.environ.setdefault("OMP_NUM_THREADS", "2")
sys.path.insert(0, args.tree)
sys.path.insert(0, str(Path(__file__).parent))
import camb  # noqa: E402
import sympy as sp  # noqa: E402

# build symbolic truth / shipped from the sympy module (import executes its prints; silence them)
import contextlib  # noqa: E402
import io  # noqa: E402

with contextlib.redirect_stdout(io.StringIO()):
    import verify_hdddot_sympy as S  # noqa: E402
syms = [S.sa, S.sH, S.sHd, S.sHdd, S.sP, S.sPd, S.sPdd, S.sL, S.sLd, S.sLdd] + list(S.so[:5])
true_f = sp.lambdify(syms, S.true3, "numpy")
ship_f = sp.lambdify(syms, S.ship3, "numpy")
for s_ in ("Hddd", "Pddd", "Lddd"):
    assert sp.Symbol(s_) not in S.true3.free_symbols

print("camb from", camb.__file__)
tmp = tempfile.mkdtemp(prefix="vh3_", dir=f"{ROOT}/work/verify/hdddot")
agrid = np.logspace(np.log10(0.01), 0, 3000)
a0 = np.array([0.05, 0.2, 0.5, 0.9])
for m in args.models.split(","):
    ini = os.path.join(tmp, f"{m}.ini")
    src = Path(m) if m.endswith(".ini") else REPO / "fortran/eftcamb_test/parameters" / f"{m}.ini"
    m = Path(m).stem
    Path(ini).write_text(f"DEFAULT({src})\noutput_root = {tmp}/{m}\nEFTCAMB_write_background = F\nfeedback_level = 0\n")
    cp = camb.read_ini(ini, no_validate=True)
    res = camb.get_background(cp)
    _, G = cp.EFTCAMB.get_eft_functions(res, agrid)
    _, C = cp.EFTCAMB.get_eft_functions(res, a0)
    tau = G["tau"]
    assert np.all(np.diff(tau) > 0), "tau field not monotonic"
    num = CubicSpline(tau, G["Hdotdot"])(C["tau"], 1)
    Ldd_num = CubicSpline(tau, G["EFTLambdadot"])(C["tau"], 1) - 2 * C["adotoa"] * C["EFTLambdadot"]
    O4_num = CubicSpline(agrid, G["EFTOmegaPPP"])(a0, 1)
    base = [a0, C["adotoa"], C["Hdot"], C["Hdotdot"], C["gpresm_t"], C["gpresdotm_t"], C["gpresdotdotm_t"],
            C["EFTLambda"], C["EFTLambdadot"]]
    om = [C["EFTOmegaV"], C["EFTOmegaP"], C["EFTOmegaPP"], C["EFTOmegaPPP"]]
    truth = true_f(*base, Ldd_num, *om, O4_num)
    ship0 = ship_f(*base, C["EFTLambdadotdot"], *om, C["EFTOmegaPPPP"])
    print(f"\n== {m}  model={cp.EFTCAMB.model_name()}  cache EFTLambdadotdot={C['EFTLambdadotdot'].tolist()} "
          f"EFTOmegaPPPP={C['EFTOmegaPPPP'].tolist()}  max|Omega|={np.max(np.abs(G['EFTOmegaV'])):.3e}")
    print(f"{'a':>5s} {'field':>12s} {'dHdd/dtau spl':>14s} {'sympy-true':>12s} {'shipped(sym)':>13s}")
    for i in range(len(a0)):
        print(f"{a0[i]:5.2f} {C['Hdotdotdot'][i]:12.4e} {num[i]:14.4e} {truth[i]:12.4e} {ship0[i]:13.4e}")
    sc = np.max(np.abs(num))
    print(f"  max|field-spline|/max|spline|     = {np.max(np.abs(C['Hdotdotdot'] - num)) / sc:.3e}")
    print(f"  max|sympytrue-spline|/max|spline| = {np.max(np.abs(truth - num)) / sc:.3e}")
    print(f"  max|shipped(sym)-field|/max|field| = "
          f"{np.max(np.abs(ship0 - C['Hdotdotdot'])) / max(np.max(np.abs(C['Hdotdotdot'])), 1e-300):.3e}")
    # control: Hdotdot vs spline derivative of Hdot
    numc = CubicSpline(tau, G["Hdot"])(C["tau"], 1)
    print(f"  control max|Hdotdot-dHdot/dtau|/max = {np.max(np.abs(C['Hdotdot'] - numc)) / np.max(np.abs(numc)):.3e}")
