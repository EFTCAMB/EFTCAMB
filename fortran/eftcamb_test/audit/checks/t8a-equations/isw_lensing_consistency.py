"""Consistency of the EFT line-of-sight source quantities filled in derivs (A0 tree, python evolution fork):

    EFTLensing = k (Phi + Psi) = sigma' + eta k        (EFT branch)
    EFTISW     = k (Phi' + Psi') = sigma'' + (eta k)'
so EFTISW must equal d EFTLensing / d tau, and Psi, Phi must satisfy k Psi = sigma' + H sigma, k Phi = eta k - H sigma.
Also PsiDot + PhiDot must equal d(Psi + Phi)/dtau.  Checked at late times (a in [0.05, 1]) on a dense grid with a
4th-order finite difference; the residual is normalised by max|EFTISW| on the grid.

For Horava (Gamma4 != 0) the pi'' bug (see repro output-pidotdot-zero.py) makes EFTISW inconsistent; the other models
are controls for the remaining terms.

Usage: [AUDIT_INTTOL=8] python isw_lensing_consistency.py [MODEL ...]
"""

import json
import os
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness  # noqa: E402

TREE = harness.AUDIT_ROOT / "trees" / "A0"
MODELS = sys.argv[1:] or ["1_EFT_GR", "2_PEFT_Omega_const_1", "4_Des_F_R_LCDM_1", "5_hdsk_TG_1", "5_hdsk_TG_2",
                          "5_Horava_eta_3"]

WORKER = r"""
import sys, json
import numpy as np
sys.path.insert(0, sys.argv[1])
import camb
pars = camb.read_ini(sys.argv[2], no_validate=True)
pars.WantTensors = False
pars.InitPower.r = 0
import os
if os.environ.get('AUDIT_INTTOL'):
    pars.Accuracy.IntTolBoost = float(os.environ['AUDIT_INTTOL'])
res = camb.get_results(pars)
ks = [0.002, 0.02, 0.1]
a = np.logspace(np.log10(0.05), 0, 2000)
fields, arr = pars.EFTCAMB.get_scale_evolution(res, ks, a)
out = {}
for ik, k in enumerate(ks):
    c = arr[ik]
    tau = c["tau"]
    sl = slice(10, -10)
    def d(f):
        return np.gradient(f, tau, edge_order=2)
    n_isw = np.max(np.abs(c["EFTISW"][sl]))
    n_psi = np.max(np.abs(c["Psi"][sl])) + np.max(np.abs(c["Phi"][sl]))
    n_dot = np.max(np.abs(c["PsiDot"][sl] + c["PhiDot"][sl]))
    out[str(k)] = {
        "isw_vs_dlensing": float(np.max(np.abs(c["EFTISW"][sl] - d(c["EFTLensing"])[sl])) / n_isw),
        "isw_vs_dlensing_p99": float(np.percentile(np.abs(c["EFTISW"][sl] - d(c["EFTLensing"])[sl]), 99) / n_isw),
        "a_at_max": float(a[sl][np.argmax(np.abs(c["EFTISW"][sl] - d(c["EFTLensing"])[sl]))]),
        "lensing_vs_k(psi+phi)": float(np.max(np.abs(c["EFTLensing"][sl] - k * (c["Psi"][sl] + c["Phi"][sl]))) /
                                       np.max(np.abs(c["EFTLensing"][sl]))),
        "psidot+phidot_vs_fd": float(np.max(np.abs(c["PsiDot"][sl] + c["PhiDot"][sl] - d(c["Psi"] + c["Phi"])[sl]))
                                     / n_dot),
        "sigmadot_vs_fd": float(np.max(np.abs(c["sigmadot"][sl] - d(c["sigma"])[sl])) /
                                np.max(np.abs(c["sigmadot"][sl]))),
        "fd_resolution(pidot_vs_fd_pi)": float(np.max(np.abs(c["pidot"][sl] - d(c["pi"])[sl])) /
                                               max(np.max(np.abs(c["pidot"][sl])), 1e-300)),
    }
print(json.dumps(out))
"""

res = {}
for m in MODELS:
    env = dict(os.environ, OMP_NUM_THREADS="2")
    env.pop("PYTHONPATH", None)
    with harness.slot(2):
        p = subprocess.run(["nice", "-n", "10", sys.executable, "-c", WORKER, str(TREE),
                            str(harness.PARAMS_DIR / f"{m}.ini")], env=env, capture_output=True, text=True)
    if p.returncode:
        res[m] = f"rc={p.returncode}: {p.stderr[-300:]}"
        print(m, res[m])
        continue
    res[m] = json.loads(p.stdout.strip().splitlines()[-1])
    for k, v in res[m].items():
        print(f"{m:22s} k={k:>6s} " + "  ".join(f"{n}={x:.1e}" for n, x in v.items()))
out = harness.AUDIT_ROOT / "work" / "t8a-equations" / f"isw_lensing_consistency{os.environ.get('AUDIT_INTTOL', '')}.json"
out.write_text(json.dumps(res, indent=1))
