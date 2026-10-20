"""Repro: EFTCAMB_evolve_metric_h = T (evolve h' with the ii Einstein equation) disagrees with the two constraint-based
formulations of the same Horndeski model by ~1% in low-ell TT and by O(10) in P(k) at k ~ 100/Mpc.

equations.f90, derivs, delta_phi branch ("EFTCAMB MOD START: equation of motion. Now we use the non-RSA expression for
dz."): with EFTCAMB_evolve_metric_h = T, z = h'/(2k) is taken from the integrated variable y(w_ix+2) whose derivative is
metric_h_ddot (ii equation); with F (the default, used by all test inis) h' comes from the 00 constraint. The same model
can also be run with the pi-field equations (EFTCAMB_evolve_delta_phi = F, 5_hdsk_TG_2.ini).
Expected: all three agree to the accuracy of the code (TG_1 vs TG_2 agree to ~2e-4).
Observed: printed (the evolve_metric_h run is the outlier).  B1 (pre-rebase EFTCAMB) behaves the same, so this is not
a rebase regression. The evolved h' drifts off the 00 constraint; sigma = (h' + 6 eta')/(2k) is a near-cancellation on
large scales, so the drift shows up in the ISW source, and at high k the evolution is not stable (CAMB itself takes z
from the constraint "for numerical stability").

Run: /Users/marco/Software/miniforge3/bin/python evolve-metric-h-inconsistent.py
"""

import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "checks"))
import harness  # noqa: E402

OUT = harness.AUDIT_ROOT / "work" / "t8a-equations" / "repro_metric_h"
OUT.mkdir(parents=True, exist_ok=True)
P = harness.PARAMS_DIR
RUNS = {
    "A0 TG_1 delta_phi, h' from 00 constraint (default)": ("A0", P / "5_hdsk_TG_1.ini", {}),
    "A0 TG_1 delta_phi, h' evolved (evolve_metric_h=T)": ("A0", P / "5_hdsk_TG_1.ini", {"EFTCAMB_evolve_metric_h": "T"}),
    "A0 TG_2 pi-field equations": ("A0", P / "5_hdsk_TG_2.ini", {}),
    "B1 TG_1 delta_phi, h' evolved (evolve_metric_h=T)": ("B1", P / "5_hdsk_TG_1.ini", {"EFTCAMB_evolve_metric_h": "T"}),
}
OV = {}  # default test settings (redshift windows on -> P(k) computed to k ~ 100/Mpc)


def run(item):
    i, (label, (tree, ini, ov)) = item
    return label, harness.run_model(tree, ini, OUT / f"run{i}.npz", overrides={**OV, **ov}, threads=1)


with ThreadPoolExecutor(4) as ex:
    st = dict(ex.map(run, enumerate(RUNS.items())))
for label, s in st.items():
    if s["status"] != "ok":
        raise SystemExit(f"{label}: {s['status']} {s['log_tail'][-5:]}")

labels = list(RUNS)
data = [harness.load(OUT / f"run{i}.npz") for i in range(len(labels))]
ref = data[0]
ells = [2, 5, 10, 30, 100, 1000]
ks = [0.1, 1.0, 10.0, 50.0, 100.0]
print("reference: " + labels[0])
print(f"{'case':55s} " + " ".join(f"dTT/TT(l={e})" for e in ells) + "  " + " ".join(f"dP/P(k={k:g})" for k in ks))
for label, d in zip(labels[1:], data[1:]):
    n = min(d["cl_unlensed_scalar"].shape[0], ref["cl_unlensed_scalar"].shape[0])
    tt = d["cl_unlensed_scalar"][:n, 0] / ref["cl_unlensed_scalar"][:n, 0] - 1
    kmax = min(ref["pk_k"].max(), d["pk_k"].max())
    pr = np.interp(np.log(ks), np.log(ref["pk_k"]), np.log(ref["pk"][0]))
    pd = np.interp(np.log(ks), np.log(d["pk_k"]), np.log(d["pk"][0]))
    pk = np.where(np.array(ks) <= kmax, np.exp(pd - pr) - 1, np.nan)
    print(f"{label:55s} " + " ".join(f"{tt[e]:+12.2e}" for e in ells) + "  " + " ".join(f"{x:+10.2e}" for x in pk))
print("expected: all rows ~1e-4 or smaller (accuracy of the code); observed: the evolve_metric_h rows are off.")
