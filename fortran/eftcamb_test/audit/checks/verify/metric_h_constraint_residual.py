"""Verifier (metric-h-ic-violates-constraint): independent check via the per-k time evolution cache.

For 5_hdsk_TG_1 with evolve_delta_phi = T, run evolve_metric_h = F and T on frozen A0 and, along the evolution of
selected k, compare the h' actually used (2 k z from the cache) with the Horndeski 00 constraint evaluated from the
cached etak, dgrho, Meff2, alphaB, adotoa, s00*, delta_phi:
    h'_c = (4 k etak + 2 dgrho/Meff2)/((2 - alphaB_BS) adotoa) + s00p dphi' + (s00 + s00k k^2) dphi,  alphaB_BS = -2 alphaB
Residual r = 2 k z / h'_c - 1. Expected ~0 on the constraint path (sanity check of the formula); a non-zero value on
the evolved path just after the turn-on measures the IC violation. Also prints GR 2kz/h'_c at the turn-on from cache.
"""

import json
import os
import sys
from pathlib import Path

os.environ["OMP_NUM_THREADS"] = "1"
REPO = Path(__file__).resolve().parents[5]
sys.path.insert(0, str(REPO / "fortran" / "eftcamb_test" / "audit" / "checks"))
import harness  # noqa: E402

import numpy as np  # noqa: E402

TREE = Path(sys.argv[2]) if len(sys.argv) > 2 else harness.AUDIT_ROOT / "trees" / "A0"
TAG = "" if len(sys.argv) <= 2 else "_" + Path(sys.argv[2]).name
OUT = harness.AUDIT_ROOT / "work" / "verify" / "metric_h"
KS = [0.01, 1.0, 10.0]
A = np.logspace(-8.0, 0.0, 161)


def main():
    mode = sys.argv[1]
    sys.path.insert(0, str(TREE))
    import camb

    with harness.slot(1):
        ini = harness.PARAMS_DIR / "5_hdsk_TG_1.ini"
        OUT.mkdir(parents=True, exist_ok=True)
        tmp = OUT / f"run_{mode}.ini"
        tmp.write_text(f"DEFAULT({ini})\nEFTCAMB_evolve_metric_h = {mode}\nEFTCAMB_write_background = F\n")
        pars = camb.read_ini(str(tmp))
        res = camb.get_results(pars)
        a = A[A > 1.0001e-8]
        fields, arr = pars.EFTCAMB.get_scale_evolution(res, KS, a)
    out = {}
    for i, k in enumerate(KS):
        c = arr[i] if arr.shape[0] == len(KS) else arr[:, i]
        alb = -2 * c["alphaB"]
        hc = (4 * k * c["etak"] + 2 * c["dgrho"] / c["Meff2"]) / ((2 - alb) * c["adotoa"]) + c["s00p"] * c[
            "delta_phidot"] + (c["s00"] + c["s00k"] * k * k) * c["delta_phi"]
        hgr = 2 * k * (0.5 * c["dgrho"] / k + c["etak"]) / c["adotoa"]
        h = 2 * k * c["z"]
        r = h / hc - 1
        out[str(k)] = {"a": a.tolist(), "r": r.tolist(), "gr_over_c": (hgr / hc - 1).tolist(),
                       "Meff2": c["Meff2"].tolist(), "alphaB": c["alphaB"].tolist()}
        sel = [0, 5, 20, 40, 80, 120, len(a) - 1]
        print(f"mode={mode} k={k}: " + "  ".join(f"a={a[j]:.1e}: r={r[j]:+.2e}" for j in sel))
        print(f"      2kz_GR/h'_c-1 at first a: {hgr[0] / hc[0] - 1:+.3e}; Meff2={c['Meff2'][0]:.4f} "
              f"alphaB_EFT={c['alphaB'][0]:+.3e}")
    (OUT / f"resid_{mode}{TAG}.json").write_text(json.dumps(out))


if __name__ == "__main__":
    main()
