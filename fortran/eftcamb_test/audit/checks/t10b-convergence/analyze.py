"""Convergence analysis (track t10b-convergence).

python analyze.py LOW HIGH [--models pat] [--json out.json]
LOW/HIGH are run labels: 'def' = reference default runs (AUDIT_ROOT/runs/A0), otherwise WORK/runs/<label>.
For every model present in both: max fractional Cl differences (ell 2..2500; TE / sqrt(TT EE)) for the lensed total,
unlensed scalar and lensing potential, max fractional P(k) difference (cubic spline in log-log on the common k
range, 1e-4 < k < 1 /Mpc) and sigma8. Prints a table normalised to the GR (1_EFT_GR) level.
"""

from __future__ import annotations

import argparse
import fnmatch
import json
import sys
from pathlib import Path

import numpy as np
from scipy.interpolate import CubicSpline

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness  # noqa: E402

WORK = harness.AUDIT_ROOT / "work" / "t10b-convergence"
LMAX = 2500


def run_dir(label: str) -> Path:
    if label == "def":
        return harness.AUDIT_ROOT / "runs" / "A0"
    if label.startswith("/"):
        return Path(label)
    return WORK / "runs" / label


def cl_diffs(A: dict, B: dict, key: str) -> dict:
    a, b = A[f"cl_{key}"], B[f"cl_{key}"]
    n = min(len(a), len(b), LMAX + 1)
    a, b = a[2:n], b[2:n]
    ell = np.arange(2, n)
    out = {}
    if key == "lens_potential":
        r = np.abs(a[:, 0] / b[:, 0] - 1)
        i = int(np.argmax(r))
        return {"PP": (float(r[i]), int(ell[i]))}
    for j, lab in enumerate(("TT", "EE", "BB")):
        if np.any(b[:, j] != 0):
            r = np.abs(a[:, j] - b[:, j]) / np.abs(b[:, j])
            i = int(np.argmax(r))
            out[lab] = (float(r[i]), int(ell[i]))
    norm = np.sqrt(np.abs(b[:, 0] * b[:, 1]))
    r = np.abs(a[:, 3] - b[:, 3]) / norm
    i = int(np.argmax(r))
    out["TE"] = (float(r[i]), int(ell[i]))
    return out


def pk_diff(A: dict, B: dict, kmin: float = 1e-4, kmax: float = 1.0) -> tuple[float, float]:
    ka, kb = A["pk_k"], B["pk_k"]
    lo, hi = max(ka[0], kb[0], kmin), min(ka[-1], kb[-1], kmax)
    k = np.logspace(np.log10(lo), np.log10(hi), 2000)
    pa = np.exp(CubicSpline(np.log(ka), np.log(A["pk"][0]))(np.log(k)))
    pb = np.exp(CubicSpline(np.log(kb), np.log(B["pk"][0]))(np.log(k)))
    r = np.abs(pa / pb - 1)
    i = int(np.argmax(r))
    return float(r[i]), float(k[i])


def diffs(A: dict, B: dict) -> dict:
    res = {}
    for key in ("total", "unlensed_scalar", "lens_potential"):
        if f"cl_{key}" in A and f"cl_{key}" in B:
            for lab, v in cl_diffs(A, B, key).items():
                res[f"{key[:3]}_{lab}"] = v
    if "pk" in A and "pk" in B:
        res["pk"] = pk_diff(A, B)
        res["sigma8"] = (float(abs(A["sigma8"][0] / B["sigma8"][0] - 1)), 0)
    return res


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("low")
    p.add_argument("high")
    p.add_argument("--models", default="*")
    p.add_argument("--gr", default="1_EFT_GR")
    p.add_argument("--json", default="")
    p.add_argument("--keys", default="tot_TT,tot_EE,tot_BB,tot_TE,unl_TT,unl_EE,unl_TE,len_PP,pk,sigma8")
    a = p.parse_args()
    dl, dh = run_dir(a.low), run_dir(a.high)
    pats = a.models.split(",")
    stems = sorted({x.stem for x in dh.glob("*.npz")} & {x.stem for x in dl.glob("*.npz")})
    stems = [s for s in stems if any(fnmatch.fnmatch(s, pt) for pt in pats) or s == a.gr]
    results = {}
    for s in stems:
        try:
            results[s] = diffs(harness.load(dl / f"{s}.npz"), harness.load(dh / f"{s}.npz"))
        except Exception as e:  # noqa: BLE001
            results[s] = {"error": repr(e)}
    keys = a.keys.split(",")
    gr = results.get(a.gr, {})
    print(f"{a.low} vs {a.high}; GR level: " + " ".join(f"{k}={gr[k][0]:.1e}" for k in keys if k in gr))
    print(f"{'model':28s} " + " ".join(f"{k:>9s}" for k in keys) + "  flags(>max(3GR,1e-3))")
    for s, r in results.items():
        if "error" in r:
            print(f"{s:28s} ERROR {r['error']}")
            continue
        row, flags = [], []
        for k in keys:
            if k not in r:
                row.append(f"{'-':>9s}")
                continue
            v = r[k][0]
            row.append(f"{v:9.1e}")
            thr = max(3 * gr.get(k, (0, 0))[0], 1e-3)
            if v > thr:
                flags.append(f"{k}@{r[k][1]:g}")
        print(f"{s:28s} " + " ".join(row) + "  " + ",".join(flags))
    if a.json:
        Path(a.json).write_text(json.dumps(results, indent=1))


if __name__ == "__main__":
    main()
