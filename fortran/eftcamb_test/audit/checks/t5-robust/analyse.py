"""Track 5: compare every variation run against the default A0 reference run of the same model.

Metrics (max over ell 2..2500 of the fractional difference): lensed total TT, EE, TE/sqrt(TT EE), BB (total, includes
the r=1 tensors), unlensed scalar TT/EE, lensing potential PP; tensor TT/EE/BB restricted to ell <= 300;
P(k, z=0) max fractional difference; sigma8 fractional difference.

Usage: python analyse.py [VAR,...] [--json out.json] [--ref DIR]
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness  # noqa: E402

WORK = harness.AUDIT_ROOT / "work" / "t5-robust"
REF = harness.AUDIT_ROOT / "runs" / "A0"


def frac(a, b, lmin=2, lmax=2500):
    n = min(len(a), len(b), lmax + 1)
    x, y = a[lmin:n], b[lmin:n]
    m = y != 0
    return float(np.max(np.abs(x[m] - y[m]) / np.abs(y[m]))) if m.any() else 0.0


def te(a, b, lmax=2500):
    n = min(len(a), len(b), lmax + 1)
    x, y = a[2:n], b[2:n]
    norm = np.sqrt(np.abs(y[:, 0] * y[:, 1]))
    m = norm > 0
    return float(np.max(np.abs(x[m, 3] - y[m, 3]) / norm[m]))


def metrics(A: dict, B: dict) -> dict:
    out = {}
    if "cl_total" in A and "cl_total" in B:
        for i, s in enumerate(("TT", "EE", "BB")):
            out[f"cl_{s}"] = frac(A["cl_total"][:, i], B["cl_total"][:, i])
        out["cl_TE"] = te(A["cl_total"], B["cl_total"])
    if "cl_unlensed_scalar" in A and "cl_unlensed_scalar" in B:
        for i, s in enumerate(("TT", "EE")):
            out[f"uls_{s}"] = frac(A["cl_unlensed_scalar"][:, i], B["cl_unlensed_scalar"][:, i])
    if "cl_lens_potential" in A and "cl_lens_potential" in B:
        out["cl_PP"] = frac(A["cl_lens_potential"][:, 0], B["cl_lens_potential"][:, 0])
    if "cl_tensor" in A and "cl_tensor" in B:
        for i, s in enumerate(("TT", "EE", "BB")):
            out[f"ten_{s}"] = frac(A["cl_tensor"][:, i], B["cl_tensor"][:, i], lmax=300)
    if "pk" in A and "pk" in B:
        if A["pk"].shape == B["pk"].shape and np.allclose(A["pk_k"], B["pk_k"], rtol=1e-10):
            out["pk"] = float(np.max(np.abs(A["pk"][0] - B["pk"][0]) / np.abs(B["pk"][0])))
        else:
            kmin, kmax = max(A["pk_k"][0], B["pk_k"][0]), min(A["pk_k"][-1], B["pk_k"][-1])
            k = np.logspace(np.log10(kmin), np.log10(kmax), 300)
            pa = np.exp(np.interp(np.log(k), np.log(A["pk_k"]), np.log(A["pk"][0])))
            pb = np.exp(np.interp(np.log(k), np.log(B["pk_k"]), np.log(B["pk"][0])))
            out["pk"] = float(np.max(np.abs(pa - pb) / pb))
        out["sigma8"] = float(abs(A["sigma8"][0] / B["sigma8"][0] - 1))
    scal = [v for k, v in out.items() if k.startswith(("cl_", "uls_"))]
    out["worst_scalar_cl"] = max(scal) if scal else float("nan")
    out["worst_tensor_cl"] = max([v for k, v in out.items() if k.startswith("ten_")] or [float("nan")])
    return out


def analyse(var: str, ref: Path = REF) -> dict:
    res = {}
    d = WORK / "runs" / var
    for js in sorted(d.glob("*.json")):
        st = json.loads(js.read_text())
        model = js.stem
        entry = {"status": st["status"]}
        npz = js.with_suffix(".npz")
        if st["status"] == "ok" and npz.exists():
            A = harness.load(npz)
            rf = ref / f"{model}.npz"
            if not rf.exists():
                rf = harness.AUDIT_ROOT / "work" / "prep" / "runs" / f"{model}.npz"
            B = harness.load(rf)
            entry.update(metrics(A, B))
            entry["errors"] = A["meta"].get("errors", {})
        else:
            entry["log_tail"] = st.get("log_tail", [])[-6:]
        res[model] = entry
    return res


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("variations")
    p.add_argument("--json", default="")
    p.add_argument("--ref", default=str(REF))
    a = p.parse_args()
    allres = {}
    for var in a.variations.split(","):
        r = analyse(var, Path(a.ref))
        allres[var] = r
        print(f"=== {var} ({len(r)} models)")
        for m, e in r.items():
            if e["status"] != "ok":
                print(f"  {m:32s} {e['status']}  {' | '.join(x.strip() for x in e.get('log_tail', [])[-2:])[:150]}")
                continue
            print(f"  {m:32s} scal={e['worst_scalar_cl']:.2e} ten={e['worst_tensor_cl']:.2e} "
                  f"pk={e.get('pk', float('nan')):.2e} s8={e.get('sigma8', float('nan')):.2e} "
                  f"TT={e.get('cl_TT', 0):.1e} EE={e.get('cl_EE', 0):.1e} PP={e.get('cl_PP', 0):.1e}"
                  + (f" ERR={e['errors']}" if e.get("errors") else ""))
    if a.json:
        Path(a.json).write_text(json.dumps(allres, indent=1))


if __name__ == "__main__":
    main()
