"""Verify finding qs-mu-nan-vanishing-denominator by an independent route.

Part 1 (no new CAMB runs): for each shipped model, take the background timestep-cache arrays stored by the harness in
the A0 reference runs (eft_<field> from get_eft_functions, a = logspace(-8,0,401)), evaluate numerator and denominator
of the QS_mu and QS_mu_gen closed forms in numpy (qsnan_formula.py), at k = 1e-3, 0.01, 0.1, 1 /Mpc and 0.02 <= a <= 1,
and classify each point: 0/0, x/0, finite.  This predicts which models get NaN without looking at the Fortran output.

Part 2 (--evol DIR): for harness runs with evol_k, compare the Fortran QS_mu / QS_mu_gen / QS_sigma output with the
numpy evaluation on the same (k, a) cache values; report NaN agreement and the max relative difference at finite points.

Usage: python verify_qsnan_scan.py [--evol DIR]
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness  # noqa: E402
import qsnan_formula as qf  # noqa: E402

REF = Path(
    "/private/tmp/claude-501/-Users-marco-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/"
    "scratchpad/audit/runs/A0"
)
PAR = Path(__file__).resolve().parents[3] / "parameters"
C_LIGHT = 2.99792458e8


def h0_of(model: str) -> float:
    p = harness.parse_ini(PAR / f"{model}.ini")
    return float(p.get("hubble", 67.0)) / C_LIGHT * 1000.0


def classify(num, den):
    num, den = np.asarray(num), np.asarray(den)
    z0 = (den == 0) & (num == 0)
    zx = (den == 0) & (num != 0)
    bad = ~np.isfinite(num) | ~np.isfinite(den)
    return int(z0.sum()), int(zx.sum()), int(bad.sum()), int(num.size)


def part1(F):
    rows = {}
    ks = np.array([1e-3, 0.01, 0.1, 1.0])[:, None]
    for f in sorted(REF.glob("[2-9]*.npz")):
        d = harness.load(f)
        a = d["eft_a_grid"]
        sel = a >= 0.02
        C = {k[4:]: d[k][sel][None, :] for k in d if k.startswith("eft_") and k != "eft_a_grid"}
        res = qf.evaluate(F, C, a[sel][None, :], ks, h0_of(f.stem))
        row = {}
        for fld, (num, den) in res.items():
            n00, nx0, nbad, n = classify(num, den)
            row[fld] = {"0/0": n00, "x/0": nx0, "nonfinite_terms": nbad, "n": n,
                        "max|num|": float(np.nanmax(np.abs(num))), "max|den|": float(np.nanmax(np.abs(den)))}
        rows[f.stem] = row
    return rows


def part2(F, evdir: Path):
    out = {}
    for f in sorted(evdir.glob("*.npz")):
        d = harness.load(f)
        if "evol_k" not in d:
            continue
        # the harness stores the per-point cache fields 'k' and 'a' as evol_k / evol_a, shape (nk, na)
        k = d["evol_k"] if d["evol_k"].ndim == 2 else d["evol_k"][:, None]
        a = d["evol_a"] if d["evol_a"].ndim == 2 else d["evol_a"][None, :]
        C = {kk[5:]: d[kk] for kk in d if kk.startswith("evol_")}
        model = f.stem.split("__")[0]
        res = qf.evaluate(F, C, a, k, h0_of(model))
        # points where Fortran actually evaluated the EFT QS closed form: not exactly 1 (GR branch) and a >= 0.012
        act = np.broadcast_to(a, C["QS_mu"].shape) > 0.012
        row = {}
        for fld in ("QS_mu", "QS_sigma", "QS_mu_gen"):
            num, den = res[fld]
            with np.errstate(all="ignore"):
                py = num / den
            fo = C[fld]
            m = act
            nan_f, nan_p = np.isnan(fo[m]), np.isnan(py[m])
            fin = m & np.isfinite(fo) & np.isfinite(py)
            rel = float(np.max(np.abs(py[fin] - fo[fin]) / np.maximum(np.abs(fo[fin]), 1e-300))) if fin.any() else None
            n00, nx0, _, _ = classify(num[m], den[m])
            row[fld] = {"n": int(m.sum()), "fortran_nan": int(nan_f.sum()), "numpy_nan": int(nan_p.sum()),
                        "nan_mismatch": int((nan_f != nan_p).sum()), "numpy_0/0": n00, "numpy_x/0": nx0,
                        "max_rel_diff_finite": rel,
                        "fortran_sample(k0,a=1)": float(fo[0, -1])}
        for fld in ("QS_sigma_gen", "mu", "gamma", "sigma_eft"):
            v = C[fld][act]
            row[fld] = {"nonfinite": int((~np.isfinite(v)).sum()), "n": int(v.size),
                        "range": [float(np.nanmin(v)), float(np.nanmax(v))] if np.isfinite(v).any() else None}
        out[f.stem] = row
    return out


def main():
    F = qf.load_formulas()
    print("QS_mu == QS_sigma textually:", F["QS_mu"] == F["QS_sigma"])
    if "--evol" in sys.argv:
        evdir = Path(sys.argv[sys.argv.index("--evol") + 1])
        res = part2(F, evdir)
        for m, row in res.items():
            print(m)
            for fld, v in row.items():
                print(f"   {fld:13s} {v}")
        (evdir / "part2.json").write_text(json.dumps(res, indent=1))
        return
    rows = part1(F)
    nanmods = []
    for m, row in rows.items():
        flags = {fld: v["0/0"] for fld, v in row.items() if v["0/0"] or v["x/0"] or v["nonfinite_terms"]}
        x0 = {fld: v["x/0"] for fld, v in row.items() if v["x/0"]}
        if flags:
            nanmods.append(m)
        print(f"{m:30s} 0/0 points (of {row['QS_mu']['n']}): {flags or 'none'}" + (f"  x/0: {x0}" if x0 else ""))
    print(f"models with predicted 0/0: {len(nanmods)}: {nanmods}")
    out = Path(
        "/private/tmp/claude-501/-Users-marco-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/"
        "scratchpad/audit/work/verify/qsnan/part1.json"
    )
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(rows, indent=1))


if __name__ == "__main__":
    main()
