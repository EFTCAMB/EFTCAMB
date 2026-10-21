"""Track 10a: A0 (code under audit) vs B1 (pre-rebase EFTCAMB) for every model, corrected for the CAMB v2 changes.

For each model m and spectrum X the inherited-corrected deviation is the double ratio

    D_X(ell) = (A0_m / B1_m) / (A0_GR / B1_GR) - 1 = (A0_m / A0_GR) / (B1_m / B1_GR) - 1

with GR = 1_EFT_GR (EFTflag=0) run in the same configuration.  For TE (zero crossings) we use the difference of the
normalised differences: (A0-B1)_m / sqrt(TT EE)_B1,m - (A0-B1)_GR / sqrt(TT EE)_B1,GR.  P(k) single-tree ratios are
built on each tree's own k grid (identical across models within a tree) and then interpolated (log k, cubic) onto
a common grid.

Usage: python double_ratio.py RUNS_A RUNS_B OUT_JSON [--gr NAME] [--lmax 2500]
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
from scipy.interpolate import CubicSpline

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness  # noqa: E402

# GR (EFTflag=0) references with CAMB's own fluid dark energy matching the model's w(a), where one exists
DE_REFS = {"2_PEFT_wDE_wCDM_1": "t10a_GR_wCDM08", "4_Des_mc5e_wCDM_1": "t10a_GR_wCDM08",
           "4_Des_mc5e_wCDM_2": "t10a_GR_wCDM09", "4_Des_F_R_wCDM_1": "t10a_GR_wCDM09",
           "4_Des_F_R_wCDM_2": "t10a_GR_wCDM09", "4_Des_mc5e_CPL_1": "t10a_GR_CPL_09_01",
           "4_Des_mc5e_CPL_2": "t10a_GR_CPL_07_03"}
DE_REFS.update({m: f"t10a_GR_w_{m}" for m in ("4_Des_mc5e_taylor_1", "4_Des_mc5e_taylor_2", "4_Des_mc5e_turning_point_1",
                                               "4_Des_mc5e_turning_point_2", "5_quint_exponential",
                                               "5_quint_inverse_power", "5_quint_axion", "5_quint_cubic",
                                               "5_quint_quadratic", "5_ADE_1", "5_scg")})

SPECS = {
    # label: (npz key, column)
    "TT": ("cl_lensed_scalar", 0),
    "EE": ("cl_lensed_scalar", 1),
    "BB_lens": ("cl_lensed_scalar", 2),
    "TE": ("cl_lensed_scalar", 3),
    "uTT": ("cl_unlensed_scalar", 0),
    "uEE": ("cl_unlensed_scalar", 1),
    "uTE": ("cl_unlensed_scalar", 3),
    "tTT": ("cl_tensor", 0),
    "tEE": ("cl_tensor", 1),
    "tBB": ("cl_tensor", 2),
    "PP": ("cl_lens_potential", 0),
}


def _cl(d: dict, key: str, col: int, lmin: int, lmax: int) -> np.ndarray | None:
    if key not in d:
        return None
    arr = d[key]
    if len(arr) <= lmax:
        return None
    return arr[lmin : lmax + 1, col]


def single(a: dict, b: dict, lmin=2, lmax=2500) -> dict:
    """Plain max fractional difference A vs B (no GR correction)."""
    res = {}
    for lab, (key, col) in SPECS.items():
        x, y = _cl(a, key, col, lmin, lmax), _cl(b, key, col, lmin, lmax)
        if x is None or y is None or not np.any(y):
            continue
        if lab.endswith("TE"):
            tt = _cl(b, key, 0, lmin, lmax)
            ee = _cl(b, key, 1, lmin, lmax)
            res[lab] = float(np.max(np.abs(x - y) / np.sqrt(np.abs(tt * ee))))
        else:
            m = y != 0
            res[lab] = float(np.max(np.abs(x[m] / y[m] - 1)))
    return res


def double_ratio(am: dict, bm: dict, ag: dict, bg: dict, lmin=2, lmax=2500) -> dict:
    res = {}
    ell = np.arange(lmin, lmax + 1)
    for lab, (key, col) in SPECS.items():
        xs = [_cl(d, key, col, lmin, lmax) for d in (am, bm, ag, bg)]
        if any(x is None for x in xs):
            continue
        xa, xb, ga, gb = xs
        if not np.any(xb) or not np.any(gb):
            continue
        if lab.endswith("TE"):
            na = np.sqrt(np.abs(_cl(bm, key, 0, lmin, lmax) * _cl(bm, key, 1, lmin, lmax)))
            ng = np.sqrt(np.abs(_cl(bg, key, 0, lmin, lmax) * _cl(bg, key, 1, lmin, lmax)))
            d = (xa - xb) / na - (ga - gb) / ng
            raw = np.abs(xa - xb) / na
        else:
            m = (xb != 0) & (gb != 0) & (ga != 0)
            d = np.full_like(xa, np.nan)
            d[m] = (xa[m] / xb[m]) / (ga[m] / gb[m]) - 1
            raw = np.full_like(xa, np.nan)
            raw[m] = np.abs(xa[m] / xb[m] - 1)
        i = int(np.nanargmax(np.abs(d)))
        res[lab] = {
            "max": float(np.abs(d[i])),
            "ell": int(ell[i]),
            "raw": float(np.nanmax(raw)),
            # coarse profile for diagnosis: max |D| in ell bands
            "bands": [float(np.nanmax(np.abs(d[(ell >= lo) & (ell < hi)])))
                      for lo, hi in ((2, 30), (30, 200), (200, 800), (800, 1500), (1500, 2501))],
        }
    # matter power spectrum: single-tree ratios model/GR, each P(k) splined in log-log on its own grid
    if all("pk" in d for d in (am, bm, ag, bg)):
        kmin = max(d["pk_k"][0] for d in (am, bm, ag, bg))
        kmax = min(d["pk_k"][-1] for d in (am, bm, ag, bg))
        k = np.logspace(np.log10(kmin), np.log10(kmax), 800)

        def lp(d):
            return CubicSpline(np.log(d["pk_k"]), np.log(d["pk"][0]))(np.log(k))

        la, lb, lga, lgb = (lp(d) for d in (am, bm, ag, bg))
        d = np.exp((la - lga) - (lb - lgb)) - 1
        i = int(np.argmax(np.abs(d)))
        kb = ((1e-4, 1e-2), (1e-2, 0.1), (0.1, 1.0), (1.0, 100.0))
        res["Pk"] = {
            "max": float(abs(d[i])),
            "k": float(k[i]),
            "raw": float(np.max(np.abs(np.exp(la - lb) - 1))),
            "bands": [float(np.max(np.abs(d[(k >= lo) & (k < hi)]))) if np.any((k >= lo) & (k < hi)) else 0.0
                      for lo, hi in kb],
        }
        sa8 = am["sigma8"][0] / ag["sigma8"][0]
        sb8 = bm["sigma8"][0] / bg["sigma8"][0]
        res["sigma8"] = {"max": float(abs(sa8 / sb8 - 1)), "raw": float(abs(am["sigma8"][0] / bm["sigma8"][0] - 1))}
    # background
    for key in ("bg_H", "bg_tau", "bg_DA"):
        if all(key in d for d in (am, bm)):
            m = bm[key] != 0
            res[key] = float(np.max(np.abs(am[key][m] / bm[key][m] - 1)))
    return res


def eft_fields(am: dict, bm: dict) -> dict:
    """max |A-B| / max|B| for each EFT background function (grid a=1e-8..1)."""
    out = {}
    for k in am:
        if k.startswith("eft_") and k in bm and k != "eft_a_grid":
            x, y = am[k], bm[k]
            if x.shape != y.shape or x.dtype.kind != "f":
                continue
            fin = np.isfinite(x) & np.isfinite(y)
            scale = np.max(np.abs(y[fin])) if fin.any() else 0.0
            if not np.array_equal(np.isfinite(x), np.isfinite(y)):
                out[k[4:]] = "nonfinite-mismatch"
                continue
            out[k[4:]] = float(np.max(np.abs(x[fin] - y[fin])) / scale) if scale > 0 else float(np.max(np.abs(x[fin])))
    return out


def gr_ref_name(model: str) -> str:
    """GR (EFTflag=0) reference matching the model's non-EFT parameters."""
    if model.startswith(("5_Kmouflage", "5_Kmimic")):
        return "t10a_GR_nu0"
    if model.startswith("5_hdsk_"):
        return "t10a_GR_hdsk"
    return DE_REFS.get(model, "1_EFT_GR")


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("runs_a")
    p.add_argument("runs_b")
    p.add_argument("out")
    p.add_argument("--gr-a", default="", help="directory holding GR refs for A (default runs_a)")
    p.add_argument("--gr-b", default="", help="directory holding GR refs for B (default runs_b)")
    p.add_argument("--single-gr", default="", help="use this GR name for every model")
    p.add_argument("--lmax", type=int, default=2500)
    p.add_argument("--models", default="")
    a = p.parse_args()
    ra, rb = Path(a.runs_a), Path(a.runs_b)
    ga, gb = Path(a.gr_a or ra), Path(a.gr_b or rb)
    names = a.models.split(",") if a.models else sorted(x.stem for x in ra.glob("*.npz"))
    names = [n for n in names if not n.startswith("t10a_GR") and n != "1_EFT_GR"] + ["1_EFT_GR"]
    cache: dict = {}

    def gr(dirpath: Path, name: str):
        key = (str(dirpath), name)
        if key not in cache:
            cache[key] = harness.load(dirpath / f"{name}.npz")
        return cache[key]

    out = {"models": {}, "gr_dirs": [str(ga), str(gb)]}
    for name in names:
        fa, fb = ra / f"{name}.npz", rb / f"{name}.npz"
        if not (fa.exists() and fb.exists()):
            out["models"][name] = {"missing": [str(f) for f in (fa, fb) if not f.exists()]}
            continue
        grn = a.single_gr or gr_ref_name(name)
        if not ((ga / f"{grn}.npz").exists() and (gb / f"{grn}.npz").exists()):
            grn = "1_EFT_GR" if not name.startswith(("5_Kmouflage", "5_Kmimic", "5_hdsk_")) else grn
        ag, bg = gr(ga, grn), gr(gb, grn)
        am, bm = harness.load(fa), harness.load(fb)
        rec = {"gr_ref": grn, "double": double_ratio(am, bm, ag, bg, lmax=a.lmax), "single": single(am, bm, lmax=a.lmax)}
        rec["eft"] = eft_fields(am, bm)
        rec["errors_a"] = am["meta"].get("errors", {})
        rec["errors_b"] = bm["meta"].get("errors", {})
        out["models"][name] = rec
    Path(a.out).write_text(json.dumps(out, indent=1))
    # summary table
    cols = ["TT", "EE", "TE", "BB_lens", "tBB", "PP", "Pk"]
    print(f"{'model':34s} " + " ".join(f"{c:>9s}" for c in cols) + "   gr_ref")
    for name, rec in out["models"].items():
        if "double" not in rec:
            print(name, "MISSING")
            continue
        dd = rec["double"]
        vals = [dd[c]["max"] if isinstance(dd.get(c), dict) and "max" in dd[c] else float("nan") for c in cols]
        print(f"{name:34s} " + " ".join(f"{v:9.2e}" for v in vals) + f"   {rec['gr_ref']}")


if __name__ == "__main__":
    main()
