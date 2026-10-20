"""Compare legacy EFTCAMB text outputs (B2, fortran/eftcamb_test/eftcamb_legacy/spectra_results) with harness runs.

Legacy formats (CAMB fortran output, CMB_outputscale = T_cmb^2 in muK^2):
  *_lensedCls.dat        L, D_TT, D_EE, D_BB, D_TE                    D = l(l+1) C_l / 2pi  [muK^2]
  *_tensCls.dat          L, D_TT, D_EE, D_BB, D_TE                    (tensor only)
  *_lenspotentialCls.dat L, TT, EE, BB, TE (unlensed total), PP, TP, EP  PP = [l(l+1)]^2 C_phi / 2pi
  *_matterpower.dat      k/h [h/Mpc], P [(Mpc/h)^3]
Harness npz: raw C_l in muK^2 (cl_lensed_scalar, cl_tensor, cl_lens_potential[:,0] = C_phi), pk_k [1/Mpc], pk [Mpc^3].

Two comparisons per model and run set X (A0 or B1):
  raw    max |X/B2 - 1| (TE: |dTE|/sqrt(TT EE))
  ratio  max |(X_m/X_GR) / (B2_m/B2_GR) - 1|  (EFT/GR response; GR = 1_EFT_GR legacy, or none for omnuh2=0 / hdsk)

Usage: python legacy_compare.py RUNS_DIR OUT_JSON [--gr-dir DIR] [--lmax 2500]
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

LEG = harness.REPO / "fortran" / "eftcamb_test" / "eftcamb_legacy" / "spectra_results"
LEGACY_NAME = {"5_scg": "5_sgc"}


def legacy(name: str) -> dict | None:
    root = LEG / LEGACY_NAME.get(name, name)
    f = Path(f"{root}_lensedCls.dat")
    if not f.exists():
        return None
    out = {}
    d = np.loadtxt(f"{root}_lensedCls.dat")
    ell = d[:, 0].astype(int)
    fac = 2 * np.pi / (ell * (ell + 1.0))
    arr = np.zeros((ell.max() + 1, 4))
    arr[ell] = d[:, 1:5] * fac[:, None]
    out["cl_lensed_scalar"] = arr
    d = np.loadtxt(f"{root}_tensCls.dat")
    ell = d[:, 0].astype(int)
    fac = 2 * np.pi / (ell * (ell + 1.0))
    arr = np.zeros((ell.max() + 1, 4))
    arr[ell] = d[:, 1:5] * fac[:, None]
    out["cl_tensor"] = arr
    d = np.loadtxt(f"{root}_lenspotentialCls.dat")
    ell = d[:, 0].astype(int)
    arr = np.zeros((ell.max() + 1, 3))
    arr[ell, 0] = d[:, 5] * 2 * np.pi / (ell * (ell + 1.0)) ** 2
    out["cl_lens_potential"] = arr
    d = np.loadtxt(f"{root}_matterpower.dat")
    out["pk_kh"], out["pk_h"] = d[:, 0], d[:, 1]
    return out


def h_of(name: str) -> float:
    ini = harness.PARAMS_DIR / f"{name}.ini"
    return float(harness.parse_ini(ini)["hubble"]) / 100.0


SPEC = {"TT": ("cl_lensed_scalar", 0), "EE": ("cl_lensed_scalar", 1), "BB": ("cl_lensed_scalar", 2),
        "TE": ("cl_lensed_scalar", 3), "tBB": ("cl_tensor", 2), "tTT": ("cl_tensor", 0), "PP": ("cl_lens_potential", 0)}


def _arr(d, key, col, lmin, lmax):
    return d[key][lmin:lmax + 1, col]


def raw_cmp(x: dict, b: dict, lmin=2, lmax=2500) -> dict:
    res = {}
    for lab, (key, col) in SPEC.items():
        xa, xb = _arr(x, key, col, lmin, lmax), _arr(b, key, col, lmin, lmax)
        if lab == "TE":
            nrm = np.sqrt(np.abs(_arr(b, key, 0, lmin, lmax) * _arr(b, key, 1, lmin, lmax)))
            res[lab] = float(np.max(np.abs(xa - xb) / nrm))
        else:
            m = xb != 0
            res[lab] = float(np.max(np.abs(xa[m] / xb[m] - 1))) if m.any() else float("nan")
    return res


def ratio_cmp(xm, xg, bm, bg, lmin=2, lmax=2500) -> dict:
    res = {}
    for lab, (key, col) in SPEC.items():
        a, ag, b, bgr = (_arr(d, key, col, lmin, lmax) for d in (xm, xg, bm, bg))
        if lab == "TE":
            na = np.sqrt(np.abs(_arr(bm, key, 0, lmin, lmax) * _arr(bm, key, 1, lmin, lmax)))
            ng = np.sqrt(np.abs(_arr(bg, key, 0, lmin, lmax) * _arr(bg, key, 1, lmin, lmax)))
            res[lab] = float(np.max(np.abs((a - b) / na - (ag - bgr) / ng)))
        else:
            m = (b != 0) & (bgr != 0) & (ag != 0)
            res[lab] = float(np.max(np.abs((a[m] / ag[m]) / (b[m] / bgr[m]) - 1))) if m.any() else float("nan")
    return res


def pk_interp(d: dict, k: np.ndarray, h: float | None) -> np.ndarray:
    if "pk_kh" in d:
        kk, pp = d["pk_kh"] * h, d["pk_h"] / h**3
    else:
        kk, pp = d["pk_k"], d["pk"][0]
    return CubicSpline(np.log(kk), np.log(pp))(np.log(k))


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("runs")
    p.add_argument("out")
    p.add_argument("--gr-dir", default="")
    p.add_argument("--lmax", type=int, default=2500)
    a = p.parse_args()
    runs = Path(a.runs)
    grdir = Path(a.gr_dir or runs)
    bgr = legacy("1_EFT_GR")
    xgr = harness.load(grdir / "1_EFT_GR.npz")
    out = {}
    for ini in harness.model_inis():
        name = ini.stem
        f = runs / f"{name}.npz"
        b = legacy(name)
        if b is None or not f.exists():
            out[name] = {"missing": True}
            continue
        x = harness.load(f)
        h = h_of(name)
        rec = {"raw": raw_cmp(x, b, lmax=a.lmax)}
        kmin = max(x["pk_k"][0], b["pk_kh"][0] * h) * 1.001
        kmax = min(x["pk_k"][-1], b["pk_kh"][-1] * h) * 0.999
        k = np.logspace(np.log10(kmin), np.log10(kmax), 600)
        lx, lb = pk_interp(x, k, h), pk_interp(b, k, h)
        rec["raw"]["Pk"] = float(np.max(np.abs(np.exp(lx - lb) - 1)))
        if not name.startswith(("5_Kmouflage", "5_Kmimic", "5_hdsk")):
            rec["ratio"] = ratio_cmp(x, xgr, b, bgr, lmax=a.lmax)
            lxg, lbg = pk_interp(xgr, k, 0.7), pk_interp(bgr, k, 0.7)
            rec["ratio"]["Pk"] = float(np.max(np.abs(np.exp((lx - lxg) - (lb - lbg)) - 1)))
        out[name] = rec
    Path(a.out).write_text(json.dumps(out, indent=1))
    cols = ["TT", "EE", "TE", "BB", "tBB", "PP", "Pk"]
    print(f"{'model':30s} kind  " + " ".join(f"{c:>9s}" for c in cols))
    for name, rec in out.items():
        for kind in ("raw", "ratio"):
            if kind in rec:
                print(f"{name:30s} {kind:5s} " + " ".join(f"{rec[kind][c]:9.2e}" for c in cols))


if __name__ == "__main__":
    main()
