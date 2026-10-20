"""Separate accuracy noise from systematic A0-B1 differences for selected models.

For each model m and spectrum X, with G the matching GR reference run in the same configuration:
  R_t(cfg)   = X_m / X_G in tree t and configuration cfg (the 'MG response'), TE as difference / sqrt(TT EE)
  D(cfg)     = R_A0(cfg) / R_B1(cfg) - 1                 (inherited-corrected A0-B1 deviation)
  conv_t     = R_t(def) / R_t(ab2) - 1                   (default-accuracy error of each tree's response)
Printed: max over ell (2..2500) / k of |D(def)|, |D(ab2)|, |conv_A0|, |conv_B1|.

Usage: python convergence.py WORKDIR [--models a,b] [--def-a runs_A0_itb2] [--def-b runs_B1_itb1]
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

SPEC = {"TT": ("cl_lensed_scalar", 0), "EE": ("cl_lensed_scalar", 1), "TE": ("cl_lensed_scalar", 3),
        "BB": ("cl_lensed_scalar", 2), "uEE": ("cl_unlensed_scalar", 1), "tBB": ("cl_tensor", 2),
        "PP": ("cl_lens_potential", 0)}
KGRID = np.logspace(np.log10(2e-5), np.log10(1.3), 500)


def gr_name(m: str) -> str:
    if m.startswith(("5_Kmouflage", "5_Kmimic")):
        return "t10a_GR_nu0"
    if m.startswith("5_hdsk_"):
        return "t10a_GR_hdsk"
    return DE_REFS.get(m, "1_EFT_GR")


def response(m: dict, g: dict, lab: str):
    if lab == "Pk":
        f = lambda d: CubicSpline(np.log(d["pk_k"]), np.log(d["pk"][0]))(np.log(KGRID))  # noqa: E731
        return f(m) - f(g)  # log response
    key, col = SPEC[lab]
    x, y = m[key][2:2501, col], g[key][2:2501, col]
    if lab == "TE":
        nm = np.sqrt(np.abs(m[key][2:2501, 0] * m[key][2:2501, 1]))
        ng = np.sqrt(np.abs(g[key][2:2501, 0] * g[key][2:2501, 1]))
        return ("TE", x / nm, y / ng)
    return np.log(x / y)


def dev(ra, rb):
    if isinstance(ra, tuple):
        return np.max(np.abs((ra[1] - ra[2]) - (rb[1] - rb[2])))
    return float(np.max(np.abs(np.expm1(ra - rb))))


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("work")
    p.add_argument("--models", default="")
    p.add_argument("--def-a", default="runs_A0_itb2")
    p.add_argument("--def-b", default="runs_B1_itb1")
    p.add_argument("--hi-a", default="runs_A0_ab2")
    p.add_argument("--hi-b", default="runs_B1_ab2")
    p.add_argument("--out", default="")
    a = p.parse_args()
    w = Path(a.work)
    da, db, ha, hb = (w / x for x in (a.def_a, a.def_b, a.hi_a, a.hi_b))
    models = a.models.split(",") if a.models else sorted(x.stem for x in ha.glob("*.npz")
                                                          if not x.stem.startswith(("t10a_", "1_EFT_GR")))
    labs = list(SPEC) + ["Pk"]
    res = {}
    print(f"{'model':28s} {'spec':5s} {'D(def)':>9s} {'D(hi)':>9s} {'convA0':>9s} {'convB1':>9s}")
    for m in models:
        fs = [d / f"{m}.npz" for d in (da, db, ha, hb)] + [d / f"{gr_name(m)}.npz" for d in (da, db, ha, hb)]
        if not all(f.exists() for f in fs):
            print(f"{m:28s} missing {[str(f) for f in fs if not f.exists()]}")
            continue
        M = [harness.load(f) for f in fs]
        res[m] = {}
        for lab in labs:
            try:
                r = [response(M[i], M[i + 4], lab) for i in range(4)]
            except (KeyError, ValueError, IndexError) as e:
                res[m][lab] = repr(e)
                continue
            vals = (dev(r[0], r[1]), dev(r[2], r[3]), dev(r[0], r[2]), dev(r[1], r[3]))
            res[m][lab] = dict(zip(("D_def", "D_hi", "conv_A0", "conv_B1"), map(float, vals)))
            print(f"{m:28s} {lab:5s} " + " ".join(f"{v:9.2e}" for v in vals))
    if a.out:
        Path(a.out).write_text(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
