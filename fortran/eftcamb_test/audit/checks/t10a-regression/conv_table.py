"""Accuracy budget of the A0-B1 inherited-corrected deviation for selected models.

For model m and spectrum X, R_t(cfg) = X_m / X_GR (GR = matching EFTflag=0 reference run with the same tree and
configuration; TE uses (TE_m/sqrt(TT EE)_m) - (TE_GR/sqrt(TT EE)_GR)).  Printed per ell band (max |.|):
  D      = R_A0(def) / R_B1(def) - 1         the double ratio of the main table
  cA0    = R_A0(def) / R_A0(hi) - 1          A0 default-accuracy error (hi = accuracy_boost 2)
  B1hi   = R_B1(def) / R_A0(hi) - 1          B1 default vs A0 high accuracy
  cB1t   = R_B1(def) / R_B1(itb2) - 1        B1 time-integration error alone (IntTolBoost 2)
  Dhi    = R_A0(hi) / R_B1(hi) - 1           (if B1 high-accuracy runs exist)
Bands: ell 2-30, 30-200, 200-800, 800-2500 (tBB: 2-30, 30-200, 200-500).

Usage: python conv_table.py WORK --models a,b --def-a DIR --def-b DIR --hi-a DIR [--hi-b DIR] [--itb2-b DIR]
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import harness  # noqa: E402
from double_ratio import gr_ref_name  # noqa: E402

SPEC = {"TT": ("cl_lensed_scalar", 0), "EE": ("cl_lensed_scalar", 1), "TE": ("cl_lensed_scalar", 3),
        "tBB": ("cl_tensor", 2), "PP": ("cl_lens_potential", 0)}
BANDS = ((2, 30), (30, 200), (200, 800), (800, 2501))
TBANDS = ((2, 30), (30, 200), (200, 501))


def resp(m, g, lab):
    key, col = SPEC[lab]
    if lab == "TE":
        f = lambda d: d[key][:2501, 3] / np.sqrt(np.abs(d[key][:2501, 0] * d[key][:2501, 1]))  # noqa: E731
        return ("d", f(m) - f(g))
    return ("r", m[key][:2501, col] / g[key][:2501, col])


def dev(r1, r2):
    if r1[0] == "d":
        return np.abs(r1[1] - r2[1])
    return np.abs(r1[1] / r2[1] - 1)


def bands(x, lab):
    bb = TBANDS if lab == "tBB" else BANDS
    return [float(np.nanmax(x[lo:hi])) for lo, hi in bb]


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("work")
    p.add_argument("--models", required=True)
    p.add_argument("--def-a", required=True)
    p.add_argument("--def-b", required=True)
    p.add_argument("--hi-a", required=True)
    p.add_argument("--hi-b", default="")
    p.add_argument("--itb2-b", default="")
    p.add_argument("--gr", default="", help="force this GR reference name")
    a = p.parse_args()
    w = Path(a.work)
    dirs = {"da": w / a.def_a, "db": w / a.def_b, "ha": w / a.hi_a}
    if a.hi_b:
        dirs["hb"] = w / a.hi_b
    if a.itb2_b:
        dirs["ib"] = w / a.itb2_b
    for m in a.models.split(","):
        g = a.gr or gr_ref_name(m)
        R = {}
        for k, d in dirs.items():
            fm, fg = d / f"{m}.npz", d / f"{g}.npz"
            if fm.exists() and fg.exists():
                R[k] = (harness.load(fm), harness.load(fg))
        if not {"da", "db", "ha"} <= set(R):
            print(f"{m}: missing runs {sorted({'da', 'db', 'ha'} - set(R))}")
            continue
        print(f"== {m}  (GR ref {g})")
        for lab in SPEC:
            rr = {k: resp(*v, lab) for k, v in R.items()}
            cols = {"D": dev(rr["da"], rr["db"]), "cA0": dev(rr["da"], rr["ha"]), "B1hi": dev(rr["db"], rr["ha"])}
            if "ib" in rr:
                cols["cB1t"] = dev(rr["db"], rr["ib"])
            if "hb" in rr:
                cols["Dhi"] = dev(rr["ha"], rr["hb"])
            print(f"  {lab:4s} " + "  ".join(f"{k}=" + ",".join(f"{v:.0e}" for v in bands(x, lab)) for k, x in cols.items()))


if __name__ == "__main__":
    main()
