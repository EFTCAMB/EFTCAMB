"""Detect isolated single-k-mode glitches in P(k) and single-ell glitches in C_ell (track t10b-convergence).

python spikes.py LABEL [LABEL ...] [--models pat] [--thr 0.05]
For each run: max |second difference of log P(k)| over 1e-3<k<1 /Mpc (z index 0) and of log D_ell^TT over
30<ell<2400, with its location. A smooth spectrum (GR, BAO and acoustic peaks) gives values well below thr.
"""

from __future__ import annotations

import argparse
import fnmatch
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness  # noqa: E402
from analyze import run_dir  # noqa: E402


def pk_spike(d: dict, kmin: float = 1e-3, kmax: float = 1.0) -> tuple[float, float]:
    k, p = d["pk_k"], d["pk"][0]
    m = (k > kmin) & (k < kmax)
    k, lp = k[m], np.log(p[m])
    lk = np.log(k)
    # second divided difference normalised to the local log-k spacing (non-uniform grid)
    h1, h2 = lk[1:-1] - lk[:-2], lk[2:] - lk[1:-1]
    pred = (lp[:-2] * h2 + lp[2:] * h1) / (h1 + h2)
    d2 = lp[1:-1] - pred
    i = int(np.argmax(np.abs(d2)))
    return float(d2[i]), float(k[i + 1])


def cl_spike(d: dict, lmin: int = 30, lmax: int = 2400) -> tuple[float, int]:
    c = d["cl_unlensed_scalar"][:, 0]
    ell = np.arange(len(c))
    lc = np.log(np.abs(c[lmin:lmax]) * ell[lmin:lmax] * (ell[lmin:lmax] + 1))
    d2 = lc[1:-1] - 0.5 * (lc[:-2] + lc[2:])
    i = int(np.argmax(np.abs(d2)))
    return float(d2[i]), int(ell[lmin + i + 1])


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("labels", nargs="+")
    p.add_argument("--models", default="*")
    p.add_argument("--thr", type=float, default=0.05)
    a = p.parse_args()
    pats = a.models.split(",")
    stems = sorted({x.stem for lab in a.labels for x in run_dir(lab).glob("*.npz")})
    stems = [s for s in stems if any(fnmatch.fnmatch(s, pt) for pt in pats)]
    print(f"{'model':28s} " + " ".join(f"{lab + ':pk':>18s} {lab + ':cl':>14s}" for lab in a.labels))
    for s in stems:
        row, flag = [], False
        for lab in a.labels:
            f = run_dir(lab) / f"{s}.npz"
            if not f.exists():
                row.append(f"{'-':>18s} {'-':>14s}")
                continue
            d = harness.load(f)
            try:
                ps, pk_at = pk_spike(d)
                cs, cl_at = cl_spike(d)
            except Exception:  # noqa: BLE001
                row.append(f"{'err':>18s} {'err':>14s}")
                continue
            flag |= abs(ps) > a.thr or abs(cs) > a.thr
            row.append(f"{ps:+8.1e}@{pk_at:<8.4f} {cs:+8.1e}@{cl_at:<5d}")
        print(f"{s:28s} " + " ".join(row) + ("  <-- glitch" if flag else ""))


if __name__ == "__main__":
    main()
