"""Quick per-ell comparison of two runs (track t10b-convergence).

python clcmp.py A.npz B.npz [--key unlensed_scalar] [--ells 2,10,...]
Prints max fractional differences per spectrum (ell 2..2500, TE/sqrt(TT EE)), location, and TT ratios at sample ells.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness  # noqa: E402
from analyze import diffs  # noqa: E402


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("a")
    p.add_argument("b")
    p.add_argument("--key", default="unlensed_scalar")
    p.add_argument("--ells", default="2,5,10,20,50,100,200,300,500,1000,1500,2000,2500")
    a = p.parse_args()
    A, B = harness.load(a.a), harness.load(a.b)
    d = diffs(A, B)
    print(" ".join(f"{k}={v[0]:.2e}@{v[1]:g}" for k, v in d.items()))
    ca, cb = A[f"cl_{a.key}"][:, 0], B[f"cl_{a.key}"][:, 0]
    print(f"{a.key} TT A/B-1:", " ".join(f"{l}:{ca[l] / cb[l] - 1:+.2e}" for l in map(int, a.ells.split(","))))


if __name__ == "__main__":
    main()
