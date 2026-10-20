"""Detail comparison of the runs made by verify_lag_conv.py: Cl differences by ell range and P(k) differences by k,
and whether the code under audit (a0) or the fixed-point fix (convfix) is closer to the accuracy_boost=2 answer.

Usage: python verify_lag_conv_detail.py [tag]   (tag = ini stem, default 008p5-beyond-horndeski)
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import harness  # noqa: E402

RUNS = Path("/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/"
            "2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit/work/verify/lag/runs")
TAG = sys.argv[1] if len(sys.argv) > 1 else "008p5-beyond-horndeski"
ELL = [(2, 30), (30, 500), (500, 2501)]


def L(n):
    return harness.load(RUNS / f"{TAG}_{n}.npz")


def cl_rel(a, b, key="cl_unlensed_scalar", col=0):
    x, y = a[key][:, col], b[key][:, col]
    out = []
    for lo, hi in ELL:
        d = np.abs(x[lo:hi] / y[lo:hi] - 1)
        j = int(np.argmax(d))
        out.append(f"l[{lo},{hi}) {d[j]:.2e}@{lo + j}")
    return "  ".join(out)


def pk_on(d, k):
    return np.exp(np.interp(np.log(k), np.log(d["pk_k"]), np.log(d["pk"][0])))


def main() -> None:
    a0, cf, a2, cf2 = L("a0"), L("convfix"), L("a0_acc2"), L("convfix_acc2")
    print(f"model {TAG}")
    for lab, (x, y) in {"LAG EFFECT acc1  convfix vs a0": (cf, a0), "LAG EFFECT acc2  convfix vs a0": (cf2, a2),
                        "accuracy acc2 vs acc1 (a0)": (a2, a0), "accuracy acc2 vs acc1 (convfix)": (cf2, cf)}.items():
        print(f"{lab}:")
        print("   unlensed TT  ", cl_rel(x, y))
        print("   unlensed EE  ", cl_rel(x, y, col=1))
        print("   lensed TT    ", cl_rel(x, y, key="cl_total"))
        print("   phiphi       ", cl_rel(x, y, key="cl_lens_potential"))
    # P(k) on the acc1 grid (k >= 1e-4 to avoid the low-k region where acc2 blows up)
    k = a0["pk_k"]
    m = k >= 1e-4
    p_a0, p_cf, p_a2, p_cf2 = a0["pk"][0], cf["pk"][0], pk_on(a2, k), pk_on(cf2, k)
    print("\nP(k, z=0) at selected k (acc1 grid); rel. diffs w.r.t. convfix_acc2 (best available):")
    print(f"{'k[1/Mpc]':>10s} {'a0-cf':>10s} {'a0-best':>10s} {'cf-best':>10s} {'a0acc2-best':>11s}")
    for kk in [1e-3, 1e-2, 0.05, 0.1, 0.3, 1.0, 2.0, 5.0, 10.0, 15.0, 20.0, k[-1]]:
        i = int(np.argmin(np.abs(k - kk)))
        print(f"{k[i]:10.4g} {p_a0[i] / p_cf[i] - 1:+10.2e} {p_a0[i] / p_cf2[i] - 1:+10.2e} {p_cf[i] / p_cf2[i] - 1:+10.2e}"
              f" {p_a2[i] / p_cf2[i] - 1:+11.2e}")
    for lo, hi in [(1e-4, 0.1), (0.1, 1.4), (1.4, 10), (10, 30)]:
        s = m & (k >= lo) & (k < hi)
        print(f"k in [{lo},{hi}): max|a0/cf-1| = {np.max(np.abs(p_a0[s] / p_cf[s] - 1)):.2e}   "
              f"max|a0/best-1| = {np.max(np.abs(p_a0[s] / p_cf2[s] - 1)):.2e}   "
              f"max|cf/best-1| = {np.max(np.abs(p_cf[s] / p_cf2[s] - 1)):.2e}")
    print(f"sigma8: a0 {float(a0['sigma8'][0]):.8f} convfix {float(cf['sigma8'][0]):.8f} "
          f"a0_acc2 {float(a2['sigma8'][0]):.8f} convfix_acc2 {float(cf2['sigma8'][0]):.8f}")


if __name__ == "__main__":
    main()
