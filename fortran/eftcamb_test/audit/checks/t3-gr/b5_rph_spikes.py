"""t3-gr (b): isolated single-k P(k) glitches of RPH near GR (alpha_K = 0.1 const, alpha_M = RPHmassP0 tiny).

RPH__s1em06 (RPHmassP0=1e-6) deviates from LCDM by 2.6e-3 at k=0.173/Mpc only (neighbours ~1e-4), s1em07 by 6.8e-4
at k=0.544, s1em09 by 5.8e-4 at k=5.39, while RPHmassP0=0 gives a smooth 1e-5 floor. Check: pre-rebase EFTCAMB (B1),
and A0 with tighter integration (IntTolBoost=4) at RPHmassP0=1e-6, 1e-7.
"""

import sys
from concurrent.futures import ThreadPoolExecutor

import numpy as np

from common import PARAMS, RUNS, harness

INI = PARAMS / "3_RPH_MK_const_1.ini"
SUB = RUNS / "b5"
REF = RUNS / "b1" / "REF_base__A0g.npz"
JOBS = {
    "B1_m1e-6": ("B1", {"RPHmassP0": "1e-06"}, {}),
    "B1_m1e-7": ("B1", {"RPHmassP0": "1e-07"}, {}),
    "B1_m0": ("B1", {"RPHmassP0": "0.0"}, {}),
    "A0_m1e-6_tol4": ("A0", {"RPHmassP0": "1e-06"}, {"Accuracy.IntTolBoost": 4.0}),
    "A0_m1e-7_tol4": ("A0", {"RPHmassP0": "1e-07"}, {"Accuracy.IntTolBoost": 4.0}),
    "A0_m1e-6_acc2": ("A0", {"RPHmassP0": "1e-06", "accuracy_boost": "2"}, {}),
}


def one(item):
    name, (tree, ov, ps) = item
    out = SUB / f"{name}.npz"
    if not out.exists():
        harness.run_model(tree, INI, out, overrides=ov, python_set=ps, threads=1, timeout=1500)
    return name, out


def main() -> None:
    ref = harness.load(REF)
    refs = {"A0": ref, "B1": None}
    with ThreadPoolExecutor(int(sys.argv[1]) if len(sys.argv) > 1 else 3) as ex:
        res = dict(ex.map(one, JOBS.items()))
    b1ref = harness.load(res["B1_m0"]) if res["B1_m0"].exists() else None
    for name, out in res.items():
        if not out.exists():
            print(name, "FAILED")
            continue
        d = harness.load(out)
        r0 = b1ref if name.startswith("B1") else ref
        if r0 is None:
            continue
        if d["pk_k"].shape != r0["pk_k"].shape or not np.allclose(d["pk_k"], r0["pk_k"]):
            # different k grid (accuracy_boost): log-log interpolate the reference (cubic) onto this grid
            from scipy.interpolate import CubicSpline
            lk = np.log(d["pk_k"])
            m = (lk > np.log(r0["pk_k"][1])) & (lk < np.log(r0["pk_k"][-2]))
            ref_i = np.exp(CubicSpline(np.log(r0["pk_k"]), np.log(r0["pk"][0]))(lk[m]))
            d = {"pk_k": d["pk_k"][m], "pk": d["pk"][:, m]}
            r = np.abs(d["pk"][0] / ref_i - 1)
        else:
            r = np.abs(d["pk"][0] / r0["pk"][0] - 1)
        # spikiness: deviation of a k point from the median of its 6 neighbours
        nb = np.array([np.median(np.r_[r[max(0, i - 3):i], r[i + 1:i + 4]]) for i in range(len(r))])
        sp = r - nb
        j = int(np.argmax(sp))
        print(f"   {name}: largest isolated spike {sp[j]:.1e} at k={d['pk_k'][j]:.3g} (neighbour median {nb[j]:.1e})")
        top = np.argsort(r)[::-1][:4]
        print(f"{name:16s} vs {'B1 m=0' if name.startswith('B1') else 'A0g'}: max {r.max():.2e} at k={d['pk_k'][top[0]]:.3g};"
              f" top: " + ", ".join(f"{d['pk_k'][j]:.3g}:{r[j]:.1e}" for j in top) + f"; median {np.median(r):.1e}")


if __name__ == "__main__":
    main()
