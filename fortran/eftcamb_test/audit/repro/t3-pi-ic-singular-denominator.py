"""Repro (t3-gr): the standard pi initial condition is singular where C + k^2 D1 + k^4 D2 crosses zero, giving isolated
single-k spikes in P(k) (and C_l) that do not vanish as the model approaches GR.

09_EFTCAMB_IC.f90 (EFT_IC_type = 1, lines ~295-305) sets at the perturbation turn-on
    pi = -h0 E/(C + k^2 D1 + k^4 D2),   pi' = h0 (E (C' + k^2 D')/(...)^2 a H - E'/(...))
guarded only by `(...) /= 0`. If C and D1 have opposite signs at the turn-on there is a k* = sqrt(-C/D1) where the
denominator vanishes and pi_ini blows up for the k modes next to it. Model: 3_RPH_MK_const_1 (alpha_K = 0.1,
alpha_M = RPHmassP0, LCDM background, turn-on a = 0.01). At a = 0.0103, C = -6.8e-14 (independent of alpha_M) and
D1 ~ alpha_M, so k* ~ alpha_M^(-1/2): the spike moves to higher k but does not shrink as alpha_M -> 0; its size
depends on how close a k-grid point falls to k* (full test-ini settings: 2.6e-3 at k=0.173 for 1e-6, 6.8e-4 at
k=0.544 for 1e-7, 5.8e-4 at k=5.39 for 1e-9; RPHmassP0 = 0 exactly: D1 = 0, no crossing, smooth 2e-5 agreement).
Same numbers in the pre-rebase B1, not changed by IntTolBoost=4, larger with accuracy_boost=2 (1.9e-2).
(Separately, EFTpiDfunction in the same file returns D1 + k^2 D1 instead of D1 + k^2 D2.)
Expected: |P/P_GR - 1| smooth in k and -> 0 as RPHmassP0 -> 0. Observed: printed.
Run: /Users/marco/Software/miniforge3/bin/python t3-pi-ic-singular-denominator.py   (6 fast CAMB runs, ~1 min)
"""

import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(REPO / "fortran" / "eftcamb_test" / "audit" / "checks"))
import harness  # noqa: E402

INI = harness.PARAMS_DIR / "3_RPH_MK_const_1.ini"
GR = harness.PARAMS_DIR / "1_EFT_GR.ini"
OUT = harness.AUDIT_ROOT / "work" / "t3-gr" / "repro_pi_ic"
PK = {"get_scalar_cls": "F", "get_tensor_cls": "F", "do_lensing": "F", "transfer_kmax": "10",
      "num_redshiftwindows": "0"}
VALS = ["0.0", "1e-09", "1e-07", "1e-06", "1e-04"]


def job(args):
    name, ini, ov, ps, evk = args
    out = OUT / f"{name}.npz"
    # get_scale_evolution needs WantScalars (upstream CAMB time evolution crashes in a transfer-only run)
    base = PK if evk is None else {"get_tensor_cls": "F", "do_lensing": "F", "l_max_scalar": "100",
                                   "num_redshiftwindows": "0"}
    st = harness.run_model("A0", ini, out, overrides={**base, **ov}, python_set=ps, evol_k=evk, threads=1, timeout=600)
    return name, st, out


def main() -> None:
    jobs = [("GR", GR, {}, {"Accuracy.IntTolBoost": 2.0}, None)]
    jobs += [(f"RPH_m{v}", INI, {"RPHmassP0": v}, {}, None) for v in VALS]
    with ThreadPoolExecutor(3) as ex:
        res = {n: (st, out) for n, st, out in ex.map(job, jobs)}
    ref = harness.load(res["GR"][1])
    k = ref["pk_k"]
    print(f"{'RPHmassP0':>10s} {'max|P/P_GR-1|':>14s} {'at k [1/Mpc]':>13s} {'median':>9s}")
    spike = {}
    for v in VALS:
        st, out = res[f"RPH_m{v}"]
        if st["status"] != "ok":
            print(f"{v:>10s} {st['status']}")
            continue
        r = np.abs(harness.load(out)["pk"][0] / ref["pk"][0] - 1)
        i = int(np.argmax(r))
        spike[v] = k[i]
        print(f"{v:>10s} {r[i]:14.2e} {k[i]:13.4g} {np.median(r):9.1e}")
    # IC denominator around the alpha_M = 1e-6 spike, at the first sampled time after the turn-on
    i0 = int(np.argmin(abs(k - spike.get("1e-06", 0.17))))
    ks = [float(x) for x in k[i0 - 2:i0 + 3]]
    _, st, out = job(("RPH_m1e-06_evol", INI, {"RPHmassP0": "1e-06"}, {}, ks))
    d = harness.load(out)
    a = np.atleast_2d(d["evol_a"])[0]
    j = int(np.searchsorted(a, 0.01))
    print(f"alpha_M = 1e-6: pi-equation coefficients at a = {a[j]:.4f} (turn-on 0.01)")
    for n, kv in enumerate(ks):
        C, D1, D2, E = (d[f"evol_{f}"][n, j] for f in ("EFTpiC", "EFTpiD1", "EFTpiD2", "EFTpiE"))
        den = C + kv**2 * D1 + kv**4 * D2
        print(f"   k={kv:.4f}  C={C:+.3e}  k^2 D1={kv * kv * D1:+.3e}  C+k^2D1+k^4D2={den:+.2e}  -E/den={-E / den:+.2e}")
    print("expected: no isolated spike; max deviation -> ~1e-5 (the RPHmassP0=0 value) as RPHmassP0 -> 0")


if __name__ == "__main__":
    main()
