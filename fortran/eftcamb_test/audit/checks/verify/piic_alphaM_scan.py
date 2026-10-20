"""Verify t3-pi-ic-singular-denominator: how often the IC pole corrupts C_l / P(k) across the near-GR region.

3_RPH_MK_const_1 (alpha_K=0.1) at default test-ini settings, RPHmassP0 log-spaced in [1e-5, 3e-3] (C<0 at turn-on,
so a pole k* exists), compared with RPHmassP0=0. Prints max |dCl/Cl| (unlensed TT, l<=2000), max |dP/P| (z=0) and the
location of the max.
"""

import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import harness  # noqa: E402

OUT = harness.AUDIT_ROOT / "work" / "verify" / "piic" / "scan"
INI = harness.PARAMS_DIR / "3_RPH_MK_const_1.ini"
VALS = ["0.0"] + [f"{x:.3g}" for x in np.logspace(-5, np.log10(3e-3), 14)]


def run(v):
    return v, harness.run_model("A0", INI, OUT / f"m{v}.npz", overrides={"RPHmassP0": v, "transfer_kmax": "5"},
                                threads=1, timeout=900)


def main() -> None:
    with ThreadPoolExecutor(4) as ex:
        res = dict(ex.map(run, VALS))
    ref = harness.load(OUT / "m0.0.npz")
    print(f"{'alpha_M':>9s} {'max dCl/Cl TT':>14s} {'at ell':>7s} {'max dP/P':>10s} {'at k':>8s} {'median dP/P':>12s}")
    for v in VALS[1:]:
        if res[v]["status"] != "ok":
            print(f"{v:>9s} {res[v]['status']}")
            continue
        d = harness.load(OUT / f"m{v}.npz")
        rc = np.abs(d["cl_unlensed_scalar"][2:2001, 0] / ref["cl_unlensed_scalar"][2:2001, 0] - 1)
        rp = np.abs(d["pk"][0] / ref["pk"][0] - 1)
        ic, ip = int(np.argmax(rc)), int(np.argmax(rp))
        print(f"{v:>9s} {rc[ic]:14.3e} {ic + 2:7d} {rp[ip]:10.3e} {d['pk_k'][ip]:8.4f} {np.median(rp):12.1e}")


if __name__ == "__main__":
    main()
