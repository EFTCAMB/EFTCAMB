"""Verify t3-pi-ic-singular-denominator: (1) same P(k) spikes in pre-rebase B1; (2) footprint in CMB C_l on A0.

Runs through harness.run_model (threads=1). Output dir: work/verify/piic/.
"""

import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import harness  # noqa: E402

OUT = harness.AUDIT_ROOT / "work" / "verify" / "piic"
INI = harness.PARAMS_DIR / "3_RPH_MK_const_1.ini"
PK = {"get_scalar_cls": "F", "get_tensor_cls": "F", "do_lensing": "F", "transfer_kmax": "10",
      "num_redshiftwindows": "0"}


def run(args):
    tree, name, ov = args
    out = OUT / f"{tree}_{name}.npz"
    return name, tree, harness.run_model(tree, INI, out, overrides=ov, threads=1, timeout=900), out


def main() -> None:
    jobs = [("B1", f"pk_m{v}", {**PK, "RPHmassP0": v}) for v in ("0.0", "1e-06", "1e-04")]
    jobs += [("A0", f"pk_m{v}", {**PK, "RPHmassP0": v}) for v in ("0.0", "1e-06", "1e-04")]
    jobs += [("A0", f"cl_m{v}", {"RPHmassP0": v}) for v in ("0.0", "1e-04", "1e-06")]
    with ThreadPoolExecutor(4) as ex:
        res = {(t, n): (st, out) for n, t, st, out in ex.map(run, jobs)}
    for key, (st, _) in res.items():
        if st["status"] != "ok":
            print(key, st["status"], st.get("log_tail", [])[-3:])
    print("P(k) z=0 vs same model with RPHmassP0=0 (alpha_K=0.1 kept):")
    for tree in ("B1", "A0"):
        ref = harness.load(res[(tree, "pk_m0.0")][1])
        for v in ("1e-06", "1e-04"):
            d = harness.load(res[(tree, f"pk_m{v}")][1])
            r = np.abs(d["pk"][0] / ref["pk"][0] - 1)
            i = int(np.argmax(r))
            nb = r[max(i - 3, 0):i + 4]
            print(f"  {tree} alpha_M={v}: max {r[i]:.3e} at k={d['pk_k'][i]:.4g}; neighbours {np.array2string(nb, precision=1)}"
                  f"; median {np.median(r):.1e} fallback={d['meta'].get('eft_init_fallback', False)}")
    print("C_l (lensed total) vs RPHmassP0=0, A0 default test-ini settings:")
    ref = harness.load(res[("A0", "cl_m0.0")][1])["cl_total"]
    for v in ("1e-04", "1e-06"):
        d = harness.load(res[("A0", f"cl_m{v}")][1])["cl_total"]
        for j, s in enumerate(("TT", "EE")):
            r = np.abs(d[2:2500, j] / ref[2:2500, j] - 1)
            i = int(np.argmax(r))
            print(f"  alpha_M={v} {s}: max {r[i]:.3e} at ell={i + 2}; median {np.median(r):.1e}; "
                  f"ell 2-50 max {r[:49].max():.1e}")
        ur = harness.load(res[("A0", f"cl_m{v}")][1])["cl_unlensed_scalar"]
        uref = harness.load(res[("A0", "cl_m0.0")][1])["cl_unlensed_scalar"]
        r = np.abs(ur[2:2500, 0] / uref[2:2500, 0] - 1)
        i = int(np.argmax(r))
        prof = ", ".join(f"{ell}:{r[ell - 2]:.1e}" for ell in (50, 100, 150, 200, 250, 300, 400, 600, 1000, 2000))
        print(f"  alpha_M={v} unlensed TT: max {r[i]:.3e} at ell={i + 2}; profile {prof}")


if __name__ == "__main__":
    main()
