"""t3-gr (a): EFTflag=0 transfer-only runs with small transfer_kmax (cmbmain q_switch_lowk 'small fix' in EFTCAMB).

B0 (plain CAMB v2) stops with 'thermo out of bounds' for transfer_kmax=0.01 when no Cls are requested; A0 runs.
Scan kmax to find where B0 fails and check that A0's P(k) at small kmax agrees with A0/B0 runs at larger kmax.
"""

import sys
from concurrent.futures import ThreadPoolExecutor

import numpy as np

from common import harness, run

SUB = "a10"
BASE = {"do_lensing": "F", "get_tensor_cls": "F", "get_scalar_cls": "F", "num_redshiftwindows": "0"}
KMAX = ["0.005", "0.01", "0.02", "0.03", "0.05", "0.2"]


def one(job):
    tree, km = job
    out = run(tree, f"{tree}_k{km}", overrides={**BASE, "transfer_kmax": km}, python_set={"Reion.max_redshift": 50.0},
              sub=SUB)
    return job, out


def main() -> None:
    jobs = [(t, k) for k in KMAX for t in ("A0", "B0")]
    res = {}
    with ThreadPoolExecutor(int(sys.argv[1]) if len(sys.argv) > 1 else 4) as ex:
        for job, out in ex.map(one, jobs):
            res[job] = out
    ref = harness.load(res[("B0", "0.2")])
    for k in KMAX:
        for t in ("A0", "B0"):
            p = res[(t, k)]
            if not p.exists():
                print(f"{t} kmax={k:6s} FAILED")
                continue
            d = harness.load(p)
            kk = d["pk_k"]
            m = (kk > ref["pk_k"][0]) & (kk < float(k))
            pr = np.exp(np.interp(np.log(kk[m]), np.log(ref["pk_k"]), np.log(ref["pk"][0])))
            dev = np.max(np.abs(d["pk"][0][m] / pr - 1)) if m.any() else float("nan")
            print(f"{t} kmax={k:6s} ok  nk={len(kk):4d} k_last={kk[-1]:.3e}  max|P/P_B0(kmax=0.2)-1| (k<kmax)={dev:.2e} "
                  f"sigma8={float(d['sigma8'][0]) if d['sigma8'].size else float('nan'):.6f}")
        pa, pb = res[("A0", k)], res[("B0", k)]
        if pa.exists() and pb.exists():
            r = harness.compare(pa, pb)
            print(f"   A0 vs B0 kmax={k}: bit={r['bit_identical']} worst_pk={r['worst_pk']:.2e}")


if __name__ == "__main__":
    main()
