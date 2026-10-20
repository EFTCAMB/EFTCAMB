"""Verify t3-fr-rgr-switch-gr-approach: designer f(R) P(k) vs B0 at several EFTCAMB_GR_threshold values.

Usage: python frrgr_threshold_scan.py TREE B0[,B0...] THR[,THR...] [timeout_s] [jobs]
Prints, per (B0, threshold): status, elapsed, RGR time (from log), P(k)/P_GR - 1 at k = 0.3, 1, 3, 7 /Mpc.
The GR reference is 1_EFT_GR.ini at the same accuracy settings as the f(R) runs (no IntTolBoost change).
"""

import re
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness  # noqa: E402

OUT = harness.AUDIT_ROOT / "work" / "verify" / "frrgr"
INI = harness.PARAMS_DIR / "4_Des_F_R_LCDM_1.ini"
GR = harness.PARAMS_DIR / "1_EFT_GR.ini"
PK = {"get_scalar_cls": "F", "get_tensor_cls": "F", "do_lensing": "F", "transfer_kmax": "10",
      "num_redshiftwindows": "0"}
KS = (0.3, 1.0, 3.0, 7.0)


def main() -> None:
    tree = sys.argv[1]
    b0s = [float(x) for x in sys.argv[2].split(",")]
    thrs = [float(x) for x in sys.argv[3].split(",")]
    timeout = float(sys.argv[4]) if len(sys.argv) > 4 else 300
    jobs = int(sys.argv[5]) if len(sys.argv) > 5 else 3
    tag = Path(tree).name
    ref_out = OUT / f"{tag}_GR.npz"
    if not ref_out.exists():
        st = harness.run_model(tree, GR, ref_out, overrides=PK, threads=1, timeout=300)
        assert st["status"] == "ok", st
    ref = harness.load(ref_out)
    k = ref["pk_k"]
    ik = [int(np.argmin(abs(k - x))) for x in KS]

    def job(args):
        b0, thr = args
        out = OUT / f"{tag}_B0{b0:.2e}_thr{thr:.0e}.npz"
        st = harness.run_model(tree, INI, out, overrides={**PK, "EFTB0": repr(b0), "EFTCAMB_GR_threshold": repr(thr)},
                               threads=1, timeout=timeout)
        return b0, thr, st, out

    combos = [(b, t) for b in b0s for t in thrs]
    with ThreadPoolExecutor(jobs) as ex:
        res = list(ex.map(job, combos))
    print(f"{'B0':>9s} {'thr':>7s} {'status':>8s} {'t[s]':>6s} {'RGR a':>8s} "
          + " ".join(f"{'k=%.2g' % k[i]:>11s}" for i in ik) + "  sigma8")
    for b0, thr, st, out in res:
        if st["status"] != "ok":
            print(f"{b0:9.2e} {thr:7.0e} {st['status']:>8s} {st['elapsed']:6.1f}")
            continue
        rgr = re.findall(r"Return to GR time: *([-+0-9.E]+)", out.with_suffix(".log").read_text())
        d = harness.load(out)
        assert np.allclose(d["pk_k"], k)
        r = d["pk"][0] / ref["pk"][0] - 1
        print(f"{b0:9.2e} {thr:7.0e} {st['status']:>8s} {st['elapsed']:6.1f} {float(rgr[-1]) if rgr else -1:8.4f} "
              + " ".join(f"{r[i]:+11.3e}" for i in ik) + f"  {d['sigma8'][0]:.6f}")
    print(f"GR sigma8 {ref['sigma8'][0]:.6f}")


if __name__ == "__main__":
    main()
