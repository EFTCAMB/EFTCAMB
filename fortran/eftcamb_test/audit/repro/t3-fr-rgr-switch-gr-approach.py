"""Repro (t3-gr): designer f(R) approaching GR: the return-to-GR switch, not the physics, sends the small-scale
P(k) deviation to zero.

09_EFTCAMB_RGR.f90 turns the EFT perturbations on at the first a where an EFT function exceeds
EFTCAMB_GR_threshold (1e-8). For f(R) (4_Des_F_R_LCDM_1, wDE = -1) the enhancement of P(k) below the Compton scale
(k_C ~ H0/sqrt(B0) ~ 1/Mpc for B0 ~ 5e-8) is O(10%) whatever the size of Omega = f_R, but the RGR time moves from
a = 0.73 (B0 = 1.2e-7) to > 1 (B0 = 5.2e-8): P(k=7/Mpc) goes 12% -> 1.4% -> exactly LCDM (7e-6) over a factor 2.3 in
B0, i.e. the result is set by when the threshold is crossed. Running the same point with the switch disabled
(EFTCAMB_GR_threshold = -1) is not possible in practice (the pi equation is stiff, m ~ H/sqrt(B0); >480 s timeout
for P(k) only), so there is no cheap way for the user to get the right answer. See also fr-rgr-threshold-high-k-pk.py
(same mechanism at B0 = 1e-4).
Expected: P(k) deviation varying smoothly (no RGR-driven collapse) as B0 -> 0. Observed: printed.
Run: /Users/marco/Software/miniforge3/bin/python t3-fr-rgr-switch-gr-approach.py   (5 fast P(k)-only runs, <1 min)
"""

import re
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(REPO / "fortran" / "eftcamb_test" / "audit" / "checks"))
import harness  # noqa: E402

INI = harness.PARAMS_DIR / "4_Des_F_R_LCDM_1.ini"
GR = harness.PARAMS_DIR / "1_EFT_GR.ini"
OUT = harness.AUDIT_ROOT / "work" / "t3-gr" / "repro_fr_rgr"
PK = {"get_scalar_cls": "F", "get_tensor_cls": "F", "do_lensing": "F", "transfer_kmax": "10",
      "num_redshiftwindows": "0"}
VALS = [1.2e-7, 7.9e-8, 6.0e-8, 5.2e-8]


def job(args):
    name, ini, ov, ps = args
    out = OUT / f"{name}.npz"
    st = harness.run_model("A0", ini, out, overrides={**PK, **ov}, python_set=ps, threads=1, timeout=300)
    return name, st, out


def main() -> None:
    jobs = [("GR", GR, {}, {"Accuracy.IntTolBoost": 2.0})] + [(f"fR_{v:.1e}", INI, {"EFTB0": repr(v)}, {}) for v in VALS]
    with ThreadPoolExecutor(3) as ex:
        res = {n: (st, out) for n, st, out in ex.map(job, jobs)}
    ref = harness.load(res["GR"][1])
    k = ref["pk_k"]
    ik = [int(np.argmin(abs(k - x))) for x in (0.3, 1.0, 3.0, 7.0)]
    print(f"{'B0':>9s} {'RGR a':>7s} " + " ".join(f"{'P(k=%.2g)-1' % k[i]:>12s}" for i in ik))
    for v in VALS:
        st, out = res[f"fR_{v:.1e}"]
        if st["status"] != "ok":
            print(f"{v:9.1e} {st['status']}")
            continue
        rgr = re.findall(r"Return to GR time: *([-+0-9.E]+)", out.with_suffix(".log").read_text())
        r = harness.load(out)["pk"][0] / ref["pk"][0] - 1
        print(f"{v:9.1e} {float(rgr[-1]):7.3f} " + " ".join(f"{r[i]:+12.3e}" for i in ik))
    print("expected: smooth dependence on B0; observed: deviation collapses to 0 as the RGR time crosses a = 1")


if __name__ == "__main__":
    main()
