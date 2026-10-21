"""Repro (t5-robust): with EFTCAMB_evolve_metric_h = T (Horndeski, delta_phi formulation) the initial h' is the GR
value 2 k z, which violates the Horndeski 00 constraint at the turn-on; the matter power spectrum is suppressed by
up to 13% at k = 10/Mpc (42% at k = 20) and the result depends strongly on EFTCAMB_turn_on_time.

09_EFTCAMB_IC.f90:290   if (CP%EFTCAMB%EFTCAMB_evolve_metric_h) y(EV%w_ix + 2) = 2._dl*k*z
with z = (0.5 dgrho/k + etak)/adotoa (GR, 09_EFTCAMB_IC.f90:199-200).  With evolve_metric_h = F the code instead uses
the constraint (equations.f90:3147)
    h' = (4 k etak + 2 dgrho/Meff2)/((2 - alphaB) adotoa) + s00p dphi' + (s00 + s00k k^2) dphi,
which differs from 2kz whenever Meff2 != 1 or alphaB != 0 at the turn-on (5_hdsk_TG_1: Meff2 = 0.9 at a = 1e-8).
The flag is the documented remedy printed by the Horndeski module when alphaB crosses 2.  (Same in pre-rebase B1.)

Runs (frozen A0): 5_hdsk_TG_1 with evolve_metric_h = F / T at the test-ini turn-on (1e-8) and at 1e-7.
If work/t5-robust/A0_icsw exists, also evolve_metric_h = T with EFT_T5_HIC=1 (h'_ini from the constraint above,
delta_phi_ini = 0; nothing else changed).

Run: /Users/marco/Software/miniforge3/bin/python metric-h-ic-violates-constraint.py   (4-6 CAMB runs, ~1.5 min)
"""

import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(REPO / "fortran" / "eftcamb_test" / "audit" / "checks"))
import harness  # noqa: E402

ICSW = harness.AUDIT_ROOT / "work" / "t5-robust" / "A0_icsw"
OUT = harness.AUDIT_ROOT / "work" / "t5-robust" / "repro_metric_h"
KS = (0.1, 1.0, 5.0, 10.0, 20.0)


def job(args):
    name, tree, ov, env = args
    out = OUT / f"{name}.npz"
    st = harness.run_model(tree, harness.PARAMS_DIR / "5_hdsk_TG_1.ini", out, overrides=ov, threads=1, timeout=600,
                           extra_env=env)
    if st["status"] != "ok":
        raise RuntimeError(f"{name}: {st['status']}")
    return name, harness.load(out)


def pk_at(d, k):
    return float(np.exp(np.interp(np.log(k), np.log(d["pk_k"]), np.log(d["pk"][0]))))


def main():
    mh = {"EFTCAMB_evolve_metric_h": "T"}
    jobs = [("h_constraint_on1e-8", "A0", {}, {}), ("h_evolved_on1e-8", "A0", mh, {}),
            ("h_constraint_on1e-7", "A0", {"EFTCAMB_turn_on_time": "1e-7"}, {}),
            ("h_evolved_on1e-7", "A0", {**mh, "EFTCAMB_turn_on_time": "1e-7"}, {})]
    fix = (ICSW / "camb" / "camblib.so").exists()
    if fix:
        jobs.append(("h_evolved_constraintIC_on1e-8", str(ICSW), mh, {"EFT_T5_HIC": "1"}))
    with ThreadPoolExecutor(3) as ex:
        res = dict(ex.map(job, jobs))
    print("5_hdsk_TG_1: P(k; evolve_metric_h=T)/P(k; F) - 1 at k = " + ", ".join(f"{k:g}" for k in KS) + " /Mpc, and TT")
    rows = [("turn-on 1e-8 (test ini)", "h_evolved_on1e-8", "h_constraint_on1e-8"),
            ("turn-on 1e-7", "h_evolved_on1e-7", "h_constraint_on1e-7")]
    if fix:
        rows.append(("turn-on 1e-8, h'_ini from constraint (A0_icsw)", "h_evolved_constraintIC_on1e-8",
                     "h_constraint_on1e-8"))
    for label, x, y in rows:
        a, b = res[x], res[y]
        tt = a["cl_total"][2:2501, 0] / b["cl_total"][2:2501, 0] - 1
        print(f"  {label:48s} " + "  ".join(f"{pk_at(a, k) / pk_at(b, k) - 1:+.2e}" for k in KS)
              + f" | dTT/TT ell=4: {tt[2]:+.2e}, max ell>=30: {np.abs(tt[28:]).max():.1e}")
    print("expected: evolving h with its second-order equation reproduces the constraint solution (< 1e-4)")


if __name__ == "__main__":
    main()
