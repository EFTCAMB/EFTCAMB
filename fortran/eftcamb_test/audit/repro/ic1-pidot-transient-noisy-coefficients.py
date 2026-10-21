"""Repro (t5-robust): the EFT_IC_type = 1 initial pi' (finite-difference derivative of the pi-equation coefficients)
puts the pi field far off its attractor for full-mapping models at an early turn-on, so the results depend on
EFTCAMB_turn_on_time and can blow up silently.

09_EFTCAMB_IC.f90:296-302 sets pi_ini = -h0 E/(C + k^2 D1 + k^4 D2) and
pi'_ini = h0 [E (C' + k^2 D') aH/(...)^2 - E'/(...)], with C' and D' from dfridr (Ridders extrapolation, h = 0.03 a)
of EFTpiCfunction / EFTpiDfunction.  For K-mimic (5_Kmimic_1, Gamma1 ~ 1e14..1e26 at a <= 1e-6) and for the
Horndeski TG model at a = 1e-8 the coefficients vary by orders of magnitude between neighbouring a (printed by
get_scale_evolution: C + k^2 D1 = -1.2e8, -4.0e9, -1.6e9 at a = 1.00, 1.12, 1.26 e-9), the "attractor" and its
derivative are meaningless, and pi' starts O(10^2) and never decays.

Part 1 (A0): 5_Kmimic_1 P(k)/P_GR - 1 for EFTCAMB_turn_on_time = 1e-9, 1e-8 (test-ini default), 1e-7.
Part 2 (A0): 5_hdsk_TG_1 (EFTCAMB_evolve_delta_phi = T, delta_phi IC = 0) vs 5_hdsk_TG_2 (same model, pi variable,
             IC type 1): the two formulations must agree; they differ by 2e-4 in lensed TT/EE at ell ~ 125
             (pre-rebase B1: 7e-7 on the same grid).
If the audit tree work/t5-robust/A0_icsw exists, the same runs are repeated there with EFT_T5_PIDOT0=1 (pi'_ini = 0,
nothing else changed; A0_icsw also contains the EFTpiD2 typo fix) to show the turn-on independence / agreement that
the IC should give.

Run: /Users/marco/Software/miniforge3/bin/python ic1-pidot-transient-noisy-coefficients.py   (~10 CAMB runs, 2-3 min)
"""

import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(REPO / "fortran" / "eftcamb_test" / "audit" / "checks"))
import harness  # noqa: E402

ICSW = harness.AUDIT_ROOT / "work" / "t5-robust" / "A0_icsw"
OUT = harness.AUDIT_ROOT / "work" / "t5-robust" / "repro_ic1_pidot_transient"
KS = (1.0, 5.0, 10.0)
TURN_ON = ("1e-9", "1e-8", "1e-7")


def job(args):
    name, tree, model, ov, env = args
    out = OUT / f"{name}.npz"
    st = harness.run_model(tree, harness.PARAMS_DIR / f"{model}.ini", out, overrides=ov, threads=1, timeout=600,
                           extra_env=env)
    return name, (harness.load(out) if st["status"] == "ok" else st["status"])


def pk_at(d, k):
    return float(np.exp(np.interp(np.log(k), np.log(d["pk_k"]), np.log(d["pk"][0]))))


def cl_diff(a, b):
    out = []
    for j in (0, 1):
        x, y = a["cl_total"][2:2501, j], b["cl_total"][2:2501, j]
        r = np.abs(x / y - 1)
        out.append((float(r.max()), int(np.argmax(r)) + 2))
    return out


def main():
    trees = [("A0", "A0", {})]
    if (ICSW / "camb" / "camblib.so").exists():
        trees.append(("pidot0", str(ICSW), {"EFT_T5_PIDOT0": "1"}))
    jobs = [("GR", "A0", "1_EFT_GR", {}, {})]
    for tag, tree, env in trees:
        for t in TURN_ON:
            jobs.append((f"{tag}_Kmimic1_on{t}", tree, "5_Kmimic_1", {"EFTCAMB_turn_on_time": t}, env))
        jobs.append((f"{tag}_TG1", tree, "5_hdsk_TG_1", {}, env))
        jobs.append((f"{tag}_TG2", tree, "5_hdsk_TG_2", {}, env))
    with ThreadPoolExecutor(4) as ex:
        res = dict(ex.map(job, jobs))
    gr = res["GR"]
    for tag, _, _ in trees:
        label = "observed (A0, release IC)" if tag == "A0" else "reference (A0_icsw, pi'_ini = 0)"
        print(f"== {label}")
        print("  5_Kmimic_1 P(k)/P_GR - 1 at k = " + ", ".join(f"{k:g}" for k in KS) + " /Mpc:")
        for t in TURN_ON:
            d = res[f"{tag}_Kmimic1_on{t}"]
            if isinstance(d, str):
                print(f"    turn-on {t}: {d}")
                continue
            print(f"    turn-on {t}: " + "  ".join(f"{pk_at(d, k) / pk_at(gr, k) - 1:+.3e}" for k in KS))
        a, b = res[f"{tag}_TG1"], res[f"{tag}_TG2"]
        if not isinstance(a, str) and not isinstance(b, str):
            (tt, ltt), (ee, lee) = cl_diff(b, a)
            print(f"  5_hdsk_TG_2 (pi) vs 5_hdsk_TG_1 (delta_phi): max|dTT/TT| = {tt:.2e} (ell {ltt}), "
                  f"max|dEE/EE| = {ee:.2e} (ell {lee})")
    print("expected: results independent of the turn-on time (the model is integrated from its attractor), no "
          "divergence, and the pi / delta_phi formulations agreeing to << 1e-4")


if __name__ == "__main__":
    main()
