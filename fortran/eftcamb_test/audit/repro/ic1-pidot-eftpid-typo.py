"""Repro (t5-robust): EFTpiDfunction returns D1 + k^2 D1 instead of D1 + k^2 D2, so the EFT_IC_type = 1 initial pi'
is wrong at high k; for the K-mimic test models P(k) at k >= 1/Mpc is wrong by 2% .. >100%.

09_EFTCAMB_IC.f90:589   EFTpiDfunction = eft_cache%EFTpiD1 + k*k*eft_cache%EFTpiD1      (should be ... k*k*eft_cache%EFTpiD2)
It feeds EFTpiDdotFunction, used at 09_EFTCAMB_IC.f90:301 in
    pi'_ini = h0 [ E (C' + k^2 D') aH / (C + k^2 D1 + k^4 D2)^2 - E'/(C + k^2 D1 + k^4 D2) ],
so the bracket contains C' + k^2 D1' + k^4 D1' instead of C' + k^2 D1' + k^4 D2': an error growing like k^4 D1'.
(Same line in the pre-rebase EFTCAMB, B1.)

The script runs 5_Kmimic_1 and 5_Kmimic_2 (default test inis) on
  A0  = frozen release build of the code under audit, and
  FIX = work/t5-robust/A0_pidfix = A0 with only that line corrected (release build),
and 1_EFT_GR for reference, through the audit harness (machine-wide CPU slots). Prints P(k)/P_GR(k) - 1 and the change
of the lensing potential C_L^phiphi and lensed TT.  Both builds are converged in accuracy_boost (checked separately:
accuracy_boost = 2 changes these numbers by < 4e-4), so the difference is the initial condition, not noise.

Run: /Users/marco/Software/miniforge3/bin/python ic1-pidot-eftpid-typo.py      (6 short CAMB runs, ~1 min)
"""

import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(REPO / "fortran" / "eftcamb_test" / "audit" / "checks"))
import harness  # noqa: E402

FIX = harness.AUDIT_ROOT / "work" / "t5-robust" / "A0_pidfix"
OUT = harness.AUDIT_ROOT / "work" / "t5-robust" / "repro_ic1_pidot_typo"
KS = (1.0, 2.0, 5.0, 10.0, 20.0)


def job(args):
    name, tree, model = args
    out = OUT / f"{name}.npz"
    st = harness.run_model(tree, harness.PARAMS_DIR / f"{model}.ini", out, threads=1, timeout=600)
    if st["status"] != "ok":
        raise RuntimeError(f"{name}: {st['status']} {st.get('log_tail', [])[-3:]}")
    return name, harness.load(out)


def pk_at(d, k):
    return float(np.exp(np.interp(np.log(k), np.log(d["pk_k"]), np.log(d["pk"][0]))))


def main():
    jobs = [("GR", "A0", "1_EFT_GR")]
    for m in ("5_Kmimic_1", "5_Kmimic_2"):
        jobs += [(f"A0_{m}", "A0", m), (f"FIX_{m}", str(FIX), m)]
    with ThreadPoolExecutor(3) as ex:
        res = dict(ex.map(job, jobs))
    gr = res["GR"]
    for m in ("5_Kmimic_1", "5_Kmimic_2"):
        a, f = res[f"A0_{m}"], res[f"FIX_{m}"]
        print(f"{m}: P(k)/P_GR(k) - 1 at k = " + ", ".join(f"{k:g}" for k in KS) + " /Mpc (z = 0)")
        print("   observed (A0, typo)    " + "  ".join(f"{pk_at(a, k) / pk_at(gr, k) - 1:+.3e}" for k in KS))
        print("   expected (D2 restored) " + "  ".join(f"{pk_at(f, k) / pk_at(gr, k) - 1:+.3e}" for k in KS))
        pp_a, pp_f = a["cl_lens_potential"][:, 0], f["cl_lens_potential"][:, 0]
        tt_a, tt_f = a["cl_total"][:, 0], f["cl_total"][:, 0]
        print("   C_L^phiphi(A0)/C_L^phiphi(fixed) - 1 at L = 500, 1000, 2000, 2500: "
              + "  ".join(f"{pp_a[L] / pp_f[L] - 1:+.2e}" for L in (500, 1000, 2000, 2500)))
        print("   lensed TT(A0)/TT(fixed) - 1 at ell = 1000, 2000, 2500:            "
              + "  ".join(f"{tt_a[L] / tt_f[L] - 1:+.2e}" for L in (1000, 2000, 2500)))
    print("expected: A0 equal to the 'D2 restored' row (the corrected line only fixes the k^4 term of pi'_ini); "
          "observed: the typo changes P(k >= 2/Mpc) by tens of percent and C_L^phiphi by ~1e-3")


if __name__ == "__main__":
    main()
