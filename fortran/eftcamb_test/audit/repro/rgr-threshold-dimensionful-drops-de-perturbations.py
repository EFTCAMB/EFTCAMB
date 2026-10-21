"""Repro (t5-robust): the return-to-GR test compares dimensionful EFT functions (in Mpc^-2 and Mpc^-3) with the
dimensionless EFTCAMB_GR_threshold; for minimally coupled dark energy (|c/a^2|, |Lambda/a^2 + grhov| ~ H0^2 ~ 5e-8
Mpc^-2) raising the threshold from 1e-8 to 1e-7 switches the dark-energy perturbations off for the whole history
(RGR time = 1.1), while the w(a) != -1 background is kept.

09_EFTCAMB_RGR.f90:153-170 scans a in [turn-on, 1] and turns the perturbations on at the first a where any of
|OmegaV|, |OmegaP|, ..., |EFTc/a^2|, |EFTLambda/a^2 + grhov|, |EFTcdot/a^2|, |EFTLambdadot/a^2|, |Gamma_i| exceeds the
threshold (09_EFTCAMB_RGR.f90:312-319); if none does it returns 1.1 and EFTCAMB_pert_turn_on is set to 1.1
(results.f90:602-604): no DE perturbations at all, no warning beyond the feedback line.  c and Lambda carry units
of H^2 (Mpc^-2), so the outcome depends on the unit system / H0, not on how GR-like the model is.

Runs (frozen A0): 4_Des_mc5e_CPL_1 (minimally coupled DE, CPL w) and 5_quint_exponential with the test-ini threshold
(1e-8) and 1e-7.  Prints the RGR time from the log, sigma8, and the CMB changes.

Run: /Users/marco/Software/miniforge3/bin/python rgr-threshold-dimensionful-drops-de-perturbations.py  (4 runs, ~1 min)
"""

import re
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(REPO / "fortran" / "eftcamb_test" / "audit" / "checks"))
import harness  # noqa: E402

OUT = harness.AUDIT_ROOT / "work" / "t5-robust" / "repro_rgr_dimensionful"
MODELS = ("4_Des_mc5e_CPL_1", "5_quint_exponential")
H0_MPC = 70.0 / 299792.458


def job(args):
    model, thr = args
    out = OUT / f"{model}_thr{thr}.npz"
    st = harness.run_model("A0", harness.PARAMS_DIR / f"{model}.ini", out, overrides={"EFTCAMB_GR_threshold": thr},
                           threads=1, timeout=600)
    if st["status"] != "ok":
        raise RuntimeError(f"{model} {thr}: {st['status']}")
    rgr = re.findall(r"Return to GR time:\s*(\S+)", out.with_suffix(".log").read_text(errors="replace"))
    return (model, thr), (harness.load(out), float(rgr[0]) if rgr else float("nan"))


def main():
    with ThreadPoolExecutor(4) as ex:
        res = dict(ex.map(job, [(m, t) for m in MODELS for t in ("1e-8", "1e-7")]))
    for m in MODELS:
        (d0, r0), (d1, r1) = res[(m, "1e-8")], res[(m, "1e-7")]
        a = d0["eft_a_grid"]
        c = np.abs(d0["eft_EFTc"] / a**2)
        print(f"{m}: max_a |EFTc/a^2| = {c.max():.2e} Mpc^-2 (H0^2 = {H0_MPC**2:.2e} Mpc^-2)")
        print(f"  threshold 1e-8: RGR time = {r0:.4g}   sigma8 = {d0['sigma8'][0]:.5f}")
        print(f"  threshold 1e-7: RGR time = {r1:.4g}   sigma8 = {d1['sigma8'][0]:.5f}  "
              f"(dsigma8/sigma8 = {d1['sigma8'][0] / d0['sigma8'][0] - 1:+.2e})")
        tt = d1["cl_total"][2:2501, 0] / d0["cl_total"][2:2501, 0] - 1
        pp = d1["cl_lens_potential"][2:2501, 0] / d0["cl_lens_potential"][2:2501, 0] - 1
        print(f"  TT change: ell=2 {tt[0]:+.2e}, max ell>=30 {np.abs(tt[28:]).max():.2e};  "
              f"C_L^phiphi change max {np.abs(pp).max():.2e}")
    print("expected: a 10x change of the GR threshold changes the observables of a model that is far from GR today "
          "(w != -1) by << 1e-3; observed: the DE perturbations are dropped for the whole history")


if __name__ == "__main__":
    main()
