"""Repro (t3-gr): EFTCAMB models ignore spatial curvature in H(a) used by the perturbations.

Pure EFT with all EFT functions zero and w=-1 is LCDM. With omk=0 it matches plain CAMB v2 (B0) to ~1e-5.
With omk=+-0.05 the same GR-equivalent model differs from B0 LCDM at the same cosmology by O(1-10%) in the CMB
Cls / P(k): compute_adotoa of the models builds adotoa from matter+radiation+DE only (no grhok), while CAMB's
background (dtauda, distances, Bessel functions) is curved. No check rejects omk/=0 with EFTflag/=0.

Runs four CAMB processes (each ~1 min at normal load) with the frozen audit trees A0 (code under audit) and B0 (CAMB).
Expected: |Delta| similar for omk=0 and omk/=0 (<~1e-4). Observed: see output.
"""

import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

REPO = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(REPO / "fortran" / "eftcamb_test" / "audit" / "checks"))
import harness  # noqa: E402

OUT = harness.AUDIT_ROOT / "work" / "t3-gr" / "repro_curv"
PEFT = harness.PARAMS_DIR / "2_PEFT_Omega_const_1.ini"
GR = harness.PARAMS_DIR / "1_EFT_GR.ini"
COMMON = {"num_redshiftwindows": "0", "get_tensor_cls": "F", "l_max_scalar": "2200"}


def job(args):
    tree, ini, ov, name = args
    return harness.run_model(tree, ini, OUT / f"{name}.npz", overrides=ov, threads=1)


def main() -> None:
    jobs = []
    for omk in ("0", "0.05", "-0.05"):
        ov = {**COMMON, "omk": omk}
        jobs.append(("A0", PEFT, {**ov, "PureEFTmodelOmega": "0"}, f"A0_PEFTzero_omk{omk}"))
        jobs.append(("B0", GR, ov, f"B0_LCDM_omk{omk}"))
    with ThreadPoolExecutor(4) as ex:
        st = list(ex.map(job, jobs))
    for s in st:
        if s["status"] != "ok":
            print(s["out"], s["status"], s.get("log_tail", [])[-5:])
    print(f"{'omk':>6s} {'TT':>9s} {'EE':>9s} {'TE':>9s} {'PP':>9s} {'P(k)':>9s} {'sigma8':>9s} {'H(z)':>9s}")
    for omk in ("0", "0.05", "-0.05"):
        r = harness.compare(OUT / f"A0_PEFTzero_omk{omk}.npz", OUT / f"B0_LCDM_omk{omk}.npz")
        t = r["cl"]["total"]
        print(f"{omk:>6s} {t['TT']:9.2e} {t['EE']:9.2e} {t['TE']:9.2e} {r['cl']['lens_potential']['PP']:9.2e} "
              f"{r['worst_pk']:9.2e} {r['sigma8']:9.2e} {r['bg_H']:9.2e}")
    print("expected: all rows <~1e-4 (GR-equivalent EFT model vs LCDM at the same cosmology)")


if __name__ == "__main__":
    main()
