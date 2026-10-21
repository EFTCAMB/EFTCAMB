"""Breadth check for curvature-silently-wrong-with-eft: other EFT models (EFT perturbations active).

Double ratio  R(omk) = [Cl_EFT(omk) / Cl_CAMB(omk)] / [Cl_EFT(0) / Cl_CAMB(0)] - 1,  CAMB = B0 upstream LCDM.
For a correct curved implementation the curvature response of a near-LCDM model is close to LCDM's, so |R| should
be small (<~1e-2 for these models); the defect shows up as O(1-50) at l=2.
Usage: python curvsil_other_models.py
"""

import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import harness  # noqa: E402

REF = harness.AUDIT_ROOT / "work" / "verify" / "curvsil"
OUT = REF / "other"
P = harness.PARAMS_DIR
COMMON = {"num_redshiftwindows": "0", "get_tensor_cls": "F", "l_max_scalar": "2200"}
MODELS = {
    "fR_B0_1": (P / "4_Des_F_R_LCDM_1.ini", {}),
    "fR_B0_1e-4": (P / "4_Des_F_R_LCDM_1.ini", {"EFTB0": "1e-4"}),
    "quint_exp": (P / "5_quint_exponential.ini", {}),
}
OMKS = ("0", "0.02", "-0.02")


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    jobs = [(m, omk) for m in MODELS for omk in OMKS]

    def one(job):
        m, omk = job
        ini, ov = MODELS[m]
        f = OUT / f"{m}_omk{omk}.npz"
        if f.exists():
            return m, omk, "cached"
        r = harness.run_model("A0", ini, f, overrides=dict(COMMON, omk=omk, **ov), threads=1, timeout=900)
        return m, omk, r["status"]

    with ThreadPoolExecutor(4) as ex:
        st = list(ex.map(one, jobs))
    for s in st:
        print(*s)
    for m in MODELS:
        e0 = harness.load(OUT / f"{m}_omk0.npz")["cl_total"]
        c0 = harness.load(REF / "B0_GR_omk0.npz")["cl_total"]
        for omk in OMKS[1:]:
            f = OUT / f"{m}_omk{omk}.npz"
            if not f.exists():
                continue
            d = harness.load(f)
            e = d["cl_total"]
            c = harness.load(REF / f"B0_GR_omk{omk}.npz")["cl_total"]
            n = min(len(e), len(c), 2201)
            R = (e[2:n, 0] / c[2:n, 0]) / (e0[2:n, 0] / c0[2:n, 0]) - 1
            Re = (e[2:n, 1] / c[2:n, 1]) / (e0[2:n, 1] / c0[2:n, 1]) - 1
            print(f"{m:12s} omk={omk:6s} TT double ratio l=2 {R[0]:+.3e} l=10 {R[8]:+.3e} l=30 {R[28]:+.3e} "
                  f"max l>=30 {abs(R[28:]).max():.2e} | EE max {abs(Re).max():.2e} | meta errors {d['meta']['errors']}")


if __name__ == "__main__":
    main()
