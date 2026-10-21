"""Verify t3-stability-dtauda-without-matter (part B): independent check of the mechanism.

Different code path from the stability check: initialise the model with EFTCAMB_skip_stability=T (background only,
no thermo), then evaluate the EFT functions via the python bridge (get_eft_functions -> compute_background_EFT_functions
+ compute_adotoa, which use grhom_t) and CAMB's own dtauda (results.hubble_parameter -> equations.f90 dtauda, which
sets grhoa2 = grho_no_de(a) = 8 pi G rho a^4). From these, rebuild what EFTCAMBFullModelComputeDtauda returns
  temp(grhoa2) = (grhoa2 + 2 c a^2 - Lambda a^2) / 3 / (1 + Omega + a Omega'),   dtauda = 1/sqrt(temp)
with grhoa2 = 0 (what EFTStabilityComputation leaves) and with grhoa2 = grhom_t a^2, and compare with
(a * adotoa)^2 and (a^2 H_CAMB)^2. Also print the grhoa2 field returned by the bridge vs grhom_t a^2.

Run: python verify_stab_dtauda_mechanism.py
"""

import json
import os
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import harness  # noqa: E402

REPO = HERE.parents[4]
PARAMS = REPO / "fortran" / "eftcamb_test" / "parameters"
INIS = REPO / "fortran" / "eftcamb_test" / "audit" / "inis"
OUT = harness.AUDIT_ROOT / "work" / "verify" / "stabdtauda" / "mechanism"
R5 = harness.AUDIT_ROOT / "work" / "t3-gr" / "R5"
A_EVAL = [1e-10, 1e-9, 1e-8, 1e-6, 1e-3, 1.0]
C_KMS = 299792.458

CASES = {
    "OLg_GR": (INIS / "007p5-OL-gamma.ini", {"OLOmegamodel": "0"}),
    "OLg_1e-6": (INIS / "007p5-OL-gamma.ini", {"OLOmega0": "1e-6"}),
    "OLg_1e-2": (INIS / "007p5-OL-gamma.ini", {"OLOmega0": "1e-2"}),
    "OLa_1e-6": (INIS / "007p6-OL-alpha.ini", {"OLmass0": "1e-6"}),
    "Kmou_GR": (PARAMS / "5_Kmouflage_1.ini", {"eps2_0": "0"}),
    "Kmim_GR": (PARAMS / "5_Kmimic_1.ini", {"eps2_0": "0"}),
}


def child(cfg: dict) -> None:
    sys.path.insert(0, cfg["tree"])
    import numpy as np

    import camb

    cp = camb.read_ini(cfg["ini"], no_validate=True)
    res = camb.get_background(cp, no_thermo=True)
    a = np.array(A_EVAL)
    _, arr = cp.EFTCAMB.get_eft_functions(res, a)
    h_kms = np.asarray(res.hubble_parameter(1.0 / a - 1.0))
    a2h_camb = a * a * h_kms / C_KMS  # a^2 H = 1/dtauda (conformal units, 1/Mpc)
    out = {"a": a.tolist(), "a2H_camb": a2h_camb.tolist()}
    for f in ("EFTc", "EFTLambda", "EFTOmegaV", "EFTOmegaP", "grhom_t", "adotoa", "grhoa2"):
        out[f] = np.asarray(arr[f]).tolist()
    print("JSON:" + json.dumps(out))


def run(tree: str, name: str, ini: Path, ov: dict) -> dict:
    OUT.mkdir(parents=True, exist_ok=True)
    wrapper = harness._wrapper_ini(ini, OUT, name, {"EFTCAMB_skip_stability": "T", **ov})
    cfg = {"tree": str(harness.tree_path(tree)), "ini": str(wrapper)}
    env = dict(os.environ, OMP_NUM_THREADS="1")
    env.pop("PYTHONPATH", None)
    with harness.slot(1):
        p = subprocess.run([sys.executable, __file__, "_child", json.dumps(cfg)], capture_output=True, text=True,
                           env=env, cwd=OUT, timeout=300)
    (OUT / f"{name}.log").write_text(p.stdout + p.stderr)
    for line in p.stdout.splitlines():
        if line.startswith("JSON:"):
            return json.loads(line[5:])
    return {"error": (p.stdout + p.stderr)[-600:]}


def main() -> None:
    import numpy as np

    allres = {}
    for tree, tag in (("A0", "A0"), (str(R5), "R5")):
        for c, (ini, ov) in CASES.items():
            if tag == "R5" and c not in ("OLg_GR",):
                continue
            r = run(tree, f"{tag}__{c}", ini, ov)
            allres[f"{tag}__{c}"] = r
            print(f"=== {tag} {c}")
            if "error" in r:
                print("  ERROR:", r["error"].strip().splitlines()[-1:])
                continue
            a = np.array(r["a"])
            c_, lam, om, omp = (np.array(r[k]) for k in ("EFTc", "EFTLambda", "EFTOmegaV", "EFTOmegaP"))
            grhom, adot, ga2 = np.array(r["grhom_t"]), np.array(r["adotoa"]), np.array(r["grhoa2"])
            den = 3.0 * (1.0 + om + a * omp)
            t0 = (2 * c_ * a**2 - lam * a**2) / den
            t1 = (grhom * a**2 + 2 * c_ * a**2 - lam * a**2) / den
            with np.errstate(all="ignore"):
                d0 = 1.0 / np.sqrt(t0)
            print(f"  {'a':>8s} {'c':>11s} {'Lambda':>11s} {'temp(grhoa2=0)':>15s} {'dtauda(0)':>10s} "
                  f"{'temp(grhom a2)':>15s} {'(a adotoa)^2':>13s} {'(a^2 H_CAMB)^2':>15s} {'bridge grhoa2':>14s} "
                  f"{'grhom_t a^2':>12s}")
            for i in range(len(a)):
                print(f"  {a[i]:8.0e} {c_[i]:11.3e} {lam[i]:11.3e} {t0[i]:15.4e} {d0[i]:10.3e} {t1[i]:15.4e} "
                      f"{(a[i] * adot[i])**2:13.4e} {r['a2H_camb'][i]**2:15.4e} {ga2[i]:14.4e} {grhom[i] * a[i]**2:12.4e}")
    (OUT / "mechanism.json").write_text(json.dumps(allres, indent=1))


if __name__ == "__main__":
    if len(sys.argv) > 2 and sys.argv[1] == "_child":
        child(json.loads(sys.argv[2]))
    else:
        main()
