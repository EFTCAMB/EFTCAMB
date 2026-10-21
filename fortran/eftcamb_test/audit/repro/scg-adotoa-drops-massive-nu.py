"""Repro: Scaling Cubic Galileon (EFTflag=4, FullMappingEFTmodel=6) drops massive neutrinos from its expansion rate.

EFTCAMBScalingCubicComputeAdotoa (fortran/eftcamb/08f_full_models/008p6_Scaling_Cubic.f90:1123) overwrites
eft_cache%grhom_t with (grhornomass+grhog)/a^2 + (grhob+grhoc)/a, i.e. without the massive-neutrino density, and
derives adotoa from it.  derivs (equations.f90, "EFTCAMB MOD START: compute adotoa") uses this adotoa for every
perturbation equation whenever EFTflag /= 0, while CAMB's dtauda (EFTCAMBFullModelComputeDtauda, used for tau(a),
distances, thermal history, results.hubble_parameter) uses grhoa2, which includes massive neutrinos.

Observed: adotoa from get_eft_functions (= the derivs value) vs a*H(z)/c from results.hubble_parameter.
Expected: agreement to ~1e-6 (as with massless neutrinos only, second case below).

Run: /Users/marco/Software/miniforge3/bin/python scg-adotoa-drops-massive-nu.py [TREE]
"""

import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

try:  # share the audit machine-wide CPU slots when the harness is available (falls back to no limit)
    sys.path.insert(1, str(Path(__file__).resolve().parents[1] / "checks"))
    from harness import slot
except Exception:  # noqa: BLE001
    from contextlib import nullcontext as slot

TREE = sys.argv[1] if len(sys.argv) > 1 else (
    "/private/tmp/claude-501/-Users-marco-IW-eftcamb-update-EFTCAMB/"
    "2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit/trees/A0")
PARAMS = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "parameters"))

WORKER = r"""
import sys, json
import numpy as np
sys.path.insert(0, sys.argv[1])
import camb
cp = camb.read_ini(sys.argv[2], no_validate=True)
res = camb.get_background(cp)
a = np.array([1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 0.1, 0.5, 1.0])
_, arr = cp.EFTCAMB.get_eft_functions(res, a)
Hcamb = np.asarray(res.hubble_parameter(1 / a - 1)) / (camb.constants.c / 1e3) * a
print(json.dumps({"a": a.tolist(), "adotoa": arr["adotoa"].tolist(), "camb": Hcamb.tolist(),
                  "grhom_t": arr["grhom_t"].tolist(), "grhoa2_a4": (arr["grhoa2"] * a**4).tolist(),
                  "omnuh2": cp.omnuh2, "num_nu_massive": int(cp.num_nu_massive)}))
"""


def run(extra: str) -> dict:
    tmp = tempfile.mkdtemp(prefix="scg_nu_")
    ini = os.path.join(tmp, "scg.ini")
    with open(ini, "w") as f:
        f.write(f"DEFAULT({PARAMS}/5_scg.ini)\n{extra}output_root = {tmp}/out\nEFTCAMB_write_background = F\n"
                "feedback_level = 0\nEFTCAMB_feedback_level = 0\n")
    w = os.path.join(tmp, "w.py")
    with open(w, "w") as f:
        f.write(WORKER)
    env = dict(os.environ, OMP_NUM_THREADS="1")
    with slot(1):
        out = subprocess.run([sys.executable, w, TREE, ini], capture_output=True, text=True, env=env, cwd=tmp,
                             preexec_fn=lambda: os.nice(10))
    line = [x for x in out.stdout.splitlines() if x.startswith("{")]
    if not line:
        raise RuntimeError(out.stdout[-2000:] + out.stderr[-2000:])
    return json.loads(line[-1])


def main() -> None:
    print("tree:", TREE)
    cases = [("5_scg.ini as shipped (omnuh2=0.01, 1 massive eigenstate)", ""),
             ("same, massless neutrinos only", "omnuh2 = 0\nmassive_neutrinos = 0\nmassless_neutrinos = 3.046\n"
              "nu_mass_eigenstates = 0\n"),
             ("same, Planck-like omnuh2=0.00064", "omnuh2 = 0.00064\n")]
    for title, extra in cases:
        r = run(extra)
        print(f"\n== {title}  (omnuh2={r['omnuh2']}, num_nu_massive={r['num_nu_massive']})")
        print("        a     adotoa(derivs)/aH_camb-1   grhom_t/(grhoa2 a^4)-1")
        for a, h, hc, g, g2 in zip(r["a"], r["adotoa"], r["camb"], r["grhom_t"], r["grhoa2_a4"]):
            print(f"   {a:8.1e}   {h / hc - 1:+.3e}               {g / g2 - 1:+.3e}")


if __name__ == "__main__":
    main()
