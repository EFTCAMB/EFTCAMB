"""Model-level independent check of the cubic spline (flag 9) defect, via EFTCAMB's own timestep cache.

Pure EFT, Omega(a) = a. Reference: PureEFTmodelOmega=2 (linear, analytic). Test: PureEFTmodelOmega=9 with nodes on
exactly linear data. Reads the EFT functions EFTOmegaV/P/PP as stored by the model (get_eft_functions path, not the
standalone Fortran driver of the repro) and compares with the analytic Omega=a, Omega'=1, Omega''=0.

Cases (tree, nodes):
  A0 N=6  on [0,1]      -- the reported configuration
  A0 N=21 on [0,1]      -- denser nodes: an interpolation-accuracy problem would shrink, an end-node defect does not
  A0 N=10 on [-0.4,1.4] -- end nodes moved outside the a-range the code evaluates: should mostly cure it
  B1 N=6  on [0,1]      -- pre-rebase EFTCAMB: is the defect pre-existing?
usage: python verify_spline_cubic_model.py
"""

import json
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness  # noqa: E402

OUT = harness.AUDIT_ROOT / "work" / "verify" / "spline_cubic"
INI = harness.PARAMS_DIR / "2_PEFT_Omega_lin_1.ini"


def spline_ov(xs):
    ov = {"PureEFTmodelOmega": "9", "EFTOmega_Spline_Pixels": str(len(xs)), "EFTOmega_null_value": "0.0"}
    for i, x in enumerate(xs):
        ov[f"EFTOmegax{i + 1}"] = repr(float(x))
        ov[f"EFTOmegav{i + 1}"] = repr(float(x))
    return ov


CASES = {
    "A0_ref": ("A0", {}),
    "A0_N6": ("A0", spline_ov(np.linspace(0, 1, 6))),
    "A0_N21": ("A0", spline_ov(np.linspace(0, 1, 21))),
    "A0_N10ext": ("A0", spline_ov(np.linspace(-0.4, 1.4, 10))),
    "B1_ref": ("B1", {}),
    "B1_N6": ("B1", spline_ov(np.linspace(0, 1, 6))),
}


def run(name):
    tree, ov = CASES[name]
    out = OUT / f"{name}.npz"
    if out.exists() and json.loads(out.with_suffix(".json").read_text())["status"] == "ok":
        return name, {"status": "ok", "cached": True}
    return name, harness.run_model(tree, INI, out, overrides=ov, threads=1)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    with ThreadPoolExecutor(4) as ex:
        for name, st in ex.map(run, CASES):
            print(f"{name:10s} {st['status']}")
    probes = [0.5, 0.94, 0.99, 1.0]
    for name in CASES:
        d = harness.load(OUT / f"{name}.npz")
        a = d["eft_a_grid"]
        meta = d["meta"] if isinstance(d["meta"], dict) else json.loads(str(d["meta"]))
        print(f"\n== {name}: model={meta.get('eft_model_name')}  init_fallback={meta.get('eft_init_fallback')}")
        for p in probes:
            i = int(np.argmin(np.abs(a - p)))
            print(f"   a={a[i]:.4f}: OmegaV={d['eft_EFTOmegaV'][i]: .6f} (exp {a[i]:.4f})  "
                  f"OmegaP={d['eft_EFTOmegaP'][i]: .6f} (exp 1)  OmegaPP={d['eft_EFTOmegaPP'][i]: .3e} (exp 0)")
        m = (a > 1e-3)
        print(f"   max|OmegaV-a| (a>1e-3) = {np.max(np.abs(d['eft_EFTOmegaV'][m] - a[m])):.3e}   "
              f"max|OmegaP-1| = {np.max(np.abs(d['eft_EFTOmegaP'][m] - 1)):.3e}   "
              f"max|OmegaPP| = {np.max(np.abs(d['eft_EFTOmegaPP'][m])):.3e}")
        ref = OUT / f"{name.split('_')[0]}_ref.npz"
        if not name.endswith("ref"):
            c = harness.compare(OUT / f"{name}.npz", ref)
            print(f"   vs {ref.stem}: worst Cl {c['worst_cl']:.4e}  worst P(k) {c['worst_pk']:.4e}  "
                  f"sigma8 {c['sigma8']:.4e}")


if __name__ == "__main__":
    main()
