"""EFTflag=0 with non-Lambda CAMB dark energy (fluid, PPF): code under audit (A0) vs plain CAMB v2 (B0).

EFTCAMB moves the dark-energy perturbation variables to the front of the scalar y vector (SetupScalarArrayIndices) and
wraps the DarkEnergy%PerturbedStressEnergy / PerturbationEvolve / diff_rhopi_Add_Term calls in derivs; this checks
that the CAMB dark-energy interface still gives CAMB's results (Reion.max_redshift aligned to 50, flat).

Usage: python de_interface_gr.py
"""

import json
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness  # noqa: E402

W = harness.AUDIT_ROOT / "work" / "t8a-equations" / "de_gr"
W.mkdir(parents=True, exist_ok=True)
base = harness.PARAMS_DIR / "base_params.ini"
CASES = {
    "fluid_w-0.9": {"dark_energy_model": "fluid", "w": "-0.9", "wa": "0", "cs2_lam": "1"},
    "fluid_w-0.9_cs0.1": {"dark_energy_model": "fluid", "w": "-0.9", "wa": "0", "cs2_lam": "0.1"},
    "ppf_w-0.9_wa0.3": {"dark_energy_model": "ppf", "w": "-0.9", "wa": "0.3"},
}
jobs = []
for name, ov in CASES.items():
    ov = dict(ov, EFTflag="0", num_redshiftwindows="0", get_tensor_cls="F")
    jobs.append(("A0", name, ov, {"Reion.max_redshift": 50.0}))
    jobs.append(("B0", name, ov, {}))


def run(j):
    tree, name, ov, ps = j
    return harness.run_model(tree, base, W / f"{tree}_{name}.npz", overrides=ov, python_set=ps, threads=1)


with ThreadPoolExecutor(4) as ex:
    st = list(ex.map(run, jobs))
res = {}
for name in CASES:
    a, b = W / f"A0_{name}.npz", W / f"B0_{name}.npz"
    if not (a.exists() and b.exists()):
        res[name] = [s["status"] for s in st]
        print(name, res[name])
        continue
    c = harness.compare(a, b)
    res[name] = {"worst_cl": c["worst_cl"], "worst_pk": c["worst_pk"], "bit": c["bit_identical"],
                 "unlensed_scalar": c["cl"].get("unlensed_scalar"), "w_A0": harness.load(a)["meta"].get("derived")}
    print(name, f"worst_cl={c['worst_cl']:.2e} worst_pk={c['worst_pk']:.2e}", c["cl"].get("unlensed_scalar"))
(W / "de_interface_gr.json").write_text(json.dumps(res, indent=1, default=str))
