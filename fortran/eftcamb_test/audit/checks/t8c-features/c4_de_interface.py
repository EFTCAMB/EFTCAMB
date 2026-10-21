"""Check 4: CAMB dark-energy model interface with EFT on (double counting?).

With EFTflag != 0 the EFT model owns the DE background and perturbations.  Setting a CAMB dark-energy model / w / wa
as well must not change anything.  Variants (all with the same EFT model):
  base      : model ini as is (fluid, w=-1)
  ini_fluid : ini override w=-0.7, wa=0.2, cs2_lam=0.5 (fluid)
  ini_ppf   : ini override dark_energy_model=ppf, w=-0.7, wa=0.2
  py_fluid  : python cp.DarkEnergy.set_params(w=-0.7, wa=0.2) after read_ini
  py_ppf    : python set_classes(dark_energy_model='ppf') + set_params(w=-0.7, wa=0.2)
  py_quint  : python set_classes(dark_energy_model='EarlyQuintessence') (num_perturb_equations=2 after Init in GR)
Observables compared to base: CMB Cls, linear P(k), H(z) -- expected bit-identical.
Also: GR (EFTflag=0) control showing that the same overrides do change GR results.
"""

from __future__ import annotations

import json
import sys

import numpy as np
from t8c_lib import WORK, cl_diffs, load, maxrel, model_ini, pk_diff, run_many

OUT = WORK / "c4"
MODELS = ["4_Des_F_R_wCDM_1", "2_PEFT_wDE_wCDM_1", "5_quint_exponential", "2_PEFT_gamma3_const_1", "1_EFT_GR"]
OV = {"get_tensor_cls": "F", "feedback_level": "0", "num_redshiftwindows": "0", "transfer_kmax": "5"}
OUTS = {"outputs": ["cmb", "pk_lin", "background"], "cl_lmax": 2500}
VARIANTS = {
    "base": ({}, {}),
    "ini_fluid": ({"w": "-0.7", "wa": "0.2", "cs2_lam": "0.5"}, {}),
    "ini_ppf": ({"dark_energy_model": "ppf", "w": "-0.7", "wa": "0.2"}, {}),
    "py_fluid": ({}, {"dark_energy": {"w": -0.7, "wa": 0.2, "cs2": 0.5}}),
    "py_ppf": ({}, {"dark_energy": {"model": "ppf", "w": -0.7, "wa": 0.2}}),
}


def jobs():
    js = []
    for m in MODELS:
        for tag, (ov, py) in VARIANTS.items():
            js.append(dict(tree="A0", ini=model_ini(m), task="generic", out=OUT / f"{m}_{tag}.npz",
                           overrides=OV | ov, args=OUTS | py))
    return js


def report():
    res = {}
    for m in MODELS:
        try:
            b = load(OUT / f"{m}_base.npz")
        except FileNotFoundError:
            res[m] = "missing"
            continue
        row = {}
        for tag in VARIANTS:
            if tag == "base":
                continue
            try:
                a = load(OUT / f"{m}_{tag}.npz")
            except FileNotFoundError:
                row[tag] = "missing"
                continue
            d = cl_diffs(a["cl_total"], b["cl_total"], lmax=2500)
            row[tag] = {"cl": max(d.values()), "pk": pk_diff(a, b), "npk": [a["pk"].shape[1], b["pk"].shape[1]], "H": maxrel(a["bg_H"], b["bg_H"]),
                        "de_rho": maxrel(a["bg_de_rho"], b["bg_de_rho"]),
                        "bit_identical": bool(np.array_equal(a["cl_total"], b["cl_total"])
                                              and np.array_equal(a["pk"], b["pk"])),
                        "DE": a["meta"]["DarkEnergy"]}
        res[m] = row
    print(json.dumps(res, indent=1))
    return res


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else "run"
    if cmd == "run":
        for r in run_many(jobs(), nproc=int(sys.argv[2]) if len(sys.argv) > 2 else 4):
            print(r["out"].split("/")[-1], r["status"], r["elapsed"])
    report()
