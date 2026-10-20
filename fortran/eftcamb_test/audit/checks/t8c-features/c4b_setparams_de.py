"""Check 4b: pure-Python camb.set_params path with EFT keys plus CAMB dark-energy keys.

designer f(R) wCDM (EFTw0=-0.9) through camb.set_params:
  base          : EFT keys only
  de_ppf_w07    : + dark_energy_model='ppf', w=-0.7, wa=0.2
  de_fluid_w07  : + w=-0.7 (default fluid)
  eftcamb_class : dark_energy_model='EFTCAMB' (set_classes branch) + w=-0.7
All must be bit-identical to base (EFT owns the DE sector).
"""

from __future__ import annotations

import json
import sys

import numpy as np
from t8c_lib import WORK, cl_diffs, load, pk_diff, run_many

OUT = WORK / "c4b"
EFT = {"EFTflag": 3, "DesignerEFTmodel": 1, "EFTwDE": 1, "EFTw0": -0.9, "EFTB0": 1.0,
       "EFT_ghost_stability": True, "EFT_gradient_stability": True, "EFTCAMB_turn_on_time": 0.01,
       "feedback_level": 0}
COSMO = {"H0": 70.0, "ombh2": 0.0226, "omch2": 0.112, "mnu": 0.06, "As": 2.1e-9, "ns": 0.96, "lmax": 2000,
         "lens_potential_accuracy": 1, "WantTransfer": True, "kmax": 2.0}
VARIANTS = {
    "base": {},
    "de_ppf_w07": {"dark_energy_model": "ppf", "w": -0.7, "wa": 0.2},
    "de_fluid_w07": {"w": -0.7},
    "eftcamb_class": {"dark_energy_model": "EFTCAMB", "w": -0.7},
}


def jobs():
    return [dict(tree="A0", ini=WORK / "dummy.ini", task="setparams", out=OUT / f"{tag}.npz",
                 args={"params": COSMO | EFT | v, "outputs": ["cmb", "pk_lin"], "cl_lmax": 2000})
            for tag, v in VARIANTS.items()]


def report():
    res = {}
    b = load(OUT / "base.npz")
    res["base_EFTflag"] = b["meta"]["EFTflag"]
    for tag in VARIANTS:
        if tag == "base":
            continue
        try:
            a = load(OUT / f"{tag}.npz")
        except FileNotFoundError:
            st = json.loads((OUT / f"{tag}.json").read_text())
            res[tag] = st["status"] + " " + str(st.get("log_tail", [])[-2:])
            continue
        res[tag] = {"cl": cl_diffs(a["cl_total"], b["cl_total"], lmax=2000), "pk": pk_diff(a, b),
                    "bit_identical": bool(np.array_equal(a["cl_total"], b["cl_total"])), "DE": a["meta"]["DarkEnergy"],
                    "EFTflag": a["meta"]["EFTflag"]}
    print(json.dumps(res, indent=1))
    return res


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else "run"
    if cmd == "run":
        (WORK / "dummy.ini").write_text("# unused\n")
        for r in run_many(jobs(), nproc=2):
            print(r["out"].split("/")[-1], r["status"], r["elapsed"])
    report()
